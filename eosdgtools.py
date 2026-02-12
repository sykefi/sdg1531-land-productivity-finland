import copy
import glob
import os
import re

import numpy as np
import pandas as pd
import psutil
import rasterio
import rioxarray
from dask.distributed import Client, LocalCluster
from scipy.stats import norm

COMPRESSION='LZW'
COMPRESSION='DEFLATE'

def get_local_cluster(args):
    """
    Create local Dask cluster from shared args
    """
    # Create a Dask cluster with specific memory and number of workers
    cluster = LocalCluster(
        n_workers=args.dask_n_workers,  #          # Number of worker processes
        threads_per_worker=args.dask_threads,  #   # Threads per worker
        memory_limit=args.dask_memlim,  # ,         # Limit memory per worker
    )
    client = Client(cluster)
    return client

def add_common_args(parser):
    """
    Shared input arguments for command line tools
    """
    # Common args for command line tools
    n_workers = max(os.cpu_count() - 1, 1)
    mem = psutil.virtual_memory()
    total_gb = mem.total / (1024 ** 3)
    worker_mem = (total_gb*0.9)/n_workers
    parser.add_argument(
        '--dask_n_workers', type=int, nargs='?', help='Dask workers', default=n_workers
    )
    parser.add_argument(
        '--dask_threads', type=int, nargs='?', help='Dask threads', default=1
    )
    parser.add_argument(
        '--dask_memlim', type=str, nargs='?', help='Dask worker memory in GB', default=f'{worker_mem}G'
    )
    parser.add_argument(
        '--blocksize', type=int, nargs='?', help='Processing block size', default=1024
    )
    parser.add_argument(
        '--nodata',
        type=float,
        nargs='?',
        help='Custom nodata value of inputs',
        default=0,
    )
    parser.add_argument(
        '--tag',
        nargs='?',
        type=str,
        help='tag for output files (z_[TAG].tif)',
        default='',
    )
    return parser

def get_tiffargs(blocksize=1024):
    """
    Tiff specs for writing the output data
    """
    writeargs = {
        'driver': 'GTiff',
        'compress': COMPRESSION,
        'tiled': True,
        'windowed': True,
        'blockxsize': blocksize,
        'blockysize': blocksize,
        'BIGTIFF': 'YES',
        'num_threads': 'all_cpus',
    }
    return writeargs


def classify_z_score(src, blocksize=1024):
    """
    Classify z score values, input can be either path to raster file or existing array
    """
    chunksz = (1, blocksize, blocksize)
    if isinstance(src, str):
        da = rioxarray.open_rasterio(src, chunks=chunksz, masked=True)
    else:
        da = src
    degraded = (da < -1.96).astype(np.uint8).rio.write_nodata(0)
    risk = np.logical_and(da >= -1.96, da < -1.28).astype(np.uint8).rio.write_nodata(0)
    nochange = np.logical_and(da >= -1.28, da <= 1.28).astype(np.uint8).rio.write_nodata(0)
    potential = np.logical_and(da > 1.28, da <= 1.96).astype(np.uint8).rio.write_nodata(0)
    improving = (da > 1.96).astype(np.uint8).rio.write_nodata(0)

    return degraded, risk, nochange, potential, improving


# Vectorized Mann Kendall test by Guangzhi XU https://numbersmithy.com/about-me/
# Modified to take into account missing (NaN) values in time series.
def countTies(x):
    """Count number of ties in rows of a 2D matrix
    Args:
        x (ndarray): 2d matrix.
    Returns:
        result (ndarray): 2d matrix with same shape as <x>. In each
            row, the number of ties are inserted at (not really) arbitary
            locations.
            The locations of tie numbers in are not important, since
            they will be subsequently put into a formula of sum(t*(t-1)*(2t+5)).

    Inspired by: https://stackoverflow.com/a/24892274/2005415.
    """
    if np.ndim(x) != 2:
        raise Exception('<x> should be 2D.')
    m, n = x.shape
    pad0 = np.zeros([m, 1]).astype('int')
    x = copy.deepcopy(x)
    x.sort(axis=1)
    diff = np.diff(x, axis=1)
    cated = np.concatenate([pad0, np.where(diff == 0, 1, 0), pad0], axis=1)
    absdiff = np.abs(np.diff(cated, axis=1))
    rows, cols = np.where(absdiff == 1)
    rows = rows.reshape(-1, 2)[:, 0]
    cols = cols.reshape(-1, 2)
    counts = np.diff(cols, axis=1) + 1
    result = np.zeros(x.shape).astype('int')
    result[rows, cols[:, 1]] = counts.flatten()
    return result


def MannKendallTrend2D(data, tails=2, axis=0, verbose=True):
    """Vectorized Mann-Kendall tests on 2D matrix rows/columns
    Args:
        data (ndarray): ndarray with shape (m, n).
    Keyword Args:
        tails (int): 1 for 1-tail, 2 for 2-tail test.
        axis (int): 0: test trend in each column. 1: test trend in each
            row.
    Returns:
        z (ndarray): If <axis> = 0, 1d array with length <n>, standard scores
            corresponding to data in each row in <x>.
            If <axis> = 1, 1d array with length <m>, standard scores
            corresponding to data in each column in <x>.
        p (ndarray): p-values corresponding to <z>.

        Original by: https://stackoverflow.com/a/24892274/2005415.
        Modified to handle missing data.

    """
    if np.ndim(data) != 2:
        raise Exception('<data> should be 2D.')
    # alway put records in rows and do M-K test on each row
    if axis == 0:
        data = data.T
    m, n = data.shape
    mask = np.triu_indices(n, 1)
    try:
        s = np.nansum(np.sign(data[:, mask[1]] - data[:, mask[0]]), axis=1)
    except MemoryError:
        print('Out of memory!')
        # data size too big, do it in batches
        i = 1
        s = np.zeros(m)
        while True:
            k = m // (n * i)
            t = m // k + 1
            try:
                for j in range(t):
                    id1 = j * k
                    id2 = min((j + 1) * k, m)
                    sj = np.nansum(
                        np.sign(data[id1:id2, mask[1]] - data[id1:id2, mask[0]]), axis=1
                    )
                    s[id1:id2] = sj
            except Exception as e:
                print(
                    f'Caught an error while processing in batches, adjusting batch size: {e}'
                )
                i += 1
            else:
                break
    # --------------------Count ties--------------------
    counts = countTies(data)
    tt = counts * (counts - 1) * (2 * counts + 5)
    tt = np.nansum(tt, axis=1)
    # -----------------Sample Gaussian-----------------
    # wroks when there are no nan's
    var = (n * (n - 1) * (2 * n + 5) - tt) / 18.0
    nn = np.count_nonzero(~np.isnan(data), axis=1)
    var = (nn * (nn - 1) * (2 * nn + 5) - tt) / 18.0
    eps = 1e-8  # avoid dividing 0
    z = (s - np.sign(s)) / (np.sqrt(var) + eps)
    p = norm.cdf(z)
    p = np.where(p > 0.5, 1 - p, p)
    if tails == 2:
        p = p * 2
    return z, p


def init_target(
    target_name, dscr, profile, dtype=rasterio.float32, nodata=np.nan, blocksize=1042
):
    """
    Initialize (Tiff) file for block-wise writing
    """
    # Open target geotiff for blockwise writing.
    dst_profile = profile.copy()
    dst_profile['blockxsize'] = blocksize
    dst_profile['blockysize'] = blocksize
    dst_profile['compress'] = COMPRESSION
    dst_profile['BIGTIFF'] = 'YES'
    dst_profile['tiled'] = True
    dst_profile['count'] = 1
    dst_profile['dtype'] = dtype
    dst_profile['nodata'] = nodata
    dst_profile['tiled'] = True,
    dst_profile['windowed'] = True,
    dst_profile['num_threads'] = 'all_cpus',
    # return handle to opened file
    dst_img = rasterio.open(target_name, 'w', **dst_profile)
    dst_img.set_band_description(1, dscr)

    return dst_img


def getndvifiles(srcpth, start_y, end_y):
    """
    Filter subset of input files from a directory containg yearly max(NDVI) geotiff's
    Assumes all files are in single dir and file names contain single 4-digit year time stamp '_YYYY'
    """

    # TODO: Add more general way for temporal filtering, now depends on file name format

    print(f'Listing files between {start_y} {end_y}')
    g = glob.glob(os.path.join(srcpth, '*.tif'))
    years = [int(re.search(r'_(\d{4})', s).groups(0)[0]) for s in g]
    df = pd.DataFrame({'year': years, 'file': g})
    df = df.sort_values(['year'])
    datatoprocess = list(df[(df.year >= start_y) & (df.year <= end_y)].file)
    print(f'Found {len(datatoprocess)} files')
    return datatoprocess
