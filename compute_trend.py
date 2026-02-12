# Productivity state indicator workflow
import argparse
import concurrent.futures
import datetime as dt
import functools
import os

import numpy as np
import rasterio
from tqdm import tqdm

from eosdgtools import (
    MannKendallTrend2D,
    add_common_args,
    classify_z_score,
    get_local_cluster,
    get_tiffargs,
    getndvifiles,
    init_target,
)


def parse_arguments(inargs=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Mann Kendall Test for NDVI raster time series data',
    )
    parser.add_argument('srcdir', type=str, help='Input directory path', default='./NDVI')
    parser.add_argument('startyear', type=int, help='Start Year')
    parser.add_argument('endyear', type=int, help='End Year')
    parser.add_argument(
        '--odir',
        nargs='?',
        type=str,
        help='Output directory path',
        default='./ProductivityTrend',
    )
    n_workers = max(os.cpu_count() - 1, 1)
    parser.add_argument(
        '--n_workers', type=int, nargs='?', help='Number of processes for multiprocessing', default=n_workers
    )

    parser = add_common_args(parser)

    if inargs is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(inargs)

    if not args.tag:
        args.tag = f'{args.startyear}_{args.endyear}'
    return args

def processall(window, imglist=[], z_dst=None, p_dst=None, nodata=0):
    # Compute MannKendall test for a block (window) of input rasters and write output directly to disk

    # Initialize an array to store one spatial block of raster data files
    raster_data = []

    # Read each file and append its data to the raster_data list
    # MKTest needs data as Float

    for file in imglist:
        with rasterio.open(file) as src:
            data = np.float32(src.read(1, window=window, masked=True))
            data[data == nodata] = np.nan
            raster_data.append((data - 100) / 100)  # re-scale NDVI back to [-1,1]

    arr = np.array(raster_data)
    shp = arr.shape
    # dimensions for MannKendallTrend2D should be (time, pixel)
    arr.resize(shp[0], shp[1] * shp[2])
    z, p = MannKendallTrend2D(arr)
    # Resize results back to 2d
    z.resize(shp[1], shp[2])
    p.resize(shp[1], shp[2])
    # Write z and p blocks to disk
    z_dst.write_band(1, z, window=window)
    p_dst.write_band(1, p, window=window)


def main():
    """
    Compute productivity trend using MannKendall test
    """

    args = parse_arguments()
    srcpth = args.srcdir
    odir = args.odir
    # Time span to use
    start_y = args.startyear
    end_y = args.endyear

    start0 = dt.datetime.now()

    # These affect performance and memory consumption
    blocksize = args.blocksize
    n_workers = args.n_workers  # Make smaller if you encounter memory problems

    datatoprocess = getndvifiles(srcpth, start_y, end_y)
    print('Files to process')
    for f in datatoprocess:
        print(f)
        
    # TODO: On the fly subsetting as with xarray/dask workflows.
    # Currently only full tiff file processing is supported
    # if args.bbox:
    #    (minx, miny,maxx, maxy) = args.bbox

    # Open a sample file to get metadata and nodata value,
    # we assume that all inputs are equal in size/georeference
    with rasterio.open(datatoprocess[0]) as src:
        inputprofile = src.profile

    os.makedirs(odir, exist_ok=True)

    # Initialize output files, BigTiff should be enabled
    zscoretarget = os.path.join(odir, f'trend_zscore_{args.tag}.tif')
    z_dst = init_target(
        zscoretarget,
        'Z',
        inputprofile,
        blocksize=blocksize,
    )
    p_dst = init_target(
        os.path.join(odir, f'trend_p_{args.tag}.tif'),
        'p-value',
        inputprofile,
        blocksize=blocksize,
    )

    # get native blocks from output file, this determines the processing block size
    windows = [window for ij, window in z_dst.block_windows()]

    # Use onebyone for testing to get exceptions and resource estimates
    onebyone = False

    with rasterio.Env() as env:
        if onebyone:
            print('Running one-by-one')
            progressbar = tqdm(total=len(windows))
            for window in windows:
                processall(window, datatoprocess, z_dst, p_dst, nodata=np.nan)
                progressbar.update(1)
        else:
            print('Using multiprocessing')
            progressbar = tqdm(total=len(windows))
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=n_workers
            ) as executor:
                futures = [
                    executor.submit(
                        functools.partial(
                            processall,
                            imglist=datatoprocess,
                            z_dst=z_dst,
                            p_dst=p_dst,
                            nodata=args.nodata,
                        ),
                        ii,
                    )
                    for ii in windows
                ]
                for future in concurrent.futures.as_completed(futures):
                    progressbar.update(1)

    print('Calculation done!')

    z_dst.close()
    p_dst.close()
    print('Files closed.')
    print("Running classification ")

    start1 = dt.datetime.now()
    writeargs = get_tiffargs(args.blocksize)
    client = get_local_cluster(args)

    degraded, risk, nochange, potential, improving = classify_z_score(zscoretarget)
    
    f2 = degraded.rio.to_raster(
        os.path.join(args.odir, f'trend_degraded_{args.tag}.tif'), **writeargs, lock=True
    )
    f3 = risk.rio.to_raster(
        os.path.join(args.odir, f'trend_risk_{args.tag}.tif'), **writeargs, lock=True
    )
    f4 = nochange.rio.to_raster(
        os.path.join(args.odir, f'trend_nochange_{args.tag}.tif'), **writeargs, lock=True
    )
    f5 = potential.rio.to_raster(
        os.path.join(args.odir, f'trend_potential_{args.tag}.tif'), **writeargs, lock=True
    )
    f6 = improving.rio.to_raster(
        os.path.join(args.odir, f'trend_improving_{args.tag}.tif'), **writeargs, lock=True
    )
    client.compute(f2, f3, f4, f5, f6)
    print(f'Classes processed in {dt.datetime.now() - start1}')

    print(f'Total processing time {dt.datetime.now() - start0}')
    return


if __name__ == '__main__':
    main()
