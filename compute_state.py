# Productivity state indicator workflow
import argparse
import datetime as dt
import os

import numpy as np
import rioxarray
import xarray as xr

from eosdgtools import (
    add_common_args,
    classify_z_score,
    get_local_cluster,
    get_tiffargs,
    getndvifiles,
)

"""
Tools to classify land productivity based on productivity state 
"""


def parse_arguments(inargs=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute and classify z-scores for productivity state indicator',
    )
    parser.add_argument(
        'srcdir',
        type=str,
        help='Input directory path, folder with yearly max NDVI geotiffs',
        default='NDVI',
    )
    parser.add_argument(
        '--timespans',
        type=int,
        nargs=4,
        help='Processing years [start, end, ref_start, ref_end] ',
        metavar=('start', 'end', 'ref_start', 'ref_end'),
        default=[2000, 2012, 2013, 2015],
    )
    parser.add_argument(
        '--bbox',
        type=float,
        nargs=4,
        help='Optional BoundingBox for processing, for example 309590 6653761 456142 6752930 ',
        metavar=('minx', 'miny', 'maxx', 'maxy'),
        default=None,
    )
    parser.add_argument(
        '--odir',
        nargs='?',
        type=str,
        help='Output directory path',
        default='./ProductivityState',
    )

    parser.add_argument(
        '--export_all',
        action=argparse.BooleanOptionalAction,
        help='Save xhat, mu and sigma too',
        default=False,
    )

    parser = add_common_args(parser)

    if inargs is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(inargs)
    if not args.tag:
        args.tag = f'{args.timespans[0]}_{args.timespans[4]}'
    return args


def get_z(file_list_a, file_list_b, blocksize=1024, bbox=None, nodata=None):
    # Compute z-values for State indicator
    # file_list_a is the first time period and file list_b the second
    #
    # Nodata handling:
    # Masks out values defines as nodata in input geotiff metadata (all inputs are assumed equal)
    # Also allows user to manually set nodata for masking certain values as inputs may have wrong metadata.
    #
    # Here we assume all rasters have same shape and alignment and
    # load files as dask-backed DataArrays

    n_b = len(file_list_b)
    # Use same chunk size here and in rio.to_raster
    chunksz = (1, blocksize, blocksize)

    # Support on-the fly subsetting
    if bbox:
        (minx, miny, maxx, maxy) = bbox
    print('Files for first period')
    for f in file_list_a:
        print(f)
    data_arrays = [
        rioxarray.open_rasterio(fp, chunks=chunksz, mask_and_scale=True)
        for fp in file_list_a
    ]
    # Get nodata from meta
    internal_nodata = data_arrays[0].rio.nodata
    print(f'Dataset internal nodata: {internal_nodata}')
    # Stack along a new axis
    stacked = xr.concat(data_arrays, dim='time')
    stacked = stacked.where(stacked != internal_nodata)
    if bbox:
        stacked = stacked.rio.clip_box(minx=minx, miny=miny, maxx=maxx, maxy=maxy)
    if nodata not in (None,):
        print(f'Using custom nodata: {nodata}')
        stacked = stacked.where(stacked != nodata)

    # Compute pixel-wise mean and std along the 'time' axis
    stacked = (stacked - 100) / 100  # Re-scale NDVI values back to [-1,1]
    mu = stacked.mean(dim='time', skipna=True).astype(np.float32)
    sigma = stacked.std(dim='time', skipna=True).astype(np.float32)

    print('Files for second period')
    for f in file_list_b:
        print(f)
    data_arrays_b = [
        rioxarray.open_rasterio(fp, chunks=chunksz, mask_and_scale=True)
        for fp in file_list_b
    ]
    # Stack along a new axis
    stacked_b = xr.concat(data_arrays_b, dim='time')
    stacked_b = stacked_b.where(stacked_b != internal_nodata)
    if bbox:
        stacked_b = stacked_b.rio.clip_box(minx=minx, miny=miny, maxx=maxx, maxy=maxy)
    if nodata not in (None,):
        stacked_b = stacked_b.where(stacked_b != nodata)

    stacked_b = (stacked_b - 100) / 100  # Re-scale NDVI values back to [-1,1]
    xhat = stacked_b.mean(dim='time', skipna=True).astype(np.float32)

    # assure output datatype, this determines the dtype of output file automatically
    z = (xhat - mu) / (sigma / np.sqrt(n_b)).astype(np.float32)

    return z, xhat, mu, sigma


def main():
    """
    Compute productivity state indicator,

    Blocksize strongly affects performance and also memory consumption, adjust according to your enviroment.
    With large block size the outputs are not optimal for GUI usage, if needed consider reformatting data with
         gdaldal_translate z_baseline.tiff -of COG -co "num_threads=ALL_CPUS"  -co "BIGTIFF=YES" z_baseline_cog.tiff
    """

    # Default time periods
    # Baseline:
    # Pixelwise mean and std for years 2000-2012
    # Pixelwise mean for 2013-2015

    # Reporting:
    # Pixelwise mean and std for years 2004-2016
    # Pixelwise mean for 2017-2019

    args = parse_arguments()

    srcpth = args.srcdir
    if args.bbox:
        print(f'Got BBOX: {args.bbox}')
    else:
        print('Processing complete area of inputs')

    # geotiff specifications
    writeargs = get_tiffargs(args.blocksize)

    client = get_local_cluster(args)

    start0 = dt.datetime.now()

    # Lazy variable, data will not be processed yet
    z_score, xhat, mu, sigma = get_z(
        getndvifiles(srcpth, args.timespans[0], args.timespans[1]),
        getndvifiles(srcpth, args.timespans[2], args.timespans[3]),
        blocksize=args.blocksize,
        nodata=args.nodata,
        bbox=args.bbox,
    )
    # Write to disk
    os.makedirs(args.odir, exist_ok=True)

    # Write z-score to disk
    zpth = os.path.join(args.odir, f'state_zscore_{args.tag}.tif')
    f1 = z_score.rio.to_raster(zpth, lock=True, **writeargs)
    tocompute = [f1]

    # Optionally export other variables for debugging
    if args.export_all:
        f12 = xhat.rio.to_raster(
            os.path.join(args.odir, f'state_xhat_{args.tag}.tif'),  lock=True, **writeargs
        )
        f13 = mu.rio.to_raster(
            os.path.join(args.odir, f'state_mu_{args.tag}.tif'),  lock=True, **writeargs
        )
        f14 = sigma.rio.to_raster(
            os.path.join(args.odir, f'state_sigma_{args.tag}.tif'),  lock=True, **writeargs
        )
        tocompute = [f1, f12, f13, f14]
    futures = client.compute(tocompute)
    client.gather(futures)
    print(f'Z-score processed in {dt.datetime.now() - start0}')

    # Store classes too
    start = dt.datetime.now()
    degraded, risk, nochange, potential, improving = classify_z_score(zpth)
    f2 = degraded.rio.to_raster(
        os.path.join(args.odir, f'state_degraded_{args.tag}.tif'),  lock=True, **writeargs
    )
    f3 = risk.rio.to_raster(
        os.path.join(args.odir, f'state_risk_{args.tag}.tif'),  lock=True, **writeargs
    )
    f4 = nochange.rio.to_raster(
        os.path.join(args.odir, f'state_nochange_{args.tag}.tif'),  lock=True, **writeargs
    )
    f5 = potential.rio.to_raster(
        os.path.join(args.odir, f'state_potential_{args.tag}.tif'),  lock=True, **writeargs
    )
    f6 = improving.rio.to_raster(
        os.path.join(args.odir, f'state_improving_{args.tag}.tif'),  lock=True, **writeargs
    )
    futures2 = client.compute(f2, f3, f4, f5, f6)
    client.gather(futures2)
    print(f'Classes processed in {dt.datetime.now() - start}')

    print(f'Total processing time {dt.datetime.now() - start0}')
    return


if __name__ == '__main__':
    main()
