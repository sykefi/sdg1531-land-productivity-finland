# Productivity state indicator workflow
import argparse
import datetime as dt
import os

import numpy as np
import rioxarray
import xarray as xr

from eosdgtools import add_common_args, get_local_cluster, get_tiffargs, getndvifiles


def parse_arguments(inargs=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Compute mean of yearly NDVI_max arrays',
    )
    parser.add_argument(
        'srcdir',
        type=str,
        help='Input directory path',
        default='./NDVI',
    )
    parser.add_argument(
        '--timespan',
        type=int,
        nargs=2,
        help='min and max years',
        metavar=('start', 'end'),
        default=[2000, 2012],
    )
    parser.add_argument(
        '--bbox',
        type=float,
        nargs=4,
        help='BoundingBox for reporting, for example 309590 6653761 456142 6752930 ',
        metavar=('minx', 'miny', 'maxx', 'maxy'),
        default=None,
    )
    parser.add_argument(
        '--odir',
        nargs='?',
        type=str,
        help='Output directory path',
        default='./MeanNDVI',
    )


    parser = add_common_args(parser)
    if inargs is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(inargs)
    if not args.tag:
        args.tag = f'{args.timespan[0]}_{args.timespan[1]}'
    return args


def get_mean(file_list, blocksize=1024, bbox=None, nodata=None):
    # Compute pixelwise mean for given datasets

    # Use same chunk size here and in rio.to_raster
    chunksz = (1, blocksize, blocksize)

    # Support on-the fly subsetting
    if bbox:
        (minx, miny, maxx, maxy) = bbox
    print('Files to process')
    for f in file_list:
        print(f)
    data_arrays = [
        rioxarray.open_rasterio(fp, chunks=chunksz, mask_and_scale=True)
        for fp in file_list
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
    mu = stacked.mean(dim='time', skipna=True).astype(np.uint8)
    return mu


def main():
    """
    Compute pixel-wise mean for a list of large rasters

    """

    args = parse_arguments()

    srcpth = args.srcdir
    if args.bbox:
        print(f'Got BBOX: {args.bbox}')
    else:
        print('Processing complete area of inputs')

    client = get_local_cluster(args)

    start0 = dt.datetime.now()
    # Lazy datasets
    mean_ndvi_max = get_mean(
        getndvifiles(srcpth, args.timespan[0], args.timespan[1]),
        blocksize=args.blocksize,
        nodata=args.nodata,
        bbox=args.bbox,
    )

    # Write to files
    os.makedirs(args.odir, exist_ok=True)
    writeargs = get_tiffargs(args.blocksize)
    outpth = os.path.join(args.odir, f'mean_ndvi_max_{args.tag}.tif')
    # unscaled data in float32
    f1 = mean_ndvi_max.rio.to_raster(outpth, lock = True, **writeargs)
    client.compute([f1])

    print(f'Processed in {dt.datetime.now() - start0}')
    print(f'Total processing time {dt.datetime.now() - start0}')
    return


if __name__ == '__main__':
    main()
