import argparse
import concurrent.futures
import functools
import os

from osgeo import gdal
from tqdm import tqdm

os.environ['GDAL_NUM_THREADS'] = 'ALL_CPUS'


def clipinputs(pjako, srcfile, polygons, outdir):
    # Prepare data, clip with zones
    print(f'Clipping input {srcfile} {pjako}')
    sqlfilter = f'SELECT * FROM BiogeographicalZones WHERE paajakonro={pjako}'

    os.makedirs(outdir, exist_ok=True)
    pth, bn = os.path.split(srcfile)
    gdal.Warp(
        os.path.join(outdir, f'{pjako}_{bn}'),
        srcfile,
        cutlineDSName=polygons,
        cutlineSQL=sqlfilter,
        cropToCutline=True, # When processing areas that are smaller than biogeographical zone, this should be False
        format='GTiff',
        creationOptions=[
            'TILED=YES',  # Enable internal tiling
            'COMPRESS=LZW',  # Apply LZW compression
            'BIGTIFF=YES',  # Allow BigTIFF if file > 4GB
            'NUM_THREADS=ALL_CPUS',  # Explicitly tell GDAL to use all cores
        ],
        multithread=True,
    )
    print(f'PääjakoNro {pjako} done!')


def main(inargs=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Subset NDVI based on biogeographical zones',
    )
    parser.add_argument('src', type=str, help='Input raster')
    parser.add_argument(
        'zones',
        type=str,
        help='Biogeographical zones',
        default='BiogeographicalZones.shp',
    )
    parser.add_argument('output', type=str, help='Target path')
    parser.add_argument(
        '--n_workers', type=int, nargs='?', help='Parallellization', default=5
    )
    parser.add_argument(
        '--geomids', type=int, nargs='*', help='Geometry id\'s', default=[1,2,3,4,5]
    )

    if inargs is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(inargs)

    pjakonro =  args.geomids
    print (pjakonro)
    print('Using multiprocessing')
    progressbar = tqdm(total=len(pjakonro))
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=args.n_workers
    ) as executor:
        futures = [
            executor.submit(
                functools.partial(
                    clipinputs, 
                    srcfile= args.src, 
                    polygons = args.zones, 
                    outdir = args.output
                ),
                ii,
            )
            for ii in pjakonro
        ]
        for future in concurrent.futures.as_completed(futures):
            progressbar.update(1)


if __name__ == '__main__':
    main()
