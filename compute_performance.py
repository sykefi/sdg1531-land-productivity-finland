## Productivity performace indicator workflow

# 1. Compute pixelwise mean of ndvi_max time series for the period for the whole region
# 2. For each applicable corine class, compute p90 of mean(ndvi_max) of the pixels within the class
# 3. classify pixels in the class as degraded if ndvi<0.5*p90

# NOTES:
# We can either process whole Finland as one region, meaning that we consider all availble pixels when computing p90
# Or we can divide Finland into different biogeographical regions and do the classification region-wise


import argparse
import json
import os

import dask.array as da
import numpy as np
import rioxarray
import xarray as xr

from eosdgtools import add_common_args, get_local_cluster, get_tiffargs


def get_p90_preclipped(ndvi, clc, classids, blocksize=1024, nodata=None):
    # Find p90 value of the input ndvi raster pixel values for each corine class separately
    # Assumes that 
    # * the input has been scaled to 0-200 and does the re-scaling back to -1,1
    # * Input arrays are masked so that pixels outside the biogeographical region are nan's
    
    ds = rioxarray.open_rasterio(ndvi, chunks={'x': blocksize, 'y': blocksize})
    clc_ds = rioxarray.open_rasterio(clc, chunks={'x': blocksize, 'y': blocksize})
    if nodata not in (None,):
        print(f'Using custom nodata: {nodata}')
        ds = ds.where(ds != nodata)

    p90 = []
    for classid in classids:
        clcmask = clc_ds == classid # May be an empty mask, then result is nan
        flat = ds.data.reshape(-1)
        valid = flat[clcmask.data.reshape(-1) & ~da.isnan(flat)]
        thisp90 = da.percentile(
            (valid - 100) / 100, 90, internal_method='tdigest'
        )  # NOTE: default method yields significantly different results
        p90.append(thisp90)
    return p90


def classify_preclipped(ndvi, clc, classids, p90s, blocksize=1024, nodata=None):
    """
    ndvi : precomputed temporal mean of NDVI_max  over the time period
    clc: path to integer array of same size as the geotiff files (land cover class)
    classid : Corine class to process

    Find the pixels in each Corine class that are smaller tha 0.5*p90 of the current class
    Assumes that the input has been scaled to 0-200 and does the re-scaling back to -1,1
    """
    meanndvi = rioxarray.open_rasterio(ndvi, chunks={'x': blocksize, 'y': blocksize})
    if nodata not in (None,):
        print(f'Using custom nodata: {nodata}')
        meanndvi = meanndvi.where(meanndvi != nodata)
    corine = rioxarray.open_rasterio(clc, chunks={'x': blocksize, 'y': blocksize})
    clc_masks = [corine == val for val in classids]
    # run everything with scaled ndvi
    meanndvi = (meanndvi - 100) / 100
    degradeds = [
        meanndvi.where(mask) < (0.5 * p90) for mask, p90 in zip(clc_masks, p90s)
    ]
    degraded = xr.concat(degradeds, dim='band').any(dim='band')
    return degraded.astype(np.uint8)

def savep90(p90results, clcids, target):
    # Write dask results to json list, 
    # empty values are replaced with nan
    prctl = {}
    p90values = [np.nan] * len(p90results)
    for ind, f in enumerate(p90results):
        res = p90results[ind].result()
        if len(res) == 1:
            p90values[ind] = res[0]
        else:
            print(f'Invalid class for region {ind}')

    for ind, v in zip(clcids, p90values):
        prctl[int(ind)] = float(v)

    json_object = json.dumps(prctl, indent=4)
    with open(target, 'w') as outfile:
        outfile.write(json_object)

    return prctl


def loadp90(src):
    # Get p90 value list from json file
    with open(src) as f:
        d = json.loads(f.read())
    return d


def compute(srcdir, outdir, client, skipP90=True, skipFinal=True, blocksize = 1025,
            clcids = [1, 2, 3, 4, 5, 6, 7, 8], clc_y_baseline = 2012, clc_y_reporting = 2018, biogeozoneids = [1, 2, 3, 4, 5], nodata=None):
    """
    
    """
    clc_template = '{jakonro}_clc{clcyear}peruutus_recoded_sdg1531.tif'

    writeargs = get_tiffargs(blocksize)

    outdirp90 = os.path.join(outdir, 'p90')

    if not skipP90:  # get p90 values for each zone & corine class and save to disk
        os.makedirs(outdirp90, exist_ok=True)
        for jakonro in biogeozoneids:
            thisbaselinecorine = os.path.join(
                srcdir, clc_template.format(jakonro=jakonro, clcyear=clc_y_baseline)
            )
            thisreportingcorine = os.path.join(
                srcdir, clc_template.format(jakonro=jakonro, clcyear=clc_y_reporting)
            )
            # baseline
            thisndvi_base = os.path.join(
                srcdir, f'{jakonro}_mean_ndvi_max_baseline.tif'
            )
            p, bn_base = os.path.split(thisndvi_base)
            p90_base = get_p90_preclipped(
                thisndvi_base, thisbaselinecorine, clcids, blocksize=blocksize, nodata=nodata
            )
            p90_base = client.compute(p90_base)
            tmp = savep90(
                p90_base, clcids, os.path.join(outdirp90, f'p90_{bn_base}.json')
            )
            # reporting
            thisndvi_rep = os.path.join(
                srcdir, f'{jakonro}_mean_ndvi_max_reporting.tif'
            )
            p, bn_rep = os.path.split(thisndvi_rep)
            p90_rep = get_p90_preclipped(
                thisndvi_rep, thisreportingcorine, clcids, blocksize=blocksize, nodata=nodata
            )
            p90_rep = client.compute(p90_rep)
            savep90(p90_rep, clcids, os.path.join(outdirp90, f'p90_{bn_rep}.json'))

    if not skipFinal: # Do the classification using the precomputed p90 values
        os.makedirs(os.path.join(outdir, 'baseline'), exist_ok=True)
        os.makedirs(os.path.join(outdir, 'reporting'), exist_ok=True)
        for jakonro in biogeozoneids:
            thisbaselinecorine = os.path.join(
                srcdir, clc_template.format(jakonro=jakonro, clcyear=clc_y_baseline)
            )
            thisreportingcorine = os.path.join(
                srcdir, clc_template.format(jakonro=jakonro, clcyear=clc_y_reporting)
            )
            thisndvi_base = os.path.join(
                srcdir, f'{jakonro}_mean_ndvi_max_baseline.tif'
            )
            p, bn_base = os.path.split(thisndvi_base)
            prctl_base = loadp90(
                os.path.join(
                    outdirp90, f'p90_{jakonro}_mean_ndvi_max_baseline.tif.json'
                )
            )
            p90_base = [prctl_base[str(ind)] for ind in clcids]
            degraded_base = classify_preclipped(
                thisndvi_base,
                thisbaselinecorine,
                clcids,
                p90_base,
                blocksize=blocksize,
                nodata=nodata,
            )
            degraded_base.rio.write_nodata(0).rio.to_raster(
                os.path.join(outdir, 'baseline', f'degraded_{bn_base}'),
                lock=True,
                **writeargs,
            )
            # reporting
            thisndvi_rep = os.path.join(
                srcdir, f'{jakonro}_mean_ndvi_max_reporting.tif'
            )
            p, bn_rep = os.path.split(thisndvi_rep)
            prctl_rep = loadp90(
                os.path.join(
                    outdirp90, f'p90_{jakonro}_mean_ndvi_max_reporting.tif.json'
                )
            )
            p90_rep = [prctl_rep[str(ind)] for ind in clcids]
            degraded_rep = classify_preclipped(
                thisndvi_rep,
                thisreportingcorine,
                clcids,
                p90_rep,
                blocksize=blocksize,
                nodata=nodata,
            )
            degraded_rep.rio.write_nodata(0).rio.to_raster(
                os.path.join(outdir, 'reporting', f'degraded_{bn_rep}'),
                lock=True,
                **writeargs,
            )


def main(inargs=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Find degraded pixels by CLC class for the performance indicator',
    )
    parser.add_argument(
        'srcdir', type=str, help='Input data path (Folder containing ndvi and corine subsets by biogeographical classes)'
    )
    parser.add_argument('outdir', type=str, help='output directory')
    parser.add_argument(
        '--skipP90',
        action=argparse.BooleanOptionalAction,
        help='Skip p90 computation',
        default=False,
    )
    parser.add_argument(
        '--skipFinal',
        action=argparse.BooleanOptionalAction,
        help='Skip final raster computation (generate p90 only)',
        default=False,
    )

    parser.add_argument(
        '--geomids', type=int, nargs='*', help='Geometry id\'s', default=[1,2,3,4,5]
    )

    parser = add_common_args(parser)

    if inargs is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(inargs)

    # Create a Dask cluster with specific memory and number of workers
    client = get_local_cluster(args)


    os.makedirs(args.outdir, exist_ok=True)
    
    clcids = [1, 2, 3, 4, 5, 6, 7, 8] # Re-sclassified Corine classes
    biogeozoneids = args.geomids # linear index for zone subsets
    clc_y_baseline = 2012 # Corine to use for baseline
    clc_y_reporting = 2018 # Corine to use for reporting

    compute(args.srcdir, args.outdir, client, args.skipP90, args.skipFinal, blocksize = args.blocksize,
            clcids = clcids, clc_y_baseline = clc_y_baseline, clc_y_reporting = clc_y_reporting, 
            biogeozoneids = biogeozoneids, nodata=args.nodata)


if __name__ == '__main__':
    main()
