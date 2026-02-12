#!/usr/bin/env python3
"""
Mask multiple data files and export them as GeoTIFFs.

This script replaces the following GDAL operations:
- Creating a mask from a baseline raster
- Rasterizing shapefiles into the mask
- Applying the mask to multiple input files
"""

import argparse
import os
import sys

import numpy as np
from osgeo import gdal

# Enable GDAL exceptions
gdal.UseExceptions()


def create_base_mask(baseline_file, output_mask, nodata_threshold=-9999):
    """
    Create a base mask from a baseline raster.
    Uses the logic: mask = (A < nodata_threshold)
    
    Args:
        baseline_file: Input baseline raster file
        output_mask: Output mask file path
        nodata_threshold: Threshold for nodata values
    """
    print(f"Creating base mask from: {baseline_file}")
    
    # Open the baseline file
    src_ds = gdal.Open(baseline_file, gdal.GA_ReadOnly)
    if src_ds is None:
        raise RuntimeError(f"Failed to open: {baseline_file}")
    
    # Read the data
    band = src_ds.GetRasterBand(1)
    data = band.ReadAsArray()
    
    # Create mask: True where A < nodata_threshold
    mask = (data < nodata_threshold).astype(np.uint8)
    
    # Create output file
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(
        output_mask,
        src_ds.RasterXSize,
        src_ds.RasterYSize,
        1,
        gdal.GDT_Byte,
        options=['TILED=YES', 'COMPRESS=DEFLATE']
    )
    
    if out_ds is None:
        raise RuntimeError(f"Failed to create: {output_mask}")
    
    # Copy geotransform and projection
    out_ds.SetGeoTransform(src_ds.GetGeoTransform())
    out_ds.SetProjection(src_ds.GetProjection())
    
    # Write mask data
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(mask)
    out_band.FlushCache()
    
    # Close datasets
    src_ds = None
    out_ds = None
    
    print(f"Base mask created: {output_mask}")


def rasterize_shapefile(shapefile, target_raster, burn_value=255, inverse=False):
    """
    Rasterize a shapefile into an existing raster.
    
    Args:
        shapefile: Input shapefile path
        target_raster: Target raster to burn values into
        burn_value: Value to burn into the raster
        inverse: If True, burn outside the geometries (inverse mode)
    """
    print(f"Rasterizing: {shapefile} {'(inverse)' if inverse else ''}")
    
    # Open the target raster for update
    target_ds = gdal.Open(target_raster, gdal.GA_Update)
    if target_ds is None:
        raise RuntimeError(f"Failed to open for update: {target_raster}")
    
    # Rasterize options
    options = {'burnValues':[burn_value]}
    if inverse:
        options = {'inverse': True}

    # Rasterize
    gdal.Rasterize(
        target_ds,
        shapefile,  # Band 1
        **options
    )
    
    # Close datasets
    target_ds = None
    
    print("Rasterized successfully")


def main():
    parser = argparse.ArgumentParser(
        description='Create template mask',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create template mask
  %(prog)s \\
    --input-template /path/to/degraded_baseline.tif \\
    --output-mask /path/to/mask.tif \\
    --water-shapefiles /path/to/meri10.shp /path/to/Jarvi10.shp /path/to/JokiAlue10.shp \\
    --border-shapefile /path/to/SuomenValtakunta_2025_10k.shp \\
        """
    )
    
    # Mask creation arguments
    parser.add_argument('--input-template', required=True,
                        help='Template raster file for creating mask')
    parser.add_argument('--output-mask', required=True,
                        help='Output mask file path')
    parser.add_argument('--water-shapefiles', nargs='+',
                        help='Water body shapefiles to burn into mask')
    parser.add_argument('--border-shapefile',
                        help='Border shapefile to burn inversely into mask')
    parser.add_argument('--nodata-threshold', type=float, default=-9999,
                        help='Threshold for nodata values in baseline (default: -9999)')
    args = parser.parse_args()
    
    try:
        # Step 1: Create base mask
        create_base_mask(args.input_template, args.output_mask, args.nodata_threshold)
        
        # Step 2: Rasterize water body shapefiles
        if args.water_shapefiles:
            for shapefile in args.water_shapefiles:
                if os.path.exists(shapefile):
                    rasterize_shapefile(shapefile, args.output_mask, burn_value=255, inverse=False)
                else:
                    print(f"Warning: Shapefile not found: {shapefile}", file=sys.stderr)
        
        # Step 3: Rasterize border shapefile (inverse)
        if args.border_shapefile:
            if os.path.exists(args.border_shapefile):
                rasterize_shapefile(args.border_shapefile, args.output_mask, burn_value=255, inverse=True)
            else:
                print(f"Warning: Border shapefile not found: {args.border_shapefile}", file=sys.stderr)
        
        print(f"\nMask creation completed: {args.output_mask}")
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
