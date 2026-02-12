#!/usr/bin/env python3
"""
Mask multiple data files and export them as GeoTIFFs.

This script replaces the following GDAL operations:
- Creating a mask from a baseline raster
- Rasterizing shapefiles into the mask
- Applying the mask to multiple input files
"""

import argparse
import glob
import os
import sys

from osgeo import gdal

# Enable GDAL exceptions
gdal.UseExceptions()


def apply_mask(mask_file, input_file, output_file, nodata_value=0):
    """
    Apply mask to an input file and save the result.
    Uses the logic: output = (mask != 255) * input
    
    Args:
        mask_file: Mask raster file
        input_file: Input raster file
        output_file: Output raster file
        nodata_value: NoData value for output
    """
    print(f"Applying mask to: {os.path.basename(input_file)}")
    
    # Open mask and input files
    mask_ds = gdal.Open(mask_file, gdal.GA_ReadOnly)
    input_ds = gdal.Open(input_file, gdal.GA_ReadOnly)
    
    if mask_ds is None or input_ds is None:
        raise RuntimeError("Failed to open mask or input file")
    
    # Read data
    mask_data = mask_ds.GetRasterBand(1).ReadAsArray()
    input_data = input_ds.GetRasterBand(1).ReadAsArray()
    
    # Apply mask: (mask != 255) * input
    output_data = ((mask_data != 255) * input_data).astype(input_data.dtype)
    
    # Create output file
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(
        output_file,
        input_ds.RasterXSize,
        input_ds.RasterYSize,
        1,
        input_ds.GetRasterBand(1).DataType,
        options=[
            'NUM_THREADS=ALL_CPUS',
            'COMPRESS=LZW',
            'TILED=YES',
            'BIGTIFF=YES'
        ]
    )
    
    if out_ds is None:
        raise RuntimeError(f"Failed to create: {output_file}")
    
    # Copy geotransform and projection
    out_ds.SetGeoTransform(input_ds.GetGeoTransform())
    out_ds.SetProjection(input_ds.GetProjection())
    
    # Write output data
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(output_data)
    out_band.SetNoDataValue(nodata_value)
    out_band.FlushCache()
    
    # Close datasets
    mask_ds = None
    input_ds = None
    out_ds = None
    
    print(f"Saved: {output_file}")


def process_directory(mask_file, input_dir, output_dir, pattern='*.tif', nodata_value=0):
    """
    Apply mask to all files in a directory matching the pattern.
    
    Args:
        mask_file: Mask raster file
        input_dir: Input directory containing files to mask
        output_dir: Output directory for masked files
        pattern: Filename pattern to match
        nodata_value: NoData value for output
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find input files
    search_pattern = os.path.join(input_dir, pattern)
    input_files = [f for f in glob.glob(search_pattern) if not ('_zscore_' in f or '_p_' in f)]
    
    if not input_files:
        print(f"Warning: No files matching '{pattern}' found in {input_dir}")
        return
    
    print(f"\nProcessing {len(input_files)} files from {input_dir}")
    
    # Process each file
    for input_file in sorted(input_files):
        filename = os.path.basename(input_file)
        output_file = os.path.join(output_dir, filename)
        
        try:
            apply_mask(mask_file, input_file, output_file, nodata_value)
        except Exception as e:
            print(f"Error processing {filename}: {e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Mask multiple data files and export as GeoTIFFs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Mask all tiff's in directory
  %(prog)s \\
    --input-mask /path/to/precomputed_mask.tif \\
    --input-dir /path/to/raw/ProdState \\
    --output-dir /path/to/ProdState \\
        """
    )
    
    parser.add_argument('--input-mask', required=True,
                        help='Precomputed mask file')
    parser.add_argument('--input-dir',
                        help='Input directory containing files to mask')
    parser.add_argument('--output-dir',
                        help='Output directory for masked files')
    parser.add_argument('--pattern', default='*.tif',
                        help='Filename pattern to match (default: *.tif)')
    parser.add_argument('--output-nodata', type=float, default=0,
                        help='NoData value for output files (default: 0)')
    
    args = parser.parse_args()
    
    try:
        if os.path.exists(args.input_dir):
            process_directory(
                args.input_mask,
                args.input_dir,
                args.output_dir,
                args.pattern,
                args.output_nodata
            )
        else:
            print(f"Warning: State directory not found: {args.input_dir}", file=sys.stderr)
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
