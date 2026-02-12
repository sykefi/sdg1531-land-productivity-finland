#!/usr/bin/env python3
"""
Merge multiple image subsets into single COG files.

This script replaces the following GDAL operations:
- Building a VRT from multiple TIF files
- Converting the VRT to a Cloud Optimized GeoTIFF (COG)
"""

import argparse
import glob
import os
import sys

from osgeo import gdal

# Enable GDAL exceptions
gdal.UseExceptions()


def find_tif_files(directory, pattern):
    """
    Recursively find all TIF files matching the pattern.
    
    Args:
        directory: Root directory to search
        pattern: Filename pattern (e.g., 'degraded*.tif')
    
    Returns:
        List of full paths to matching TIF files
    """
    search_pattern = os.path.join(directory, '**', pattern)
    files = glob.glob(search_pattern, recursive=True)
    return sorted(files)


def build_vrt(input_files, output_vrt, src_nodata=0, vrt_nodata=0):
    """
    Build a VRT from multiple input files.
    
    Args:
        input_files: List of input file paths
        output_vrt: Output VRT file path
        src_nodata: Source nodata value
        vrt_nodata: VRT nodata value
    """
    print(f"Building VRT from {len(input_files)} files...")
    
    vrt_options = gdal.BuildVRTOptions(
        srcNodata=src_nodata,
        VRTNodata=vrt_nodata
    )
    
    vrt_ds = gdal.BuildVRT(output_vrt, input_files, options=vrt_options)
    
    if vrt_ds is None:
        raise RuntimeError(f"Failed to build VRT: {output_vrt}")
    
    vrt_ds = None  # Close the dataset
    print(f"VRT created: {output_vrt}")


def convert_to_cog(input_vrt, output_tif):
    """
    Convert VRT to Cloud Optimized GeoTIFF (COG).
    
    Args:
        input_vrt: Input VRT file path
        output_tif: Output COG file path
    """
    print(f"Converting to COG: {output_tif}...")
    
    translate_options = gdal.TranslateOptions(
        format='COG',
        outputType=gdal.GDT_Byte,
        creationOptions=['COMPRESS=DEFLATE', 'NUM_THREADS=ALL_CPUS']
    )
    
    ds = gdal.Translate(output_tif, input_vrt, options=translate_options)
    
    if ds is None:
        raise RuntimeError(f"Failed to create COG: {output_tif}")
    
    ds = None  # Close the dataset
    print(f"COG created: {output_tif}")


def process_subset(subset_dir, output_file, pattern='degraded*.tif', 
                   src_nodata=0, vrt_nodata=0, cleanup_vrt=True):
    """
    Process a subset directory: find files, build VRT, convert to COG.
    
    Args:
        subset_dir: Directory containing subset TIF files
        output_file: Output COG file path
        pattern: Filename pattern to match
        src_nodata: Source nodata value
        vrt_nodata: VRT nodata value
        cleanup_vrt: Whether to delete temporary VRT file
    """
    # Find input files
    input_files = find_tif_files(subset_dir, pattern)
    
    if not input_files:
        raise ValueError(f"No files matching '{pattern}' found in {subset_dir}")
    
    print(f"Found {len(input_files)} files matching '{pattern}'")
    
    # Create temporary VRT
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    
    vrt_file = os.path.join(output_dir, 
                            os.path.splitext(os.path.basename(output_file))[0] + '.vrt')
    
    try:
        # Build VRT
        build_vrt(input_files, vrt_file, src_nodata, vrt_nodata)
        
        # Convert to COG
        convert_to_cog(vrt_file, output_file)
        
    finally:
        # Cleanup VRT if requested
        if cleanup_vrt and os.path.exists(vrt_file):
            os.remove(vrt_file)
            print(f"Cleaned up temporary VRT: {vrt_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Merge multiple image subsets into single COG files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process baseline subset
  %(prog)s -i /path/to/subset_results/baseline -o /path/to/output/prodperf_degraded_baseline.tif
  
  # Process reporting subset
  %(prog)s -i /path/to/subset_results/reporting -o /path/to/output/prodperf_degraded_reporting.tif
  
        """
    )
    
    parser.add_argument('-i', '--input-dir', required=True,
                        help='Input directory containing subset TIF files')
    parser.add_argument('-o', '--output', required=True,
                        help='Output COG file path')
    parser.add_argument('-p', '--pattern', default='degraded*.tif',
                        help='Filename pattern to match (default: degraded*.tif)')
    parser.add_argument('--src-nodata', type=float, default=0,
                        help='Source nodata value (default: 0)')
    parser.add_argument('--vrt-nodata', type=float, default=0,
                        help='VRT nodata value (default: 0)')
    parser.add_argument('--keep-vrt', action='store_true',
                        help='Keep temporary VRT file')
    
    args = parser.parse_args()
    
    try:
        process_subset(
            subset_dir=args.input_dir,
            output_file=args.output,
            pattern=args.pattern,
            src_nodata=args.src_nodata,
            vrt_nodata=args.vrt_nodata,
            cleanup_vrt=not args.keep_vrt
        )
        
        print("Processing completed successfully!")
        return 0
        
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
