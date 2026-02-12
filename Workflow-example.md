# Example workflow for generating SDG 15.3.1 subindicators

This example AOI does not cover the whole biogeographical zone, which means that the Performance indicator will not be comparable to full-area processing. Trend and State are pixel-wise algorithms and therefore should be identical to full-area results.

The example dataset covers a small subset of the northernmost Lapland.

# First download the example data bundle 

Dataset can be obtained from [sdg1531_example_data.zip](../../releases/download/v0.1-alpha/sdg1531_example_data.zip) 

Extract to, for example, E:\tmp\sdg1531_example_data

# Create and activate the correct environment
    
    mamba create -n eosdg-env python=3.12 dask=2025.9.1 xarray rioxarray rasterio gdal crick tqdm numpy scipy pandas
    
or

    mamba env create -n eosdg -f eosdg-env.yaml
    mamba activate eosdg

You can monitor processing of state and performance computations from the Dask dashboard: http://localhost:8787/status

Trend tool updates progress bar in shell

# Input data, modify root paths as needed

    set dataroot=D:\temp\sdg1531_example_data
    set outputroot=D:\temp\sdg1531_example_output


    mkdir %outputroot%
    set ndvisrcpth=%dataroot%\NDVI
    set clcsrcpth=%dataroot%\CLC
    set biogeozones=%dataroot%\BiogeographicalZones\BiogeographicalZones.shp
    set watermasksrcpth=%dataroot%\auxdata\ranta10
    set countryborders=%dataroot%\auxdata\TietoaKuntajaosta_2025_10k\SuomenValtakunta_2025_10k.shp

# Compute Productivity State

    python compute_state.py %ndvisrcpth% --timespans 2000 2012 2013 2015 --tag baseline --odir %outputroot%\raw\ProdState\ --nodata 0 --blocksize 2048
    
    python compute_state.py %ndvisrcpth% --timespans 2004 2016 2017 2019  --tag reporting --odir %outputroot%\raw\ProdState\ --nodata  0 --blocksize 2048

# Compute Productivity Trend

Use smaller blocksize here, multiprocessing approach uses more memory

    python compute_trend.py  %ndvisrcpth% 2000 2015 --odir %outputroot%\raw\ProdTrend\ --tag baseline --nodata 0 --blocksize 1024

    python compute_trend.py  %ndvisrcpth% 2004 2019 --odir %outputroot%\raw\ProdTrend\  --tag reporting --nodata 0   --blocksize 1024

# Compute Productivity Performance

## First meanNDVI 

    python compute_mean_ndvi.py %ndvisrcpth% --timespan 2000 2015  --tag baseline --odir %outputroot%\raw\ProdPerf\MeanNDVI --nodata 0 --blocksize 2048 

    python compute_mean_ndvi.py %ndvisrcpth% --timespan 2004 2019  --tag reporting  --odir %outputroot%\raw\ProdPerf\MeanNDVI --nodata 0 --blocksize 2048 

## Then split meanNDVI and CLC to biogeographical zones, or only zone 5 in this demo case.

Due to the zone geometry, result extent will be expanded and filled with no-data

    set PROJ_NETWORK=OFF
    python prepare_data_for_prodperf.py %outputroot%\raw\ProdPerf\MeanNDVI\mean_ndvi_max_baseline.tif  %biogeozones% %outputroot%\raw\ProdPerf\subsets --geomids 5

    python prepare_data_for_prodperf.py %outputroot%\raw\ProdPerf\MeanNDVI\mean_ndvi_max_reporting.tif  %biogeozones% %outputroot%\raw\ProdPerf\subsets --geomids 5

    python prepare_data_for_prodperf.py %clcsrcpth%\clc2012peruutus_recoded_sdg1531.tif  %biogeozones% %outputroot%\raw\ProdPerf\subsets --geomids 5

    python prepare_data_for_prodperf.py %clcsrcpth%\clc2018peruutus_recoded_sdg1531.tif  %biogeozones% %outputroot%\raw\ProdPerf\subsets --geomids 5

## Get p90 and then classify

    python compute_performance.py %outputroot%\raw\ProdPerf\subsets %outputroot%\raw\ProdPerf\subsets_results --blocksize 2048 --geomids 5

## Finally merge subsets

    python merge_subsets.py -i  %outputroot%\raw\ProdPerf\subsets_results\baseline\ -o %outputroot%\raw\ProdPerf\prodperf_degraded_baseline.tif -p degraded*.tif
    
    python merge_subsets.py -i  %outputroot%\raw\ProdPerf\subsets_results\reporting\ -o %outputroot%\raw\ProdPerf\prodperf_degraded_reporting.tif -p degraded*.tif
    
# As a last step mask out all water pixels and pixels outside Finland for all subindicators

    set maskingdir=D:\temp\sdg1531_example_output\masking
    mkdir %maskingdir%
   
Get one template file to the current directory and rasterize masks to it.

    python create_base_mask.py --input-template %outputroot%\raw\ProdState\state_degraded_baseline.tif  --output-mask %maskingdir%\mask.tif  --border-shapefile %countryborders%  --water-shapefiles %watermasksrcpth%\ranta10jarvet\Jarvi10.shp %watermasksrcpth%\ranta10joet\JokiAlue10.shp

Example AOI does not include sea areas; in the general case, also include the sea mask in the `--water-shapefiles` argument list

     %watermasksrcpth%\ranta10meret\meri10.shp  

# Then use this template to mask all subindicators
        
    python mask_raw_products.py --input-mask %maskingdir%\mask.tif --input-dir %outputroot%\raw\ProdState  --output-dir %maskingdir%\ProdState

    python mask_raw_products.py --input-mask %maskingdir%\mask.tif --input-dir %outputroot%\raw\ProdTrend  --output-dir %maskingdir%\ProdTrend

# Then do the same for the performance which has slightly different AOI due to subsetting and biogeographical zone processing

    python create_base_mask.py --input-template  %outputroot%\raw\ProdPerf\prodperf_degraded_baseline.tif --output-mask %maskingdir%\mask_perf.tif --border-shapefile %countryborders%  --water-shapefiles %watermasksrcpth%\ranta10jarvet\Jarvi10.shp %watermasksrcpth%\ranta10joet\JokiAlue10.shp

    python mask_raw_products.py --input-mask %maskingdir%\mask_perf.tif --input-dir %outputroot%\raw\ProdPerf  --output-dir %maskingdir%\ProdPerf

# Final masked results are now stored in %maskingdir% 
