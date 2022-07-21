# HydroGraphs

This repository contains the code for the manuscript "A Graph Formulation for Tracing Hydrological Pollutant Transport in Surface Waters." There are three main folders containing code and data, and these are outlined below. We call the framework for building a graph of these hydrological systems `HydroGraphs`. This code was run using Python versions 3.8 and 3.9.

### NOTE: because of the size of some of the data, the files are not currently uploaded to Github. This will be updated in the future, either with the files uploaded using Git LFS or with a link to another location where the larger files can be downloaded. 

## 1. graph_construction
This folder contains the data and code for building a graph of the watershed-river-waterbody hydrological system. It uses data from the Watershed Boundary Dataset (WBD; link [here](https://apps.nationalmap.gov/downloader/#/)) and the National Hydrography Dataset (NHDPlusV2; link [here](https://www.epa.gov/waterdata/get-nhdplus-national-hydrography-dataset-plus-data)) as a basis and builds a list of directed edges. We use `NetworkX` to build and visualize the list as a graph. This folder contains the following subfolders:

* CAFOS_shp - contains the shape files for concentrated animal feeding operations (CAFOs) in Wisconsin
* Counties - contains the shape files for all the counties in Wisconsin, obtained from [here](https://data-wi-dnr.opendata.arcgis.com/datasets/wi-dnr::county-boundaries-24k/about). This shapefile was dissolved to get a polygon of the state of Wisconsin. 
* agland - contains a shapefile of the agricultural land for the state of Wisconsin obtained from [here](https://doi.org/10.15482/USDA.ADC/1520625). This file was altered to include the HUC8, HUC10, and HUC12 codes for each set of land, and the resulting GeoDataFrame is agland.df of this same folder
* lakes_rivers - contains the basic dataframes obtained from the NHDPlusV2 and WBD. These are contained within the Watershed4 and Watershed7 subfolders. These subfolders also contain the PlusFlow.dbf file from the NHDPlusV2, and the PluFlow.csv file that is a CSV format of the PlusFlow.dbf. In addition, there is a subfolder, WI, that contains the shapefile of Wisconsin that was constructed from the Wisconsin counties. 
* WIgeodataframes - contains the shapefiles of the lakes, rivers, and HUCs for the state of Wisconsin. These were obtained from the NHDPlusV2 and WBD data in lakes_rivers, but with some minor attribute additions and only over the state of Wisconsin. In addition, this folder contains the CSVs of TOCOMID and FROMCOMID list for different graphs. These lists are manipulations of the PlusFlow.csv where:
    * WItofroms.csv is the PlusFlow.csv applied to just the state of Wisconsin. This can be used to form the river graph
    * WItofroms_lakes.csv is the WItofroms.csv with waterbodies added to the list, replacing some river COMIDs
    * WItofroms_agg.csv is the aggregated form of WItofroms_lakes.csv

This folder also contains several files and scripts. These are as follows:

* WILakes.df - This is a set of the waterbodies from the NHDPlusV2 for the state of Wisconsin. It includes columns for the HUC8, HUC10, and HUC12 codes and a column for intersecting rivers
* WIRivers.df - This is a set of the rivers from the NHDPlusV2 for the state of Wisconsin. It includes columns for the HUC8, HUC10, and HUC12 codes and a column for intersecting waterbodies. 
* agland.df - This is the agricultural land with the HUC8, HUC10, and HUC12 codes added
* WI_graph_functions.py - This script contains several functions for working with `HydroGraphs`. It also contains several functions that are called to form the desired graphs from the NHDPlusV2 and WBD.
* run_WI_graph_code.py - This script runs the following other scripts, in the following order, to take the data from lakes_rivers and construct the files in WIgeodataframes and the files WILakes.df and WIRivers.df
    * build_base_dataframes.py - takes the NHDPlusV2 data and removes the areas outside of Wisconsin.
    * add_hucs_to_lakes_rivers.py - adds the HUC8, HUC10, and HUC12 codes to the river and waterbody GeoDataFrames
    * add_river_lake_nodes.py - adds the column of intersecting waterbodies to the river GeoDataFrame and adds the column of intersecting rivers to the waterbody GeoDataFrame. 
    * add_to_comid_list.py - builds the first two tofrom csvs in the WIgeodataframes folder
    * aggregate_graph.py - aggregates the graph and forms the WItofroms_agg.csv
* Visualizations.ipynb - This script builds some of the visualizations formed in the manuscript mentioned above, and it contains a test of connectivity to compare the aggregated graph to the original graph. 

## 2. case_studies

This folder contains three .ipynb files for three separate case studies. These three case studies focus on how "Hydrology Graphs" can be used to analyze pollutant impacts in surface waters. Details of these case studies can be found inthe manuscript above.

## 3. DNR_data

This folder contains data from the Wisconsin Department of Natural Resources (DNR) on water quality in several Wisconsin lakes. The data was obtained from [here](https://dnr.wi.gov/lakes/waterquality/) using the file "Web_scraping_script.py". The original downloaded reports are found in the folder "original_lake"reports." These reports were then cleaned and reformatted using the script "DNR_data_filter.ipynb." The resulting, cleaned reports are found in the "Lakes" folder. Each subfolder of the "Lakes" folder contains data for a single lake. Teh two .csvs "lake_index_WBIC.csv" contain an index for what lake each numbered subfolder corresponds. In addition, we added the corresponding COMID in "lake_index_WBIC_COMID.csv" by matching the NHDPlusV2 data to the HYDROLakes [data](https://doi.org/10.1038/ncomms13603). The DNR's reported data only matches lakes to a waterbody identification code (WBIC), so we use HYDROLakes (indexed by WBIC) to match to the COMID. This is done in the "DNR_data_filter.ipynb" script as well. 
