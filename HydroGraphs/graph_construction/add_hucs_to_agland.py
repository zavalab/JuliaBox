import geopandas as gpd
from WI_graph_functions import *

# Read in the agricultural land file
agland = gpd.GeoDataFrame.from_file("agland/WisconsinFieldBoundaries2019.shp")

# Read in HUC8, HUC10, and HUC12 files
HUC8   = gpd.GeoDataFrame.from_file("WIgeodataframes/HUC8/HUC8.shp")
HUC10  = gpd.GeoDataFrame.from_file("WIgeodataframes/HUC10/HUC10.shp")
HUC12  = gpd.GeoDataFrame.from_file("WIgeodataframes/HUC12/HUC12.shp")

# Set HUC8, HUC10, and HUC12 codes for agricultural land geodataframe
agland = add_HUC8(agland, HUC8)
agland = add_HUC10(agland, HUC10)
agland = add_HUC12(agland, HUC12)

# Save the agricultural land file
agland.to_pickle("agland.df")