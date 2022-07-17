import warnings

import time
import os
import pandas as pd
import numpy as np
import json
from shapely.ops import nearest_points
from shapely.geometry import MultiPoint
import matplotlib.dates as mdates
import math
pd.options.display.max_columns = 50
pd.options.display.max_rows=50
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors as colors
import rtree
import geopandas as gpd
import pandas as pd
from scipy.interpolate import interp1d
import datetime
import scipy.stats as st
from scipy.optimize import curve_fit
from tqdm import tqdm
from WI_graph_functions import *
warnings.filterwarnings('ignore')

# Read in the Wisconsin dataframes for lakes and rivers
WILakes  = pd.read_pickle("WILakes.df")
WIRivers = pd.read_pickle("WIRivers.df")

# Read in the list of tofroms that includes waterbodies
WItofroms_lakes = pd.read_csv("WIgeodataframes/WItofroms_lakes.csv")

# Aggregate the intermediate river nodes
WItofroms_agg = aggregate_river_nodes(WItofroms_lakes, WIRivers, WILakes)

# Save the aggregated list of tofroms
WItofroms_agg.to_csv("WIgeodataframes/WItofroms_agg.csv")
