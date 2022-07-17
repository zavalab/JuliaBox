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


import build_base_dataframes
import add_hucs_to_lakes_rivers
import add_river_lake_nodes
import add_to_comid_list
import aggregate_graph