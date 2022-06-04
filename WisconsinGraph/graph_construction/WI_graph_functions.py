import pandas as pd
import numpy as np
from shapely.ops import nearest_points
from shapely.geometry import MultiPoint
import matplotlib.dates as mdates

pd.options.display.max_columns = 250
pd.options.display.max_rows=150

import matplotlib.pyplot as plt
import matplotlib.colors as colors


import geopandas as gpd
from geopandas.tools import sjoin
import pandas as pd
from shapely.geometry import Point
from scipy.interpolate import interp1d
from tqdm import tqdm
import networkx as nx


def add_huc_col_to_hucs(HUC8, HUC10, HUC12):
    """
    Add columns to the HUC10 and HUC12 dataframes to give what HUC8 and HUC10 watersheds
    they lie in. The HUC8, HUC10, and HUC12 shapefiles have their identifying codes given by 
    strings in the dataframe. This allows us to identify the HUC8 and HUC10 (for HUC12 watersheds) 
    watersheds and the HUC8 (for HUC10 watersheds) watersheds. These identifiers are useful for 
    working with other shape files that lie in their same watersheds
    """
    
    # Add watershed columns
    HUC10["huc8"]  = ''
    HUC12["huc8"]  = ''
    HUC12["huc10"] = ''

    # Fill these columns based on their HUC10 or HUC12 codes
    for i in range(len(HUC10)):
        HUC10.huc8.iloc[i]     = HUC10.huc10.iloc[i][0:8]
    for i in range(len(HUC12)):
        HUC12.huc8.iloc[i]     = HUC12.huc12.iloc[i][0:8]
        HUC12.huc10.iloc[i]    = HUC12.huc12.iloc[i][0:10]
        
    
    # Convert the strings of the above dataframes to integers
    HUC8['huc8']   = HUC8['huc8'].astype(int) 
    HUC10['huc8']  = HUC10['huc8'].astype(int)
    HUC10['huc10'] = HUC10['huc10'].astype(int)
    HUC12['huc8']  = pd.to_numeric(HUC12['huc8'])
    HUC12['huc10'] = pd.to_numeric(HUC12['huc10'])
    HUC12['huc12'] = pd.to_numeric(HUC12['huc12'])
    
    # Return the HUC8, HUC10, and HUC12 watersheds with their idnetifiers
    return HUC8, HUC10, HUC12

def add_HUC8(shp, huc8):
    """
    Add HUC8 location to shp file. If centroid of object is in shp file is in a HUC8,
    then that HUC8 code is added to the object
    """
    
    print("Running add_HUC8 over GeoPandas dataframe")
    
    shp["huc8"] = 0
    polys  = huc8.geometry
    points = shp.centroid
    
    for i in tqdm(range(len(shp))):
        p = points.iloc[i]
        for j in range(len(huc8)):
            if p.within(polys.iloc[j]):
                shp['huc8'].iloc[i] = huc8.huc8.iloc[j]

    return shp

def add_HUC10(shp, huc10):
    """
    Add HUC12 location to shp file. If centroid of object is in shp file is in a HUC10,
    then that HUC10 code is added to the object

    Requires that the shp file has a column called "huc8"
    """
    
    print("Running add_HUC10 over GeoPandas dataframe")
    if 'huc8' not in shp.columns:
        raise ValueError("GeoPandas dataframe does not have a 'huc8' columnn")
    
    shp['huc10'] = 0
    points = shp.centroid

    for i in tqdm(range(len(shp))):
        p     = points.iloc[i]
        
        hucs  = huc10[huc10.huc8 == shp.huc8.iloc[i]]
        polys = hucs.geometry
        for j in range(len(polys)):
            if p.within(polys.iloc[j]):
                shp['huc10'].iloc[i] = hucs.huc10.iloc[j]
                break
                
    return shp

def add_HUC12(shp, huc12):

    """
    Add HUC12 location to shp file. If centroid of object is in shp file is in a HUC12,
    then that HUC12 code is added to the object

    Requires that the shp file has a column called "huc10"
    """
    
    print("Running add_HUC12 over GeoPandas dataframe")
    if 'huc10' not in shp.columns:
        raise ValueError("GeoPandas dataframe does not have a 'huc10' column")
    
    shp['huc12'] = 0
    points = shp.centroid
    
    for i in tqdm(range(len(shp))):
        p     = points.iloc[i] 
        hucs = huc12[huc12.huc10 == shp.huc10.iloc[i]]
        polys = hucs.geometry
        for j in range(len(polys)):
            if p.within(polys.iloc[j]):
                shp['huc12'].iloc[i] = hucs.huc12.iloc[j]
                break
                
    return shp

def add_lake_river_node_column(river_gdf,lake_gdf):
    """
    Adds a column to each gdf to identify which rivers or lakes the objects
    in that gdf intersect. Returns a river and lake gdf containing a new column
    listing the intersecting objects
    """
    
    lines = river_gdf.geometry
    polys = lake_gdf.geometry
    
    river_gdf["lake_nodes"] = ""
    lake_gdf["river_nodes"] = ""
    river_gdf = river_gdf.reset_index(drop=True)
    
    # Add empty list to river and lake gdfs
    for i in tqdm(range(len(polys))):
        lake_gdf.river_nodes.iloc[i] = []
    for i in tqdm(range(len(lines))):
        river_gdf.lake_nodes.iloc[i] = []

    # If a line intersects a lake polygon, add that lake polygon to the river's list of lakes
    # and add that river segment to the lake's list of rivers
    for i in tqdm(range(len(lines))):
        line = lines.iloc[i]
        
        for j in range(len(polys)):
            if line.intersects(polys.iloc[j]):
                river_gdf.lake_nodes.iloc[i].append(lake_gdf.COMID.iloc[j])
                lake_gdf.river_nodes.iloc[j].append(river_gdf.COMID.iloc[i])
                

    return river_gdf, lake_gdf

def add_to_tofroms(tofrom, river_gdf,lake_gdf):
    """
    From a given list of to-from COMIDs for the river network, add the lakes/waterbodies to the network.
    This function returns a new to-from list that includes COMIDs for waterbodies
    """
    
    # lake_list_applied is an optional list that keeps track of how many lakes
    # are impacted by the second for loop below
    lake_list_applied = []

    # Copy the to froms list so that the original dataframe is not impacted
    tofroms = tofrom.copy(deep=True)

    # Iterate through all rivers; if a river segment only intersects a single waterbody, 
    # then replace that river in the to-from list with the waterbody COMID
    for i in tqdm(range(len(river_gdf))):

        # Test to make sure the river segment only intersects one waterbody
        if len(river_gdf.lake_nodes.iloc[i]) == 1:
            # Identify the COMID of the river segment
            old_COMID = river_gdf.COMID.values[i]
            # Identify the COMID of the waterbody that is going to replace the river segment
            new_COMID = river_gdf.lake_nodes.iloc[i][0]
            
            new_lake_gdf =  lake_gdf[lake_gdf.COMID == new_COMID].copy(deep=True)
            
            # Replace the river COMID with the new lake COMID
            if len(new_lake_gdf) > 0:
                tofroms["TOCOMID"]   = tofroms["TOCOMID"].replace(to_replace=old_COMID, value=new_COMID)
                tofroms["FROMCOMID"] = tofroms["FROMCOMID"].replace(to_replace=old_COMID, value=new_COMID)

    # Iterate through all rivers; if a river segment intersects multiple waterbodies, it will add waterbodies 
    # to the network only if that waterbody wasn't already added to the network by the first for loop
    for i in tqdm(range(len(river_gdf))):
        

        # Test to make sure the river segment intersects multiple waterbodies
        if len(river_gdf.lake_nodes.iloc[i]) > 1:
            
            # Get the set of all waterbodies that are intersected by the river segment
            all_lakes = river_gdf.lake_nodes.iloc[i]
            
            # Get the gdf of the waterbodies intersected by the river segment
            new_lake_gdf = lake_gdf[lake_gdf.COMID.isin(all_lakes)].copy(deep=True)
            
            # Identify the COMID of the river segment
            old_COMID = river_gdf.COMID.values[i]
            
            # Create a list for the COMIDs that need to be added
            # This list will omit any waterbodies already added through the first for loop
            lake_comids_to_add = []
            
            # Loop through the gdf of intersected waterbodies to build the lake_comids_to_add list
            for k in range(len(new_lake_gdf)):
                # Get the COMID of the waterbody to test
                new_lake_comid = new_lake_gdf.COMID.iloc[k]

                # Test if the waterbody has already been added to the to-from list; if it hasn't,
                # then add it to the lake_comids_to_add list
                if len(tofroms[(tofroms.TOCOMID == new_lake_comid) | (tofroms.FROMCOMID == new_lake_comid)]) == 0:
                    lake_comids_to_add = np.append(lake_comids_to_add, new_lake_comid)
                    lake_list_applied  = np.append(lake_list_applied, new_lake_comid)
                    
            # Identify the length of the lake_comids to add
            length = len(lake_comids_to_add)
            
            # Loop through the lake_comids_to_add list; since the river segment intersects multiple lakes, the river 
            # segment is going to be replaced by every waterbody in this list; these waterbodies will not be directly
            # connected to each other, but upstream and downstream connections outside of thew waterbodies will be consistent
            for j in range(length-1):
                
                # Create a copy of the to-froms list, and replace the river segment COMID with the waterbody COMID
                tofroms_copy = tofroms[(tofroms.FROMCOMID == old_COMID) | (tofroms.TOCOMID == old_COMID)].copy(deep=True)
                new_COMID = lake_comids_to_add[j]
                
                tofroms_copy["TOCOMID"]   = tofroms_copy["TOCOMID"].replace(to_replace=old_COMID, value=new_COMID)
                tofroms_copy["FROMCOMID"] = tofroms_copy["FROMCOMID"].replace(to_replace=old_COMID, value=new_COMID)
                
                tofroms = tofroms.append(tofroms_copy)
            
            if length>0:

                # Replace the river segment COMID with the final waterbody COMID
                new_COMID = lake_comids_to_add[length-1]
            
                tofroms["TOCOMID"]   = tofroms["TOCOMID"].replace(to_replace=old_COMID, value=new_COMID)
                tofroms["FROMCOMID"] = tofroms["FROMCOMID"].replace(to_replace=old_COMID, value=new_COMID)

    # Reset the index; indices are now off because we have added new lines to the to-from list
    tofroms = tofroms.reset_index(drop=True)
    
    # Uncomment the two lines below to get the number of waterbodies that are impacted by the second for loop
    #print("The number of affected lakes is")
    #print(len(np.unique(lake_list_applied)))

    # Return the new to-from COMID list
    return tofroms

def remove_river_nodes(tofrom, river_gdf, lake_gdf):
    """
    Aggregate river nodes together if they fulfill certain criteria. 
    """
    
    # Make a copy of the to-froms list
    tofroms = tofrom.copy(deep=True)
    
    # Make a list of all waterbody COMIDs
    lake_list = [i for i in lake_gdf.COMID.values]
    # Iterate through the to-from list
    for i in tqdm(range(len(tofroms))):
        # Identify the to and from nodes based on TOCOMID and FROMCOMID
        to_node    = tofroms.TOCOMID.values[i]
        from_node  = tofroms.FROMCOMID.values[i]

        # Skip iteration if the to node and from node are identical
        if to_node == from_node:
            continue
    
        # Test if both the to and from nodes are both lakes; if so, skip the iteration
        if np.isin(to_node,lake_list) and np.isin(from_node,lake_list):
            continue
            
        # If both the to and from nodes are not in the lake list, test if the from node is in the lake list
        # If this statement is true, then the from_node is a waterbody and the to_node is a river
        elif np.isin(from_node, lake_list):

            # Get the list of all immediate upstream nodes to the to_node river segment
            list_of_froms = list(np.unique(tofroms.FROMCOMID.values[tofroms.TOCOMID == to_node]))

            # If the to_node is in the list_of_froms, remove it
            if to_node in list_of_froms:
                list_of_froms.remove(to_node)

            # Get the list of all immediate downstream nodes to the to_node river segment
            list_of_tos = list(np.unique(tofroms.TOCOMID.values[tofroms.FROMCOMID == to_node]))

            # If the to_node is in the list_of_tos, remove it
            if to_node in list_of_tos:
                list_of_tos.remove(to_node)
            
            # If there is only one value in the list_of_froms, then that upstream node is the from_node waterbody
            # We then will replace the to_node in the to-froms list with the from_node (i.e., we aggregate the 
            # river node with the waterbody from_node). This only occurs if the two nodes are in the same watershed
            if len(list_of_froms) == 1 and list_of_froms[0] == from_node:

                # Identify the HUC12 code for both the to and from node
                to_HUC   = river_gdf.huc12.values[river_gdf.COMID == to_node]
                from_HUC = lake_gdf.huc12.values[lake_gdf.COMID == from_node]
                
                # If the to and from nodes are in the same HUC12, replace the to_node river segment with the waterbody node
                if to_HUC == from_HUC:
                    tofroms["FROMCOMID"] = tofroms["FROMCOMID"].replace(to_replace=to_node, value=from_node)
                    tofroms["TOCOMID"]   = tofroms["TOCOMID"].replace(to_replace=to_node, value=from_node)
                    

        # If both the to and from nodes are not in the lake list, test if the to node is in the lake list
        # If this statement is true, then the to_node is a waterbody and the from_node is a river 
        elif np.isin(to_node, lake_list):

            # Get the list of all immediate upstream nodes to the from_node river segment
            list_of_froms = list(np.unique(tofroms.FROMCOMID.values[tofroms.TOCOMID == from_node]))

            # If the from_node is in the list_of_froms, remove it
            if from_node in list_of_froms:
                list_of_froms.remove(from_node)
            
            # Get the list of all immediate downstream nodes to the from_node river segment
            list_of_tos = list(np.unique(tofroms.TOCOMID.values[tofroms.FROMCOMID == from_node]))

            # If the from_node is in the list_of_tos, remove it
            if from_node in list_of_tos:
                list_of_tos.remove(from_node)

            # If there is only one value in the list_of_tos, then that downstream node is the to_node waterbody
            # We then will replace the from_node in the to-froms list with the to_node (i.e., we aggregate the 
            # river node with the waterbody to_node). This only occurs if the two nodes are in the same watershed
            if len(list_of_tos) == 1 and list_of_tos[0] == to_node:
                
                # Identify the HUC12 code for both the to and from node
                to_HUC   = lake_gdf.huc12.values[lake_gdf.COMID == to_node][0]
                from_HUC = river_gdf.huc12.values[river_gdf.COMID == from_node][0]
                
                # If the to and from nodes are in the same HUC12, replace the from_node river segment with the waterbody node
                if to_HUC == from_HUC:
                    tofroms["FROMCOMID"] = tofroms["FROMCOMID"].replace(to_replace=from_node, value=to_node)
                    tofroms["TOCOMID"]   = tofroms["TOCOMID"].replace(to_replace=from_node, value=to_node)
    
        # If both the to and from nodes are not lakes, then they must both be river nodes
        # We will then aggregate the river nodes if certain criteria are met
        else:
            
            # Get the list of nodes connected to the to_node river segment
            list_of_froms = list(np.unique(tofroms.FROMCOMID.values[tofroms.TOCOMID == to_node]))
            
            # If the list_of_froms has more than one value in it, then we will run further tests
            if len(list_of_froms) > 1:
                
                # Define an indicator. If certain criteria are met, we set this indicator to 1 and exit the loop
                indicator = 0

                # Iterate through the list of froms
                for k in list_of_froms:

                    # If the from nodes flow into more than one to_node, then we will not do any aggregation
                    # by setting the indicator to 1. This prevents us from aggregating a downstream node farther upstream
                    # and thus creating connections downstream that don't exist if we run this funciton again
                    if len(list(np.unique(tofroms.TOCOMID.values[tofroms.FROMCOMID == k]))) > 1:
                        indicator = 1
                # If the to_node is connected has an immediate upstream connection to a lake, set indicator to 1.
                # This prevents us from moving a lake downstream and creating connections that don't exist if we call the function again
                if  len(lake_gdf[lake_gdf.COMID.isin(list_of_froms)]) > 0:
                    indicator = 1
                
                # If we have set the indicator to 1, skip this iteration. Otherwise, continue to the else statement
                if indicator == 1:
                    continue

                # Test to make sure that all nodes to be aggregated are in the same HUC12. If not, we will not aggregate
                else: 

                    # Identify the HUC12 of the to_node
                    to_huc = river_gdf.huc12.values[river_gdf.COMID == to_node][0]

                    # Make a list to put all HUC12 values into
                    all_hucs = [to_huc]

                    # Iterate through all the list_of_froms and add their HUC12 to the list
                    for k in list_of_froms:
                        from_huc = river_gdf.huc12.values[river_gdf.COMID == k][0]
                        all_hucs = np.append(all_hucs, from_huc)
                    
                    # If all the potentially aggregated nodes are in the same HUC12, aggregate
                    # Note that we are aggregating all immediate upstream nodes with the to_node Because of our tests above,
                    # none of these 
                    if len(np.unique(all_hucs)) == 1:
                        for k in list_of_froms:                    
                            tofroms["FROMCOMID"] = tofroms["FROMCOMID"].replace(to_replace=k, value=to_node)
                            tofroms["TOCOMID"]   = tofroms["TOCOMID"].replace(to_replace=k, value=to_node)
                
            # If the list_of_froms is only length 1, then we will aggregate the from node into the two node if they are in the same HUC12
            else:

                # Identify HUC12 codes for each river node
                to_HUC   = river_gdf.huc12.values[river_gdf.COMID == to_node][0]
                from_HUC = river_gdf.huc12.values[river_gdf.COMID == from_node][0]
                
                # If the river nodes are in the same HUC12, then aggregate
                if to_HUC == from_HUC:
                    tofroms["FROMCOMID"] = tofroms["FROMCOMID"].replace(to_replace=from_node, value=to_node)
                    tofroms["TOCOMID"]   = tofroms["TOCOMID"].replace(to_replace=from_node, value=to_node)
                
                
    tofroms = tofroms.reset_index(drop=True)
    return tofroms.copy(deep=True)

def remove_tofrom_duplicates(tofroms):
    """
    Removes lines of the to-from list that contain the same TOCOMID and FROMCOMID value
    """

    # Make a copy of the tofroms list
    tofroms_copy = tofroms.copy(deep=True)
    
    # Iterate through the list
    # If the TOCOMID is the same as the FROMCOMID, drop that line
    for i in tqdm(range(len(tofroms))):
        if tofroms.TOCOMID.iloc[i] == tofroms.FROMCOMID.iloc[i]:
            tofroms_copy = tofroms_copy.drop(tofroms.index[i])
            
    # Return shortened to-from list
    return tofroms_copy.copy(deep=True)

def aggregate_river_nodes(tofroms, river_gdf, lake_gdf):
    """
    Aggregate river nodes from the to-from list. This function will remove several river nodes
    from the graph. However, it does not remove any waterbody nodes, and it does not merge any 
    nodes that lie in different watersheds
    """

    # Run the function remove_river_nodes and then remove the duplicates from the resulting to-from list
    tofroms = remove_river_nodes(tofroms.copy(), river_gdf.copy(), lake_gdf.copy())
    tofroms = remove_tofrom_duplicates(tofroms.copy())
    
    
    # Build a graph from the aggregated to from list
    G = build_graph(tofroms)
    
    # Define variables to determine if the aggregation change the graph at all. We will aggregate
    # nodes until there is no longer any change happening between calls of the aggregate function
    x = 0
    y = len([i for i in G.nodes])
    
    # While there is a change between calls to the aggregate (remove_river_nodes) function, then 
    # keep running the loop
    while y - x != 0:
        x = y

        # Run the function remove_river_nodes and then remove the duplicates from teh resulting to-from list
        # This function aggregates river nodes together
        tofroms = remove_river_nodes(tofroms.copy(), river_gdf.copy(), lake_gdf.copy())
        tofroms = remove_tofrom_duplicates(tofroms.copy())

        # Build the new graph to determine the number of nodes it now has
        G = build_graph(tofroms)
        
        # If the number of nodes before and after aggregation changes, rerun this loop
        y = len([i for i in G.nodes])
        print("new number of nodes = ", y)
        print("old number of nodes = ", x)
        if y - x != 0:
            print("Running another iteration")
            print("number of nodes changed")
            print()
        
    # Return aggregate to-from list
    return tofroms

def build_graph(tofroms):
    """
    Builds a directed graph network from the to-froms list
    """

    # Define the directed graph using NetworkX
    G = nx.DiGraph()

    # Iterate through to-froms
    for j in range(len(tofroms)):
        
        # Define the from node and to node based on to-from list
        from_node = tofroms.FROMCOMID.iloc[j]
        to_node   = tofroms.TOCOMID.iloc[j]
        
        # If the TOCOMID and FROMCOMID are different, add a new edge to the graph
        # Adding the edge also adds nodes to the graph that were not previously defined
        if from_node != to_node:
            G.add_edge(from_node, to_node)

    # Return the directed graph
    return G

def get_pos_dict(G, lake_gdf, riv_gdf, n_size = 10,lake_color="blue", river_color="red"):
    """
    Return a dictionary containing geographic locations of nodes for plotting.
    Also return a list of node colors and node_sizes. These will be used to plot
    the graph network and give different colors to different node types

    Takes inputs of a networkX directed graph (G), a waterbody GeoDataFrame (lake_gdf),
    and a river GeoDataFrame (riv_gdf)
    """
    
    # Define a dictionary to fill with geographic positions for the nodes
    G_pos = dict()
    
    # Define empty lists for the node colors and node sizes
    node_colors = []
    node_size   = []
    
    # Iterate through the nodes of G
    for j in tqdm((G.nodes)):
        
        # If the given node from G is in the river gdf, then take the geographic
        # location from riv_gdf; also set the color and size accordingly
        if np.isin(j, riv_gdf.COMID.values):
            poly_val = riv_gdf.geometry[riv_gdf.COMID == j]
            cent_val = poly_val.centroid
            node_colors.append(river_color)
            node_size.append(n_size)
        
        # If the given node from G is not a river gdf, then it must be a waterbody gdf.
        # Take the geographic location from lake_gdf, and set the color and size
        else:
            poly_val = lake_gdf.geometry[lake_gdf.COMID == j]
            cent_val = poly_val.centroid
            node_colors.append(lake_color)
            node_size.append(n_size)
            
        # Set the position as the x and y coordinate of the centroid of the river or waterbody object
        G_pos[j] = (np.array([cent_val.x.values[0], cent_val.y.values[0]]))

    # Return the dictionary of positions and the lists of node colors and sizes
    return G_pos, node_colors, node_size

def count_cycles(G):

    """
    Count the number of cycles in the graph, G. Print the result
    """

    counter = 0
    for i in nx.simple_cycles(G):
        counter += 1
    print(counter)

def get_downstream_graph(G, node, lake_gdf, river_gdf):

    """
    Construct the downstream graph of a given node
    Return the downstream graph and the set of node colors based on whether the 
    nodes are waterbody or river nodes
    """

    # Build a list of all nodes that lie downstream of a given node
    downstream = [n for n in nx.traversal.bfs_tree(G, node, reverse=False)]
    
    # Define a directed graph
    subG = nx.DiGraph()

    subG.add_node(node)
    
    for j in downstream:
        edge_list = G.edges(j)
        
        for i in edge_list:
            if i[0] in downstream and i[1] in downstream:
                subG.add_edge(i[0], i[1])
    
    node_colors = []
    
    for i in subG.nodes:
        if i in river_gdf.COMID.values:
            node_colors.append("red")
        elif i in lake_gdf.COMID.values:
            node_colors.append("blue")
        else:
            node_colors.append("orange")
            
    return subG, node_colors

def get_upstream_graph(G, node, lake_gdf, river_gdf):
    downstream = [n for n in nx.traversal.bfs_tree(G, node, reverse=True)]
    
    subG = nx.DiGraph()

    subG.add_node(node)
    
    for j in downstream:
        edge_list = G.edges(j)
        
        for i in edge_list:
            if i[0] in downstream and i[1] in downstream:
                subG.add_edge(i[0], i[1])
    
    node_colors = []

    for i in subG.nodes:
        if i in river_gdf.COMID.values:
            node_colors.append("red")
        elif i in lake_gdf.COMID.values:
            node_colors.append("blue")
        else:
            node_colors.append("orange")
            
    return subG, node_colors


def add_CAFOS_to_graph(G_old, lake_gdf, river_gdf, CAFOS):
    G = G_old.copy()
    
    all_nodes = [i for i in G.nodes]
    node_df = pd.concat([lake_gdf[["COMID", 'huc12', 'geometry']][lake_gdf.COMID.isin(all_nodes)].copy(), river_gdf[['COMID', 'huc12','geometry']][river_gdf.COMID.isin(all_nodes)].copy()])
    
    for i in range(len(CAFOS)):
        huc = CAFOS.huc12.iloc[i]
        point = CAFOS.geometry.iloc[i]
        cafo_name = CAFOS.Node.iloc[i]
        
        df_in_huc = node_df[node_df.huc12 == huc].copy()
        df_in_huc = df_in_huc.reset_index(drop=True)
        
        distances = df_in_huc.distance(point)
        
        if len(df_in_huc) != 0:
            comid_to_connect = df_in_huc.COMID.values[distances == min(distances)][0]
        
            G.add_edge(cafo_name, comid_to_connect)
    return G
        
def build_graph_with_pollutants(tofroms,lake_gdf, river_gdf, source):
    node_df = pd.concat([lake_gdf[["COMID", 'huc12']], river_gdf[['COMID', 'huc12']]])
    
    G = nx.DiGraph()
    
    for j in range(len(tofroms)):
        from_node = tofroms.FROMCOMID.iloc[j]
        to_node   = tofroms.TOCOMID.iloc[j]
        
        if from_node != to_node:
            G.add_edge(from_node, to_node)
    
    all_nodes = [i for i in G.nodes]
    
    for i in range(len(source)):
        huc = source.huc12.iloc[i]
        new_df = node_df[(node_df.huc12 == huc) & (node_df.COMID.isin(all_nodes))]
        
        for j in range(len(new_df)):
            from_node = source.Node.iloc[i]
            to_node = new_df.COMID.iloc[j]
            G.add_edge(from_node, to_node)
    return G

def get_pos_dict_with_pollutant(G, lake_gdf, riv_gdf, source):
    G_pos = dict()
    
    node_colors= []
    node_size = []
    
    
    for j in tqdm((G.nodes)):
        
        if np.isin(j, riv_gdf.COMID.values):
            poly_val = riv_gdf.geometry[riv_gdf.COMID == j]
            cent_val = poly_val.centroid
            node_colors.append("red")
            node_size.append(10)
        elif np.isin(j, source.Node.values):
            cent_val = source.geometry[source.Node == j]
            node_colors.append("orange")
            node_size.append(30)
        else:
            poly_val = lake_gdf.geometry[lake_gdf.COMID == j]
            cent_val = poly_val.centroid
            node_colors.append("blue")
            node_size.append(10)
            
        G_pos[j] = (np.array([cent_val.x.values[0], cent_val.y.values[0]]))

    return G_pos, node_colors, node_size


def add_source_to_graph(G, lake_gdf, river_gdf, source):
    all_nodes = [i for i in G.nodes]
    node_df = pd.concat([lake_gdf[["COMID", 'huc12', 'geometry']][lake_gdf.COMID.isin(all_nodes)].copy(), river_gdf[['COMID', 'huc12','geometry']][river_gdf.COMID.isin(all_nodes)].copy()])
    
    for i in range(len(source)):
        huc = source.huc12.iloc[i]
        point = source.geometry.iloc[i]
        cafo_name = source.Node.iloc[i]
        
        df_in_huc = node_df[node_df.huc12 == huc].copy()
        df_in_huc = df_in_huc.reset_index(drop=True)
        
        distances = df_in_huc.distance(point)
        
        if len(df_in_huc) != 0:
            comid_to_connect = df_in_huc.COMID.values[distances == min(distances)][0]
        
            G.add_edge(cafo_name, comid_to_connect)
    return G






