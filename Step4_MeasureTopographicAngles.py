########################################################################
# TTools
# Step 4: Measure Topographic Angles - v 0.97
# Ryan Michie

# Measure_Topographic_Angles will take an input point feature 
# (from Step 1) and calculate the maximum topographic elevation
# and the the slope angle from each node in different directions.

# INPUTS
# 0: Input TTools point feature class(nodes_fc)
# 1: input the directions to 
#    sample (topo_directions) 1. [W,S,E], 2. [NE,E,SE,S,SW,W,NW,N]
# 2: input the maximum km distance to search (searchDistance_max_km)
# 3: input elevation raster (z_raster)
# 4: input elevation raster z units (z_units) 
#     "Feet", "Meters", or "Other"
# 5: output sample point file name/path (topo_fc)
# 6: input flag if existing data can be over
#     written (overwrite_data) True or False

# OUTPUTS
# 0. point feature class (edit nodes_fc) - Added fields with topographic 
#     shade angles for each direction at each node
# 1. point feature class (new) - point for each x/y location 
#     of the maximum elevation.

# Future Updates
# find the bug that generates off array samples (see s. willamette)
# update fc after each topo line sample so data isn't lost if the 
# script crashes
# eliminate arcpy and use gdal for reading/writing feature class data

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

########################################################################

# Import system modules
from __future__ import division, print_function
import sys
import os
import gc
import time
import traceback
from datetime import timedelta
import arcpy
from arcpy import env
from math import radians, sin, cos, hypot, ceil
from collections import defaultdict
import numpy as np

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
topo_directions = 1
searchDistance_max_km = 10
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
z_units = "Meters"
topo_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\topo_samples"
block_size = "#" # OPTIONAL defualt to 5
overwrite_data = True
# End Fill in Data

# Used for debugging. Currently turned off. 
topo_line_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\topo_line"
block_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\blocks"
plot_dir = r"D:\Projects\TTools_9\plots"
# ----------------------------------------------------------------------

# Parameter fields for python toolbox
#nodes_fc = parameters[0].valueAsText
#topo_directions = parameters[1].valueAsText # Needs to be a long
#searchDistance_max_km = parameters[2].valueAsText
#z_raster = parameters[3].valueAsText
#z_units = parameters[4].valueAsText
#topo_fc = parameters[5].valueAsText
#overwrite_data = parameters[6].valueAsText True/False

def nested_dict(): 
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the STREAM_ID,
    NODE_ID, and X/Y coordinates as a nested dictionary"""
    
    print("Reading nodes feature class")
    
    nodeDict = nested_dict()
    incursorFields = ["NODE_ID", "STREAM_ID", "Z_NODE", "SHAPE@X","SHAPE@Y"]

    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)     

    # Check to see if the 1st field exists if yes add it.
    if overwrite_data is False and (addFields[0] in existingFields) is True:
        incursorFields.append(addFields[0])
    else:
        overwrite_data = True

    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc,incursorFields,"",proj) as Inrows:
        if overwrite_data:
            for row in Inrows:
                nodeDict[row[0]]["STREAM_ID"] = row[1]
                nodeDict[row[0]]["Z_NODE"] = row[2]
                nodeDict[row[0]]["POINT_X"] = row[3]
                nodeDict[row[0]]["POINT_Y"] = row[4]
                
        else:
            for row in Inrows:
                # if the data is null or zero (0 = default for shapefile),
                # it is retreived and will be overwritten.
                if row[5] is None or row[5] == 0 or row[5] < -9998:
                    nodeDict[row[0]]["STREAM_ID"] = row[1]
                    nodeDict[row[0]]["Z_NODE"] = row[2]
                    nodeDict[row[0]]["POINT_X"] = row[3]
                    nodeDict[row[0]]["POINT_Y"] = row[4]
    if len(nodeDict) == 0:
        sys.exit("The fields checked in the input point feature class "+
                 "have existing data. There is nothing to process. Exiting")
            
    return(nodeDict)

def create_topo_line_fc(topo_line, streamID, nodeID, a, topo_line_fc, proj):
    
    poly_array = arcpy.Array()
    pnt = arcpy.Point()
    cursorfields = ["STREAM_ID", "NODE_ID", "AZIMUTH"]   
    
    # Check to see if the block fc exists, if not create it
    if not arcpy.Exists(topo_line_fc):  
        #Create an empty output with the same projection as the input polyline
        arcpy.CreateFeatureclass_management(os.path.dirname(topo_line_fc),
                                            os.path.basename(topo_line_fc),
                                            "POLYLINE","","DISABLED","DISABLED",
                                            proj)
        
        # Determine Stream ID field properties
        sid_type = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].type
        sid_precision = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].precision
        sid_scale = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].scale
        sid_length = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].length    
    
        # Add attribute fields # TODO add dictionary 
        # of field types so they aren't all double
        for f in cursorfields:
            if f == "STREAM_ID":
                arcpy.AddField_management(topo_line_fc, f, sid_type,
                                          sid_precision, sid_scale,
                                          sid_length, "", "NULLABLE",
                                          "NON_REQUIRED")
            else:
                arcpy.AddField_management(topo_line_fc, f, "DOUBLE", "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(topo_line_fc, ["SHAPE@"] + cursorfields) as cursor:
        for pnt_x, pnt_y in topo_line:
            pnt.X = pnt_x
            pnt.Y = pnt_y
            poly_array.add(pnt)
        poly = arcpy.Polyline(poly_array)
        cursor.insertRow([poly, streamID, nodeID, a])
        poly_array.removeAll() 

def update_topo_fc(topo_list, topo_fc, nodes_fc, nodes_to_update, overwrite_data, proj):
    """Creates/updates the output topo point feature
    class using the data from the topo list"""
    
    #Create an empty output with the same projection as the input polyline
    cursorfields = ["POINT_X","POINT_Y","STREAM_ID","NODE_ID",
                    "AZIMUTH","TOPOANGLE","TOPO_ELE","NODE_ELE",
                    "ELE_CHANGE","TOPODIS","SEARCHDIS","NA_SAMPLES"]    
    
    # Check if the output exists and create if not
    if not arcpy.Exists(topo_fc):
        arcpy.CreateFeatureclass_management(os.path.dirname(topo_fc),
                                            os.path.basename(topo_fc),
                                            "POINT","","DISABLED","DISABLED",proj)
        
        # Determine Stream ID field properties
        sid_type = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].type
        sid_precision = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].precision
        sid_scale = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].scale
        sid_length = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].length    
    
        # Add attribute fields # TODO add dictionary of field types 
        # so they aren't all double
        for f in cursorfields:
            if f == "STREAM_ID":
                arcpy.AddField_management(topo_fc, f, sid_type,
                                          sid_precision, sid_scale, sid_length,
                                          "", "NULLABLE", "NON_REQUIRED")
                
            else:
                arcpy.AddField_management(topo_fc, f, "DOUBLE", "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")
        
    if not overwrite_data:
        # Build a query to retreive existing rows from the nodes 
        # that need updating
        whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_update))        
        
        # delete those rows
        with arcpy.da.UpdateCursor(topo_fc,["NODE_ID"], whereclause) as cursor:  
            for row in cursor:
                cursor.deleteRow()    

    with arcpy.da.InsertCursor(topo_fc, ["SHAPE@X","SHAPE@Y"] +
                               cursorfields) as cursor:
        for row in topo_list:
            cursor.insertRow(row)

def update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_query):
    """Updates the input point feature class with data from the nodes dictionary"""
    
    # Build a query to retreive just the nodes that needs updating
    whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))

    with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:  
        for row in cursor:
            for f, field in enumerate(addFields):
                nodeID =row[0]
                row[f+1] = nodeDict[nodeID][field]
                cursor.updateRow(row)

def to_meters_con(inFeature):
    """Returns the conversion factor to get
    from the input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that is".format(inFeature)+
                       "not projected or not recognized. Use a projected "+
                       "coordinate system preferably in linear units of "+
                       "feet or meters.")
        sys.exit("Coordinate system is not projected or not recognized. "+
                 "Use a projected coordinate system, preferably in"+
                 "linear units of feet or meters.")   
    return con_to_m
    
def from_meters_con(inFeature):
    """Returns the conversion factor to get from
    meters to the spatial units of the input feature class"""
    try:
        con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that ".format(inFeature)+
                       "is not projected or not recognized. Use a "+
                       "projected coordinate system preferably in "+
                       "linear units of feet or meters.")
        sys.exit("Coordinate system is not projected or not recognized. "+
                 "Use a projected coordinate system, preferably in "+
                 "linear units of feet or meters.")   
    return con_from_m

def from_z_units_to_meters_con(zUnits):
    """Returns the converstion factor to
    get from the input z units to meters"""
        
    try:
        con_z_to_m = float(zunits)
    except:
        if zUnits == "Meters":
            con_z_to_m = 1.0 
        elif zUnits == "Feet":
            con_z_to_m = 0.3048
        else: con_z_to_m = None # The conversion factor will not be used
    
    return con_z_to_m

def build_search_array(searchDistance_min, searchDistance_max, cellsize, use_skippy):
    """Build a numpy array from the minimum to the max
    search distance by increments of the cellsize."""
    
    # use next cell over to avoid divide by zero errors
    if searchDistance_min <= 0: searchDistance_min = cellsize
    
    if use_skippy:
        # This is a modified version of the skippy 
        # algorithm from Greg Pelletier. It is not being used but I'm 
        # keeping it in here in case someone wants to turn it on.
        # It needs to be fixed so the ncells are adjusted based on distance
        searchDistance = searchDistance_min
        distanceList = [searchDistance_min]
        
        #ncells = max([int(ceil((block_x_max - block_x_min)/ x_cellsize)), 1])
        
        ncells = 1
        while not searchDistance > searchDistance_max:
            if ncells <= 10:
                searchDistance_m = searchDistance + (cellsize)
            if 10 < ncells <= 20:
                searchDistance = searchDistance + (cellsize * 3)
            if 20 < ncells <= 40:
                searchDistance = searchDistance + (cellsize * 6)
            if 40 < ncells <= 50:
                searchDistance = searchDistance + (cellsize * 12)
            if 50 < ncells <= 60:
                searchDistance = searchDistance + (cellsize * 25)
            if ncells > 60:
                searchDistance = searchDistance + (cellsize * 50)
            distanceList.append(searchDistance)
            ncells = ncells + 1
        distance_array = np.array(distanceList)
    else:
        if (searchDistance_max - searchDistance_min >= cellsize):
            distance_array = np.arange(searchDistance_min, searchDistance_max, cellsize)
        else:
            distance_array = np.array([searchDistance_min])
    return distance_array

def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append(int((easting - block_x_min) / x_cellsize))  # col, x
    xy.append(int((northing - block_y_max) / y_cellsize * -1))  # row, y 
    return xy

def plot_it(pts1, pts2, nodeID, a, b, b0, plot_dir):
    """plots the block and topo line"""
    
    import matplotlib.pyplot as plt
    
    x1 = [i[0] for i in pts1]
    y1 = [i[1] for i in pts1]
    
    x2 = [i[0] for i in pts2]
    y2 = [i[1] for i in pts2]    

    plt.plot(x1, y1, 'b--', x2, y2, 'r--')
    plt.savefig(r"{0}\node{1}_a{2}_{3}b_{4}b0.png".format(plot_dir, nodeID, a, b, b0))

def create_block_fc(block_segments, b, block_fc, proj):
    """Creates a poly line feature class of the block segments"""
    
    poly_array = arcpy.Array()
    pnt = arcpy.Point()
    cursorfields = ["BLOCK", "SEGMENT"]   
    
    # Check to see if the block fc exists, if not create it
    if not arcpy.Exists(block_fc):  
        #Create an empty output with the same projection as the input polyline
        arcpy.CreateFeatureclass_management(os.path.dirname(block_fc),
                                            os.path.basename(block_fc),
                                            "POLYLINE","","DISABLED","DISABLED",
                                            proj)
        
        # Add attribute fields
        for f in cursorfields:
            arcpy.AddField_management(block_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(block_fc, ["SHAPE@"] + cursorfields) as cursor:
        for s, segment in enumerate(block_segments):
            for pnt_x, pnt_y in segment:
                pnt.X = pnt_x
                pnt.Y = pnt_y
                poly_array.add(pnt)
            poly = arcpy.Polyline(poly_array)
            cursor.insertRow([poly, b, s])
            poly_array.removeAll() 

def create_blocks(NodeDict, block_size, last_azimuth, searchDistance_max):
    """Returns two lists, one containting the coordinate extent
    for each block that will be itterativly extracted to an array
    and the other containing the start and stop distances for each
    topo line that is within the block extent."""
    
    print("Preparing blocks and topo sampling")
    
    # Create a dictionary to lookup which nodesIDs can be updated after
    # the block has been sampled.
    blockDict = nested_dict()
        
    # Get a list of the nodes, sort them
    nodes = nodeDict.keys()
    nodes.sort()    
    
    topo_list = []
    x_coord_list = []
    y_coord_list = []
    
    for nodeID in nodes:
        node_x = nodeDict[nodeID]["POINT_X"]
        node_y = nodeDict[nodeID]["POINT_Y"]
        streamID = nodeDict[nodeID]["STREAM_ID"]
        z_node = nodeDict[nodeID]["Z_NODE"]
    
        for a in azimuths:
            # calculate x/y coordinates at max search distance
            end_x = ((searchDistance_max * sin(radians(a))) + node_x)
            end_y = ((searchDistance_max * cos(radians(a))) + node_y)
            
            x_coord_list.append(end_x)
            y_coord_list.append(end_y)
            
            topo_list.append([nodeID, streamID, a, z_node, node_x, node_y, end_x, end_y])
    
    # calculate bounding box extent for samples
    x_min = min(x_coord_list)
    x_max = max(x_coord_list)
    y_min = min(y_coord_list)
    y_max = max(y_coord_list)
    
    x_width = int(x_max - x_min + 1)
    y_width = int(y_max - y_min + 1)
    
    block_extents = []
    block_samples = []
    b = 0
     
    # Build blocks
    for x in range(0, x_width, block_size):
        for y in range(0, y_width, block_size):
            
            # Lower left coordinate of block (in map units)
            block_x_min = min([x_min + x, x_max])
            block_y_min = min([y_min + y, y_max])
            # Upper right coordinate of block (in map units)
            block_x_max = min([block_x_min + block_size, x_max])
            block_y_max = min([block_y_min + block_size, y_max])
            
            nodes_to_update = []
            topo_in_block = []
            
            block_segments = (((block_x_min, block_y_max), (block_x_min, block_y_min)),
                              ((block_x_min, block_y_min), (block_x_max, block_y_min)),
                              ((block_x_max, block_y_min), (block_x_max, block_y_max)),
                              ((block_x_max, block_y_max), (block_x_min, block_y_max)))
            
            #--------------------------------------------------------
            # This is an argument for plot_it(). I used it for debugging.
            # It's a tuple of the x/y coords for each block segment.
            block_for_plot = ((block_x_min, block_y_max),
                             (block_x_min, block_y_min),
                             (block_x_min, block_y_min),
                             (block_x_max, block_y_min),
                             (block_x_max, block_y_min),
                             (block_x_max, block_y_max),
                             (block_x_max, block_y_max),
                             (block_x_min, block_y_max))
            #--------------------------------------------------------
             
            # Now start itterating through the topo list to evaluate
            # if any part of the topo line is in the block extent
            for nodeID, streamID, a, z_node, node_x, node_y, end_x, end_y in topo_list:
                
                #--------------------------------------------------------
                # This was used for debugging. 
                # It slows the script WAY down.
                #topo_line = ((node_x, node_y),(end_x, end_y))
                #plot_it(block_for_plot, topo_line, nodeID, a, b, b0, plot_dir)
                #create_topo_line_fc(topo_line, streamID, nodeID, a,
                #                    topo_line_fc, proj)
                #--------------------------------------------------------
                
                contains_node = False
                contains_end = False
                last_sample = False
                
                # check if the node is inside the block
                if (block_x_min <= node_x <= block_x_max and
                    block_y_min <= node_y <= block_y_max):
                    contains_node = True
                
                # check if the topo line end point is inside the block
                if (block_x_min <= end_x <= block_x_max and
                    block_y_min <= end_y <= block_y_max):
                    contains_end = True
                    
                # check if this is the last block to process for this node
                if a == last_azimuth:
                    if a == 45:
                        # we have to get the coordinate at 
                        # the far upper right corner
                        searchDistance_last = (hypot(searchDistance_max, searchDistance_max))
                        last_x = ((searchDistance_last * sin(radians(a))) + node_x)
                        last_y = ((searchDistance_last * cos(radians(a))) + node_y)
                    else:
                        last_x = end_x
                        lasy_y = end_y
                        
                    if (block_x_min <= last_x <= block_x_max and
                        block_y_min <= lasy_y <= block_y_max):
                        
                        last_sample = True      
            
                # check if the entire topo line segment is inside the block
                if contains_node and contains_end:
                    block_search_start = 0
                    block_search_end = searchDistance_max
                    
                    topo_in_block.append([nodeID, streamID, a,
                                          z_node, 
                                          node_x, node_y,
                                          end_x, end_y,
                                          block_search_start,
                                          block_search_end])
                    
                    if last_sample: nodes_to_update.append(nodeID)
                
                # check if and where the topo segment cross the block    
                else:
                    distance = []
                    
                    # check each block segment for an intersection
                    for i, block_segment in enumerate(block_segments):
                        
                        intersects, inter1_x, inter1_y, inter2_x, inter2_y = find_intersection(block_segment[0],
                                                            block_segment[1],
                                                            (node_x, node_y),
                                                            (end_x, end_y), True)
                        
                        # if there is an intersection calculate the 
                        # distance in fc units from node to 
                        # that intersection
                        if intersects:
                            if a in [0, 180]:
                                # need to use the y directon for these
                                distance.append((inter1_y - node_y) / cos(radians(a)))
                            else:
                                distance.append((inter1_x - node_x) / sin(radians(a)))

                    if (len(distance) == 1 and
                        (0 < distance[0] < searchDistance_max) and
                        (contains_node or contains_end)):
                        # one intersection
                        if contains_node:
                            # This will be changed from zero to 
                            # the cell size in build_search_array()
                            # Should probably just change it here
                            block_search_start = 0 
                            block_search_end = distance[0]
                        elif contains_end:
                            # end of the topo line
                            block_search_start = distance[0]
                            block_search_end = searchDistance_max
                        
                        # part of the topo line is in the block, add it
                        topo_in_block.append([nodeID, streamID, a,
                                              z_node, 
                                              node_x, node_y,
                                              end_x, end_y,
                                              block_search_start,
                                              block_search_end])
                        
                        if last_sample and not nodeDict[nodeID]["updated"]:
                            nodes_to_update.append(nodeID)
                            nodeDict[nodeID]["updated"] = True
                    
                    elif len(distance) > 1:
                        # two intersections, crosses the block
                        # two intersections, end or start on block line
                        # three intersections, collienar over length of block
                        # three intersections, collienar w/ end or start on block line
                        block_search_start = min(i for i in distance if i is not None)
                        block_search_end = max(distance)
                    
                        topo_in_block.append([nodeID, streamID, a,
                                              z_node, 
                                              node_x, node_y,
                                              end_x, end_y,
                                              block_search_start,
                                              block_search_end])
                    
                        if last_sample and contains_end and not nodeDict[nodeID]["updated"]:
                            nodes_to_update.append(nodeID)
                            nodeDict[nodeID]["updated"] = True
                    
                    elif last_sample and not contains_end and not nodeDict[nodeID]["updated"]:
                        nodes_to_update.append(nodeID)
                        nodeDict[nodeID]["updated"] = True
                        
                    del distance[:]
                    
            if topo_in_block:
                
                # order 0 left,      1 bottom,    2 right,     3 top
                blockDict[b]["extent"] = (block_x_min, block_y_min,
                                          block_x_max, block_y_max)
                blockDict[b]["samples"] = topo_in_block
                blockDict[b]["nodes_to_update"] = nodes_to_update
                
                #--------------------------------------------------
                # Creates a feature class of the blocks.  
                #create_block_fc(block_segments, b, block_fc, proj)
                #--------------------------------------------------
                
            b = b + 1
    return blockDict

def find_intersection(a, b, c, d, check_collinear=True):
    """Calculates 2D intersection coordinates of segments a-b and c-d.
    a,b,c,d are tuples in the form of (x,y). If checking for collinearity
    then segments that overlap or touch but do not cross will be
    considered an interesection. Returns the min and max points of
    intersection for overlap and the same points for intersections."""
    
    # a-b = block segment
    # c-d = topo line
    
    Dx_Cx = d[0] - c[0]
    Ay_Cy = a[1] - c[1]
    Dy_Cy = d[1] - c[1]
    Ax_Cx = a[0] - c[0]
    Bx_Ax = b[0] - a[0]
    By_Ay = b[1] - a[1]
    
    Cx_Bx = c[0] - b[0]
    Dx_Ax = d[0] - a[0]
    Cy_By = c[1] - b[1]
    Dy_Ay = d[1] - a[1]    
    
    numerator_a = Dx_Cx * Ay_Cy - Dy_Cy * Ax_Cx
    numerator_b = Bx_Ax * Ay_Cy - By_Ay * Ax_Cx
    denominator = Dy_Cy * Bx_Ax - Dx_Cx * By_Ay
        
    if (check_collinear and numerator_a == 0 and numerator_b == 0 and denominator == 0):
        # Lines are collinear
        
        # if the signs are different the 
        # segments have some overlap
        overlap_x = (Cx_Bx < 0) != (Dx_Ax < 0)
        overlap_y = (Cy_By < 0) != (Dy_Ay < 0)
        
        point_overlap = (a == d or b == c)
        
        # if any are True there is an intersection
        if (overlap_x or overlap_y or point_overlap):
            # There is overlap 
            x = sorted((a[0], b[0], c[0], d[0]))
            y = sorted((a[1], b[1], c[1], d[1]))
            
            # min          
            ixa = x[1]
            iya = y[1]
            
            # max
            ixb = x[2]
            iyb = y[2]
            
            return True, ixa, iya, ixb, iyb
        return False, None, None, None, None
    
    if denominator == 0:
        # Lines are parallel or 
        # maybe collinear if check_collinear=False
        
        return False, None, None, None, None
    
    u_a = numerator_a / denominator
    u_b = numerator_b / denominator
    
    if (u_a >= 0) and (u_a <= 1) and (u_b >= 0) and (u_b <= 1):
        # segments intersect
        ixa = a[0] + Bx_Ax * u_a
        iya = a[1] + By_Ay * u_a
        
        ixb = c[0] + Dx_Cx * u_b
        iyb = c[1] + Dy_Cy * u_b
        
        return True, ixa, iyb, ixa, iyb
    return False, None, None, None, None

def get_topo_angles(nodeDict, block_extent, block_samples, z_raster, azimuthdisdict, searchDistance_max_m, con_z_to_m):
    """This gets the maximum topographic angle and other informaiton for
    each topo line within the block. The data is saved to the nodeDict
    as a list."""
    
    nodata_to_value = -9999 / con_z_to_m
    
    x_cellsize = arcpy.Describe(z_raster).meanCellWidth
    y_cellsize = arcpy.Describe(z_raster).meanCellHeight
    
    # localize the block extent values
    block_x_min = block_extent[0]
    block_y_min = block_extent[1]
    block_x_max = block_extent[2]
    block_y_max = block_extent[3]    

    # Get the coordinates extent of the input raster
    # this could be in main so it doesn't have to run each time
    raster_x_min = float(arcpy.GetRasterProperties_management(z_raster, "LEFT").getOutput(0))
    raster_y_min = float(arcpy.GetRasterProperties_management(z_raster, "BOTTOM").getOutput(0))
    raster_x_max = float(arcpy.GetRasterProperties_management(z_raster, "RIGHT").getOutput(0))
    raster_y_max = float(arcpy.GetRasterProperties_management(z_raster, "TOP").getOutput(0))
      
    # Calculate the X and Y offset from the upper left node 
    # coordinates bounding box    
    x_minoffset = (block_x_min - raster_x_min)%x_cellsize
    y_minoffset = (block_y_min - raster_y_min)%y_cellsize
    x_maxoffset = (raster_x_max - block_x_max)%x_cellsize 
    y_maxoffset = (raster_y_max - block_y_max)%y_cellsize
    
    # adjust so the coordinates are at the raster cell corners
    block_x_min = block_x_min - x_minoffset 
    block_y_min = block_y_min - y_minoffset
    block_x_max = block_x_max + x_maxoffset
    block_y_max = block_y_max + y_maxoffset
    
    # Get the lower left cell center coordinate. This is for ESRI's
    # RastertoNumpyArray function which defaults to the adjacent 
    # lower left cell
    block_x_min_center = block_x_min + (x_cellsize / 2)
    block_y_min_center = block_y_min + (y_cellsize / 2)    
    
    # calculate the number or cols/ros from the lower left
    ncols = max([int(ceil((block_x_max - block_x_min)/ x_cellsize)), 1])
    nrows = max([int(ceil((block_y_max - block_y_min)/ y_cellsize)), 1])
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        z_array = arcpy.RasterToNumPyArray(z_raster, arcpy.Point(block_x_min_center, block_y_min_center),
                                       ncols, nrows, nodata_to_value)
    except:
        tbinfo = traceback.format_exc()
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
        sys.exit(pymsg)    
    
    # convert array values to meters if needed
    z_array = z_array * con_z_to_m
    
    topo_samples = []    
    if z_array.max() > -9999:
        # There is at least one pixel of data
        for (nodeID, streamID, a, z_node,
             node_x, node_y, end_x, end_y,
             block_search_start, block_search_end) in block_samples:
            
            # create list of distance movements along the topo line 
            # from the block edges
            # this is in units of the fc
            cellsize = azimuthdisdict[a]
            distance_array = build_search_array(block_search_start,
                                                  block_search_end,
                                                  cellsize,
                                                  use_skippy=False)
            z_topo_list = []
            for distance in distance_array:
                pt_x = ((distance * sin(radians(a))) + node_x)
                pt_y = ((distance * cos(radians(a))) + node_y)
               
                xy = coord_to_array(pt_x, pt_y, block_x_min, block_y_max, x_cellsize, y_cellsize)
                z_topo_list.append(z_array[xy[1], xy[0]])
    
            # convert distances to meters
            distance_array_m = distance_array * con_to_m            
            z_topo_array = np.array(z_topo_list)
            # Calculate the topo angles along the topo line
            angle_array = np.degrees(np.arctan((z_topo_array - z_node) / distance_array_m))
            # remove the off raster samples
            naindex = np.where(z_topo_array < -9998)
            for x in naindex[0]: angle_array[x] = -9999            
            # Find the max topo angle
            topoAngle = angle_array.max()
            # array index at the max topo angle 
            arryindex = np.where(angle_array==topoAngle)[0][0]
            z_topo = z_topo_array[arryindex]
            # elevation change between topo angle location and node elevation
            z_change = z_topo - z_node
            # distance from the node to topo angle location in units of fc
            topoAngleDistance = distance_array[arryindex]
            topoAngle_x = (topoAngleDistance * sin(radians(a))) + node_x
            topoAngle_y = (topoAngleDistance * cos(radians(a))) + node_y
            topoAngleDistance_m = topoAngleDistance * con_to_m
            off_rastersamples = (z_topo_array < -9998).sum()
            
            topo_samples.append([topoAngle_x, topoAngle_y,
                                  topoAngle_x, topoAngle_y,
                                  streamID, nodeID, a,
                                  topoAngle, z_topo,
                                  z_node, z_change,
                                  topoAngleDistance_m,
                                  searchDistance_max_m,
                                  off_rastersamples])
        
    return topo_samples

#enable garbage collection
gc.enable()

try:
    print("Step 4: Measure Topographic Angles") 
    
    #keeping track of time
    startTime= time.time()
    
    # Check if the node fc exists
    if not arcpy.Exists(nodes_fc):
        arcpy.AddError("\nThis output does not exist: \n" +
                       "{0}\n".format(nodes_fc))
        sys.exit("This output does not exist: \n" +
                 "{0}\n".format(nodes_fc))     
        
    # Check if the topo output exists and delete if needed
    if arcpy.Exists(topo_fc) and overwrite_data:
        arcpy.Delete_management(topo_fc)
    
    # Check if the block fc exists and delete or throw an error       
    if arcpy.Exists(block_fc):
        if overwrite_data:
            arcpy.Delete_management(block_fc)
        else:
            arcpy.AddMessage("\nThis output already exists: \n" +
                           "{0}\n".format(block_fc) + 
                           "overwrite data = False. New data will be " +
                           "appended to the existing feature class.")
            print("This output already exists: \n" +
                           "{0}\n".format(block_fc) + 
                           "overwrite data = False. New data will be " +
                           "appended to the existing feature class.")
    
    if overwrite_data: 
        env.overwriteOutput = True
    else:
        env.overwriteOutput = False    
    
    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference
    proj_ele = arcpy.Describe(z_raster).spatialReference

    # Check to make sure the raster and input 
    # points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not have ".format(nodes_fc,z_raster)+
                       "the same projection. Please reproject your data.")
        sys.exit("Input points and elevation raster do not have the "+
                 "same projection. Please reproject your data.")

    # Get the units conversion factors
    con_z_to_m = from_z_units_to_meters_con(z_units)    
    con_to_m = to_meters_con(nodes_fc)
    con_from_m = from_meters_con(nodes_fc)
    searchDistance_max_m = searchDistance_max_km * 1000 # in meters
    
    # search distance in units of the fc
    searchDistance_max = int(con_from_m * searchDistance_max_m)
    
    # convert block size from km to meters to units of the node fc
    # in the future block size should be estimated based on availiable memory
    # memorysize = datatypeinbytes*nobands*block_size^2
    # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
    if block_size in ["#", ""]:
        block_size = int(con_from_m * 5000)
    else:
        block_size = int(con_from_m * block_size * 1000)    

    # Get the elevation raster cell size in units of the raster
    x_cellsize = arcpy.Describe(z_raster).meanCellWidth
    y_cellsize = arcpy.Describe(z_raster).meanCellHeight
    
    if topo_directions == 2: # All directions
        azimuths = [45,90,135,180,225,270,315,365]
        last_azimuth = 45
    else:        
        azimuths = [270,180,90]
        last_azimuth = 90

    azimuthdict = {45:"TOPO_NE",90:"TOPO_E",135:"TOPO_SE",
                   180:"TOPO_S",225:"TOPO_SW",270:"TOPO_W",
                   315:"TOPO_NW",365:"TOPO_N"}
        
    # build a dictionary of the cell size in azimuth directions
    disXY = (hypot(x_cellsize,y_cellsize))
    disX = (x_cellsize)
    disY = (y_cellsize)
    
    azimuthdisdict = {45:disXY,90:disX,135:disXY,180:disY,225:disXY,270:disX,315:disXY,365:disY}
    del(disX,disXY,disY)
    
    # Build the topo field names
    addFields = [azimuthdict[a] for a in azimuths]    
    
    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)     

    # Check to see if the field exists and add it if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")
        
    # Read the feature class data into a nested dictionary
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
    
    # Build the blockDict
    blockDict = create_blocks(nodeDict, block_size, last_azimuth,
                              searchDistance_max)    

    # Itterate through each block
    total_samples = 0
    blockIDs = blockDict.keys()
    blockIDs.sort()

    for p, blockID in enumerate(blockIDs):
        block_extent = blockDict[blockID]["extent"]
        block_samples = blockDict[blockID]["samples"]
        block_samples.sort()
        
        print("Processing block {0} of {1}".format(p + 1, len(blockIDs)))
        
        # calculate coordinates along the 
        # portion of the topo line in the block, 
        # convert raster to array, sample the raster
        # calculate the topo angles and other info
        topo_samples = get_topo_angles(nodeDict, block_extent , block_samples,
                                   z_raster, azimuthdisdict,
                                   searchDistance_max, con_z_to_m)
        if topo_samples:
            # Update the nodeDict
            for sample in topo_samples:
                nodeID = sample[5]
                a = sample[6]
                topoAngle = sample[7]
                
                # Create a key to hold the topo list info for this block
                topo_key = azimuthdict[a] + "_list"
                
                if azimuthdict[a] in nodeDict[nodeID]:
                    if nodeDict[nodeID][azimuthdict[a]] < topoAngle:
                        nodeDict[nodeID][azimuthdict[a]] = topoAngle
                        nodeDict[nodeID][topo_key] = sample
                        
                else:
                    nodeDict[nodeID][azimuthdict[a]] = topoAngle
                    nodeDict[nodeID][topo_key] = sample
                    
            del topo_samples
        
        # Check if any nodes can be updated in the node and topo fc
        if blockDict[blockID]["nodes_to_update"]:
            nodes_to_update = blockDict[blockID]["nodes_to_update"]
            
            # Write the topo data to the TTools point feature class
            update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_update)
            
            # Build/add to the output topo feature class
            topo_list = []
            for nodeID in nodes_to_update:
                for field in addFields:
                    topo_key = field + "_list"
                    topo_list.append(nodeDict[nodeID][topo_key])
                    # delete some data
                    nodeDict[nodeID].pop(field, None)
                    nodeDict[nodeID].pop(topo_key, None)
                    nodeDict[nodeID].pop("updated", None)
            update_topo_fc(topo_list, topo_fc, nodes_fc,
                           nodes_to_update, overwrite_data, proj)
            
            del topo_list
            gc.collect()
    
    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / len(nodeDict.keys())).microseconds
    print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))
    #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))

# For arctool errors
except arcpy.ExecuteError:
    msgs = arcpy.GetMessages(2)
    #arcpy.AddError(msgs)
    print(msgs)

# For other errors
except:
    tbinfo = traceback.format_exc()

    pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" +str(sys.exc_info()[1])
    msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

    #arcpy.AddError(pymsg)
    #arcpy.AddError(msgs)

    print(pymsg)
    print(msgs)