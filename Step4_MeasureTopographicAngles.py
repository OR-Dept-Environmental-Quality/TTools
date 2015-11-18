########################################################################
# TTools
# Step 4: Measure Topographic Angles - v 0.961
# Ryan Michie

# Measure_Topographic_Angles will take an input point feature 
# (from Step 1) and calculate the maximum topographic elevation
# and the the slope angle from each node in different directions.

# INPUTS
# 0: Input TTools point feature class(nodes_fc)
# 1: input the directions to 
#    sample (directions) 1. [W,S,E], 2. [NE,E,SE,S,SW,W,NW,N]
# 2: input the maximum km distance to search (searchDistance_max_km)
# 3: input elevation raster (z_raster)
# 4: input elevation raster z units (z_units) 
#     1. "Feet", 2. "Meters" 3. "Other"
# 5: output sample point file name/path (topo_fc)
# 6: input flag if existing data can be over
#     written (overwrite_data) 1. True, 2. False

# OUTPUTS
# 0. point feature class (edit nodes_fc) - Added fields with topographic 
#     shade angles for each direction at each node
# 1. point feature class (new) - point for each x/y location 
#     of the maximum elevation.

# Future Updates
# partial runs don't correctly update the topo sample fc
# update after each sample so data isn't lost after larger runs
# check to see if the outputs exist before running process
# faster for larger runs
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
directions = 1
searchDistance_max_km = 10
skippy = False
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
z_units = "Meters"
topo_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\topo_samples"
overwrite_data = True
# End Fill in Data
# ----------------------------------------------------------------------

# Parameter fields for python toolbox
#nodes_fc = parameters[0].valueAsText
#directions = parameters[1].valueAsText # Needs to be a long
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
    nodeDict = nested_dict()
    incursorFields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]

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
        if overwrite_data is True:
            for row in Inrows:
                nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                nodeDict[row[0]][row[1]]["POINT_X"] = row[3]
                nodeDict[row[0]][row[1]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # if the data is null or zero (0 = default for shapefile),
                # it is retreived and will be overwritten.
                if row[5] is None or row[5] == 0 or row[5] < -9998:
                    nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                    nodeDict[row[0]][row[1]]["POINT_X"] = row[3]
                    nodeDict[row[0]][row[1]]["POINT_Y"] = row[4]
    if len(nodeDict) == 0:
        sys.exit("The fields checked in the input point feature class "+
                 "have existing data. There is nothing to process. Exiting")
            
    return(nodeDict)

def create_topo_fc(pointList, topo_fc, nodes_fc, proj):
    """Creates the output topo point feature class
    using the data from the nodes list"""
    #arcpy.AddMessage("Exporting Data")
    print("Exporting data to topo sample feature class")
        
    #Create an empty output with the same projection as the input polyline
    cursorfields = ["POINT_X","POINT_Y","STREAM_ID","NODE_ID",
                    "AZIMUTH","TOPOANGLE","TOPO_ELE","NODE_ELE",
                    "ELE_CHANGE","TOPODIS","SEARCHDIS","NA_SAMPLES"]
    arcpy.CreateFeatureclass_management(os.path.dirname(topo_fc),
                                        os.path.basename(topo_fc),
                                        "POINT","","DISABLED","DISABLED",
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
            arcpy.AddField_management(topo_fc, f, sid_type,
                                      sid_precision, sid_scale,
                                      sid_length, "", "NULLABLE",
                                      "NON_REQUIRED")
        else:
            arcpy.AddField_management(topo_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(topo_fc, ["SHAPE@X","SHAPE@Y"] +
                               cursorfields) as cursor:
        for row in pointList:
            cursor.insertRow(row)
            
def update_topo_fc(pointList, topo_fc, nodes_fc, proj):
    """Creates and updates the output topo point feature class
    using the data from the nodes list"""
    #arcpy.AddMessage("Exporting Data")
    #print("Exporting data to topo sample feature class")
    
    # Check to see if the topo fc exists, if not create it
    if not arcpy.Exists(topo_fc):  
        #Create an empty output with the same projection as the input polyline
        cursorfields = ["POINT_X","POINT_Y","STREAM_ID","NODE_ID",
                        "AZIMUTH","TOPOANGLE","Z_TOPO","Z_NODE",
                        "Z_CHANGE","TOPODIS","SEARCHDIS","NA_SAMPLES"]
        arcpy.CreateFeatureclass_management(os.path.dirname(topo_fc),
                                            os.path.basename(topo_fc),
                                            "POINT","","DISABLED","DISABLED",
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
                arcpy.AddField_management(topo_fc, f, sid_type,
                                          sid_precision, sid_scale,
                                          sid_length, "", "NULLABLE",
                                          "NON_REQUIRED")
            else:
                arcpy.AddField_management(topo_fc, f, "DOUBLE", "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(topo_fc, ["SHAPE@X","SHAPE@Y"] +
                               cursorfields) as cursor:
        for row in pointList:
            cursor.insertRow(row)
            
def update_nodes_fc(nodeDict, nodes_fc, addFields): 
    """Updates the input point feature class
    with data from the nodes dictionary"""
    print("Updating input point feature class")

    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)     

    # Check to see if the field exists and add it if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "",
                                      "", "", "NULLABLE", "NON_REQUIRED")   

    with arcpy.da.UpdateCursor(nodes_fc,["STREAM_ID","NODE_ID"] +
                               addFields) as cursor:
        for row in cursor:
            for f, field in enumerate(addFields):
                streamID = row[0]
                nodeID =row[1]
                row[f+2] = nodeDict[streamID][nodeID][field]
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

def build_search_array(searchDistance_min_m, searchDistance_max_m, cellsize, skippy):
    """Build a numpy array from the minimum to the max
    search distance by increments of the cellsize or
    via the skippy algorithm"""
    if skippy is True:
        # This is a modified version of the skippy 
        # algorithm from Greg Pelletier
        searchDistance_m = searchDistance_min_m
        distanceList = [searchDistance_min_m]
        i = 1
        while not searchDistance_m > searchDistance_max_m:
            if i <= 10:
                searchDistance_m = searchDistance_m + (cellsize)
            if 10 < i <= 20:
                searchDistance_m = searchDistance_m + (cellsize * 3)
            if 20 < i <= 40:
                searchDistance_m = searchDistance_m + (cellsize * 6)
            if 40 < i <= 50:
                searchDistance_m = searchDistance_m + (cellsize * 12)
            if 50 < i <= 60:
                searchDistance_m = searchDistance_m + (cellsize * 25)
            if i > 60:
                searchDistance_m = searchDistance_m + (cellsize * 50)
            distanceList.append(searchDistance_m)
            i = i + 1
        distance_array = np.array(distanceList)
    else:
        distance_array = np.arange(searchDistance_min_m, searchDistance_max_m, cellsize)
    return distance_array

def coord_to_array(easting, northing, bbox_upper_left):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append((easting - bbox_upper_left[0]) / bbox_upper_left[2])  # col, x
    xy.append((northing - bbox_upper_left[1]) / bbox_upper_left[3] * -1)  # row, y    
    return xy

def create_coord_list(origin_x, origin_y, azimuth, azimuthdisdict, con_from_m):
    
    coordList = []
    for searchdistance in azimuthdisdict[a]:
        # Calculate the x and y coordinate of the landcover 
        # sample location in the spatial units of the inputFC                
        topo_x = (searchdistance * con_from_m * sin(radians(a))) + origin_x
        topo_y = (searchdistance * con_from_m * cos(radians(a))) + origin_y

        # Add the all the coordinates to the list
        coordList.append([topo_x, topo_y])
    return coordList

def get_topo_angles(coordList, z_raster, a, azimuthdisdict, con_z_to_m):
    """This gets the maximum topographic angle for a set of
    sample coordinates in an azimuth direction. The data is
    output as a list. This version is slower but will not crash
    as the arrays are small"""
    
    x_cellsize = arcpy.Describe(z_raster).meanCellWidth
    y_cellsize = arcpy.Describe(z_raster).meanCellHeight
    
    # Get the coordinates of the upper left cell corner of the input raster
    top_y = float(arcpy.GetRasterProperties_management(z_raster, "TOP").getOutput(0))
    left_x = float(arcpy.GetRasterProperties_management(z_raster, "LEFT").getOutput(0))
    
    # calculate lower left corner and nrows/cols for the bounding box
    # first transpose the list so x and y coordinates are in the same list
    tlist = [list(x) for x in zip(*coordList)]
    
    # calculate the buffer distance (in raster spatial units) to
    # add to the raster bounding box when extracting to an array
    buffer = x_cellsize * 1
        
    x_min = min(tlist[0]) - buffer
    y_min = min(tlist[1]) - buffer
    y_max = max(tlist[1]) + buffer
    
    # Calculate the X and Y offset from the upper left 
    # node coordinates bounding box
    x_minoffset = ((left_x - x_min)%x_cellsize) - x_cellsize
    y_maxoffset = (top_y - y_max)%y_cellsize
     
    ncols = int(ceil((max(tlist[0]) + buffer - x_min) / x_cellsize)) + 1
    nrows = int(ceil((y_max - y_min) / y_cellsize)) + 1
    bbox_lower_left = arcpy.Point(x_min, y_min) # must be in raster map units
    bbox_upper_left = [x_min + x_minoffset, y_max + y_maxoffset, x_cellsize, y_cellsize]
    nodata_to_value = -9999 / con_z_to_m
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        z_array = arcpy.RasterToNumPyArray(z_raster, bbox_lower_left,
                                           ncols, nrows, nodata_to_value)
    except:
        import traceback
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the search distance"       
        sys.exit(pymsg)    
    
    # convert array values to meters if needed
    z_array = z_array * con_z_to_m             
    
    z_topo_list = []
    
    #print("Extracting raster values")
    for i in range(0,len(coordList)):
        xy = coord_to_array(coordList[i][0], coordList[i][1], bbox_upper_left)
        z_topo_list.append(z_array[xy[1], xy[0]])
        #topoList[i].append(z_array[xy[1], xy[0]])
    
    z_topo_array = np.array(z_topo_list)
    # in meters
    distance_array = azimuthdisdict[a]
    angle_array = np.degrees(np.arctan((z_topo_array - z_topo_array[0]) / distance_array))
    # remove the off raster samples
    naindex = np.where(z_topo_array < -9998)
    for x in naindex[0]: angle_array[x] = -9999            
    # Find the max topo angle
    topoAngle = angle_array.max()
    # array index at the max topo angle 
    arryindex = np.where(angle_array==topoAngle)[0][0]
    z_topo = z_topo_array[arryindex]
    # elevation change between topo angle location and node elevation
    z_change = z_topo - z_topo_array[0]
    # distance from the node to topo angle location (meters)
    topoAngleDistance = azimuthdisdict[a][arryindex]
    searchDistance_m = azimuthdisdict[a].max()
    topoAngle_x = (topoAngleDistance * con_from_m * sin(radians(a))) + coordList[0][0]
    topoAngle_y = (topoAngleDistance * con_from_m * cos(radians(a))) + coordList[0][1]
    off_rastersamples = (z_topo_array < -9998).sum()
    
    topo_data_list = [topoAngle_x, topoAngle_y, topoAngle_x, topoAngle_y, 
                    streamID, nodeID, a, topoAngle, z_topo, z_topo_array[0], 
                    z_change, topoAngleDistance, searchDistance_m, off_rastersamples]
    return topo_data_list

#enable garbage collection
gc.enable()

try:
    print("Step 4: Measure Topographic Angles") 
    
    #keeping track of time
    startTime= time.time()
    
    # Check if the output exists
    if not arcpy.Exists(nodes_fc):
        arcpy.AddError("\nThis output does not exist: \n" +
                       "{0}\n".format(nodes_fc))
        sys.exit("This output does not exist: \n" +
                 "{0}\n".format(nodes_fc))     
    
    if overwrite_data is True: 
        env.overwriteOutput = True
    else:
        env.overwriteOutput = False    
    
    # Check if the topo output exists and delete or throw an error
    if arcpy.Exists(topo_fc):
        if overwrite_data is True:
            arcpy.Delete_management(topo_fc)
        else:
            arcpy.AddError("\nThis output already exists: \n" +
                           "{0}\n".format(topo_fc) + 
                           "Please rename your output or"
                           "choose to overwrite your data.")
            sys.exit("This output already exists: \n" +
                           "{0}\n".format(topo_fc) + 
                           "Please rename your output or choose"
                           "to overwrite your data.")            
            

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

    # Get the elevation raster cell size in meters
    x_cellsize = arcpy.Describe(z_raster).meanCellWidth * con_to_m
    y_cellsize = arcpy.Describe(z_raster).meanCellHeight * con_to_m
    
    if directions == 2: # All directions
        azimuths = [45,90,135,180,225,270,315,365]
    else:        
        azimuths = [270,180,90]

    azimuthdict = {45:"TOPO_NE",90:"TOPO_E",135:"TOPO_SE",
                   180:"TOPO_S",225:"TOPO_SW",270:"TOPO_W",
                   315:"TOPO_NW",365:"TOPO_N"}
        
    # Build a numpy array from the minimum to the max 
    # search distance by increments of the cellsize
    # using 0.000001 to avoid divide by zero errors
    disX = build_search_array(0.000001, searchDistance_max_m, x_cellsize, skippy)
    disXY = build_search_array(0.000001, searchDistance_max_m, hypot(x_cellsize,y_cellsize), skippy)
    disY = build_search_array(0.000001, searchDistance_max_m, y_cellsize, skippy)
    azimuthdisdict = {45:disXY,90:disX,135:disXY,180:disY,225:disXY,270:disX,315:disXY,365:disY}
    del(disX,disXY,disY)

    # Build the topo field names
    addFields = [azimuthdict[a] for a in azimuths]
        
    # Read the feature class data into a nested dictionary
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
    
    topoList = []
    n = 1
    for streamID in nodeDict:
        
        nodes = nodeDict[streamID].keys()
        nodes.sort()

        for i, nodeID in enumerate(nodes):
            print("Processing stream {0} of {1}. {2:.0%} complete".format(n, len(nodeDict), i / len(nodes)))
            origin_x = nodeDict[streamID][nodeID]["POINT_X"]
            origin_y = nodeDict[streamID][nodeID]["POINT_Y"]
            for a in azimuths:
                coordList = create_coord_list(origin_x, origin_y, a, azimuthdisdict, con_from_m)
                topo_data_list = get_topo_angles(coordList, z_raster, a, azimuthdisdict, con_z_to_m)
                topoList.append(topo_data_list)   
    
    # Update the nodeDict
    for row in topoList:
        k = azimuthdict[row[6]]
        nodeDict[row[4]][row[5]][k] = row[7]    
    
    update_nodes_fc(nodeDict, nodes_fc, addFields)
    create_topo_fc(topoList, topo_fc, nodes_fc, proj)
    
    gc.collect()

    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / len(topoList) / len(azimuths)).microseconds
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