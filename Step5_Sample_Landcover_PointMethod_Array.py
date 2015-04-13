########################################################################
# TTools
# Step 5: Sample Landcover - Star Pattern, Point Method v 0.99
# Ryan Michie

# Sample_Landcover_PointMethod will take an input point 
# feature (from Step 1) and sample input landcover rasters in a user 
# specificed number of cardinal directions with point samples spaced 
# at a user defined distance moving away from the stream.

# General steps include:
# 1. read nodes feature class
# 2. build a lcsample output feature
# 3. calculate all sample easting/northing coordinates, iterate by stream
# 4. calculate the bounding box of the coordinates + extra
# 5. import elevation/veght rasters and convert to arrays for the 
#     stream bounding box
# 6. sample array
# 7. save to dictionary
# 8. output to lc sample feature class
# 9. save to nodes feature class

# INPUTS
# 0: Input TTools point feature class (nodes_fc)
# 1: input number of transects per node (trans_count)
# 2: input number of samples per transect (transsample_count)
# 3: input The distance between transect samples (transsample_distance)
# 4: use heatsource 8 methods 1. True 2. False (heatsource8)
# 5: input landcover code or height raster (lc_raster)
# 6: input (optional) landcover height z units 
#     1. "Feet", 2. "Meters" 3. "None" Float (lc_units)
# 7: input (optional) landcover data type. 
#     1."CanopyCover", or 2."LAI" (canopy_data_type)
# 8: input (optional) canopy cover or LAI raster (canopy_raster)
# 9: input (optional) k coeffcient raster (k_raster)
# 10: input (optional) overhang raster (oh_raster)
# 11: input elevation raster (z_raster)
# 12: input elvation raster z units (z_units) 
#      1. "Feet", 2. "Meters" 3. "Other"
# 13: output sample point file name/path (lc_point_fc)
# 14: input stream km distance to process within each array (block_size)
# 15: input flag if existing data can be over 
#     written (overwrite_data) 1. True, 2. False

# OUTPUTS
# 0. point feature class (edit nodes_fc) - added fields with 
#     Landcover and elevation data for each azimuth direction at each node
# 1. point feature class (new) - point at each x/y sample and 
#     the sample raster values

# Future Updates
# include stream sample in transect count (True/False)
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
from datetime import timedelta
from math import radians, sin, cos, ceil
from collections import defaultdict
import numpy
import arcpy
from arcpy import env

env.overwriteOutput = True

# Parameter fields for python toolbox
#nodes_fc = parameters[0].valueAsText
#trans_count = parameters[1].valueAsText # LONG
#transsample_count = parameters[2].valueAsText # LONG
#transsample_distance = parameters[3].valueAsText # LONG
#heatsource8 = parameters[4].valueAsText # True/False
#lc_raster = parameters[5].valueAsText # This is either landcover height or codes
#lc_units = parameters[6].valueAsText # OPTIONAL One of these: 1. "Feet", 2. "Meters" 3. Float
#canopy_data_type = parameters[7].valueAsText # OPTIONAL One of these: 1."CanopyCover", or 2."LAI"
#canopy_raster = parameters[8].valueAsText # OPTIONAL This is either canopy cover or LAI raster
#k_raster = parameters[9].valueAsText # OPTIONAL The k value raster for LAI
#oh_raster = = parameters[10].valueAsText # OPTIONAL
#z_raster = parameters[11].valueAsText
#z_units = parameters[12].valueAsText One of these: 1. "Feet", 2. "Meters" 3. Float
#lc_point_fc = parameters[13].valueAsText
#block_size = parameters[14].valueAsText
#overwrite_data = parameters[15].valueAsText # True/False

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes_ogic"
trans_count = 8 
transsample_count = 4 # does not include a sample at the stream node
transsample_distance = 8
heatsource8 = False
lc_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_vght_m_mosaic"
lc_units = "Meters"
canopy_data_type = "#" # OPTIONAL This is either 1. "CanopyCover", or 2."LAI"
canopy_raster = "#" # OPTIONAL This is either canopy cover or a LAI raster
k_raster = "#" # OPTIONAL This is the k value for LAI
oh_raster = "#" # OPTIONAL This is the overhang raster
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
z_units = "Meters"
lc_point_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\LC_samplepoint_ogic"
block_size = "#" # OPTIONAL defualt to 5
overwrite_data = True
# End Fill in Data
# ----------------------------------------------------------------------

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

    # Check to see if the last field exists if yes add it. 
    # Grabs last field becuase often the first field, emergent, is zero
    if overwrite_data is False and (addFields[len(addFields)-1] in existingFields) is True:
        incursorFields.append(addFields[len(addFields)-1])
    else:
        overwrite_data = True

    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc, incursorFields,"",proj) as Inrows:
        if overwrite_data is True:
            for row in Inrows:
                nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                nodeDict[row[0]][row[1]]["POINT_X"] = row[3]
                nodeDict[row[0]][row[1]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[5] is None or row[5] == 0 or row[5] < -9998:
                    nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                    nodeDict[row[0]][row[1]]["POINT_X"] = row[3]
                    nodeDict[row[0]][row[1]]["POINT_Y"] = row[4]
    
    if len(nodeDict) == 0:
        sys.exit("The fields checked in the input point feature class " +
                 "have existing data. There is nothing to process. Exiting")
              
    return nodeDict

def create_lc_point_fc(pointList, LCFields, lc_point_fc, nodes_fc, proj):
    """Creates the output landcover sample point feature
    class using the data from the point list"""
    print("Exporting data to land cover sample feature class")

    arcpy.CreateFeatureclass_management(os.path.dirname(lc_point_fc),
                                        os.path.basename(lc_point_fc),
                                        "POINT","","DISABLED","DISABLED",proj)
    
    # Determine Stream ID field properties
    sid_type = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].type
    sid_precision = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].precision
    sid_scale = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].scale
    sid_length = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].length    

    cursorfields = ["POINT_X","POINT_Y"] +["STREAM_ID","NODE_ID",
                                            "AZIMUTH","TRANSNUM",
                                            "SAMPLENUM"] +LCFields

    # Add attribute fields # TODO add dictionary of field types 
    # so they aren't all double
    for f in cursorfields:
        if f == "STREAM_ID":
            arcpy.AddField_management(lc_point_fc, f, sid_type,
                                      sid_precision, sid_scale, sid_length,
                                      "", "NULLABLE", "NON_REQUIRED")
        else:
            arcpy.AddField_management(lc_point_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(lc_point_fc, ["SHAPE@X","SHAPE@Y"] +
                               cursorfields) as cursor:
        for row in pointList:
            cursor.insertRow(row)

def setup_LC_data_headers(transsample_count, trans_count,
                          canopy_data_type, stream_sample, heatsource8):
    """Generates a list of the landcover data file
    column header names and data types"""
    
    type = ["LC","ELE"]
     
    #Use LAI methods   
    if canopy_data_type == "LAI":
        type = type + ["LAI","k","OH"]
    
    #Use Canopy Cover methods    
    if canopy_data_type == "CanopyCover":  
        type = type + ["CAN","OH"]
    
    lcdataheaders =[]
    # a flag indicating the model should use the heat source 8 methods 
    # (same as 8 directions but no north)
    if heatsource8 is True:
        dir = ['T' + str(x) for x in range(1, 8)]
    else:        
        dir = ['T' + str(x) for x in range(1, trans_count + 1)]

    zone = range(1,int(transsample_count)+1)
    
    # Concatenate the type, dir, and zone and order in the correct way
    for t in type:
        for d in range(0,len(dir)):
            for z in range(0,len(zone)):
                if stream_sample is True and t !="ELE" and d==0 and z==0:
                    #lcdataheaders.append(t+"_EMERGENT") # add emergent
                    lcdataheaders.append(t+"_T0_S0") # add emergent
                    lcdataheaders.append(t+"_"+dir[d]+"_S"+str(zone[z]))
                else:
                    lcdataheaders.append(t+"_"+dir[d]+"_S"+str(zone[z]))
    
    return lcdataheaders, type

def coord_to_array(easting, northing, bbox_upper_left):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append((easting - bbox_upper_left[0]) / bbox_upper_left[2])  # col, x
    xy.append((northing - bbox_upper_left[1]) / bbox_upper_left[3] * -1)  # row, y 
    return xy

def create_lc_point_list(nodeDict, streamID, block_size, dir, zone, transsample_distance):
    """This builds a unique long form list of information for all the
    landcover samples. This list is used to create the output point
    feature class. The outer list holds all the nodes within a specified
    km extent (block_size). This is done for memory managment when the
    raster is converted to an array. Really large arrays will use up all
    the memory and cause a crash."""   
    
    lc_pointList = []
    nodeBlocks = []
    # Build a list of km every 10 km from 10 to 6700. 
    # if there is a stream longer than 6700 km we are not on Earth
    km_blocks = [x for x in range(block_size, 6700, block_size)]     
    i = 0
    
    nodes = nodeDict.keys()
    nodes.sort()

    for nodeID in nodes:
        origin_x = nodeDict[nodeID]["POINT_X"]
        origin_y = nodeDict[nodeID]["POINT_Y"]
        stream_km = nodeDict[nodeID]["STREAM_KM"]
        
        if stream_km < km_blocks[i]:
            # This is the emergent/stream sample
            nodeBlocks.append([origin_x, origin_y, origin_x, origin_y,
                               streamID, nodeID, 0, 0, 0])
                
            for d in range(0,len(dir)):
                for z in zone:
                    # Calculate the x and y coordinate of the 
                    # landcover sample location
                    lc_x = (z * transsample_distance * con_from_m *
                            sin(radians(dir[d]))) + origin_x
                    lc_y = (z * transsample_distance * con_from_m *
                            cos(radians(dir[d]))) + origin_y
    
                    # Add the all the data to the list
                    nodeBlocks.append([lc_x, lc_y, lc_x, lc_y, streamID,
                                       nodeID, dir[d], d+1, z])
        else: # New block
            lc_pointList.append([list(x) for x in zip(*nodeBlocks)]) 
            nodeBlocks = []
            # This is the emergent/stream sample
            nodeBlocks.append([origin_x, origin_y, origin_x, origin_y,
                               streamID, nodeID, 0, 0, 0])
                
            for d in range(0,len(dir)):
                for z in zone:
                    # Calculate the x and y coordinate of the 
                    # landcover sample location
                    lc_x = (z * transsample_distance * con_from_m *
                            sin(radians(dir[d]))) + origin_x
                    lc_y = (z * transsample_distance * con_from_m *
                            cos(radians(dir[d]))) + origin_y
    
                    # Add the all the data to the list
                    nodeBlocks.append([lc_x, lc_y, lc_x, lc_y, streamID,
                                       nodeID, dir[d], d+1, z])
            i = i + 1
    lc_pointList.append([list(x) for x in zip(*nodeBlocks)])            
    return lc_pointList

def sample_raster(x_coordList, y_coordList, raster, con):
    
    
    x_cellsize = arcpy.Describe(raster).meanCellWidth
    y_cellsize = arcpy.Describe(raster).meanCellHeight    
    
    # Get the coordinates of the upper left cell corner of the input raster
    top_y = float(arcpy.GetRasterProperties_management(raster, "TOP").getOutput(0))
    left_x = float(arcpy.GetRasterProperties_management(raster, "LEFT").getOutput(0))
    
    # calculate the buffer distance (in raster spatial units) to add to 
    # the raster bounding box when extracting to an array
    buffer = x_cellsize * 2
    
    # calculate lower left corner and nrows/cols for the bounding box    
    x_min = min(x_coordList) - buffer
    y_min = min(y_coordList) - buffer
    y_max = max(y_coordList) + buffer
    
    # Calculate the X and Y offset from the upper left node 
    # coordinates bounding box
    x_minoffset = ((left_x - x_min)%x_cellsize) - x_cellsize
    y_maxoffset = (top_y - y_max)%y_cellsize
     
    ncols = int(ceil((max(x_coordList) + buffer - x_min) / x_cellsize)) + 1
    nrows = int(ceil((y_max - y_min) / y_cellsize)) + 1
    bbox_lower_left = arcpy.Point(x_min, y_min) # must be in raster map units
    bbox_upper_left = [x_min + x_minoffset, y_max + y_maxoffset, x_cellsize, y_cellsize]
    nodata_to_value = -9999 / con_z_to_m
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        raster_array = arcpy.RasterToNumPyArray(raster, bbox_lower_left,
                                                ncols, nrows, nodata_to_value)
    except:
        import traceback
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
        sys.exit(pymsg)    
    
    # convert array values to meters if needed
    if con is not None:
        raster_array = raster_array * con
    
    lc_list = []
    
    #print("Extracting raster values")
    for i in range(0,len(x_coordList)):
        xy = coord_to_array(x_coordList[i], y_coordList[i], bbox_upper_left)
        lc_list.append(raster_array[xy[1], xy[0]])
    return lc_list

def update_nodes_fc(nodeDict, nodes_fc, addFields): 
    """Updates the input point feature class with data from the nodes dictionary"""
    print("Updating input point feature class")

    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)     

    # Check to see if the field exists and add it if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")   

    with arcpy.da.UpdateCursor(nodes_fc,["STREAM_ID","NODE_ID"] + addFields) as cursor:
        for row in cursor:
            for f in xrange(0,len(addFields)):
                streamID = row[0]
                nodeID =row[1]
                row[f+2] = nodeDict[streamID][nodeID][addFields[f]]
                cursor.updateRow(row)

def from_meters_con(inFeature):
    """Returns the conversion factor to get from meters to the
    spatial units of the input feature class"""
    try:
        con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that".format(inFeature)+
        "is not projected or not recognized. Use a projected"+
        "coordinate system preferably"+
        "in linear units of feet or meters.")
        sys.exit("Coordinate system is not projected or not recognized. "+
                 "Use a projected coordinate system, preferably in linear "+
                 "units of feet or meters.")   
    return con_from_m

def from_z_units_to_meters_con(zUnits):
    """Returns the converstion factor to get from
    the input z units to meters"""
        
    try:
        con_z_to_m = float(zunits)
    except:
        if zUnits == "Meters":
            con_z_to_m = 1.0 
        elif zUnits == "Feet":
            con_z_to_m = 0.3048
        else: con_z_to_m = None # The conversion factor will not be used
    
    return con_z_to_m

#enable garbage collection
gc.enable()

try:
    print("Step 5: Sample Landcover - Star Pattern, Point Method")
    
    #keeping track of time
    startTime= time.time()
    
    # Check if the output exists and delete or throw an error
    if arcpy.Exists(lc_point_fc):
        if overwrite_data is True:
            arcpy.Delete_management(lc_point_fc)
        else:
            arcpy.AddError("\nThis output already exists: \n" +
                           "{0}\n".format(lc_point_fc) + 
                           "Please rename your output or choose to"+
                           "overwrite your data.")
            sys.exit("This output already exists: \n" +
                           "{0}\n".format(lc_point_fc) + 
                           "Please rename your output or choose to" +
                           "overwrite your data.")    
    
    # Determine input spatial units and set unit conversion factors
    proj = arcpy.Describe(nodes_fc).SpatialReference
    proj_ele = arcpy.Describe(z_raster).spatialReference
    proj_lc = arcpy.Describe(lc_raster).spatialReference
    
    con_from_m = from_meters_con(nodes_fc)
    con_lc_to_m = from_z_units_to_meters_con(lc_units)
    con_z_to_m = from_z_units_to_meters_con(z_units)
    
    # Check to make sure the raster and input 
    # points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not ".format(nodes_fc,z_raster)+
                       "have the same projection."+
                       "Please reproject your data.")
        sys.exit("Input points and elevation raster do not have the "+
                 "same projection. Please reproject your data.")    
    
    if proj_lc.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not ".format(proj_lc,z_raster)+
                       "have the same projection."+
                       "Please reproject your data.")
        sys.exit("The landcover and elevation rasters do not have the "+
                 "same projection. Please reproject your data.")    
    
    if block_size is "#": block_size = 5
    
    # Setup the raster list
    typeraster = [lc_raster, z_raster]
    if canopy_data_type == "LAI":  #Use LAI methods
        if canopy_raster is not "#": typeraster.append(canopy_raster)
        else: typeraster.append(None)
        if k_raster is not "#": typeraster.append(k_raster)
        else: typeraster.append(None)
        if oh_raster is not "#": typeraster.append(oh_raster)
        else: typeraster.append(None)
        
    if canopy_data_type == "CanopyCover":  #Use Canopy Cover methods
        if canopy_raster is not "#": typeraster.append(canopy_raster)
        else: typeraster.append(None)
        if oh_raster is not "#": typeraster.append(oh_raster)
        else: typeraster.append(None)   
    
    stream_sample = True
    # flag indicating the model should use the heat source 8 methods 
    # (same as 8 directions but no north)
    if heatsource8 is True:
        dir = [45,90,135,180,225,270,315]
    else:        
        dir = [x * 360.0 / trans_count for x in range(1,trans_count+ 1)]

    zone = range(1,int(transsample_count+1))
    
    # TODO 
    # This is a future function that may replace the emergent methods.
    # If True there is a regular landcover sample at the stream node
    # for each azimuth direction vs a single emergent sample at the 
    # stream node.
    #if stream_sample == "TRUE":
        #zone = range(0,int(transsample_count))
    #else:
        #zone = range(1,int(transsample_count+1))
        
    addFields, type = setup_LC_data_headers(transsample_count, trans_count,
                                            canopy_data_type, stream_sample,
                                            heatsource8)
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
       
    lc_pointList = []
    lc_pointList2 = []
    n = 1 
    for streamID in nodeDict:
        print("Processing stream %s of %s" % (n, len(nodeDict)))
        lc_pointList = create_lc_point_list(nodeDict[streamID], streamID,
                                            block_size, dir, zone,
                                            transsample_distance)
        for NodeBlock in lc_pointList:
            for raster in typeraster:
                if raster is None:
                    lc_list = [-9999 for i in NodeBlock[0]]
                else: 
                    if raster == z_raster:
                        con = con_z_to_m
                    elif raster == lc_raster:
                        con = con_lc_to_m
                    else:
                        con = None
                    lc_list = sample_raster(NodeBlock[0],NodeBlock[1], raster, con)
                NodeBlock.append(lc_list)
            # Transpose the list and append to ouput
            lc_pointList2 = lc_pointList2 + [list(x) for x in zip(*NodeBlock)]
        n = n + 1       
    
    # Update the nodeDict
    for row in lc_pointList2:
        for t in range(0,len(type)):
            lc_key = type[t]+'_T'+str(row[7])+'_S'+str(row[8])
            nodeDict[row[4]][row[5]][lc_key] = row[9 + t]

    endTime = time.time()
    arcpy.ResetProgressor()		
    gc.collect()

    # Create the landcover headers to be added to the TTools point 
    # feature class
    #addFields = addFields + ["NUM_DIR","NUM_ZONES","SAMPLE_DIS"]    
    
    # Write the landcover data to the TTools point feature class
    update_nodes_fc(nodeDict, nodes_fc, addFields)    
    
    # Build the output point feature class using the data 
    # from the LCPointList
    create_lc_point_fc(lc_pointList2, type, lc_point_fc, nodes_fc, proj)

    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspersample = timedelta(seconds=(endTime - startTime) /
                            len(lc_pointList2)).microseconds
    print("Process Complete in %s minutes. %s microseconds per sample" % (elapsedmin, mspersample))    
    #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per sample" % (elapsedmin, mspersample))


# For arctool errors
except arcpy.ExecuteError:
    msgs = arcpy.GetMessages(2)
    #arcpy.AddError(msgs)
    print(msgs)

# For other errors
except:
    import traceback, sys
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]

    pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" +str(sys.exc_info()[1])
    msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

    #arcpy.AddError(pymsg)
    #arcpy.AddError(msgs)

    print(pymsg)
    print(msgs)