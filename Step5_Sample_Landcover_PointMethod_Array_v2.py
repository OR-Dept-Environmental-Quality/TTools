########################################################################
# TTools
# Step 5: Sample Landcover - Star Pattern, Point Method v 0.995
# Ryan Michie

# Sample_Landcover_PointMethod will take an input point 
# feature (from Step 1) and sample input landcover rasters in a user 
# specificed number of cardinal directions with point samples spaced 
# at a user defined distance moving away from the stream.

# General scripts steps include:
# 1. open nodes fc. iterate and read info from node fc into a dict

# 2. create a list with the x/y and related info for each lc sample
#    calculate the extent bounding box for the entire dataset

# 3. create a list holding the bounding box coords for each block itteration

# 4. loop through each block
    #- sample the raster for all samples in the block
    #- update the node dict
    #- update the node fc
    #- update the lc sample fc
    #- continue to the next block

# INPUTS
# 0: TTools point feature class (nodes_fc)
# 1: Number of transects per node (trans_count)
# 2: Number of samples per transect (transsample_count)
# 3: The distance between transect samples (transsample_distance)
# 4: True/False flag if using heatsource 8 methods (heatsource8)
# 5: Land cover code or height raster (lc_raster)
# 6: input (optional) landcover height z units 
#     1. "Feet", 2. "Meters" 3. "None" Float (lc_units)
# 7: OPTIONAL - landcover data type:
#     1."CanopyCover", or 2."LAI" (canopy_data_type)
# 8: OPTIONAL - canopy cover or LAI raster (canopy_raster)
# 9: OPTIONAL - k coeffcient raster (k_raster)
# 10: OPTIONAL - overhang raster (oh_raster)
# 11: Elevation raster (z_raster)
# 12: Elvation raster z units (z_units) 
#      1. "Feet", 2. "Meters" 3. "Other"
# 13: Path/name of output sample point file (lc_point_fc)
# 14: OPTIONAL - km distance to process within each array (block_size)
# 15: True/False flag if existing data can be over written (overwrite_data)

# OUTPUTS
# 0. point feature class (edit nodes_fc) - added fields with 
#     Landcover and elevation data for each azimuth direction at each node
# 1. point feature class (new) - point at each x/y sample and 
#     the sample raster values

# Future Updates
# -Change the node dict so the node is the primary key
# -Build the block list based on the nodes and then build point list itterativly
# instead of building them into one huge list. The huge list results
# in a memory error for large areas
# -Include stream sample in transect count (True/False)
# -Eliminate arcpy and use gdal for reading/writing feature class data

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
from math import radians, sin, cos, ceil, sqrt
from collections import defaultdict
import numpy
import arcpy
from arcpy import env

env.overwriteOutput = True

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
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
lc_point_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\LC_samplepoint"
block_size = 5 # OPTIONAL defualt to 5
overwrite_data = True
# End Fill in Data
# ----------------------------------------------------------------------

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
    
    cursorfields = ["POINT_X","POINT_Y"] +["STREAM_ID","NODE_ID",
                                           "TRANS_AZIMUTH","TRANSNUM",
                                           "SAMPLENUM"] +LCFields    
    
    # Check if the output exists and create if not
    if not arcpy.Exists(lc_point_fc):
        arcpy.CreateFeatureclass_management(os.path.dirname(lc_point_fc),
                                            os.path.basename(lc_point_fc),
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

def create_lc_point_list(origin_x, origin_y, streamID,
                          nodeID, dir, zone, transsample_distance):
    """This builds a unique long form list of information for all the
    landcover samples. This list is used to create the output point
    feature class. The outer list holds all the nodes within a specified
    km extent (block_size). This is done for memory managment when the
    raster is converted to an array. Really large arrays will use up all
    the memory and cause a crash."""
    
    point_list = []
    
    # This is the emergent/stream sample
    point_list.append([origin_x, origin_y, origin_x, origin_y,
                       streamID, nodeID, 0, 0, 0])

    for d in range(0,len(dir)):
        for z in zone:
            # Calculate the x and y coordinate of the 
            # landcover sample location
            pt_x = (z * transsample_distance * con_from_m *
                    sin(radians(dir[d]))) + origin_x
            pt_y = (z * transsample_distance * con_from_m *
                    cos(radians(dir[d]))) + origin_y

            # Add to the list          
            point_list.append([pt_x, pt_y, pt_x, pt_y,
                               streamID, nodeID, dir[d], d+1, z])
     
    return point_list

def create_block_list(lc_pointList, block_size):
        
    x_coord_list = [i[0] for i in lc_pointList]
    y_coord_list = [i[1] for i in lc_pointList]
    
    # calculate the buffer distance (in raster spatial units) to add to 
    # the raster bounding box when extracting to an array
    buffer = 0     
    
    # calculate bounding box extent for samples
    x_min = min(x_coord_list)
    x_max = max(x_coord_list)
    y_min = min(y_coord_list) 
    y_max = max(y_coord_list)
    
    x_width = int(x_max - x_min + 1)
    y_width = int(y_max - y_min + 1)
    
    block_extents = []
    block_points = []
      
    # Build data blocks
    for x in range(0, x_width, block_size):
        for y in range(0, y_width, block_size):

            # Lower left coordinate of block (in map units)
            block_x_min = min([x_min + x, x_max])
            block_y_min = min([y_min + y, y_max])
            # Upper right coordinate of block (in map units)
            block_x_max = min([block_x_min + block_size, x_max])
            block_y_max = min([block_y_min + block_size, y_max])
            
            samples_in_block = []
            for sample in lc_pointList:
                if block_x_min <= sample[0] <= block_x_max and block_y_min <= sample[1] <= block_y_max:
                    samples_in_block.append(sample)
            
            if samples_in_block:      
                # order is left, bottom, right, top
                block_extents.append([block_x_min - buffer, block_y_min - buffer,
                                      block_x_max + buffer, block_y_max - + buffer])
                # 0 left,      1 bottom,    2 right,     3 top
                # block_x_min, block_y_min, block_x_max, block_y_max, ncols, nrows             
                block_points.append(samples_in_block)
    
    return block_extents, block_points
    
def sample_raster(block, block_point_list, raster, con):
        
    nodata_to_value = -9999 / con
    
    x_cellsize = arcpy.Describe(raster).meanCellWidth
    y_cellsize = arcpy.Describe(raster).meanCellHeight    

    # Get the coordinates of the upper left cell corner of the input raster
    raster_y_max = float(arcpy.GetRasterProperties_management(raster, "TOP").getOutput(0))
    raster_x_min = float(arcpy.GetRasterProperties_management(raster, "LEFT").getOutput(0))
    
    # block list order:
    # 0 left,      1 bottom,    2 right,     3 top
    # block_x_min, block_y_min, block_x_max, block_y_max, ncols, nrows
    ncols = max([int(ceil((block[2] - block[0]) / x_cellsize)), 1])
    nrows = max([int(ceil((block[3] - block[1]) / y_cellsize)), 1])
    
    # Calculate the X and Y offset from the upper left node 
    # coordinates bounding box
    x_minoffset = ((raster_x_min - block[0])%x_cellsize) - x_cellsize
    y_maxoffset = ((raster_y_max - block[3])%y_cellsize) - y_cellsize
    
    block_lower_left = arcpy.Point(block[0], block[1])
    bbox_upper_left = [block[0] + x_minoffset, block[3] +
               y_maxoffset, x_cellsize, y_cellsize]

    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        raster_array = arcpy.RasterToNumPyArray(raster, block_lower_left,
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
    
    point_list = []
    for point in block_point_list:
        xy = coord_to_array(point[0], point[1], bbox_upper_left)
        point.append(raster_array[xy[1], xy[0]])
        point_list.append(point)
    return point_list            

def update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_query):
    """Updates the input point feature class with data from the nodes dictionary"""
    print("Updating input point feature class")
    
    # Build a query to retreive just the nodes that need updating    
    whereclause = """%s IN %s""" % (arcpy.AddFieldDelimiters(nodes_fc, "NODE_ID"), tuple(nodes_to_query))

    with arcpy.da.UpdateCursor(nodes_fc,["STREAM_ID","NODE_ID"] + addFields, whereclause) as cursor:  
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
            arcpy.AddMessage("\nThis output already exists: \n" +
                           "{0}\n".format(lc_point_fc) + 
                           "overwrite data = False. New data will be " +
                           "appended to the existing feature class.")
            print("This output already exists: \n" +
                           "{0}\n".format(lc_point_fc) + 
                           "overwrite data = False. New data will be " +
                           "appended to the existing feature class.")    
    
    # Determine input spatial units and set unit conversion factors
    proj = arcpy.Describe(nodes_fc).SpatialReference
    proj_ele = arcpy.Describe(z_raster).spatialReference
    proj_lc = arcpy.Describe(lc_raster).spatialReference
    
    con_from_m = from_meters_con(nodes_fc)
    con_lc_to_m = from_z_units_to_meters_con(lc_units)
    con_z_to_m = from_z_units_to_meters_con(z_units)
    
    # convert block size from km to meters into units of the node fc
    # in the future block size should be estimated based on availiable memory
    # memorysize = datatypeinbytes*nobands*block_size^2
    # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
    if block_size in ["#", ""]:
        block_size = 5
    else:
        block_size = int(con_from_m * block_size / 1000)
    
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
    
    # Add the new fields to nodes_fc
    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)     
    
    # Check to see if the field exists and add it if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")
            
    lc_pointList = []
    n = 1
    for streamID in nodeDict:
        
        nodes = nodeDict[streamID].keys()
        nodes.sort()
        
        for nodeID in nodes:
            origin_x = nodeDict[streamID][nodeID]["POINT_X"]
            origin_y = nodeDict[streamID][nodeID]["POINT_Y"]
            
            lc_pointList += create_lc_point_list(origin_x, origin_y,
                                                  streamID, nodeID,
                                                  dir, zone,
                                                  transsample_distance)
    block_extents, block_points = create_block_list(lc_pointList, block_size)
    
    p = 0
    for block in block_extents:
        print("Processing block %s of %s" % (p + 1, len(block_extents)))
        #point_list = block_points[p]
        for raster in typeraster:
            if raster is None:
                for i in range(0, len(point_list)):
                    point_list[i].append(-9999)
                    
            else: 
                if raster == z_raster:
                    con = con_z_to_m
                elif raster == lc_raster:
                    con = con_lc_to_m
                else:
                    con = None
                point_list = sample_raster(block, block_points[p], raster, con)    
            
        # Update the node fc
        nodes_to_query = []
        for row in point_list:
            nodes_to_query.append(row[5])
            for t in range(0,len(type)):                
                lc_key = type[t]+'_T'+str(row[7])+'_S'+str(row[8])
                nodeDict[row[4]][row[5]][lc_key] = row[9 + t]
        
        # Write the landcover data to the TTools point feature class 
        update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_query)
        
        # Build the output point feature class using the data         
        create_lc_point_fc(point_list, type, lc_point_fc, nodes_fc, proj)
        p = p + 1
        del point_list
    
    endTime = time.time()
    gc.collect()
    
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspersample = timedelta(seconds=(endTime - startTime) /
                            len(lc_pointList)).microseconds
    print("Process Complete in %s minutes. %s microseconds per sample" % (elapsedmin, mspersample))    
    #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per sample" % (elapsedmin, mspersample))


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