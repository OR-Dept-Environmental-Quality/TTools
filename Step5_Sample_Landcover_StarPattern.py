#!/usr/bin/python

"""
TTools Step 5: Sample Landcover - Star Pattern

The star pattern sampling method will take an input point feature (from Step 1) and sample input landcover rasters
along transects in any number of azimuth directions oriented outward from the stream node. The numer of transects,
sample points along each transect, and the distance between samples is user defined. The measured distance along each
transect starts at the stream node.

The star pattern sampling method can be used to develop landcover inputs for Heat Source 7 - 9.

REQUIREMENTS
TTools steps 1 - 3 must be run before Step 5.
ESRI ArcGIS Pro w/ Spatial Analyst extension
Python 3.7+

INPUT VARIABLES
0: nodes_fc:
Path to the TTools point feature class.

1: trans_count:
Number of transects per node. The degrees of separation between each transect is equal to 360 / trans_count. If
trans_count = 8, the heading of the first transect is 45 degrees (Northeast), the second transect is 90 degrees
(West), and the eighth is 360 degrees (North). Transects are oriented clockwise from north. trans_count is ignored if
heatsource8 = True.

2: transsample_count:
Number of samples per transect. This number DOES NOT include the sample at the stream node.

3: transsample_distance:
The distance between transect samples measured in meters.

4: zone_sample:
Boolean (True/False) flag to indicate if the sample should represent a zone and be centered relative to
the transsample_distance as measured from the stream node. If this is True the distance from the stream
node to the sample location along each transect is equal to the transsample_distance * (sample number - 0.5).
This should be True if using heat source 7, heat source 8.0.1 - 8.0.5, or heat source 8.0.7 -  8.0.8 when the model is
configured to use the zone method. This should be False if using heat source 9 or when using heat source 8.0.7 -  8.0.8
and the model is configured to use the point method.

5. heatsource8:
Boolean (True/False) flag to indicate if the star pattern with 7 transects should be used. Heat source
version 7 and 8 use 7 transects around each stream node in the following directions: Northeast, East, Southeast, South,
Southwest, West, Northwest.

6: sampleID_for_code
Boolean (True/False) flag to indicate if the sample id should be used as the landcover value. If True the lc_raster
will be ignored and only the z_raster will be sampled. Generally this should be set to False.

7: lc_raster:
Path and name of the land cover code, height, or elevation raster.

8: lc_units:
z units of the lc_raster (aka units of height or elevation). Use "Feet", "Meters", or "None" if the lc_raster values are
codes and do not represent elevation or height units.

9: canopy_data:
The input canopy data type being sampled in canopy_raster. Can be one of "CanopyCover", "LAI", or "None" if not sampling
a canopy raster. If canopy_data = None, the variables canopy_raster, k_raster, and oh_raster are ignored.

10: canopy_raster
Path and name of the canopy or LAI raster. Ignored if canopy_data = "None".

11: k_raster
Path and name of the k coefficient raster. There must be a k_raster if canopy_data = "LAI". Ignored if
canopy_data = "None".

12: oh_raster
Path and name of the vegetation overhang raster. Ignored if canopy_data = "None".

13: z_raster:
Path and name of the ground elevation raster.

14: z_units:
z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is not in feet or meters the elevation
values must be converted.

15: lc_point_fc:
Path and name of output sample point feature file.

16: block_size:
The x and y size in kilometers for each raster block pulled into an array. Start with 5 if you aren't sure and reduce
if there is an error. To increase processing speed rasters are subdivided iteratively into smaller blocks and pulled
into arrays for processing. Very large block sizes may use a lot of memory and result in an error.

17: overwrite_data:
True/False flag if existing data in nodes_fc and lc_point_fc can be overwritten.

OUTPUTS
0. nodes_fc:
New fields are added into nodes_fc with the Landcover and elevation values for each transect sample

1. lc_point_fc:
New point feature class created with a point at each x/y sample and the sample raster values

"""
# Import system modules
import sys
import os
import gc
import time
import traceback
from datetime import timedelta
from math import radians, sin, cos, ceil
from collections import defaultdict, OrderedDict
import arcpy
from arcpy import env

# ----------------------------------------------------------------------
# Start input variables
nodes_fc = r"C:\workspace\ttools_tests\TTools_py39\jc_test_py39.gdb\jc_stream_nodes_py39"
trans_count = 8
transsample_count = 5
transsample_distance = 8
zone_sample = False
heatsource8 = False
sampleID_for_code = False
lc_raster = r"C:\workspace\ttools_tests\JohnsonCreek.gdb\jcw_vght_m_mosaic"
lc_units = "Meters"
canopy_data = None
canopy_raster = None
k_raster = None
oh_raster = None
z_raster = r"C:\workspace\ttools_tests\JohnsonCreek.gdb\jcw_be_m_mosaic"
z_units = "Meters"
lc_point_fc = r"C:\workspace\ttools_tests\TTools_py39\jc_test_py39.gdb\jc_LC_samples"
block_size = 5
overwrite_data = True
correctTTools = False

# End input variables
# ----------------------------------------------------------------------

# General scripts steps include:
# 1. open nodes fc. iterate and read info from node fc into a dict

# 2. create a list with the x/y and related info for each lc sample
#    calculate the extent bounding box for the entire dataset

# 3. create a list holding the bounding box coords for each block iteration

# 4. loop through each block
    #- sample the raster for all samples in the block
    #- update the node dict
    #- update the node fc
    #- update the lc sample fc
    #- continue to the next block

# Future Updates
# -Change the node dict so the node is the primary key
# -Build the block list based on the nodes and then build point list iteratively
# instead of building them into one huge list. The huge list results
# in a memory error for large areas
# -Include stream sample in transect count (True/False)
# -Eliminate arcpy and use gdal for reading/writing feature class data

env.overwriteOutput = True

def nested_dict(): 
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the STREAM_ID,
    NODE_ID, and X/Y coordinates as a nested dictionary"""
    nodeDict = nested_dict()
    incursorFields = ["NODE_ID", "STREAM_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]

    # Get a list of existing fields
    existingFields = [f.name for f in arcpy.ListFields(nodes_fc)] 

    # Check to see if the last field exists if yes add it. 
    # Grabs last field because often the first field, emergent, is zero
    if overwrite_data is False and (addFields[len(addFields)-1] in existingFields) is True:
        incursorFields.append(addFields[len(addFields)-1])
    else:
        overwrite_data = True

    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc, incursorFields,"",proj) as Inrows:
        if overwrite_data:
            for row in Inrows:
                nodeDict[row[0]]["STREAM_ID"] = row[1]
                nodeDict[row[0]]["STREAM_KM"] = row[2] 
                nodeDict[row[0]]["POINT_X"] = row[3]
                nodeDict[row[0]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[5] is None or row[5] == 0 or row[5] < -9998:
                    nodeDict[row[0]]["STREAM_ID"] = row[1]
                    nodeDict[row[0]]["STREAM_KM"] = row[2] 
                    nodeDict[row[0]]["POINT_X"] = row[3]
                    nodeDict[row[0]]["POINT_Y"] = row[4]
    
    if len(nodeDict) == 0:
        sys.exit("The fields checked in the input point feature class " +
                 "have existing data. There is nothing to process. Exiting")
              
    return nodeDict

def update_lc_point_fc(lc_point_list, type, lc_point_fc, nodes_fc,
                       nodes_in_block, overwrite_data, proj):
    """Creates/updates the output landcover sample point feature
    class using the data from the landcover point list"""

    cursorfields = ["POINT_X","POINT_Y"] + ["STREAM_ID","NODE_ID", "SAMPLE_ID",
                                            "TRANS_DIR", "TRANSECT",
                                            "SAMPLE", "KEY"] + type
    
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
        
        typeDict = {"POINT_X": "DOUBLE", 
                    "POINT_Y": "DOUBLE", 
                    "NODE_ID": "LONG",
                    "SAMPLE_ID": "LONG",
                    "TRANS_DIR": "DOUBLE",
                    "TRANSECT": "SHORT",
                    "SAMPLE": "SHORT",
                    "KEY": "TEXT",}
                
        for t in type:
            typeDict[t] = "DOUBLE"
    
        # Add attribute fields
        for f in cursorfields:
            if f == "STREAM_ID":
                arcpy.AddField_management(lc_point_fc, f, sid_type,
                                          sid_precision, sid_scale, sid_length,
                                          "", "NULLABLE", "NON_REQUIRED")
            else:
                arcpy.AddField_management(lc_point_fc, f, typeDict[f], "",
                                          "", "", "", "NULLABLE", "NON_REQUIRED")                 
    
    if not overwrite_data:
        # Build a query to retrieve existing rows from the nodes
        # that need updating
        whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_in_block))        
        
        # delete those rows
        with arcpy.da.UpdateCursor(lc_point_fc,["NODE_ID"], whereclause) as cursor:  
            for row in cursor:
                cursor.deleteRow()     

    with arcpy.da.InsertCursor(lc_point_fc, ["SHAPE@X","SHAPE@Y"] +
                               cursorfields) as cursor:
        for row in lc_point_list:
            cursor.insertRow(row)

def setup_lcdata_headers(transsample_count, trans_count, canopy_data):
    """Generates a list of the landcover data file
    column header names and data types"""
    
    type = ["ELE"]
     
    #Use LAI methods   
    if canopy_data == "LAI":
        type = type + ["LAI","k","OH"]
    
    #Use Canopy Cover methods    
    if canopy_data == "CanopyCover":
        type = type + ["CAN","OH"]
    
    lcheaders = []
    otherheaders = []
     
    dirs = ["T{0}".format(x) for x in range(1, trans_count + 1)]

    zones = range(1,int(transsample_count)+1)
    
    # Concatenate the type, dir, and zone and order in the correct way
    
    for d, dir in enumerate(dirs):
        for z, zone in enumerate(zones):
            if d==0 and z==0:
                lcheaders.append("LC_T0_S0") # add emergent
                lcheaders.append("LC_{0}_S{1}".format(dir, zone))
            else:
                lcheaders.append("LC_{0}_S{1}".format(dir, zone))
    
    for t in type:
        for d, dir in enumerate(dirs):
            for z, zone in enumerate(zones):
                if t !="ELE" and d==0 and z==0:
                    otherheaders.append(t+"_T0_S0") # add emergent
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))
                else:
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))    

    return lcheaders, otherheaders

def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
    """converts x/y coordinates to col and row of the array"""
    col_x = int((easting - block_x_min) / x_cellsize)  # col, x
    row_y = int((block_y_max - northing) / y_cellsize)  # row, y
    return [col_x, row_y]

def create_lc_point_list(nodeDict, nodes_in_block, dirs, zones, transsample_distance, zone_sample):
    """This builds a unique long form list of information for all the
    landcover samples in the block. This list is used to
    create/update the output feature class."""
    
    lc_point_list = []
    numDirs = len(dirs)
    numZones = len(zones)
    zonesPerNode = (numDirs * numZones) + 1

    if zone_sample:
        adjust = 0.5
    else:
        adjust = 0.0

    for nodeID in nodes_in_block:
        origin_x = nodeDict[nodeID]["POINT_X"]
        origin_y = nodeDict[nodeID]["POINT_Y"]
        streamID = nodeDict[nodeID]["STREAM_ID"]
        sampleID = nodeID * zonesPerNode
        
        # This is the emergent/stream sample
        lc_point_list.append([origin_x, origin_y, origin_x, origin_y,
                              streamID, nodeID, sampleID,
                              0, 0, 0, "T0_S0"])
    
        for d, dir in enumerate(dirs):
            for zone in zones:
                # Calculate the x and y coordinate of the 
                # landcover sample location
                pt_x = ((zone - adjust) * transsample_distance * con_from_m *
                        sin(radians(dir))) + origin_x
                pt_y = ((zone - adjust) * transsample_distance * con_from_m *
                        cos(radians(dir))) + origin_y
                
                key = 'T{0}_S{1}'.format(d+1, zone)
                tran_num = '{:{}{}}'.format(d+1, 0, 3)
                samp_num = '{:{}{}}'.format(zone, 0, 2)
                sampleID = (nodeID * zonesPerNode) + (d * numZones) + zone
        
                # Add to the list          
                lc_point_list.append([pt_x, pt_y, pt_x, pt_y,
                                      streamID, nodeID, sampleID,
                                      dir, d+1, zone, key])    
     
    return lc_point_list

def create_block_list(nodes, block_size):
    """Returns two lists, one containing the coordinate extent
    for each block that will be iteratively extracted to an array
    and the other containing node IDs within each block extent."""
    
    print("Calculating block extents")    
    x_coord_list = [nodeDict[nodeID]["POINT_X"] for nodeID in nodes]
    y_coord_list = [nodeDict[nodeID]["POINT_Y"] for nodeID in nodes]
    
    # calculate the buffer distance (in raster spatial units) to add to 
    # the base bounding box when extracting to an array. The buffer is 
    # equal to the sample distance + 1 to make sure the block includes 
    # all the landcover samples for each node.
    buffer = int((transsample_count + 1) * transsample_distance * con_from_m)
    
    # calculate bounding box extent for samples
    x_min = min(x_coord_list)
    x_max = max(x_coord_list)
    y_min = min(y_coord_list) 
    y_max = max(y_coord_list)
    
    x_width = int(x_max - x_min + 1)
    y_width = int(y_max - y_min + 1)
    
    block_extents = []
    block_nodes = []
      
    # Build data blocks
    for x in range(0, x_width, block_size):
        for y in range(0, y_width, block_size):

            # Lower left coordinate of block (in map units)
            block0_x_min = min([x_min + x, x_max])
            block0_y_min = min([y_min + y, y_max])
            # Upper right coordinate of block (in map units)
            block0_x_max = min([block0_x_min + block_size, x_max])
            block0_y_max = min([block0_y_min + block_size, y_max])
            
            block_x_min = block0_x_max
            block_x_max = block0_x_min
            block_y_min = block0_y_max
            block_y_max = block0_y_min
            
            nodes_in_block = []
            for nodeID in nodes:
                node_x = nodeDict[nodeID]["POINT_X"]
                node_y = nodeDict[nodeID]["POINT_Y"]
                if (block0_x_min <= node_x <= block0_x_max and
                    block0_y_min <= node_y <= block0_y_max):
                    
                    nodes_in_block.append(nodeID)
                    
                    # Minimize the size of the block0 by the true 
                    # extent of the nodes in the block
                    if block_x_min > node_x: block_x_min = node_x
                    if block_x_max < node_x: block_x_max = node_x
                    if block_y_min > node_y: block_y_min = node_y
                    if block_y_max < node_y: block_y_max = node_y
            
            if nodes_in_block:
                # add the block extent for processing
                # order 0 left,      1 bottom,    2 right,     3 top
                block_extents.append((block_x_min - buffer, block_y_min - buffer,
                                      block_x_max + buffer, block_y_max + buffer))           
                block_nodes.append(nodes_in_block)
    
    return block_extents, block_nodes
    
def sample_raster(block, lc_point_list, raster, con):
    
    if con is not None:
        nodata_to_value = -9999 / con
    else:
        nodata_to_value = -9999
        
    # localize the block extent values
    block_x_min = block[0]
    block_y_min = block[1]
    block_x_max = block[2]
    block_y_max = block[3]
    
    x_cellsize = float(arcpy.GetRasterProperties_management(raster, "CELLSIZEX").getOutput(0))
    y_cellsize = float(arcpy.GetRasterProperties_management(raster, "CELLSIZEY").getOutput(0)) 

    # Get the coordinates extent of the input raster
    raster_x_min = float(arcpy.GetRasterProperties_management(raster, "LEFT").getOutput(0))
    raster_y_min = float(arcpy.GetRasterProperties_management(raster, "BOTTOM").getOutput(0))
    raster_x_max = float(arcpy.GetRasterProperties_management(raster, "RIGHT").getOutput(0))
    raster_y_max = float(arcpy.GetRasterProperties_management(raster, "TOP").getOutput(0))
      
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
        raster_array = arcpy.RasterToNumPyArray(raster, arcpy.Point(block_x_min_center, block_y_min_center),
                                                ncols, nrows, nodata_to_value)
    except:
        tbinfo = traceback.format_exc()
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
        sys.exit(pymsg)    
    
    # convert array values to meters if needed
    if con is not None:
        raster_array = raster_array * con
    
    lc_point_list_new = []
    if raster_array.max() > -9999:
        # There is at least one pixel of data
        for point in lc_point_list:
            xy = coord_to_array(point[0], point[1], block_x_min, block_y_max, x_cellsize, y_cellsize)
            point.append(raster_array[xy[1], xy[0]])
            lc_point_list_new.append(point)
    else:
        # No data, add -9999
        for point in lc_point_list:
            point.append(-9999)
            lc_point_list_new.append(point)
    return lc_point_list_new            

def update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_query):
    """Updates the input point feature class with data from the
    nodes dictionary"""

    # Build a query to retrieve just the nodes that needs updating
    whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))

    with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:  
        for row in cursor:
            for f, field in enumerate(addFields):
                nodeID =row[0]
                row[f+1] = nodeDict[nodeID][field]
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
    """Returns the conversion factor to get from
    the input z units to meters"""
        
    try:
        con_z_to_m = float(zUnits)
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
    
    # Check if the node fc exists
    if not arcpy.Exists(nodes_fc):
        arcpy.AddError("\nThis output does not exist: \n" +
                       "{0}\n".format(nodes_fc))
        sys.exit("This output does not exist: \n" +
                 "{0}\n".format(nodes_fc))
    
    # Check if the lc point fc exists and delete if needed
    if arcpy.Exists(lc_point_fc) and overwrite_data:
        arcpy.Delete_management(lc_point_fc)
    
    
    # Determine input spatial units and set unit conversion factors
    proj = arcpy.Describe(nodes_fc).SpatialReference
    proj_ele = arcpy.Describe(z_raster).spatialReference
    proj_lc = arcpy.Describe(lc_raster).spatialReference
    
    con_from_m = from_meters_con(nodes_fc)
    con_lc_to_m = from_z_units_to_meters_con(lc_units)
    con_z_to_m = from_z_units_to_meters_con(z_units)
    
    # convert block size from km to meters to units of the node fc
    # in the future block size should be estimated based on available memory
    # memorysize = datatypeinbytes*nobands*block_size^2
    # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
    if block_size in ["#", "", None]:
        block_size = int(con_from_m * 5000)
    else:
        block_size = int(con_from_m * block_size * 1000)
    
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

    # Setup the raster dictionary. It is ordered because
    # key list needs to correspond to the order of the attribute fields
    rasterDict = OrderedDict({"LC": lc_raster,
                              "ELE": z_raster})
    if sampleID_for_code:
        rasterDict["LC"] = None

    if canopy_data == "LAI":
        # Use LAI methods
        rasterDict["LAI"] = canopy_raster
        rasterDict["k"] = k_raster
        rasterDict["OH"] = oh_raster

    if canopy_data == "CanopyCover":
        # Use Canopy Cover methods
        rasterDict["CAN"] = canopy_raster
        rasterDict["OH"] = oh_raster

    # flag indicating the model should use the heat source 8 methods 
    # (same as 8 directions but no north)
    if heatsource8:
        dirs = [45,90,135,180,225,270,315]
        trans_count = 7
    else:        
        dirs = [x * 360.0 / trans_count for x in range(1,trans_count+ 1)]

    zones = range(1,int(transsample_count+1))

    lcheaders, otherheaders = setup_lcdata_headers(transsample_count, trans_count, canopy_data)
    
    addFields = lcheaders + otherheaders
    
    # Get a list of existing fields
    existingFields = [f.name for f in arcpy.ListFields(nodes_fc)] 
        
    # Check to see if the field exists and add it if not
    for f in lcheaders:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "TEXT", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")    

    # Check to see if the field exists and add it if not
    for f in otherheaders:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")    
    
    # read the node data into the node dictionary
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
    
    # Get a list of the nodes, sort them
    nodes = list(nodeDict.keys())
    nodes.sort()
   
    # Build the block list
    block_extents, block_nodes = create_block_list(nodes, block_size)
    
    # Iterate through each block, calculate sample coordinates,
    # convert raster to array, sample the raster
    total_samples = 0
    for p, block in enumerate(block_extents):
        nodes_in_block = block_nodes[p]
        nodes_in_block.sort()
        print("Processing block {0} of {1}".format(p + 1, len(block_extents)))
        
        # build the landcover sample list
        lc_point_list = create_lc_point_list(nodeDict, nodes_in_block,
                                            dirs, zones, transsample_distance, zone_sample)
        
        for t, (type, raster) in enumerate(rasterDict.items()):
            if raster is None:
                for i in range(0, len(lc_point_list)):
                    lc_point_list[i].append(-9999)
            else: 
                if raster == z_raster:
                    con = con_z_to_m
                elif raster == lc_raster:
                    con = con_lc_to_m
                else:
                    con = None  
            
                lc_point_list = sample_raster(block, lc_point_list, raster, con)
        
            # Update the node dict
            if (sampleID_for_code and type == "LC"):
                for row in lc_point_list:
                    key = "{0}_{1}".format(type, row[10])
                    nodeDict[row[5]][key] = row[6]
            
            else:
                for row in lc_point_list:
                    key = "{0}_{1}".format(type, row[10])
                    nodeDict[row[5]][key] = row[11 + t]
        
        # Write the landcover data to the TTools point feature class 
        update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_in_block)
        
        # Build the output point feature class using the data         
        update_lc_point_fc(lc_point_list=lc_point_list, type=list(rasterDict.keys()), lc_point_fc=lc_point_fc,
                           nodes_fc=nodes_fc, nodes_in_block=nodes_in_block, overwrite_data=overwrite_data, proj=proj)
    
        total_samples = total_samples + len(lc_point_list)
        del lc_point_list
        gc.collect()
   
    endTime = time.time()
    
    elapsedmin = ceil(((endTime - startTime) / 60)* 10)/10
    mspersample = timedelta(seconds=(endTime - startTime) /
                            total_samples).microseconds
    print("Process Complete in {0} minutes. {1} microseconds per sample".format(elapsedmin, mspersample))
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