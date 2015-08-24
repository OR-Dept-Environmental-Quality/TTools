########################################################################
# TTools
# Step 3: Sample Stream Elevations/ Gradient - v 0.961
# Ryan Michie

# Sample_ElevationsGradient will take an input point feature 
# (from Step 1) and sample the input raster elevation to find the 
# lowest elevation in a user defined search radius and calculate the 
# gradient for each node in the downstream direction.

# INPUTS
# 0: Input TTools point feature class (nodes_fc)
# 1: input the number of cells to search out around the node 
#     for the lowest elevation (searchCells)
# 2: input flag for smoothing if gradient is zero 
#     or negative (smooth_flag) 1. True, 2. False
# 3: input elevation raster (z_raster)
# 4: input elevation raster z units (z_units) 1. "Feet", 2. "Meters"
# 5: OPTIONAL - km distance to process within each array (block_size)
# 6: input flag if existing data can 
#     be over written (overwrite_data) 1. True, 2. False

# OUTPUTS
# 0. point feature class (edit nodes_fc) - Added fields 
#     with ELEVATION and GRADIENT for node

# Future Updates
# eliminate arcpy and use gdal for reading/writing feature class data

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

########################################################################

# Import system modules
from __future__ import division, print_function
import sys
import gc
import time
import traceback
from datetime import timedelta
import arcpy
import itertools
from arcpy import env
from math import ceil
from collections import defaultdict

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
searchCells = 0
smooth_flag = True
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
z_units = "Meters"
block_size = 5 # OPTIONAL defualt to 5
overwrite_data = True
# End Fill in Data
# ----------------------------------------------------------------------

# Parameter fields for python toolbox
#nodes_fc = parameters[0].valueAsText
#searchCells = parameters[1].valueAsText # Needs to be a int
#smooth_flag = parameters[2].valueAsText # Needs to be a int
#z_raster = parameters[3].valueAsText
#z_units = parameters[4].valueAsText
#block_size =  parameters[5].valueAsText
#overwrite_data = parameters[6].valueAsText

def nested_dict(): 
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc1(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the
    NODE_ID and X/Y coordinates as a nested dictionary"""
    nodeDict = nested_dict()
    incursorFields = ["NODE_ID", "SHAPE@X","SHAPE@Y"]

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
        
    # Check to see if all the new fields exist and add them if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")    

    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc, incursorFields,"",proj) as Inrows:
        if overwrite_data is True:
            for row in Inrows:
                nodeDict[row[0]]["POINT_X"] = row[1]
                nodeDict[row[0]]["POINT_Y"] = row[2]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[3] is None or row[3] == 0 or row[3] < -9998:
                    nodeDict[row[0]]["POINT_X"] = row[1]
                    nodeDict[row[0]]["POINT_Y"] = row[2]
              
    return nodeDict

def read_nodes_fc2(nodes_fc, overwrite_data, addFields):
    """Reads the input point file, adds new fields, and returns the
    STREAM_ID, NODE_ID, STREAM_KM, LENGTH, ELEVATION, and X/Y coordinates
    as a nested dictionary"""
    
    nodeDict = nested_dict()
    incursorFields = ["STREAM_ID", "STREAM_KM", "NODE_ID", "LENGTH", "ELEVATION", "SHAPE@X","SHAPE@Y"]
    
    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)    

    # Check to see if the 1st field exists if yes add it to 
    # the cursorfields to be retreived.
    if overwrite_data is False and (addFields[0] in existingFields) is True:
        incursorFields.append(addFields[0])
    else:
        overwrite_data = True

    # Check to see if all the new fields exist and add them if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")    

    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc,incursorFields,"",proj) as Inrows:
        if overwrite_data is True:
            for row in Inrows:
                nodeDict[row[0]][row[1]]["NODE_ID"] = row[2]
                nodeDict[row[0]][row[1]]["LENGTH"] = row[3]
                nodeDict[row[0]][row[1]]["ELEVATION"] = row[4]
                nodeDict[row[0]][row[1]]["POINT_X"] = row[5]
                nodeDict[row[0]][row[1]]["POINT_Y"] = row[6]
        else:
            for row in Inrows:
                # if the data is null or zero (0 = default for shapefile),
                # it is retreived and will be overwritten.
                if row[7] is None or row[7] < -9998:
                    nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                    nodeDict[row[0]][row[1]]["LENGTH"] = row[3]
                    nodeDict[row[0]][row[1]]["ELEVATION"] = row[4]
                    nodeDict[row[0]][row[1]]["POINT_X"] = row[5]
                    nodeDict[row[0]][row[1]]["POINT_Y"] = row[6]
    if len(nodeDict) == 0:
        sys.exit("The gradient field checked in the input point feature class "+
                 "have existing data. There is nothing to process. Exiting")    
    return nodeDict

def calculate_gradient(zList, len_list, smooth_flag):
    
    skipupNodes = [0]
    gradientList = [0 for i in zList]
    
    for i in range(1,len(zList)):
        z = zList[i]
        zUp = zList[i - 1 - max(skipupNodes)]
    
        # Check if the gradient is <= 0, if yes keep going until
        # it is positive again. Then recalculate the gradient over 
        # the longer distance		
        if z > zUp and smooth_flag is True:
            skipupNodes.append(max(skipupNodes) + 1)
        else:
            dx_meters = sum(len_list[i:i+max(skipupNodes)+1])
            #dx_meters = float(kmList[i] - kmList[i- max(skipdownNodes)]) * 1000
            gradient = (zUp - z) / dx_meters
            for Skip in skipupNodes:
                gradientList[i-Skip] = gradient
            skipupNodes = [0]
            
    return (gradientList)

def update_nodes_fc1(nodeDict, nodes_fc, addFields, nodes_to_query):
    """Updates the input point feature class with data from the nodes dictionary"""
    #print("Updating input point feature class")
    
    # Build a query to retreive just the nodes that needs updating
    if len(nodes_to_query) == 1:
        whereclause = """%s = %s""" % (arcpy.AddFieldDelimiters(nodes_fc, "NODE_ID"), nodes_to_query[0])
    else:
        whereclause = """%s IN %s""" % (arcpy.AddFieldDelimiters(nodes_fc, "NODE_ID"), tuple(nodes_to_query))

    with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:  
        for row in cursor:
            for f in xrange(0,len(addFields)):
                nodeID =row[0]
                row[f+1] = nodeDict[nodeID][addFields[f]]
                cursor.updateRow(row)

def update_nodes_fc2(nodeDict, nodes_fc, addFields, nodes_to_query):
    """Updates the input point feature class with data from
    the nodes dictionary"""
    print("Updating input point feature class")

    # Build a query to retreive just the nodes that needs updating
    if len(nodes_to_query) == 1:
        whereclause = """%s = %s""" % (arcpy.AddFieldDelimiters(nodes_fc, "NODE_ID"), nodes_to_query[0])
    else:
        whereclause = """%s IN %s""" % (arcpy.AddFieldDelimiters(nodes_fc, "NODE_ID"), tuple(nodes_to_query)) 

    with arcpy.da.UpdateCursor(nodes_fc,["STREAM_ID","NODE_ID","STREAM_KM"] + addFields, whereclause) as cursor:  
        for row in cursor:
            for f in xrange(0,len(addFields)):
                streamID = row[0]
                stream_km =row[2]
                row[f+3] = nodeDict[streamID][stream_km][addFields[f]]
                cursor.updateRow(row)

def create_block_list(nodeDict, nodes, buffer, block_size):
    """Returns two lists, one containting the coordinate extent
    for each block that will be itterativly extracted to an array
    and the other containing the stream and node IDs within each block extent."""
    
    print("Calculating block extents")    
    x_coord_list = [nodeDict[nodeID]["POINT_X"] for nodeID in nodes]
    y_coord_list = [nodeDict[nodeID]["POINT_Y"] for nodeID in nodes]    
    
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
                    nodes_in_block.append([nodeID, node_x, node_y])
                    
                    # Minimize the size of the block0 by the true 
                    # extent of the nodes in the block
                    if block_x_min > node_x: block_x_min = node_x
                    if block_x_max < node_x: block_x_max = node_x
                    if block_y_min > node_y: block_y_min = node_y
                    if block_y_max < node_y: block_y_max = node_y
            
            if nodes_in_block:
                # add the block extent for processing
                # order 0 left,      1 bottom,    2 right,     3 top
                block_extents.append([block_x_min - buffer, block_y_min - buffer,
                                      block_x_max + buffer, block_y_max + buffer])           
                block_nodes.append(nodes_in_block)
    
    return block_extents, block_nodes

def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append(int((easting - block_x_min) / x_cellsize))  # col, x
    xy.append(int((northing - block_y_max) / y_cellsize * -1))  # row, y 
    return xy

def sample_raster(block, nodes_in_block, z_raster, cellcoords, con_z_to_m):
    
    if con_z_to_m is not None:
        nodata_to_value = -9999 / con_z_to_m
    else:
        nodata_to_value = -9999
        
    # localize the block extent values
    block_x_min = block[0]
    block_y_min = block[1]
    block_x_max = block[2]
    block_y_max = block[3]
    
    x_cellsize = arcpy.Describe(z_raster).meanCellWidth
    y_cellsize = arcpy.Describe(z_raster).meanCellHeight    

    # Get the coordinates extent of the input raster
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
        raster_array = arcpy.RasterToNumPyArray(z_raster, arcpy.Point(block_x_min_center, block_y_min_center),
                                                ncols, nrows, nodata_to_value)
    except:
        tbinfo = traceback.format_exc()
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
        sys.exit(pymsg)    
    
    # convert array values to meters if needed
    if con_z_to_m is not None:
        raster_array = raster_array * con_z_to_m
    
    z_list = []
    z_node = []
    if raster_array.max() > -9999:
        # There is at least one pixel of data
        for node in nodes_in_block:
            xy = coord_to_array(node[1], node[2], block_x_min, block_y_max, x_cellsize, y_cellsize)
            z_sampleList = []
            for coord in cellcoords:
                # Calculate the cell X/Y based on the base coordinate movement
                cell_x = xy[0] + coord[0]
                cell_y = xy[1] + coord[1]
                z_sampleList.append(raster_array[cell_y,cell_x])
                
                # sample at node:
                if coord[0] == 0 and coord[1] == 0:
                    z_node = raster_array[cell_y,cell_x]
                
            # Remove no data values (-9999) unless they are all no data
            if not max(z_sampleList) < -9998:
                z_sampleList = [z for z in z_sampleList if z > -9999]
            # Get the lowest elevation        
            node.append(min(z_sampleList))
            node.append(z_node)
            z_list.append(node)
    
    else:
        # No data, add -9999 for elevation and z_node
        for node in nodes_in_block:
            node.append(-9999)
            node.append(-9999)
            z_list.append(node)
        
    return z_list

def from_z_units_to_meters_con(zUnits):
    """Returns the converstion factor to get from the input z
    units to meters"""
        
    try:
        con_z_to_m = float(zunits)
    except:
        if zUnits == "Meters":
            con_z_to_m = 1.0 
        elif zUnits == "Feet":
            con_z_to_m = 0.3048
        else: con_z_to_m = None # The conversion factor will not be used
    
    return con_z_to_m

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

#enable garbage collection
gc.enable()

try:
    print("Step 3: Sample Stream Elevations/Gradient")
    
    #keeping track of time
    startTime= time.time()

    # Check if the node fc exists
    if not arcpy.Exists(nodes_fc):
        arcpy.AddError("\nThis output does not exist: \n" +
                       "{0}\n".format(nodes_fc))
        sys.exit("This output does not exist: \n" +
                 "{0}\n".format(nodes_fc))     
    
    if overwrite_data is True: 
        env.overwriteOutput = True
    else:
        env.overwriteOutput = False

    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference
    proj_ele = arcpy.Describe(z_raster).spatialReference

    # Check to make sure the raster and input 
    # points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not ".format(nodes_fc,z_raster)+
                       "have the same projection."+
                       "Please reproject your data.")
        sys.exit("Input points and elevation raster do not have the "+
                 "same projection. Please reproject your data.")
    
    if block_size == "#": block_size = 5

    # Get the units conversion factor
    con_z_to_m = from_z_units_to_meters_con(z_units)
    con_from_m = from_meters_con(nodes_fc)
    
    # convert block size from km to meters to units of the node fc
    # in the future block size should be estimated based on availiable memory
    # memorysize = datatypeinbytes*nobands*block_size^2
    # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
    if block_size in ["#", ""]:
        block_size = int(con_from_m * 5000)
    else:
        block_size = int(con_from_m * block_size * 1000)
    
    # Get the elevation raster cell size
    cellsize = arcpy.Describe(z_raster).meanCellWidth
    
    # calculate the buffer distance (in raster spatial units) to add to 
    # the base bounding box when extracting to an array. The buffer is 
    # equal to the cellsize * the searchcells to make sure the block includes 
    # the surrounding cells at each corner
    buffer = int((searchCells + 1)* cellsize)    

    # Make a list of the base x/y coordinate movments 
    # from the node origin. These values will be 
    # multipled by the cell size.
    # searchCells = 0 samples at the node
    # searchCells = 1 cell width around node = 9 cells
    # searchCells = 2 cell widths around node = 25 cells ...
    cell_moves = [i for i in range(searchCells*-1, searchCells+1, 1)]
    cellcoords = list(itertools.product(cell_moves, cell_moves))
    
    # read the data into a nested dictionary
    addFields = ["ELEVATION", "Z_NODE"]
    nodeDict = read_nodes_fc1(nodes_fc, overwrite_data, addFields)
    if len(nodeDict) != 0:
        
        # Get a list of the nodes, sort them
        nodes = nodeDict.keys()
        nodes.sort()
        n_nodes = len(nodes)
        
        # Build the block list
        block_extents, block_nodes = create_block_list(nodeDict, nodes, buffer, block_size)
        
        # Itterate through each block, calculate sample coordinates,
        # convert raster to array, sample the raster
        total_samples = 0
        
        for p, block in enumerate(block_extents):
            nodes_in_block = block_nodes[p]
            nodes_in_block.sort()
            
            print("Processing block {0} of {1}".format(p + 1, len(block_extents)))
        
            z_list = sample_raster(block, nodes_in_block, z_raster, cellcoords, con_z_to_m)
            
            # Update the node fc
            for row in z_list:
                nodeDict[row[0]]["ELEVATION"] = row[3]
                nodeDict[row[0]]["Z_NODE"] = row[4]
                
            nodes_to_query = [row[0] for row in z_list]
            
            # Write the elevation data to the TTools point feature class 
            update_nodes_fc1(nodeDict, nodes_fc, addFields, nodes_to_query)
        
            total_samples = total_samples + len(z_list)
            del z_list
            gc.collect()
        
    else:
        print("The elevation field checked in the input point feature class " +
              "have existing data. Andvancing to gradient processing")
    del(nodeDict)
    
    # Start on gradients
    
    # read the data into a nested dictionary
    addFields = ["GRADIENT"]
    nodeDict = read_nodes_fc2(nodes_fc, overwrite_data, addFields)    
    nodes_to_query = []

    for n, streamID in enumerate(nodeDict):
        print("Calculating gradients stream {0} of {1}".format(n + 1, len(nodeDict)))
            
        stream_kms = nodeDict[streamID].keys()
        stream_kms.sort(reverse=True)        
    
        z_list = [nodeDict[streamID][km]["ELEVATION"] for km in stream_kms]
        len_list = [nodeDict[streamID][km]["LENGTH"] for km in stream_kms]
        
        # Calculate Gradient
        gradientList = calculate_gradient(z_list, len_list, smooth_flag)
        
        for i, km in enumerate(stream_kms):
            nodeDict[streamID][km]["GRADIENT"] = gradientList[i]
            nodes_to_query.append(nodeDict[streamID][km]["NODE_ID"])

    update_nodes_fc2(nodeDict, nodes_fc, addFields, nodes_to_query)
    
    endTime = time.time()
    gc.collect()  
    
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / n_nodes).microseconds
    print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))    

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