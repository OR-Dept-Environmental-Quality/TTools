#!/usr/bin/python

"""
TTools Step 3: Sample Elevation and Gradient

This script will take an input point feature (from Step 1) and sample the input raster elevation to find the
lowest elevation in a user defined search radius and calculate the gradient for each node in the downstream direction.

REQUIREMENTS
TTools steps 1 must be run before Step 3.
ESRI ArcGIS Pro
Python 3.7+

INPUT VARIABLES
0: nodes_fc:
Path to the TTools point feature class.

1: searchCells
The number of cells to search around the node for the lowest elevation.

2: smooth_flag
Boolean (True/False) flag to indicate if smoothing should occur when the gradient is <= 0. When gradients are <= 0
the algorithm moves downstream until it a positive gradient is found and then recalculate the gradient over the longer
distance and applies the resulting value to all nodes over that distance.

3: z_raster:
Path and name of the ground elevation raster.

4: z_units:
z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is not in feet or meters the elevation
values must be converted.

6: block_size:
The x and y size in kilometers for the z_raster blocks pulled into an array. Start with 5 if you aren't sure and reduce
if there is an error. To increase processing the z_raster is subdivided iteratively into smaller blocks and pulled
into arrays for processing. Very large block sizes may use a lot of memory and result in an error.

7: overwrite_data:
True/False flag if existing data in nodes_fc can be overwritten.

OUTPUTS
0. nodes_fc:
New fields ELEVATION and GRADIENT are added into nodes_fc.

"""
# Import system modules
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
# Start input variables
nodes_fc = r"C:\workspace\ttools_tests\TTools_py39\jc_test_py39.gdb\jc_stream_nodes_py39"
searchCells = 2
smooth_flag = True
z_raster = r"C:\workspace\ttools_tests\JohnsonCreek.gdb\jcw_be_m_mosaic"
z_units = "Meters"
block_size = 5
overwrite_data = True
# End input variables
# ----------------------------------------------------------------------

# Future Updates
# eliminate arcpy and use gdal for reading/writing feature class data

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
    # Grabs last field because often the first field, emergent, is zero
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
                nodeID =row[0]
                nodeDict[nodeID]["POINT_X"] = row[1]
                nodeDict[nodeID]["POINT_Y"] = row[2]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[3] is None or row[3] == 0 or row[3] < -9998:
                    nodeID =row[0]
                    nodeDict[nodeID]["POINT_X"] = row[1]
                    nodeDict[nodeID]["POINT_Y"] = row[2]
              
    return nodeDict

def read_nodes_fc2(nodes_fc, overwrite_data, addFields):
    """Reads the input point file, adds new fields, and returns the
    STREAM_ID, STREAM_KM, NODE_ID, LENGTH, ELEVATION, and X/Y coordinates
    as a nested dictionary"""
    
    nodeDict = nested_dict()
    incursorFields = ["STREAM_ID", "STREAM_KM", "NODE_ID", "LENGTH", "ELEVATION", "SHAPE@X","SHAPE@Y"]
    
    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)    

    # Check to see if the 1st field exists if yes add it to 
    # the cursorfields to be retrieved.
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
                streamID = row[0]
                stream_km = row[1]
                nodeDict[streamID][stream_km]["NODE_ID"] = row[2]
                nodeDict[streamID][stream_km]["LENGTH"] = row[3]
                nodeDict[streamID][stream_km]["ELEVATION"] = row[4]
                nodeDict[streamID][stream_km]["POINT_X"] = row[5]
                nodeDict[streamID][stream_km]["POINT_Y"] = row[6]
        else:
            for row in Inrows:
                # if the data is null or zero (0 = default for shapefile),
                # it is retrieved and will be overwritten.
                if row[7] is None or row[7] < -9998:
                    streamID = row[0]
                    stream_km = row[1]                    
                    nodeDict[streamID][stream_km]["NODE_ID"] = row[2] 
                    nodeDict[streamID][stream_km]["LENGTH"] = row[3]
                    nodeDict[streamID][stream_km]["ELEVATION"] = row[4]
                    nodeDict[streamID][stream_km]["POINT_X"] = row[5]
                    nodeDict[streamID][stream_km]["POINT_Y"] = row[6]
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
    """Updates the input point feature class with data from
    the nodes dictionary with node_id as the primary key"""
    #print("Updating input point feature class")
    
    # Build a query to retrieve just the nodes that needs updating
    whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))

    with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:
        for row in cursor:
            for f, field in enumerate(addFields):
                nodeID =row[0]
                row[f+1] = nodeDict[nodeID][field]
                cursor.updateRow(row)

def update_nodes_fc2(nodeDict, nodes_fc, addFields, nodes_to_query):
    """Updates the input point feature class with data from
    the nodes dictionary with stream IDs as the primary key"""
    print("Updating input point feature class")

    # Build a query to retreive just the nodes that needs updating
    whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))

    with arcpy.da.UpdateCursor(nodes_fc,["STREAM_ID","NODE_ID","STREAM_KM"] + addFields, whereclause) as cursor:  
        for row in cursor:
            for f, field in enumerate(addFields):
                streamID = row[0]
                stream_km =row[2]
                row[f+3] = nodeDict[streamID][stream_km][field]
                cursor.updateRow(row)

def create_block_list(nodeDict, nodes, block_size, buffer):
    """Returns two lists, one containing the coordinate extent
    for each block that will be iteratively extracted to an array
    and the other containing the stream and node IDs within each
    block extent."""
    
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

            block_x_max = block0_x_max
            block_x_min = block0_x_min
            block_y_max = block0_y_max
            block_y_min = block0_y_min
            
            nodes_in_block = []
            for nodeID in nodes:
                node_x = nodeDict[nodeID]["POINT_X"]
                node_y = nodeDict[nodeID]["POINT_Y"]
                if (block0_x_min <= node_x <= block0_x_max and
                    block0_y_min <= node_y <= block0_y_max):

                    nodes_in_block.append([nodeID, node_x, node_y])

            if nodes_in_block:

                # Minimize the size of the block0 by the true
                # extent of the nodes in the block
                node_x_min = min([nodeDict[nodeID]["POINT_X"] for nodeID in nodes_in_block])
                node_y_min = min([nodeDict[nodeID]["POINT_Y"] for nodeID in nodes_in_block])
                node_x_max = max([nodeDict[nodeID]["POINT_X"] for nodeID in nodes_in_block])
                node_y_max = max([nodeDict[nodeID]["POINT_Y"] for nodeID in nodes_in_block])

                if block0_x_min < node_x_min: block_x_min = node_x_min
                if block0_x_max > node_x_max: block_x_max = node_x_max
                if block0_y_min < node_y_min: block_y_min = node_y_min
                if block0_y_max > node_y_max: block_y_max = node_y_max

                # Add the block extent for processing and the buffer distance.
                # Because the node coordinates are used to determine if the node is within the block extent,
                # a buffer distance is added to ensure each block extent will include all samples around nodes
                # located at the edge of the block.

                # order 0 left,      1 bottom,    2 right,     3 top
                block_extents.append((block_x_min - buffer, block_y_min - buffer,
                                      block_x_max + buffer, block_y_max + buffer))
                block_nodes.append(nodes_in_block)
    
    return block_extents, block_nodes

def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
    """converts x/y coordinates to col and row of the array"""
    col_x = int(((easting - block_x_min) - ((easting - block_x_min) % x_cellsize)) / x_cellsize)
    row_y = int(((block_y_max - northing) - ((block_y_max - northing) % y_cellsize)) / y_cellsize)
    return [col_x, row_y]

def sample_raster(block, nodes_in_block, z_raster, cellcoords, con_z_to_m):
    
    if con_z_to_m is not None:
        nodata_to_value = -9999 / con_z_to_m
    else:
        nodata_to_value = -9999

    x_cellsize = float(arcpy.GetRasterProperties_management(raster, "CELLSIZEX").getOutput(0))
    y_cellsize = float(arcpy.GetRasterProperties_management(raster, "CELLSIZEY").getOutput(0))

    # localize the block extent values and add one cell distance to the size to ensure all cells that need to be
    # sampled are included in the array.
    block_x_min = block[0] - x_cellsize
    block_y_min = block[1] - y_cellsize
    block_x_max = block[2] + x_cellsize
    block_y_max = block[3] + y_cellsize

    # Get the coordinate extent of the input raster
    raster_x_min = float(arcpy.GetRasterProperties_management(raster, "LEFT").getOutput(0))
    raster_y_min = float(arcpy.GetRasterProperties_management(raster, "BOTTOM").getOutput(0))
    raster_x_max = float(arcpy.GetRasterProperties_management(raster, "RIGHT").getOutput(0))
    raster_y_max = float(arcpy.GetRasterProperties_management(raster, "TOP").getOutput(0))

    # Calculate the block x and y offset from the raster and adjust
    # the block coordinates so they are at the raster cell corners.
    # This is for ESRI's RastertoNumpyArray function which defaults to the adjacent
    # lower left cell
    block_x_min_corner = block_x_min - ((block_x_min - raster_x_min) % x_cellsize)
    block_y_min_corner = block_y_min - ((block_y_min - raster_y_min) % y_cellsize)
    block_x_max_corner = block_x_max + ((raster_x_max - block_x_max) % x_cellsize)
    block_y_max_corner = block_y_max + ((raster_y_max - block_y_max) % y_cellsize)

    # calculate the number of cols/rows from the lower left
    ncols = int((block_x_max_corner - block_x_min_corner) / x_cellsize)
    nrows = int((block_y_max_corner - block_y_min_corner) / y_cellsize)
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        raster_array = arcpy.RasterToNumPyArray(z_raster, arcpy.Point(block_x_min_corner, block_y_min_corner),
                                                ncols, nrows, nodata_to_value)
    except:
        tbinfo = traceback.format_exc()
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
        sys.exit(pymsg)    
    
    # convert array values to meters if needed
    if con_z_to_m is not None:
        raster_array = raster_array * con_z_to_m
    
    z_list = []
    if raster_array.max() > -9999:
        # There is at least one pixel of data
        for node in nodes_in_block:
            xy = coord_to_array(node[1], node[2], block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize)
            z_sampleList = []
            for coord in cellcoords:
                # Calculate the cell X/Y based on the base coordinate movement
                cell_x = xy[0] + coord[0]
                cell_y = xy[1] + coord[1]
                z_sampleList.append(raster_array[cell_y,cell_x])
                
            # sample at node:
            z_node = raster_array[xy[1], xy[0]]
                
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
    """Returns the conversion factor to get from the input z
    units to meters"""
        
    try:
        con_z_to_m = float(zUnits)
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
    # in the future block size should be estimated based on available memory
    # memorysize = datatypeinbytes*nobands*block_size^2
    # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
    if block_size in ["#", ""]:
        block_size = int(con_from_m * 5000)
    else:
        block_size = int(con_from_m * block_size * 1000)
    
    # Get the elevation raster cell size
    cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEX").getOutput(0))
    
    # calculate the buffer distance (in raster spatial units) to add to 
    # the base bounding box when extracting to an array. The buffer is 
    # equal to the cellsize * the searchcells to make sure the block includes 
    # the surrounding cells at each corner
    buffer = int((searchCells + 1)* cellsize)    

    # Make a list of the base x/y coordinate movements
    # from the node origin. These values will be 
    # multiplied by the cell size.
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
        nodes = list(nodeDict.keys())
        nodes.sort()
        n_nodes = len(nodes)
        
        # Build the block list
        block_extents, block_nodes = create_block_list(nodeDict, nodes, block_size, buffer)
        
        # Iterate through each block, calculate sample coordinates,
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
              "have existing data. Advancing to gradient processing")
    del(nodeDict)
    
    # Start on gradients
    
    # read the data into a nested dictionary
    addFields = ["GRADIENT"]
    nodeDict = read_nodes_fc2(nodes_fc, overwrite_data, addFields)    
    nodes_to_query = []

    for n, streamID in enumerate(nodeDict):
        print("Calculating gradients stream {0} of {1}".format(n + 1, len(nodeDict)))
            
        stream_kms = list(nodeDict[streamID].keys())
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