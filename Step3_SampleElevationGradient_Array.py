########################################################################
# TTools
# Step 3: Sample Stream Elevations/ Gradient - v 0.9
# Ryan Michie

# Sample_ElevationsGradient will take an input point feature 
# (from Step 1) and sample the input raster elevation to find the 
# lowest elevation in a user defined search radius and calculate the 
# gradient for each node in the downstream direction.

# INPUTS
# 0: Input TTools point feature class (nodes_fc)
# 1: input the number of samples around the 
#     node (searchCells) 1. [0],  2. [9], 3. [25]
# 2: input flag for smoothing if gradient is zero 
#     or negative (smooth_flag) 1. True, 2. False
# 3: input elevation raster (z_raster)
# 4: input elevation raster z units (z_units) 1. "Feet", 2. "Meters"
# 5: input stream km distance to process within each array (block_size)
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

# Parameter fields for python toolbox
#nodes_fc = parameters[0].valueAsText
#searchCells = parameters[1].valueAsText # Needs to be a int
#smooth_flag = parameters[2].valueAsText # Needs to be a int
#z_raster = parameters[3].valueAsText
#z_units = parameters[4].valueAsText
#block_size =  parameters[5].valueAsText
#overwrite_data = parameters[6].valueAsText

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
searchCells = 0
smooth_flag = True
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
z_units = "Meters"
block_size = 5 # OPTIONAL defualt to 5
overwrite_data = False
# End Fill in Data
# ----------------------------------------------------------------------

def nested_dict(): 
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point file, adds new fields, and returns the
    STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    nodeDict = nested_dict()
    incursorFields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]

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
    return nodeDict

def calculate_gradient(zList, kmList, smooth_flag):
    
    skipdownNodes = [1]
    gradientList = [0 for i in zList]
    for i in range(1,len(zList)):
        z = zList[i]
        zDown = zList[i- max(skipdownNodes)]
        
        # Check if the gradient is <= 0, if yes keep going until
        # it is positive again. Then recalculate the gradient over 
        # the longer distance		
        if z < zDown and smooth_flag is True:
            skipdownNodes.append(max(skipdownNodes) + 1)
        else:
            dx_meters = float(kmList[i] - kmList[i- max(skipdownNodes)]) * 1000
            gradientDown = (z - zDown) / dx_meters
            for Skip in skipdownNodes:
                gradientList[i-Skip] = gradientDown
            skipdownNodes = [1]     
    return (gradientList)

def update_nodes_fc(nodeDict, nodes_fc, addFields): 
    """Updates the input point feature class with data from
    the nodes dictionary"""
    print("Updating input point feature class")

    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)     

    # Check to see if the field exists and add it if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")   

    with arcpy.da.UpdateCursor(nodes_fc,["STREAM_ID","NODE_ID"] + addFields) as cursor:
        for row in cursor:
            for f in xrange(0,len(addFields)):
                streamID = row[0]
                nodeID =row[1]
                row[f+2] = nodeDict[streamID][nodeID][addFields[f]]
                cursor.updateRow(row)

def create_node_list(nodeDict, streamID, block_size):
    """This builds a nested list of node information. The outer
    list holds all the nodes within a specified km extent (block_size).
    This is done for memory managmentwhen the raster is converted to an
    array. Really large arrays will use up all the memory and cause
    a crash"""
    
    nodeList = []
    nodeBlocks = []
    # Build a km list at intervals identifed by the block size 
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
            nodeBlocks.append([origin_x, origin_y, streamID, nodeID, stream_km])
        else:
            # New block
            nodeList.append([list(x) for x in zip(*nodeBlocks)])
            nodeBlocks = []
            nodeBlocks.append([origin_x, origin_y, streamID, nodeID, stream_km])
            i = i + 1
    nodeList.append([list(x) for x in zip(*nodeBlocks)])
    
    return(nodeList)

def coord_to_array(easting, northing, bbox_upper_left):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append((easting - bbox_upper_left[0]) / bbox_upper_left[2])  # col, x
    xy.append((northing - bbox_upper_left[1]) / bbox_upper_left[3] * -1)  # row, y     
    
    return xy

def get_elevations(x_coordList, y_coordList, z_raster, cellcoords, con_z_to_m):
    
    x_cellsize = arcpy.Describe(z_raster).meanCellWidth
    y_cellsize = arcpy.Describe(z_raster).meanCellHeight
    
    # Get the coordinates of the upper left cell corner of the input z_raster
    top_y = float(arcpy.GetRasterProperties_management(z_raster, "TOP").getOutput(0))
    left_x = float(arcpy.GetRasterProperties_management(z_raster, "LEFT").getOutput(0))

    # calculate the buffer distance (in z_raster spatial units) to 
    # add to the z_raster bounding box when extracting to an array
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
    bbox_lower_left = arcpy.Point(x_min, y_min) # must be in z_raster map units
    bbox_upper_left = [x_min + x_minoffset, y_max + y_maxoffset, x_cellsize, y_cellsize]
    nodata_to_value = -9999 / con_z_to_m
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        z_array = arcpy.RasterToNumPyArray(z_raster, bbox_lower_left, ncols, nrows, nodata_to_value)
    except:
        import traceback
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
        sys.exit(pymsg)
        
    # convert array values to meters if needed
    z_array = z_array * con_z_to_m
    
    zList = []
    
    #print("Extracting raster values")
    for i in range(0,len(x_coordList)):
        xy = coord_to_array(x_coordList[i], y_coordList[i], bbox_upper_left)
        z_sampleList = []
        for coord in cellcoords:
            # Calculate the cell X/Y based on the base coordinate movement
            cell_x = xy[0] + coord[0]
            cell_y = xy[1] + coord[1]
            z_sampleList.append(z_array[cell_y,cell_x])
        
        # Remove no data values (-9999) unless they are all no data
        if not max(z_sampleList) < -9998:
            z_sampleList = [z for z in z_sampleList if z > -9999]
        # Get the lowest elevation
        zList.append(min(z_sampleList))
    return zList

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

#enable garbage collection
gc.enable()

try:
    print("Step 3: Sample Stream Elevations/Gradient")
    
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

    # Get the elevation raster cell size
    cellsizeResult = arcpy.GetRasterProperties_management(z_raster, "CELLSIZEX")
    cellsize = float(cellsizeResult.getOutput(0))

    # Make a list of the base x/y coordinate movments. 
    # These values will be multipled by the cell size.
    # search only at the sample point
    if searchCells == 0: 
        cellcoords = list(itertools.product([0],[0]))

    # search 1 cell width around point = 9 cells
    if searchCells == 1:
        cellcoords = list(itertools.product([-1,0,1],[-1,0,1]))

    # search 2 cell widths around point = 25 cells
    if searchCells == 2:
        cellcoords = list(itertools.product([-2,-1,0,1,2],[-2,-1,0,1,2]))

    # read the data into a nested dictionary
    addFields = ["ELEVATION","GRADIENT"]
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
    nodeList = []
    n_nodes = 0
    n = 1
    for streamID in nodeDict:
        print("Processing stream %s of %s" % (n, len(nodeDict)))
        nodeList = create_node_list(nodeDict[streamID],
                                    streamID, block_size)
        
        # Get the Elevations
        for NodeBlock in nodeList:
            zList = get_elevations(NodeBlock[0],NodeBlock[1], z_raster,
                                   cellcoords, con_z_to_m)
            NodeBlock.append(zList)
        
            # Calculate Gradient
            gradientList = calculate_gradient(zList, NodeBlock[4],
                                              smooth_flag)
            NodeBlock.append(gradientList)
            
            # Transpose the list
            NodeBlock = [list(x) for x in zip(*NodeBlock)]
        
            # Update the nodeDict
            for row in NodeBlock:
                nodeDict[row[2]][row[3]]["ELEVATION"] = row[5]
                nodeDict[row[2]][row[3]]["GRADIENT"] = row[6]
            n_nodes = n_nodes + len(NodeBlock)
        n = n + 1
    
    endTime = time.time()
    gc.collect()
       
    # Write the Elevation and Gradient to the TTools point feature class
    update_nodes_fc(nodeDict, nodes_fc, addFields)    
    
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / n_nodes).microseconds
    print("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))    

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