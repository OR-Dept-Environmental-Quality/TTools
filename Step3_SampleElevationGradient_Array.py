#######################################################################################
# TTools
# Step 3: Sample Stream Elevations/ Gradient - v 0.1
# Ryan Michie

# Sample_ElevationsGradient will take an input point feature (from Step 1) and sample the input raster elevation
# to find the lowest elevation in a user defined search radius and calculate the gradient for each node in the downstream direction.

# INPUTS
# 0: Input TTools point feature class (NodesFC)
# 1: input the number of samples around the node (searchCells) 1. [0],  2. [9], 3. [25]
# 2: input flag for smoothing if gradient is zero or negative (SmoothFlag) 1. True, 2. False
# 2: input elevation raster (EleRaster)
# 3: input elevation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 5: input flag if existing data can be over written (OverwriteData) 1. True, 2. False

# OUTPUTS
# point feature class (edit NodesFC) - Added fields with ELEVATION and GRADIENT for node

# Future Updates

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys
import os
import string
import gc
import shutil
import time
from datetime import timedelta
import arcpy
import itertools
from arcpy import env
from math import sqrt, pow, ceil
from collections import defaultdict

# Check out the ArcGIS Spatial Analyst extension license
#arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#NodesFC = parameters[0].valueAsText
#searchCells = parameters[1].valueAsText # Needs to be a int
#SmoothFlag = parameters[2].valueAsText # Needs to be a int
#EleRaster = parameters[3].valueAsText
#EleUnits = parameters[4].valueAsText
#OverwriteData = parameters[5].valueAsText

# Start Fill in Data
#NodesFC = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
#searchCells = 9
#SmoothFlag = True
#EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
#EleUnits = "Feet"
#OverwriteData = True

NodesFC = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
searchCells = 0
SmoothFlag = True
#EleRaster = r"\\DEQWQNAS01\Lidar07\Willamette.gdb\BE"
EleRaster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
EleUnits = "Meters"
BlockSize = 5 # OPTIONAL defualt to 5
OverwriteData = True
# End Fill in Data

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def ReadNodesFC(NodesFC, OverwriteData, AddFields):
    """Reads the input point file, adds new fields, and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    NodeDict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]

    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(NodesFC):
        ExistingFields.append(f.name)    

    # Check to see if the 1st field exists if yes add it to the cursorfields to be retreived.
    if OverwriteData == False and (AddFields[0] in ExistingFields) == True:
        Incursorfields.append(AddFields[0])
    else:
        OverwriteData = True

    # Check to see if all the new fields exist and add them if not
    for f in AddFields:
        if (f in ExistingFields) == False:
            arcpy.AddField_management(NodesFC, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")    

    # Determine input point spatial units
    proj = arcpy.Describe(NodesFC).spatialReference

    with arcpy.da.SearchCursor(NodesFC,Incursorfields,"",proj) as Inrows:
        if OverwriteData == True:
            for row in Inrows:
                NodeDict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                NodeDict[row[0]][row[1]]["POINT_X"] = row[3]
                NodeDict[row[0]][row[1]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # if the data is null or zero (0 = default for shapefile), it is retreived and will be overwritten.
                if row[5] == None or row[5] == 0:
                    NodeDict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                    NodeDict[row[0]][row[1]]["POINT_X"] = row[3]
                    NodeDict[row[0]][row[1]]["POINT_Y"] = row[4]
    return NodeDict

def CalculateGradient(NodeList, SmoothFlag):
    
    SkipDownNodes = [1]
    for i in range(1,len(NodeList)):
        Elev = NodeList[i][5]
        ElevDown = NodeList[i-1][5]
        
        # Check if the gradient is <= 0, if yes keep going until
        # it is positive again. Then recalculate the gradient over the longer distance		
        if Elev <= ElevDown and SmoothFlag == True:
            ElevSkip = ElevDown
            SkipDownNodes.append(max(SkipDownNodes) + 1)
        else:
            EleDownSkip = NodeList[i- max(SkipDownNodes)][5]
            dx_meters = float(NodeList[i][4] - NodeList[i- max(SkipDownNodes)][4]) * 1000
            GradientDown = (Elev - EleDownSkip) / dx_meters
            for Skip in SkipDownNodes:
                NodeList[i-Skip][6] = GradientDown
            SkipDownNodes = [1]     
    return (NodeList)

def UpdateNodesFC(pointDict, NodesFC, AddFields): 
    """Updates the input point feature class with data from the nodes dictionary"""
    print("Updating input point feature class")

    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(NodesFC):
        ExistingFields.append(f.name)     

    # Check to see if the field exists and add it if not
    for f in AddFields:
        if (f in ExistingFields) == False:
            arcpy.AddField_management(NodesFC, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")   

    with arcpy.da.UpdateCursor(NodesFC,["STREAM_ID","NODE_ID"] + AddFields) as cursor:
        for row in cursor:
            for f in xrange(0,len(AddFields)):
                streamID = row[0]
                nodeID =row[1]
                row[f+2] = pointDict[streamID][nodeID][AddFields[f]]
                cursor.updateRow(row)

def CreateNodeList(NodeDict, streamID, BlockSize):
    """This builds a nested list of node information. The outer list holds all 
    the nodes on each 5 km length stream. This is done for memory managment
    when the raster is converted to an array. We can't convert the whole raster
    so instead we minimize the raster area down to the extent of the smaller block of nodes"""
    
    NodeList = []
    NodeBlocks = []
    # Build a list of km at intervals identifed by the block size 
    # if there is a stream longer than 6700 km we are not on Earth
    km_blocks = [x for x in range(BlockSize, 6700, BlockSize)] 
    i = 0
    
    Nodes = NodeDict.keys()
    Nodes.sort()
    
    for nodeID in Nodes:  
        origin_x = NodeDict[nodeID]["POINT_X"]
        origin_y = NodeDict[nodeID]["POINT_Y"]
        stream_km = NodeDict[nodeID]["STREAM_KM"]
        if stream_km < km_blocks[i]:
            NodeBlocks.append([origin_x, origin_y, streamID, nodeID, stream_km, 0, 0])
        else: # New block
            NodeList.append(NodeBlocks)
            NodeBlocks = []
            NodeBlocks.append([origin_x, origin_y, streamID, nodeID, stream_km, 0, 0])
            i = i + 1
    NodeList.append(NodeBlocks)
    return(NodeList)

def CoordToArray(easting, northing, bbox_upper_left):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append((easting - bbox_upper_left[0]) / bbox_upper_left[2])  # col, x
    xy.append((bbox_upper_left[1] - northing) / bbox_upper_left[3])  # row, y
    return xy

def GetElevations(NodeList, raster, cellcoords, con_z_to_m):
    
    cellsizeX = arcpy.Describe(raster).meanCellWidth
    cellsizeY = arcpy.Describe(raster).meanCellHeight    
    
    # calculate the buffer distance (in raster spatial units) to add to the raster bounding box when extracting to an array
    buffer = cellsizeX * 2    
    
    # calculate lower left corner and nrows/cols for the bounding box
    # first transpose the list so x and y coordinates are in the same list
    tlist = map(lambda *i: list(i), *NodeList)
    
    Xmin = min(tlist[0]) - buffer
    Ymin = min(tlist[1]) - buffer
    Ymax = max(tlist[1]) + buffer            
    ncols = (max(tlist[0]) + buffer - Xmin) / cellsizeX + 1
    nrows = (Ymax - Ymin) / cellsizeY + 1
    bbox_lower_left = arcpy.Point(Xmin, Ymin) # must be in raster map units
    bbox_upper_left = [Xmin, Ymax, cellsizeX, cellsizeY]
    nodata_to_value = -9999 / con_z_to_m
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        arry = arcpy.RasterToNumPyArray(raster, bbox_lower_left, ncols, nrows, nodata_to_value)
    except:
        import traceback
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
        sys.exit(pymsg)
        
    # convert array values to meters if needed
    arry = arry * con_z_to_m
    
    #print("Extracting raster values")
    for i in range(0,len(NodeList)):
        xy = CoordToArray(NodeList[i][0], NodeList[i][1], bbox_upper_left)
        SampleElev = []
        for coord in cellcoords:
            # Calculate the cell X/Y based on the base coordinate movement
            cell_x = xy[0] + coord[0]
            cell_y = xy[1] + coord[1]
            SampleElev.append(arry[cell_y,cell_x])
        
        # Remove no data values (-9999) unless they are all no data
        if not max(SampleElev) < -9998:
            SampleElev = [z for z in SampleElev if z > -9999]
        # Get the lowest elevation
        NodeList[i][5] = min(SampleElev)
    return NodeList

#enable garbage collection
gc.enable()

try:
    print("Step 3: Sample Stream Elevations/Gradient")
    
    #keeping track of time
    startTime= time.time()    

    # Determine input point spatial units
    proj = arcpy.Describe(NodesFC).spatialReference
    proj_ele = arcpy.Describe(EleRaster).spatialReference

    # Check to make sure the raster and input points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(NodesFC,EleRaster))
        sys.exit("Input points and elevation raster do not have the same projection. Please reproject your data.")
    
    if BlockSize == "#": BlockSize = 5

    # Set the converstion factor to get from the input elevation z units to meters
    if EleUnits == "Meters": #International meter
        con_z_to_m = 1 
    if EleUnits == "Feet": #International foot
        con_z_to_m = 0.3048
    if EleUnits == "Other": #Some other units
        sys.exit("Please modify your raster elevation z units so they are either in meters or feet")	

    # Get the elevation raster cell size
    CellSizeResult = arcpy.GetRasterProperties_management(EleRaster, "CELLSIZEX")
    CellSize = float(CellSizeResult.getOutput(0))

    # Make a list of the base x/y coordinate movments. These values will be multipled by the cell size.
    if searchCells == 0: 
        cellcoords = list(itertools.product([0],[0]))

    if searchCells == 9:
        cellcoords = list(itertools.product([-1,0,1],[-1,0,1]))

    if searchCells == 25:
        cellcoords = list(itertools.product([-2,-1,0,1,2],[-2,-1,0,1,2]))

    # read the data into a nested dictionary
    AddFields = ["ELEVATION","GRADIENT"]
    NodeDict = ReadNodesFC(NodesFC, OverwriteData, AddFields)
    NodeList = []
    n_nodes = 0
    n = 1
    for streamID in NodeDict:
        print("Processing stream %s of %s" % (n, len(NodeDict)))
        NodeList = CreateNodeList(NodeDict[streamID], streamID, BlockSize)
        
        # Get the Elevations
        for NodeBlock in NodeList:
            NodeBlock = GetElevations(NodeBlock, EleRaster, cellcoords, con_z_to_m)
        
            # Calculate Gradient
            NodeBlock = CalculateGradient(NodeBlock, SmoothFlag)
        
            # Update the NodeDict
            for row in NodeBlock:
                NodeDict[row[2]][row[3]]["ELEVATION"] = row[5]
                NodeDict[row[2]][row[3]]["GRADIENT"] = row[6]
            n_nodes = n_nodes + len(NodeBlock)
        n = n + 1
    
    endTime = time.time()
    gc.collect()
       
    # Write the Elevation and Gradient to the TTools point feature class
    UpdateNodesFC(NodeDict, NodesFC, AddFields)    
    
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
    import traceback
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]

    pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
    msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

    #arcpy.AddError(pymsg)
    #arcpy.AddError(msgs)

    print(pymsg)
    print(msgs)