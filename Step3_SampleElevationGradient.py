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
arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#NodesFC = parameters[0].valueAsText
#searchCells = parameters[1].valueAsText # Needs to be a int
#SmoothFlag = parameters[2].valueAsText # Needs to be a int
#EleRaster = parameters[3].valueAsText
#EleUnits = parameters[4].valueAsText
#OverwriteData = parameters[5].valueAsText

# Start Fill in Data
NodesFC = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
searchCells = 9
SmoothFlag = False
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = "Feet"
OverwriteData = True
# End Fill in Data

def DistanceBetweenPoints(Xa,Ya,Xb,Yb):
    """Determines the distance between two points using the pythagorean theorem and returns distance between them in map units"""
    from math import sqrt, pow
    dist = math.sqrt(math.pow((Xa - Xb),2) + math.pow((Ya - Yb),2))
    return dist

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def ReadNodesFC(NodesFC, OverwriteData, AddFields):
    """Reads the input point file, adds new fields, and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    NodeDict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y",]

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

def ToMetersUnitConversion(inFeature):
    """Returns the conversion factor to get from the input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that is not projected or not recognized. Use a projected coordinate system preferably in linear units of feet or meters.".format(inFeature))
        sys.exit("Coordinate system is not projected or not recognized. Use a projected coordinate system, preferably in linear units of feet or meters.")   
    return con_to_m

def GetElevation(samplexy, EleRaster, LowElev, con_z_to_m):
    """Retreives the elevation value from an elevation raster at given x and y coordinate"""
    thevalue = arcpy.GetCellValue_management(EleRaster, samplexy)

    if str(thevalue.getOutput(0)) != "NoData":
        SampleElev= float(thevalue.getOutput(0)) *  con_z_to_m

    # See if the sample elevation is the lowest for this stream node
    if SampleElev < LowElev:
        LowElev = SampleElev
    return(LowElev)

def CalculateGradient(LowElev,LowElevUp):
    return (grad)

def UpdateNodesFC(NodeDict, NodesFC, AddFields): 
    """Updates the input point feature class with all new data from the nodes dictionary"""

    with arcpy.da.UpdateCursor(NodesFC,["STREAM_ID","NODE_ID"] + AddFields) as cursor:
        for row in cursor:
            for f in xrange(0,len(AddFields)):
                streamID = row[0]
                nodeID =row[1]
                row[f+2] = NodeDict[streamID][nodeID][AddFields[f]]
                cursor.updateRow(row)

def UpdateNodesFC1(NodeDict, streamID, nodeID, NodesFC, AddFields): 
    """Updates a single node in the input point feature class with data from the nodes dictionary"""

    # Build a query to retreive just the node row
    whereclause = """%s = %s""" % (arcpy.AddFieldDelimiters(NodesFC, "NODE_ID"), nodeID)

    with arcpy.da.UpdateCursor(NodesFC,["NODE_ID"] + AddFields, whereclause) as cursor:
        for row in cursor:   
            for f in xrange(0,len(AddFields)):
                row[f+1] = NodeDict[streamID][nodeID][AddFields[f]]
                cursor.updateRow(row)

def UpdateNodesFC2(UpdateValue, nodeID, NodesFC, UpdateField): 
    """Updates the input point feature class with the new value"""

    whereclause = """%s = %s""" % (arcpy.AddFieldDelimiters(NodesFC, "NODE_ID"), nodeID)

    with arcpy.da.UpdateCursor(NodesFC,["NODE_ID"] + UpdateField, whereclause) as cursor:
        for row in cursor:   
            row[1] = UpdateValue
            cursor.updateRow(row)

#enable garbage collection
gc.enable()

try:
    print("Step 3: Sample Stream Elevations/ Gradient")
    
    #keeping track of time
    startTime= time.time()    

    # Determine input point spatial units
    proj = arcpy.Describe(NodesFC).spatialReference
    proj_ele = arcpy.Describe(EleRaster).spatialReference

    # Check to make sure the raster and input points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(NodesFC,EleRaster))
        sys.exit("Input points and elevation raster do not have the same projection. Please reproject your data.")
        # reproject the inpoints
        #NodesFC_rp = os.path.dirname(NodesFC) + "\\rp_" + os.path.basename(NodesFC)
        #arcpy.Project_management(NodesFC,NodesFC_rp,proj_ele)
        #NodesFC = NodesFC_rp	

    nodexy_to_m = ToMetersUnitConversion(NodesFC)
    Elexy_to_m = ToMetersUnitConversion(EleRaster)
    units_con=  nodexy_to_m / Elexy_to_m

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
    n = 1
    for streamID in NodeDict:
        SkipDownNodes = [1]

        for nodeID in NodeDict[streamID]:

            # 1. Start with Elevation
            node_x = float(NodeDict[streamID][nodeID]["POINT_X"])
            node_y = float(NodeDict[streamID][nodeID]["POINT_Y"])
            LowElev = 99999
            UpdateField = ["GRADIENT"]

            for coord in cellcoords:
                # Calculate the cell X/Y based on the base coordinate movement
                cell_x = node_x + (coord[0] * CellSize)
                cell_y = node_y + (coord[1] * CellSize)

                samplexy = str(cell_x) + " " + str(cell_y) # xy string requirement for arcpy.GetCellValue

                # Get the lowest elevation
                LowElev = GetElevation(samplexy, EleRaster, LowElev, con_z_to_m)

            # Save the elevation back into the dictionary and update the point file
            NodeDict[streamID][nodeID]["ELEVATION"] = LowElev
            #UpdateNodesFC1(NodeDict,streamID, nodeID, NodesFC, AddFields)
            UpdateNodesFC2(LowElev, nodeID, NodesFC, UpdateField)

            # 2. CalculateGradient
            if nodeID > 0:
                UpdateField = ["GRADIENT"]
                # Calculate gradient between the current node and downstream node
                Ele = LowElev
                EleDown = NodeDict[streamID][nodeID - 1]["ELEVATION"]
                dx_meters = float(NodeDict[streamID][nodeID]["STREAM_KM"] - NodeDict[streamID][nodeID - 1]["STREAM_KM"]) * 1000

                # Check if the gradient is <= 0, if yes keep going until
                # it is positive again. Then recalculate the gradient over the longer distance		
                if Ele <= EleDown and SmoothFlag == True:
                    EleDownSkip = EleDown
                    #GradientDown = (Ele - EleDown) / (dx_meters)
                    SkipDownNodes.append(max(SkipDownNodes) + 1)
                    #NodeDict[streamID][nodeID -1]["GRADIENT"] = GradientDown
                    #UpdateNodesFC2(GradientDown, nodeID-1, NodesFC, UpdateField)
                else:
                    EleDownSkip = NodeDict[streamID][nodeID - max(SkipDownNodes)]["ELEVATION"]
                    dx_meters = float(NodeDict[streamID][nodeID]["STREAM_KM"] - NodeDict[streamID][nodeID - max(SkipDownNodes)]["STREAM_KM"]) * 1000
                    GradientDown = (Ele - EleDownSkip) / dx_meters
                    for Skip in SkipDownNodes:
                        NodeDict[streamID][nodeID - Skip]["GRADIENT"] = GradientDown
                        UpdateNodesFC2(GradientDown, nodeID - Skip, NodesFC, UpdateField)
                    SkipDownNodes = [1]
        n = n + 1
    
    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / n).microseconds
    print("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))    

    # Write the Elevation and Gradient to the TTools point feature class
    #UpdateNodesFC(NodeDict, NodesFC, AddFields)	    

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

    pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
    msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

    #arcpy.AddError(pymsg)
    #arcpy.AddError(msgs)

    print(pymsg)
    print(msgs)