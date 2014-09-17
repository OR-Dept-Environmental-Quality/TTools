#######################################################################################
# TTools
# Step 3: Sample Stream Elevations/ Gradient - v 0.1
# Ryan Michie

# Sample_ElevationsGradient will take an input point feature (from Step 1) and sample the input raster elevation
# to find the lowest elevation in a user defined search radius and calculate the gradient for each node in the downstream direction.

# INPUTS
# 0: Input TTools point feature class (PointFile)
# 1: input the number of samples around the node (searchCells) 1. [0],  2. [9], 3. [25]
# 2: input flag for smoothing if gradient is zero or negative (SmoothFlag) 1. True, 2. False
# 2: input elevation raster (EleRaster)
# 3: input elevation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 5: input flag if existing data can be over written (OverwriteData) 1. True, 2. False

# OUTPUTS
# point feature class (edit PointFile) - Added fields with ELEVATION and GRADIENT for node

# Future Updates

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
import itertools
from arcpy import env
from math import sqrt, pow, ceil
from collections import defaultdict

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#PointFile = parameters[0].valueAsText
#searchCells = parameters[1].valueAsText # Needs to be a int
#SmoothFlag = parameters[2].valueAsText # Needs to be a int
#EleRaster = parameters[3].valueAsText
#EleUnits = parameters[4].valueAsText
#OverwriteData = parameters[5].valueAsText

# Start Fill in Data
PointFile = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
searchCells = 9
SmoothFlag = False
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = 1
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

def ReadPointFile(PointFile, OverwriteData, AddFields):
    """Reads the input point file, adds new fields, and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    NodeDict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y",]

    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(PointFile):
        ExistingFields.append(f.name)    

    # Check to see if the 1st field exists if yes add it to the cursorfields to be retreived.
    if OverwriteData == False and (AddFields[0] in ExistingFields) == True:
        Incursorfields.append(AddFields[0])
    else:
        OverwriteData = True

    # Check to see if all the new fields exist and add them if not
    for f in AddFields:
        if (f in ExistingFields) == False:
            arcpy.AddField_management(PointFile, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")    

    # Determine input point spatial units
    proj = arcpy.Describe(PointFile).spatialReference

    with arcpy.da.SearchCursor(PointFile,Incursorfields,"",proj) as Inrows:
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
    unitCode = arcpy.Describe(inFeature).SpatialReference.linearUnitCode
    if unitCode == 9001: #International meter
        units_con = 1 
    if unitCode == 9002: #International foot
        units_con = 0.3048
    if unitCode == 9003: #US Survey foot
        units_con = 1200/3937
    if unitCode == 9005: #Clarke's foot
        units_con =  0.3047972654 
    if unitCode not in [9001,9002,9003,9005]:
        arcpy.AddError("{0} has an unrecognized spatial reference. Use projection with units of feet or meters.".format(inFeature))
        system.exit("Unrecognized spatial reference. Use projection with units of feet or meters.")
    return units_con

def GetElevation(samplexy, EleRaster, LowElev, eleZ_to_m):
    """Retreives the elevation value from an elevation raster at given x and y coordinate"""
    thevalue = arcpy.GetCellValue_management(EleRaster, samplexy)

    if str(thevalue.getOutput(0)) != "NoData":
        SampleElev= float(thevalue.getOutput(0)) *  eleZ_to_m

    #See if the sample elevation is the lowest for this stream node
    if SampleElev < LowElev:
        LowElev = SampleElev
    return(LowElev)

def CalculateGradient(LowElev,LowElevUp):
    return (grad)

def UpdatePointFile(pointDict, pointfile, AddFields): 
    """Updates the input point feature class with all new data from the nodes dictionary"""

    with arcpy.da.UpdateCursor(pointfile,["STREAM_ID","NODE_ID"] + AddFields) as cursor:
        for row in cursor:
            for f in xrange(0,len(AddFields)):
                streamID = row[0]
                nodeID =row[1]
                row[f+2] = pointDict[streamID][nodeID][AddFields[f]]
                cursor.updateRow(row)

def UpdatePointFile1(NodeDict,streamID, nodeID, PointFile, AddFields): 
    """Updates a single node in the input point feature class with data from the nodes dictionary"""

    # Build a query to retreive just the node row
    whereclause = """%s = %s""" % (arcpy.AddFieldDelimiters(PointFile, "NODE_ID"), nodeID)

    with arcpy.da.UpdateCursor(PointFile,["NODE_ID"] + AddFields, whereclause) as cursor:
        for row in cursor:   
            for f in xrange(0,len(AddFields)):
                row[f+1] = NodeDict[streamID][nodeID][AddFields[f]]
                cursor.updateRow(row)

def UpdatePointFile2(UpdateValue, nodeID, PointFile, UpdateField): 
    """Updates the input point feature class with the new value"""

    whereclause = """%s = %s""" % (arcpy.AddFieldDelimiters(PointFile, "NODE_ID"), nodeID)

    with arcpy.da.UpdateCursor(PointFile,["NODE_ID"] + UpdateField, whereclause) as cursor:
        for row in cursor:   
            row[1] = UpdateValue
            cursor.updateRow(row)

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()    

    # Determine input point spatial units
    proj = arcpy.Describe(PointFile).spatialReference
    proj_ele = arcpy.Describe(EleRaster).spatialReference

    # Check to make sure the raster and input points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(PointFile,EleRaster))
        sys.exit("Input points and elevation raster do not have the same projection. Please reproject your data.")
        # reproject the inpoints
        #PointFile_rp = os.path.dirname(PointFile) + "\\rp_" + os.path.basename(PointFile)
        #arcpy.Project_management(PointFile,PointFile_rp,proj_ele)
        #PointFile = PointFile_rp	

    nodexy_to_m = ToMetersUnitConversion(PointFile)
    Elexy_to_m = ToMetersUnitConversion(EleRaster)
    units_con=  nodexy_to_m / Elexy_to_m

    # Determine the elevation Z units conversion into meters
    if EleUnits == 1: # Feet
        eleZ_to_m = 0.3048
    if EleUnits == 2: # Meters
        eleZ_to_m = 1
    if EleUnits == 3: # Other
        sys.exit("Please modify your raster elevtion units to feet or meters.") 	

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
    NodeDict = ReadPointFile(PointFile, OverwriteData, AddFields)
    n = 0
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
                LowElev = GetElevation(samplexy, EleRaster, LowElev, eleZ_to_m)

            # Save the elevation back into the dictionary and update the point file
            NodeDict[streamID][nodeID]["ELEVATION"] = LowElev
            #UpdatePointFile1(NodeDict,streamID, nodeID, PointFile, AddFields)
            UpdatePointFile2(LowElev, nodeID, PointFile, UpdateField)

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
                    #UpdatePointFile2(GradientDown, nodeID-1, PointFile, UpdateField)
                else:
                    EleDownSkip = NodeDict[streamID][nodeID - max(SkipDownNodes)]["ELEVATION"]
                    dx_meters = float(NodeDict[streamID][nodeID]["STREAM_KM"] - NodeDict[streamID][nodeID - max(SkipDownNodes)]["STREAM_KM"]) * 1000
                    GradientDown = (Ele - EleDownSkip) / dx_meters
                    for Skip in SkipDownNodes:
                        NodeDict[streamID][nodeID - Skip]["GRADIENT"] = GradientDown
                        UpdatePointFile2(GradientDown, nodeID - Skip, PointFile, UpdateField)
                    SkipDownNodes = [1]

    # Write the Elevation and Gradient to the TTools point feature class
    #UpdatePointFile(NodeDict, PointFile, AddFields)	    

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