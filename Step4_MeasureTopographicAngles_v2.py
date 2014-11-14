#######################################################################################
# TTools
# Step 4: Measure Topographic Angles - v 0.95
# Ryan Michie

# Measure_Topographic_Angles will take an input point feature (from Step 1) and calculate the maximum
# topographic elevation and the the slope angle from each node in different directions.

# INPUTS
# 0: Input TTools point feature class(NodesFC)
# 1: input the directions to sample (Directions) 1. [W,S,E], 2. [NE,E,SE,S,SW,W,NW,N]
# 2: input the maximum km distance to search (MaxSearchDistance_km)
# 3: input elevation raster (EleRaster)
# 4: input elevation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 5: output sample point file name/path (outpoint_final)
# 6: input flag if existing data can be over written (OverwriteData) 1. True, 2. False

# OUTPUTS
# point feature class (edit NodesFC) - Added fields with topographic shade angles for each direction at each node
# point feature class (new) - point for each x/y location of the maximum elevation.

# Future Updates
# eliminate arcpy and use gdal for reading/writing feature class data

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from arcpy import env
from math import degrees, radians, sin, cos, atan, hypot, ceil
from collections import defaultdict
import numpy as np

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#NodesFC = parameters[0].valueAsText
#Directions = parameters[1].valueAsText # Needs to be a long
#MaxSearchDistance_km = parameters[2].valueAsText
#EleRaster = parameters[3].valueAsText
#EleUnits = parameters[4].valueAsText
#outpoint_final = parameters[5].valueAsText
#OverwriteData = parameters[6].valueAsText True/False

# Start Fill in Data
NodesFC = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
Directions = 2
MaxSearchDistance_km = 1
skippy = False
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = "Feet"
outpoint_final = r"D:\Projects\TTools_9\Example_data.gdb\topo_samples"
OverwriteData = True
# End Fill in Data

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def ReadNodesFC(NodesFC, OverwriteData, AddFields):
    """Reads the input point feature class and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y",]

    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(NodesFC):
        ExistingFields.append(f.name)     

    # Check to see if the 1st field exists if yes add it.
    if OverwriteData == False and (AddFields[0] in ExistingFields) == True:
        Incursorfields.append(AddFields[0])
    else:
        OverwriteData = True

    # Determine input point spatial units
    proj = arcpy.Describe(NodesFC).spatialReference

    with arcpy.da.SearchCursor(NodesFC,Incursorfields,"",proj) as Inrows:
        if OverwriteData == True:
            for row in Inrows:
                pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
                pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # if the data is null or zero (0 = default for shapefile), it is retreived and will be overwritten.
                if row[5] == None or row[5] == 0:
                    pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                    pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
                    pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
    if len(pnt_dict) == 0:
        sys.exit("The fields checked in the input point feature class have existing data. There is nothing to process. Exiting")
            
    return(pnt_dict)

def CreateTopoFC(pointList, NodesFC, proj):
    """Creates the output topo point feature class using the data from the nodes list"""
    #arcpy.AddMessage("Exporting Data")
    print("Exporting Data")
    
    #Create an empty output with the same projection as the input polyline
    cursorfields = ["POINT_X","POINT_Y","STREAM_ID","NODE_ID","AZIMUTH","TOPOANGLE","TOPO_ELE","NODE_ELE","ELE_CHANGE","TOPODIS","SEARCHDIS","NA_SAMPLES"]
    arcpy.CreateFeatureclass_management(os.path.dirname(NodesFC),os.path.basename(NodesFC), "POINT","","DISABLED","DISABLED",proj)

    # Add attribute fields # TODO add dictionary of field types so they aren't all double
    for f in cursorfields:
        arcpy.AddField_management(NodesFC, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(NodesFC, ["SHAPE@X","SHAPE@Y"] + cursorfields) as cursor:
        for row in pointList:
            cursor.insertRow(row)
            
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

def ToMetersUnitConversion(inFeature):
    """Returns the conversion factor to get from the input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that is not projected or not recognized. Use a projected coordinate system preferably in linear units of feet or meters.".format(inFeature))
        sys.exit("Coordinate system is not projected or not recognized. Use a projected coordinate system, preferably in linear units of feet or meters.")   
    return con_to_m
    
def FromMetersUnitConversion(inFeature):
    """Returns the conversion factor to get from meters to the spatial units of the input feature class"""
    try:
        con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that is not projected or not recognized. Use a projected coordinate system preferably in linear units of feet or meters.".format(inFeature))
        sys.exit("Coordinate system is not projected or not recognized. Use a projected coordinate system, preferably in linear units of feet or meters.")   
    return con_from_m

def BuildSearchArray(MinSearchDistance, MaxSearchDistance, cellsize, skippy):
    """Build a numpy array from the minimum to the max search distance by increments of the cellsize or via the skippy algorithm"""
    if skippy == True:
        # This is the skippy algorithm from Greg Pelletier
        SearchDistance = MinSearchDistance
        DistList = [MinSearchDistance]
        i = 1
        while not SearchDistance > MaxSearchDistance:
            if i <= 10:
                SearchDistance = SearchDistance + (cellsize)
            if 10 < i <= 20:
                SearchDistance = SearchDistance + (cellsize * 3)
            if 20 < i <= 40:
                SearchDistance = SearchDistance + (cellsize * 6)
            if 40 < i <= 50:
                SearchDistance = SearchDistance + (cellsize * 12)
            if 50 < i <= 60:
                SearchDistance = SearchDistance + (cellsize * 25)
            if i > 60:
                SearchDistance = SearchDistance + (cellsize * 50)
            DistList.append(SearchDistance)
            i = i + 1
        Distarry = np.array(DistList)
    else:
        Distarry = np.arange(MinSearchDistance, MaxSearchDistance, cellsize)
    return Distarry

def CoordToArray(easting, northing, bbox_upper_left):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append((easting - bbox_upper_left[0]) / bbox_upper_left[2])  # col, x
    xy.append((bbox_upper_left[1] - northing) / bbox_upper_left[3])  # row, y
    return xy

def GetShadeAngles(NodeDict, streamID, EleRaster, azimuths, azimuthdisdict, buffer, con_z_to_m):
    """This gets the maximum shade angle for each stream and outputs the data as a list"""
    TopoList = []
    
    cellsizeX = arcpy.Describe(EleRaster).meanCellWidth
    cellsizeY = arcpy.Describe(EleRaster).meanCellHeight       
    n = 1
    for nodeID in NodeDict:
        #print("Processing node %s of %s" % (n, len(NodeDict)))
        origin_x = NodeDict[nodeID]["POINT_X"]
        origin_y = NodeDict[nodeID]["POINT_Y"]
                 
        for a in azimuths:
            CoordList = []
            
            for searchdistance in azimuthdisdict[a]:
                # Calculate the x and y coordinate of the landcover sample location in the spatial units of the inputFC                
                topoX = (searchdistance * con_from_m * sin(radians(a))) + origin_x
                topoY = (searchdistance * con_from_m * cos(radians(a))) + origin_y

                # Add the all the coordinates to the list
                CoordList.append([topoX, topoY])    
            
            # calculate lower left corner coordinate and nrows/cols for the bounding box
            # first transpose the list so x and y coordinates are in the same vector
            tlist = map(lambda *i: list(i), *CoordList)
                
            Xmin = min(tlist[0]) - buffer
            Ymin = min(tlist[1]) - buffer
            Ymax = max(tlist[1]) + buffer            
            ncols = (max(tlist[0]) + buffer - Xmin) / cellsizeX
            nrows = (Ymax - Ymin) / cellsizeY
            bbox_lower_left = arcpy.Point(Xmin, Ymin) # must be in raster map units
            bbox_upper_left = [Xmin, Ymax, cellsizeX, cellsizeY]
            nodata_to_value = -9999 / con_z_to_m
            
            # Construct the array. Note returned array is (row, col) so (y, x)
            Zarry = arcpy.RasterToNumPyArray(EleRaster, bbox_lower_left, ncols, nrows, nodata_to_value)
            
            # convert array values to meters if needed
            Zarry = Zarry * con_z_to_m             
            
            TopoZList = []
            
            #print("Extracting raster values")
            for i in range(0,len(CoordList)):
                xy = CoordToArray(CoordList[i][0], CoordList[i][1], bbox_upper_left)
                TopoZList.append(Zarry[xy[1], xy[0]])
                #TopoList[i].append(Zarry[xy[1], xy[0]])
            
            TopoZarry = np.array(TopoZList)
            Disarry = azimuthdisdict[a] * con_to_m
            Shadearry = np.degrees(np.arctan((TopoZarry - TopoZarry[0]) / Disarry))
            # Take out the off raster samples
            naindex = np.where(TopoZarry < -9998)
            for x in naindex[0]: Shadearry[x] = -9999            
            # Find the max shade angle
            ShadeAngle = Shadearry.max()
            # array index at the max shade angle 
            arryindex = np.where(Shadearry==ShadeAngle)[0][0]
            ShadeZ = TopoZarry[arryindex]
            ZChange = ShadeZ - TopoZarry[0]
            ShadeDistance = azimuthdisdict[a][arryindex]
            SearchDistance = azimuthdisdict[a].max()
            ShadeAngle_X = (ShadeDistance * con_from_m * sin(radians(a))) + origin_x #CoordList[0][0]
            ShadeAngle_Y = (ShadeDistance * con_from_m * cos(radians(a))) + origin_y #CoordList[0][1]
            offRasterSamples = (TopoZarry > -9998).sum()

    
            TopoList.append([ShadeAngle_X, ShadeAngle_Y, ShadeAngle_X, ShadeAngle_Y, streamID, nodeID, a, ShadeAngle, ShadeZ, TopoZarry[0], ZChange, ShadeDistance, SearchDistance, offRasterSamples])
        n = n + 1
    return TopoList

def GetShadeAngles2(NodeDict, streamID, EleRaster, azimuths, azimuthdisdict, buffer, con_z_to_m):
    """This gets the maximum shade angle for each stream and outputs the data as a list"""
    TopoList = []
    CoordList = []
    
    cellsizeX = arcpy.Describe(EleRaster).meanCellWidth
    cellsizeY = arcpy.Describe(EleRaster).meanCellHeight       
    
    for nodeID in NodeDict:
        origin_x = NodeDict[nodeID]["POINT_X"]
        origin_y = NodeDict[nodeID]["POINT_Y"]      
        for a in azimuths:
            for searchdistance in azimuthdisdict[a]:
                # Calculate the x and y coordinate of the landcover sample location in the spatial units of the inputFC                
                topoX = (searchdistance * con_from_m * sin(radians(a))) + origin_x
                topoY = (searchdistance * con_from_m * cos(radians(a))) + origin_y

                # Add the all the coordinates to the list
                CoordList.append([topoX, topoY, nodeID, a])
            
    # calculate lower left corner coordinate and nrows/cols for the bounding box
    # first transpose the list so x and y coordinates are in the same vector
    tlist = map(lambda *i: list(i), *CoordList)
      
    Xmin = min(tlist[0]) - buffer
    Ymin = min(tlist[1]) - buffer
    Ymax = max(tlist[1]) + buffer            
    ncols = (max(tlist[0]) + buffer - Xmin) / cellsizeX
    nrows = (Ymax - Ymin) / cellsizeY
    bbox_lower_left = arcpy.Point(Xmin, Ymin) # must be in raster map units
    bbox_upper_left = [Xmin, Ymax, cellsizeX, cellsizeY]
    nodata_to_value = -9999 / con_z_to_m
    
    del(tlist)
    
    # Construct the array. Note returned array is (row, col) so (y, x)
    Zarry = arcpy.RasterToNumPyArray(EleRaster, bbox_lower_left, ncols, nrows, nodata_to_value)
    
    # convert array values to meters if needed
    Zarry = Zarry * con_z_to_m             
    
    # iterate through each node and azimuth
    for nodeID in NodeDict:
        for a in azimuths:
            CoordList2 = [row for row in CoordList if row[2] == nodeID and row[3] == a]
            TopoZList = []
    
            for i in range(0,len(CoordList2)):
                xy = CoordToArray(CoordList2[i][0], CoordList2[i][1], bbox_upper_left)
                TopoZList.append(Zarry[xy[1], xy[0]])
            
            TopoZarry = np.array(TopoZList)
            Disarry = azimuthdisdict[a] * con_to_m
            Shadearry = np.degrees(np.arctan((TopoZarry - TopoZarry[0]) / Disarry))
            # Take out the off raster samples
            naindex = np.where(TopoZarry < -9998)
            for x in naindex[0]: Shadearry[x] = -9999            
            # Find the max shade angle
            ShadeAngle = Shadearry.max()
            # array index at the max shade angle 
            arryindex = np.where(Shadearry==ShadeAngle)[0][0]
            ShadeZ = TopoZarry[arryindex]
            ZChange = ShadeZ - TopoZarry[0]
            ShadeDistance = azimuthdisdict[a][arryindex]
            SearchDistance = azimuthdisdict[a].max()
            ShadeAngle_X = (ShadeDistance * con_from_m * sin(radians(a))) + CoordList2[0][0]
            ShadeAngle_Y = (ShadeDistance * con_from_m * cos(radians(a))) + CoordList2[0][1]
            offRasterSamples = (TopoZarry > -9998).sum()
        
            TopoList.append([ShadeAngle_X, ShadeAngle_Y, ShadeAngle_X, ShadeAngle_Y, streamID, nodeID, a, ShadeAngle, ShadeZ, TopoZarry[0], ZChange, ShadeDistance, SearchDistance, offRasterSamples])

    return TopoList

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()    

    # Determine input point spatial units
    proj = arcpy.Describe(NodesFC).spatialReference
    proj_ele = arcpy.Describe(EleRaster).spatialReference

    # Check to make sure the raster and input points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(NodesFC,EleRaster))
        sys.exit("Input points and elevation raster do not have the same projection. Please reproject your data.")

    # Determine the elevation Z units conversion into meters
    if EleUnits == "Feet":
        con_z_to_m = 0.3048
    if EleUnits == "Meters":
        con_z_to_m = 1
    if EleUnits == "Other":
        sys.exit("Please modify your raster elevtion units to feet or meters.")
        
    con_to_m = ToMetersUnitConversion(NodesFC)
    con_from_m = FromMetersUnitConversion(NodesFC)
    MaxSearchDistance = MaxSearchDistance_km * 1000 # in meters

    # Get the elevation raster cell size in meters
    cellsizeX = arcpy.Describe(EleRaster).meanCellWidth * con_to_m
    cellsizeY = arcpy.Describe(EleRaster).meanCellHeight * con_to_m
    
    if Directions == 2: # All directions
        azimuths = [45,90,135,180,225,270,315,365]
    else:        
        azimuths = [270,180,90]

    azimuthdict = {45:"NE",90:"E",135:"SE",180:"S",225:"SW",270:"W",315:"NW",365:"N"}
        
    # Build a numpy array from the minimum to the max search distance by increments of the cellsize
    # using 0.000001 to avoid divide by zero errors
    disX = BuildSearchArray(0.000001, MaxSearchDistance, cellsizeX, skippy)
    disXY = BuildSearchArray(0.000001, MaxSearchDistance, hypot(cellsizeX,cellsizeY), skippy)
    disY = BuildSearchArray(0.000001, MaxSearchDistance, cellsizeY, skippy)
    azimuthdisdict = {45:disXY,90:disX,135:disXY,180:disY,225:disXY,270:disX,315:disXY,365:disY}
    del(disX,disXY,disY)

    # Build the topo field names
    AddFields = ["TOPO_"+ azimuthdict[a] for a in azimuths]
        
    # Read the feature class data into a nested dictionary
    NodeDict = ReadNodesFC(NodesFC, OverwriteData, AddFields)

    # calculate the buffer distance (in EleRaster spatial units) to add to the raster bounding box when extracting to an array
    buffer = cellsizeX * 3
    
    TopoList = []
    n = 1
    for streamID in NodeDict:
        print("Processing stream %s of %s" % (n, len(NodeDict)))
        #TopoList = CreateTopoList(NodeDict[streamID], streamID, a, azimuthdisdict, con_from_m)
        TopoList = TopoList + GetShadeAngles2(NodeDict[streamID], streamID, EleRaster, azimuths, azimuthdisdict, buffer, con_z_to_m)
        n = n + 1       
    
    # Update the NodeDict
    for row in TopoList:
        for key in AddFields:
            NodeDict[row[4]][row[5]][key] = row[7]    
    
    CreateTopoFC(TopoList, outpoint_final, proj)
    UpdateNodesFC(NodeDict, NodesFC, AddFields)

    gc.collect()

    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10   
    print("Process Complete in %s minutes" % (elapsedmin))
    arcpy.AddMessage("Process Complete in %s minutes" % (elapsedmin))


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