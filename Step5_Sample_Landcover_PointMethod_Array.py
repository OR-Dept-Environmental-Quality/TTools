#######################################################################################
# TTools
# Step 5: Sample Landcover - Star Pattern, Point Method v 0.98
# Ryan Michie

# Sample_Landcover_PointMethod will take an input point feature (from Step 1) and sample input landcover rasters
# in a user specificed number of cardinal directions with point samples spaced at a user defined distance
# moving away from the stream.

# General steps include:
# 1. read nodes feature class
# 2. build a lcsample output feature
# 3. calculate all sample easting/northing coordinates, iterate by stream
# 4. calculate the bounding box of the coordinates + extra
# 5. import elevation/veght rasters and convert to arrays for the stream bounding box
# 6. sample array
# 7. save to dictionary
# 8. output to lc sample feature class
# 9. save to nodes feature class

# INPUTS
# 0: Input TTools point feature class (NodesFC)
# 1: input number of transects per node (trans_count)
# 2: input number of samples per transect (transsample_count)
# 3: include stream sample in transect count (True/False)
# 4: input The distance between transect samples (transsample_distance)
# 5: input landcover data type. 1."Codes" 2. "CanopyCover", or 3."LAI" (CanopyData)
# 6. input landcover data type z units (HTUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 7: use heatsource 8 methods 1. True 2. False
# 8: input landcover code or height raster (LCRaster)
# 9: input (optional) canopy cover or LAI raster (CanopyRaster)
# 10: input (optional) k coeffcient raster (kRaster)
# 11: input elevation raster (EleRaster)
# 12: input elvation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 13: output sample point file name/path (outpoint_final)
# 14: input flag if existing data can be over written (OverwriteData) 1. True, 2. False

# OUTPUTS
# point feature class (edit NodesFC) - added fields with Landcover and elevation data for each azimuth direction at each node
# point feature class (new) - point at each x/y sample and the sample raster values

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
import arcpy
from arcpy import env
from math import radians, sin, cos
from collections import defaultdict
from operator import itemgetter
import numpy

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
from arcpy.management import *

env.overwriteOutput = True

# Parameter fields for python toolbox
#NodesFC = parameters[0].valueAsText
#trans_count = parameters[1].valueAsText # LONG
#transsample_count = parameters[2].valueAsText # LONG
#StreamSample = parameters[3].valueAsText True/False
#transsample_distance = parameters[4].valueAsText # LONG
#CanopyDataType = parameters[5].valueAsText One of these: 1."CanopyCover", or 2."LAI"
#LCRaster = parameters[6].valueAsText # This is either landcover height or codes
#CanopyRaster = parameters[7].valueAsText # OPTIONAL This is either canopy cover or LAI raster
#kRaster = parameters[8].valueAsText # OPTIONAL The k value raster for LAI
#OHRaster = = parameters[9].valueAsText
#EleRaster = parameters[10].valueAsText
#EleUnits = parameters[11].valueAsText
#outpoint_final = parameters[12].valueAsText
#OverwriteData = parameters[13].valueAsText True/False

# Start Fill in Data
NodesFC = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
trans_count = 8 
transsample_count = 4 # does not include a sample at the stream node
StreamSample = True # include a sample at the stream node (emergent sample)? (True/False)
transsample_distance = 8
LCRaster = r"D:\Projects\TTools_9\Example_data.gdb\veght_lidar_ft" # This is either landcover height or codes
CanopyDataType = "Codes"
heatsource8 = False
CanopyRaster = "" # OPTIONAL This is either canopy cover or a LAI raster
kRaster = "" # OPTIONAL This is the k value for LAI
OHRaster = "" # OPTIONAL This is the overhang raster
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = "Feet"
outpoint_final = r"D:\Projects\TTools_9\Example_data.gdb\LC_samplepoint"
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

    # Check to see if the last field exists if yes add it. 
    # Grabs last field becuase often the first field, emergent, is zero
    if OverwriteData == False and (AddFields[len(AddFields)-1] in ExistingFields) == True:
        Incursorfields.append(AddFields[len(AddFields)-1])
    else:
        OverwriteData = True

    # Determine input point spatial units
    proj = arcpy.Describe(NodesFC).spatialReference

    with arcpy.da.SearchCursor(NodesFC, Incursorfields,"",proj) as Inrows:
        if OverwriteData == True:
            for row in Inrows:
                pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
                pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[5] == None or row[5] == 0:
                    pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                    pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
                    pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
    
    if len(pnt_dict) == 0:
        sys.exit("The fields checked in the input point feature class have existing data. There is nothing to process. Exiting")
              
    return pnt_dict

def CreateLCPointFC(pointList, LCFields, NodesFC, proj):
    """Creates the output landcover sample point feature class using the data from the point list"""
    print("Exporting data to land cover sample feature class")

    arcpy.CreateFeatureclass_management(os.path.dirname(NodesFC),os.path.basename(NodesFC), "POINT","","DISABLED","DISABLED",proj)

    AddFields = ["POINT_X","POINT_Y"] + ["STREAM_ID","NODE_ID","AZIMUTH","TRANSNUM","SAMNUM"] + LCFields

    # Add attribute fields # TODO add dictionary of field types so they aren't all double
    for f in AddFields:
        arcpy.AddField_management(NodesFC, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(NodesFC, ["SHAPE@X","SHAPE@Y"] + AddFields) as cursor:
        for row in pointList:
            cursor.insertRow(row)

def SetupLCDataHeaders(transsample_count, trans_count, CanopyDataType, StreamSample, heatsource8):
    """Generates a list of the landcover data file column header names and data types"""
        
    if CanopyDataType == "Codes":        
        type = ["LC","ELE"]
        typer = [LCRaster, EleRaster]
        
    if CanopyDataType == "LAI":  #Use LAI methods
        type = ["HT","ELE","LAI","k", "OH"]
        typer = [LCRaster, EleRaster, CanopyRaster, kRaster, OHRaster]
        
    if CanopyDataType == "CanopyCover":  #Use Canopy Cover methods
        type = ["HT","ELE","CAN", "OH"]
        typer = [LCRaster, EleRaster, CanopyRaster, OHRaster]
    
    lcdataheaders =[] 
    if heatsource8 == True: # a flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
        dir = ['NE','E','SE','S','SW','W','NW']
    else:        
        dir = ['T' + str(x) for x in range(1, trans_count + 1)]

    zone = range(1,int(transsample_count)+1)
    
    # Concatenate the type, dir, and zone and order in the correct way
    for t in type:
        for d in range(0,len(dir)):
            for z in range(0,len(zone)):
                if StreamSample == True and t !="ELE" and d==0 and z==0:
                    #lcdataheaders.append(t+"_EMERGENT") # add emergent
                    lcdataheaders.append(t+"_T0_S0") # add emergent
                    lcdataheaders.append(t+"_"+dir[d]+"_S"+str(zone[z]))
                else:
                    lcdataheaders.append(t+"_"+dir[d]+"_S"+str(zone[z]))
    
    return lcdataheaders, type, typer

def CoordToArray(easting, northing, origin_a):
    """converts x/y coordinates to col and row of the array"""
    xy = []
    xy.append((easting - origin_a[0])/ origin_a[2])  # col
    xy.append((northing - origin_a[1])/ origin_a[2])  # row
    return xy

def CreateLCPointList(NodeDict, streamID, dir, zone, transsample_distance):
    """This builds a unique long form list of all nodes, transect directions, zone values.
    and xy coordinates for all the landcover samples. Azimuth direction 0 and zone 0 are the emergent samples.
    This list is used to create the output point feature class"""
    
    LCpointlist = []
    
    for nodeID in NodeDict:
        origin_x = NodeDict[nodeID]["POINT_X"]
        origin_y = NodeDict[nodeID]["POINT_Y"]
        # This is the emergent/stream sample
        LCpointlist.append([origin_x, origin_y, origin_x, origin_y, streamID, nodeID, 0, 0, 0])
            
        for d in range(0,len(dir)):
            for z in zone:
                # Calculate the x and y coordinate of the landcover sample location
                _X_ = (z * transsample_distance * con_from_m * sin(radians(dir[d]))) + origin_x
                _Y_ = (z * transsample_distance * con_from_m * cos(radians(dir[d]))) + origin_y

                # Add the all the data to the list
                LCpointlist.append([_X_, _Y_, _X_, _Y_, streamID, nodeID, dir[d], d+1, z])
                    
    return LCpointlist

def RasterToArray(raster, lower_left_corner, nrows, ncols, nodata_to_value):
    """Converts a raster to a numpy array, returns the geotransform data and number of rows and cols"""

    arry = arcpy.RasterToNumPyArray(raster, lower_left_corner, ncols, nrows, nodata_to_value)
    arry_origin = []
    arry_origin.append(arcpy.Describe(raster).extent.XMin)
    arry_origin.append(arcpy.Describe(raster).extent.YMin)
    return arry

def SampleRaster(LCpointlist, raster, buffer, con):
    
    cellsizeX = arcpy.Describe(EleRaster).meanCellWidth
    cellsizeY = arcpy.Describe(EleRaster).meanCellHeight    
    
    # calculate lower left corner and nrows/cols for the bounding box
    # first transpose the list so x and y coordinates are in the same list
    tlist = map(lambda *i: list(i), *LCpointlist)
    
    Xmin = min(tlist[0]) - buffer
    Ymin = min(tlist[1]) - buffer
    ncols = (max(tlist[0]) + buffer - Xmin) / cellsizeX
    nrows = (max(tlist[1]) + buffer - Ymin) / cellsizeY
    lower_left_corner = arcpy.Point(Xmin, Ymin)
    origin_a = [Xmin, Ymin, cellsizeY]
    nodata_to_value = -9999 / con
    
    print("Building Array")
    arry = arcpy.RasterToNumPyArray(raster, lower_left_corner, ncols, nrows, nodata_to_value)
    
    # convert array values to meters if needed
    arry = arry * con
    
    print("Extracting raster values")
    for i in range(0,len(LCpointlist)):
        xy = CoordToArray(LCpointlist[i][0], LCpointlist[i][1], origin_a)
        LCpointlist[i].append(arry[xy[1], xy[0]])
    return LCpointlist

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

def FromMetersUnitConversion(inFeature):
    """Returns the conversion factor to get from meters to the spatial units of the input feature class"""
    try:
        con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that is not projected or not recognized. Use a projected coordinate system preferably in linear units of feet or meters.".format(inFeature))
        sys.exit("Coordinate system is not projected or not recognized. Use a projected coordinate system, preferably in linear units of feet or meters.")   
    return con_from_m

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()
    
    # Determine input spatial units and set conversion factor to get from meters to the input spatial units
    proj = arcpy.Describe(NodesFC).SpatialReference
    con_from_m = FromMetersUnitConversion(NodesFC)

    # Set the converstion factor to get from the input elevation z units to meters
    if EleUnits == "Meters": #International meter
        ele_con = 1 
    if EleUnits == "Feet": #International foot
        ele_con = 0.3048
    if EleUnits == "Other": #Some other units
        sys.exit("Please modify your raster elevation z units so they are either in meters or feet")	

    if heatsource8 == True: # flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
        dir = [45,90,135,180,225,270,315]
    else:        
        dir = [x * 360.0 / trans_count for x in range(1,trans_count+ 1)]

    zone = range(1,int(transsample_count+1))
    
    # TODO This is a future function that may replace the emergent methods.
    # If True there is a regular landcover sample at the stream node
    # for each azimuth direction vs a single emergent sample at the stream node.
    #if StreamSample == "TRUE":
        #zone = range(0,int(transsample_count))
    #else:
        #zone = range(1,int(transsample_count+1))    
    
    AddFields, type, typer = SetupLCDataHeaders(transsample_count, trans_count, CanopyDataType, StreamSample, heatsource8)
    NodeDict = ReadNodesFC(NodesFC, OverwriteData, AddFields)
    
    # calculate the buffer distance to add to the raster bounding box when extracting to an array
    buffer = (transsample_count + 1) * transsample_distance * con_from_m    
    LCpointlist = []
    n = 1 
    for streamID in NodeDict:
        print("Processing stream %s of %s" % (n, len(NodeDict)))
        LCpointlist = LCpointlist + CreateLCPointList(NodeDict[streamID], streamID, dir, zone, transsample_distance)
        for raster in typer:
            if raster == EleRaster:
                con = ele_con
            else:
                con = 1
            LCpointlist = SampleRaster(LCpointlist, raster, buffer, con)
        n = n + 1       
    
    # Update the NodeDict
    for row in LCpointlist:
        for t in range(0,len(type)):
            LCkey = type[t]+'_T'+str(row[7])+'_S'+str(row[8])
            NodeDict[row[4]][row[5]][LCkey] = row[9 + t]

    arcpy.ResetProgressor()		
    gc.collect()

    ####################################################################################################### 
    # Build the output point feature class using the data from the LCPointList
    CreateLCPointFC(LCpointlist, type, outpoint_final, proj)

    # Create the landcover headers to be added to the TTools point feature class
    #AddFields = AddFields + ["NUM_DIR","NUM_ZONES","SAMPLE_DIS"]

    # Write the landcover data to the TTools point feature class
    UpdateNodesFC(NodeDict, NodesFC, AddFields)

    endTime = time.time()
    elapsedmin= (endTime - startTime) / 60	
    print("Process Complete in %s minutes" % (elapsedmin))
    arcpy.AddMessage("Process Complete at %s, %s minutes" % (endTime, elapsedmin))


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