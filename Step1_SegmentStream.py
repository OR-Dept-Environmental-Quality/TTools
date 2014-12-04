#######################################################################################
# TTools
# Step 1: Create Stream Nodes  version 0.92
# Ryan Michie

# This script will take an input polyline feature with unique stream IDs and generate 
# evenly spaced points along each unique stream ID line at a user defined spacing 
# measured from the downstream endpoint.

# INPUTS
# 0: input stream centerline polyline (inLine)
# 1: unique StreamID field (StreamIDfield)
# 2: spacing between nodes (node_dx)
# 3: output point file name/path (outpoint_final)

# OUTPUTS
# point feature class

# The output point feature class has the following fields: 
# 0: "NODE_ID" - unique node ID
# 1: "STREAM_ID" - field matching a unique polyline ID field identifed by the user,
# 2: "STREAM_KM" - double measured from the downstream end of the stream for each STREAM ID
# 3: "POINT_X" - X coordinate of the node in units of the input point projection.
# 4: "POINT_Y" - Y coordinate of the node in units of the input point projection.

# Future Updates
# implement a method to flip the direction of the stream km if the stram is digitized in the wrong direction.

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
from math import ceil
from collections import defaultdict
from operator import itemgetter
import arcpy
from arcpy import env

env.overwriteOutput = True

# Parameter fields for python toolbox
#inLine = parameters[0].valueAsText
#StreamIDfield = parameters[1].valueAsText
#node_dx = parameters[2].valueAsText
#outpoint_final = parameters[3].valueAsText

# Start Fill in Data
#inLine = r"D:\Projects\TTools_9\Example_data.gdb\McFee_flowline"
#SIDname = "GNIS_ID"
#node_dx = 50
#outpoint_final = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"

inLine = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_streams_single"
SIDname = "NAME"
node_dx = 50
outpoint_final = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
# End Fill in Data

def CreateNodeList(inLine):
    """Reads an input stream centerline file and returns the NODE ID, STREAM ID, and X/Y coordinates as a list"""
    NodeList = []
    Incursorfields = ["SHAPE@","SHAPE@LENGTH",SIDname]
    NID = 0
    # Determine input projection and spatial units
    proj = arcpy.Describe(inLine).spatialReference
    con_from_m = FromMetersUnitConversion(inLine)
    con_to_m = ToMetersUnitConversion(inLine)
    
    with arcpy.da.SearchCursor(inLine,Incursorfields,"",proj) as Inrows:
        print("Creating Nodes")	
        for row in Inrows:
            LineLength = row[1]  # These units are in the units of projection
            numNodes = int(LineLength * con_to_m / node_dx)
            nodes = range(0,numNodes+1)
            flip = CheckStreamDirection(inLine, EleRaster=None)
            arcpy.SetProgressor("step", "Creating Nodes", 0, numNodes+1, 1)
            positions = [n * node_dx * con_from_m / LineLength for n in nodes] # list of percentage of feature length to traverse
            for position in positions:
                node = row[0].positionAlongLine(abs(flip - position),True).centroid
                # list of "NODE_ID",STREAM_ID,"STREAM_KM","POINT_X","POINT_Y","SHAPE@X","SHAPE@Y"
                NodeList.append((NID, row[2], float(position * LineLength * con_to_m /1000), node.X, node.Y, node.X, node.Y ))
                NID = NID + 1
            arcpy.SetProgressorPosition()
    arcpy.ResetProgressor()
    return(NodeList)

def CreateNodesFC(NodeList, NodesFC, SIDname, proj):
    """Create the output point feature class using the data from the nodes list"""
    #arcpy.AddMessage("Exporting Data")
    print("Exporting Data")

    # Determine Stream ID field properties
    SIDtype = arcpy.ListFields(inLine,SIDname)[0].type
    SIDprecision = arcpy.ListFields(inLine,SIDname)[0].precision
    SIDscale = arcpy.ListFields(inLine,SIDname)[0].scale
    SIDlength = arcpy.ListFields(inLine,SIDname)[0].length    

    #Create an empty output with the same projection as the input polyline
    cursorfields = ["NODE_ID","STREAM_ID","STREAM_KM","LONGITUDE","LATITUDE","SHAPE@X","SHAPE@Y"]
    arcpy.CreateFeatureclass_management(os.path.dirname(NodesFC),os.path.basename(NodesFC), "POINT","","DISABLED","DISABLED",proj)

    # Add attribute fields
    arcpy.AddField_management(NodesFC, "NODE_ID", "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(NodesFC, "STREAM_ID", SIDtype, SIDprecision, SIDscale, SIDlength, "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(NodesFC, "STREAM_KM", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(NodesFC, "LONGITUDE", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(NodesFC, "LATITUDE", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(NodesFC, cursorfields) as cursor:
        for row in NodeList:
            cursor.insertRow(row)

    #Change X/Y from input spatial units to decimal degrees
    proj_dd = arcpy.SpatialReference(4269) #GCS_North_American_1983 
    with arcpy.da.UpdateCursor(NodesFC,["SHAPE@X","SHAPE@Y","LONGITUDE","LATITUDE"],"",proj_dd) as cursor:
        for row in cursor:
            row[2] = row[0] # LONGITUDE
            row[3] = row[1] # LATITUDE
            cursor.updateRow(row)

def CheckStreamDirection(inLine, EleRaster):
    """Samples the elevation raster at both ends of the stream polyline to see which is the downstream end 
    and flips the stream km if needed"""
    
    # TBA 
    # 1= flip, 0 = no flip
    flip = 0
    
    return flip
    

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

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()   

    # Create the stream nodes and return them as a list
    NodeList = CreateNodeList(inLine)

    #sort the list by stream ID and then stream km
    NodeList = sorted(NodeList, key=itemgetter(1,2), reverse=True)

    # Get the spatial projecton of the input stream lines
    proj = arcpy.Describe(inLine).SpatialReference

    # Create the output point feature class with the nodes list
    CreateNodesFC(NodeList, outpoint_final, SIDname, proj)

    gc.collect()

    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / len(NodeList)).microseconds
    print("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))
    #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))	

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