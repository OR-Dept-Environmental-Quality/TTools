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
# 2: spacing between nodes (NodeDistance)
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
# make many of these tasks into methods.

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from arcpy import env
from math import radians, sin, cos
from collections import defaultdict
from operator import itemgetter


env.overwriteOutput = True

# Parameter fields for python toolbox
#inLine = parameters[0].valueAsText
#StreamIDfield = parameters[1].valueAsText
#NodeDistance = parameters[2].valueAsText
#outpoint_final = parameters[3].valueAsText

# Start Fill in Data
inLine = r"D:\Projects\TTools_9\Example_data.gdb\McFee_flowline"
SIDname = "GNIS_ID"
NodeDistance = 50
outpoint_final = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
# End Fill in Data

def CreateNodes(inLine):
    """Reads an input stream centerline file and returns the NODE ID, STREAM ID, and X/Y coordinates as a list"""
    NODES = []
    Incursorfields = ["SHAPE@",SIDname]
    NID = 0
    # Determine input point spatial units
    proj = arcpy.Describe(inLine).spatialReference
    with arcpy.da.SearchCursor(inLine,Incursorfields,"",proj) as Inrows:
        print("Creating Nodes")	
        for row in Inrows:
            LineLength = row[0].getLength("PRESERVE_SHAPE")
            numNodes = int(LineLength / NodeDistance)
            nodes = range(0,numNodes+1)
            arcpy.SetProgressor("step", "Creating Nodes", 0, numNodes+1, 1)
            positions = [n * NodeDistance / LineLength for n in nodes] # list of Lengths in meters
            for position in positions:
                node = row[0].positionAlongLine(position,True).centroid
                # list of "NODE_ID",STREAM_ID,"STREAM_KM","POINT_X","POINT_Y","SHAPE@X","SHAPE@Y"
                NODES.append((NID, row[1], float(position * LineLength /1000), node.X, node.Y, node.X, node.Y ))
                NID = NID + 1
            arcpy.SetProgressorPosition()
    arcpy.ResetProgressor()
    return(NODES)

def CreatePointFile(pointList,pointfile, SIDname, proj):
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
    arcpy.CreateFeatureclass_management(os.path.dirname(pointfile),os.path.basename(pointfile), "POINT","","DISABLED","DISABLED",proj)

    # Add attribute fields
    arcpy.AddField_management(pointfile, "NODE_ID", "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(pointfile, "STREAM_ID", SIDtype, SIDprecision, SIDscale, SIDlength, "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(pointfile, "STREAM_KM", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(pointfile, "LONGITUDE", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
    arcpy.AddField_management(pointfile, "LATITUDE", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(pointfile, cursorfields) as cursor:
        for row in pointList:
            cursor.insertRow(row)

    #Change X/Y from input spatial units to decimal degrees
    proj_dd = arcpy.SpatialReference(4269) #GCS_North_American_1983 
    with arcpy.da.UpdateCursor(pointfile,["SHAPE@X","SHAPE@Y","LONGITUDE","LATITUDE"],"",proj_dd) as cursor:
        for row in cursor:
            row[2] = row[0] # LONGITUDE
            row[3] = row[1] # LATITUDE
            cursor.updateRow(row)

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()   

    # Create the stream nodes and return them as a list
    NODES = CreateNodes(inLine)

    #sort the list by stream ID and then stream km
    NODES = sorted(NODES, key=itemgetter(1,2), reverse=True)

    # Get the spatial projecton of the input stream lines
    proj = arcpy.Describe(inLine).SpatialReference

    # Create the output point feature class with the nodes list
    CreatePointFile(NODES,outpoint_final, SIDname, proj)

    gc.collect()

    endTime = time.time()
    elapsedmin= (endTime - startTime) / 60	
    print("Process Complete in %s minutes" % (elapsedmin))
    #arcpy.AddMessage("Process Complete at %s, %s minutes" % (endTime, elapsedmin))	

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