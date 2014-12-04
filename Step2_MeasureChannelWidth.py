
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
from math import ceil
from collections import defaultdict
#######################################################################################
# TTools
# Step 2: Measure Channel Widths / Aspects - v 0.1
# Ryan Michie

# INPUTS
# 0: Input TTools point feature class (NodesFC)
# 1: input Right Bank feature class(rb)
# 2: input Left Bank feature class(lb)
# 3: input flag if existing data can be over written (OverwriteData) 1. True, 2. False

# OUTPUTS
# point feature class (edit NodesFC) - Added fields with ASPECT and CHANWIDTH for each node

# Future Updates

# This version is for manual starts from within python.
# This script requires Python 2.6, comptypes package, and ArcGIS 10.1 or higher to run.

#######################################################################################

#Check out the ArcGIS Spatial Analyst extension license
#arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#NodesFC = parameters[0].valueAsText
#rb = parameters[1].valueAsText
#lb = parameters[2].valueAsText
#OverwriteData = parameters[3].valueAsText

# Start Fill in Data
NodesFC = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
rb = r"D:\Projects\TTools_9\Example_data.gdb\McFee_RB"
lb = r"D:\Projects\TTools_9\Example_data.gdb\McFee_LB"
OverwriteData = True
# End Fill in Data

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def ReadNodesFC(NodesFC, OverwriteData, AddFields):
    """Reads the input point feature class and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]

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

def ToMetersUnitConversion(inFeature):
    """Returns the conversion factor to get from the input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that is not projected or not recognized. Use a projected coordinate system preferably in linear units of feet or meters.".format(inFeature))
        sys.exit("Coordinate system is not projected or not recognized. Use a projected coordinate system, preferably in linear units of feet or meters.")   
    return con_to_m

def ReadPolylineGeometry(polylineFC, proj):
    """Reads an input polyline into a polyline geometry object"""
    poly_list = []
    # Get the x and y of each vertex in the polyline and save it as a list.
    for row in arcpy.da.SearchCursor(polylineFC, ["SHAPE@"]):
        for part in row[0]:
            for pnt in part:
                poly_list.append(arcpy.Point(pnt.X, pnt.Y))
    poly_array = arcpy.Array(poly_list)
    # put it into a geometry object.
    poly_geom = arcpy.Polyline(poly_array, proj)
    return(poly_geom)
        
def CalcChannelWidthAdvancedLic(nodexy, rb, lb):
    rb_near_result = arcpy.analysis.GenerateNearTable(
        nodexy, rb, r'in_memory\neartable', '', 'LOCATION','NO_ANGLE', 'CLOSEST')
    
    with arcpy.da.SearchCursor(rb_near_result, ['NEAR_DIST']) as rows:
            row = rows.next()
            rb_distance = row[0]    
    
    lb_near_result = arcpy.analysis.GenerateNearTable(
        nodexy, lb, r'in_memory\neartable', '', 'LOCATION','NO_ANGLE', 'CLOSEST')
    
    with arcpy.da.SearchCursor(lb_near_result, ['NEAR_DIST']) as rows:
                row = rows.next()
                lb_distance = row[0]    
    return(rb_distance + lb_distance)
    

def CalcChannelWidth(node_geom, rb_geom, lb_geom):
    
    rb_distance = node_geom.distanceTo(rb_geom)
    lb_distance = node_geom.distanceTo(lb_geom)
    # http://gis.stackexchange.com/questions/75605/converting-a-distance-to-a-map-distance
    
    return(rb_distance + lb_distance)

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
    

#enable garbage collection
gc.enable()

try:
    print("Step 2: Measure Channel Width") 
    
    #keeping track of time
    startTime= time.time()    

    # Determine input spatial units
    proj = arcpy.Describe(NodesFC).spatialReference
    proj_rb = arcpy.Describe(rb).spatialReference
    proj_lb = arcpy.Describe(lb).spatialReference
    
    # Check to make sure the rb/lb and input points are in the same projection.
    if proj.name != proj_rb.name:
        arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(NodesFC, rb))
        sys.exit("Input points and right bank feature class do not have the same projection. Please reproject your data.")
    
    if proj.name != proj_lb.name:
            arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(NodesFC, lb))
            sys.exit("Input points and left bank feature class do not have the same projection. Please reproject your data.")     
    
    nodexy_to_m = ToMetersUnitConversion(NodesFC)
    
    AddFields = ["CHANWIDTH"]

    # Read the feature class data into a nested dictionary
    NodeDict = ReadNodesFC(NodesFC, OverwriteData, AddFields)
    
    # Read each of the bank polylines into geometry objects
    rb_geom = ReadPolylineGeometry(rb, proj_rb)
    lb_geom = ReadPolylineGeometry(lb, proj_lb)
    
    n = 1
    for streamID in NodeDict:
        print("Processing stream %s of %s" % (n, len(NodeDict)))
        i =0 
        for nodeID in NodeDict[streamID]:
            
            node_x = float(NodeDict[streamID][nodeID]["POINT_X"])
            node_y = float(NodeDict[streamID][nodeID]["POINT_Y"])
            node_geom = arcpy.PointGeometry(arcpy.Point(node_x, node_y), proj)
            
            #NodeDict[streamID][nodeID]["ASPECT"] = CalcAspect(NodeDict[streamID], streamID)
            NodeDict[streamID][nodeID]["CHANWIDTH"] = CalcChannelWidth(node_geom, rb_geom, lb_geom) * nodexy_to_m
        n = n + 1         
        
    UpdateNodesFC(NodeDict, NodesFC, AddFields)

    gc.collect()

    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / n).microseconds
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