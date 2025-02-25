#!/usr/bin/python

"""
TTools Step 2: Measure Channel Widths

This script will measure the channel width and distance to the closest edge on the right and left banks at each
stream node. Output distances are in meters.

REQUIREMENTS
TTools steps 1 must be run before Step 2.
ESRI ArcPro
Python 3.7+

INPUT VARIABLES
0: nodes_fc:
Path to the TTools point feature class.

1: rb_fc
The right bank feature class. The right bank is on the right looking downstream.

2: lb_fc
The left bank feature class. The left bank is on the left looking downstream.

3: overwrite_data:
True/False flag if existing data in nodes_fc can be overwritten.

OUTPUTS
0. nodes_fc:
New fields listed below are added into nodes_fc.
CHANWIDTH: distance in meters between left and right banks.
LEFT: distance in meters from the stream node to the closest edge of the left bank feature.
RIGHT: distance in meters from the stream ode to the closest edge of the right bank feature.

"""
# Import system modules
import sys
import gc
import time
import traceback
from datetime import timedelta
import arcpy
from arcpy import env
from math import ceil
from collections import defaultdict

# ----------------------------------------------------------------------
# Start input variables
nodes_fc = r"C:\workspace\ttools_tests\TTools_py39\jc_test_py39.gdb\jc_stream_nodes_py39"
rb_fc = r"C:\workspace\ttools_tests\JohnsonCreek.gdb\jc_thalweg"
lb_fc = r"C:\workspace\ttools_tests\JohnsonCreek.gdb\jc_thalweg"
overwrite_data = True
# End input variables
# ----------------------------------------------------------------------

# Future Updates
# eliminate arcpy and use gdal for reading/writing feature class data

def nested_dict():
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    nodeDict = nested_dict()
    incursorFields = ["STREAM_ID", "NODE_ID", "STREAM_KM", "ASPECT", "SHAPE@X", "SHAPE@Y"]

    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)

    # Check to see if the 1st field exists if yes add it.
    if overwrite_data is False and (addFields[0] in existingFields) is True:
        incursorFields.append(addFields[0])
    else:
        overwrite_data = True

    # Determine input point spatial units
    proj_nodes = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc, incursorFields, "", proj_nodes) as Inrows:
        if overwrite_data is True:
            for row in Inrows:
                nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2]
                nodeDict[row[0]][row[1]]["ASPECT"] = row[3]
                nodeDict[row[0]][row[1]]["POINT_X"] = row[4]
                nodeDict[row[0]][row[1]]["POINT_Y"] = row[5]
        else:
            for row in Inrows:
                # if the data is null or zero (0 = default for shapefile),
                # it is retrieved and will be overwritten.
                if row[6] is None or row[6] == 0 or row[6] < -9998:
                    nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2]
                    nodeDict[row[0]][row[1]]["ASPECT"] = row[3]
                    nodeDict[row[0]][row[1]]["POINT_X"] = row[4]
                    nodeDict[row[0]][row[1]]["POINT_Y"] = row[5]
    if len(nodeDict) == 0:
        sys.exit("The fields checked in the input point feature class " +
                 "have existing data. There is nothing to process. Exiting")

    return (nodeDict)

def to_meters_con(inFeature):
    """Returns the conversion factor to get from the
    input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that ".format(inFeature) +
                       "is not projected or not recognized. Use a " +
                       "projected coordinate system preferably in linear " +
                       "units of feet or meters.")
        sys.exit("Coordinate system is not projected or not recognized. " +
                 "Use a projected coordinate system, preferably in " +
                 "linear units of feet or meters.")
    return con_to_m

def read_polyline_geometry(polyline_fc, proj_polyline):
    """Reads an input polyline into an arcpy polyline geometry object"""
    poly_list = []
    # Get the x and y of each vertex in the polyline and save
    # it as a list.
    for row in arcpy.da.SearchCursor(polyline_fc, ["SHAPE@"]):
        for part in row[0]:
            for pnt in part:
                poly_list.append(arcpy.Point(pnt.X, pnt.Y))
    poly_array = arcpy.Array(poly_list)
    # put it into a geometry object.
    poly_geom = arcpy.Polyline(poly_array, proj_polyline)
    del row
    return (poly_geom)

def calc_channel_width(node_geom, bank_geom, aspect, line_dis, proj_nodes):
    pt1 = node_geom.pointFromAngleAndDistance(aspect, line_dis, "PLANAR")
    line = arcpy.Polyline(arcpy.Array([node_geom.centroid, pt1.centroid]), proj_nodes)
    pt2 = line.intersect(bank_geom, 1)

   # allows us to rotate the line without touching the original aspect
    aspect_attempt = aspect

   # helps us choose which distance to return at the end of the function
    aspect_switch = False

    # if there's no intersection enter this loop
    if not pt2.centroid:

        # turn on the rotation switch to return the correct distance
        aspect_switch = True

        # generated line never intersected the right bank or left bank
        while not pt2.centroid:
            # ...first we try decrementing and getting that distance
            aspect_attempt = aspect_attempt - 1
            if aspect_attempt < 0:
                aspect_attempt = aspect_attempt + 360
            pt1 = node_geom.pointFromAngleAndDistance(aspect_attempt, line_dis, "PLANAR")
            line = arcpy.Polyline(arcpy.Array([node_geom.centroid, pt1.centroid]), proj_nodes)
            pt2 = line.intersect(bank_geom, 1)
        to_bank_distance_1 = node_geom.distanceTo(pt2)
        del pt2
        pt1 = node_geom.pointFromAngleAndDistance(aspect, line_dis, "PLANAR")
        line = arcpy.Polyline(arcpy.Array([node_geom.centroid, pt1.centroid]), proj_nodes)
        pt2 = line.intersect(bank_geom, 1)
        aspect_attempt = aspect
        while not pt2.centroid:
            # reset and try incrementing the angle (moving clockwise)
            aspect_attempt = aspect_attempt + 1
            if aspect_attempt < 360:
                aspect_attempt = aspect_attempt - 360
            pt1 = node_geom.pointFromAngleAndDistance(aspect_attempt, line_dis, "PLANAR")
            line = arcpy.Polyline(arcpy.Array([node_geom.centroid, pt1.centroid]), proj_nodes)
            pt2 = line.intersect(bank_geom, 1)
        to_bank_distance_2 = node_geom.distanceTo(pt2)

    if not aspect_switch:
        # if we didn't have to rotate and increment the aspect, just report out the distance
        to_bank_distance = node_geom.distanceTo(pt2)
    else:
        # if we had to rotate, compare the two distances and take the lesser
        to_bank_distance = min(to_bank_distance_1,
                               to_bank_distance_2)

    return (to_bank_distance)

def update_nodes_fc(nodeDict, nodes_fc, addFields):
    """Updates the input point feature class with
    data from the nodes dictionary"""
    print("Updating input point feature class")

    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)

    # Check to see if the field exists and add it if not
    for f in addFields:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.UpdateCursor(nodes_fc, ["STREAM_ID", "NODE_ID"] +
                                         addFields) as cursor:
        for row in cursor:
            for f, field in enumerate(addFields):
                streamID = row[0]
                nodeID = row[1]
                row[f + 2] = nodeDict[streamID][nodeID][field]
                cursor.updateRow(row)

# enable garbage collection
gc.enable()

try:
    print("Step 2: Measure Channel Width")

    # keeping track of time
    startTime = time.time()

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

    # Determine input spatial units
    proj_nodes = arcpy.Describe(nodes_fc).spatialReference
    proj_rb = arcpy.Describe(rb_fc).spatialReference
    proj_lb = arcpy.Describe(lb_fc).spatialReference

    # Check to make sure the rb_fc/lb_fc and input points are
    # in the same projection.
    if proj_nodes.name != proj_rb.name:
        arcpy.AddError("{0} and {1} do not have ".format(nodes_fc, rb_fc) +
                       "the same projection. Please reproject your data.")
        sys.exit("Input points and right bank feature class do not have " +
                 "the same projection. Please reproject your data.")

    if proj_nodes.name != proj_lb.name:
        arcpy.AddError("{0} and {1} do not have ".format(nodes_fc, lb_fc) +
                       "the same projection. Please reproject your data.")
        sys.exit("Input points and left bank feature class do not have " +
                 "the same projection. Please reproject your data.")

    con_to_m = to_meters_con(nodes_fc)

    # distance in the spatial units of nodes_fc to look for the right or left bank
    # 5000 meters should be long enough
    out_dis = 5000 * 1 / con_to_m

    addFields = ["CHANWIDTH", "LEFT", "RIGHT"]

    # Read the feature class data into a nested dictionary
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)

    # Read each of the bank polylines into geometry objects
    rb_geom = read_polyline_geometry(rb_fc, proj_rb)
    lb_geom = read_polyline_geometry(lb_fc, proj_lb)

    for n, streamID in enumerate(nodeDict):
        print("Processing stream {0} of {1}".format(n + 1, len(nodeDict)))

        for nodeID in nodeDict[streamID]:
            node_x = float(nodeDict[streamID][nodeID]["POINT_X"])
            node_y = float(nodeDict[streamID][nodeID]["POINT_Y"])
            aspect = float(nodeDict[streamID][nodeID]["ASPECT"])
            node_geom = arcpy.PointGeometry(arcpy.Point(node_x, node_y), proj_nodes)

            # calculate right and left transect directions in degrees
            dir_lb = aspect - 90
            dir_rb = aspect + 90

            if dir_lb < 0:
                dir_lb = dir_lb + 360

            if dir_rb > 360:
                dir_rb = dir_rb - 360

            lb_distance = calc_channel_width(node_geom, lb_geom, dir_lb, out_dis, proj_nodes)
            rb_distance = calc_channel_width(node_geom, rb_geom, dir_rb, out_dis, proj_nodes)

            nodeDict[streamID][nodeID]["CHANWIDTH"] = (lb_distance + rb_distance) * con_to_m
            nodeDict[streamID][nodeID]["LEFT"] = lb_distance * con_to_m
            nodeDict[streamID][nodeID]["RIGHT"] = rb_distance * con_to_m
            
    update_nodes_fc(nodeDict, nodes_fc, addFields)
    
    gc.collect()

    endTime = time.time()
    elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
    mspernode = timedelta(seconds=(endTime - startTime) / (n + 1)).microseconds
    print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))
    # arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))

# For arctool errors
except arcpy.ExecuteError:
    msgs = arcpy.GetMessages(2)
    # arcpy.AddError(msgs)
    print(msgs)

# For other errors
except:
    tbinfo = traceback.format_exc()

    pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
    msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

    # arcpy.AddError(pymsg)
    # arcpy.AddError(msgs)

    print(pymsg)
    print(msgs)