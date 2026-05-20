"""TTools Step 2: Measure Channel Widths"""

import sys
import gc
import time
import traceback
from datetime import timedelta
from math import ceil

import arcpy

from ttools.utils import to_meters_con, message, warning, read_fc, update_fc


def read_polyline_geometry(polyline_fc, proj_polyline):
    """Reads an input polyline into an arcpy polyline geometry object."""
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
    """Calculate the distance from a node to a bank along a
    perpendicular transect."""
    pt1 = node_geom.pointFromAngleAndDistance(aspect, line_dis, "PLANAR")
    line = arcpy.Polyline(arcpy.Array([node_geom.centroid, pt1.centroid]), proj_nodes)
    pt2 = line.intersect(bank_geom, 1)

    if pt2.centroid:
        to_bank_distance = node_geom.distanceTo(pt2)
    else:
        # No intersection = no bank 90 deg from aspect
        # Find the minimum distance to bank, regardless of the angle
        near_distance = bank_geom.queryPointAndDistance(node_geom.centroid)
        to_bank_distance = near_distance[2]

    return (to_bank_distance)


def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the NODE_ID,
    STREAM_ID, and X/Y coordinates as a dictionary"""

    nodeDict = read_fc(nodes_fc)

    # Check to see if the 1st field exists if yes add it.
    check_field = addFields[0]
    first_node = next(iter(nodeDict.values()))

    if not overwrite_data and check_field in first_node:
        # if the data is null or zero (0 = default for shapefile),
        # it is retrieved and will be overwritten.
        filtered = {}
        for nodeID in nodeDict:
            val = nodeDict[nodeID].get(check_field)
            if val is None or val == 0 or val < -9998:
                filtered[nodeID] = nodeDict[nodeID]
        nodeDict = filtered

    if len(nodeDict) == 0:
        raise ValueError("The fields checked in the input point feature class " +
                 "have existing data. There is nothing to process. Exiting")

    return nodeDict


def step2(nodes_fc, rb_fc, lb_fc, overwrite_data=True):
    """TTools Step 2: Measure Channel Widths

    This script will measure the channel width and distance to the right and left banks 90 degrees perpendicular
    to the stream aspect at each stream node. If there isn't a channel bank polyline 90 degrees perpendicular to the
    stream aspect, the script will return the distance to the nearest edge. Output distances are in meters.

    TTools step 1 must be run before Step 2.

    Parameters:
        nodes_fc (str): Path to the TTools point feature class.
        rb_fc (str): The right bank feature class. The right bank is on the right looking downstream.
        lb_fc (str): The left bank feature class. The left bank is on the left looking downstream.
        overwrite_data (bool): True/False flag if existing data in nodes_fc can be overwritten.

    Outputs:
        nodes_fc: New fields listed below are added into nodes_fc.
            CHANWIDTH - Distance in meters between left and right banks.
            LEFT - Distance in meters from the stream node to the closest edge of the left bank
                feature 90 degrees perpendicular to the stream aspect.
            RIGHT - Distance in meters from the stream node to the closest edge of the right bank
                feature 90 degrees perpendicular to the stream aspect.
    """

    # enable garbage collection
    gc.enable()

    try:
        message("Step 2: Measure Channel Width")

        # keeping track of time
        startTime = time.time()

        # Check if the output exists
        if not arcpy.Exists(nodes_fc):
            raise ValueError("This output does not exist: \n" +
                     "{0}\n".format(nodes_fc))

        # Determine input spatial units
        proj_nodes = arcpy.Describe(nodes_fc).spatialReference
        proj_rb = arcpy.Describe(rb_fc).spatialReference
        proj_lb = arcpy.Describe(lb_fc).spatialReference

        # Check to make sure the rb_fc/lb_fc and input points are
        # in the same projection.
        if proj_nodes.name != proj_rb.name:
            raise ValueError("Input points and right bank feature class do not have " +
                     "the same projection. Please reproject your data.")

        if proj_nodes.name != proj_lb.name:
            raise ValueError("Input points and left bank feature class do not have " +
                     "the same projection. Please reproject your data.")

        con_to_m = to_meters_con(nodes_fc)

        # distance in the spatial units of nodes_fc to look for the right or left bank
        # 5000 meters should be long enough
        out_dis = 5000 * 1 / con_to_m

        addFields = ["CHANWIDTH", "LEFT", "RIGHT"]

        # Read the feature class data into a dictionary
        nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)

        # Read each of the bank polylines into geometry objects
        rb_geom = read_polyline_geometry(rb_fc, proj_rb)
        lb_geom = read_polyline_geometry(lb_fc, proj_lb)

        nodes = list(nodeDict.keys())
        nodes.sort()

        ten_percent = max(1, len(nodes) // 10)
        for n, nodeID in enumerate(nodes):
            if n > 0 and n % ten_percent == 0:
                message("Processed {0}% of nodes.".format(int(n / len(nodes) * 100)))
            node_x = float(nodeDict[nodeID]["POINT_X"])
            node_y = float(nodeDict[nodeID]["POINT_Y"])
            aspect = float(nodeDict[nodeID]["ASPECT"])
            node_geom = arcpy.PointGeometry(arcpy.Point(node_x, node_y), proj_nodes)

            # calculate right and left transect directions in degrees
            dir_lb = aspect - 90
            dir_rb = aspect + 90

            if dir_lb < 0:
                dir_lb = dir_lb + 360

            if dir_rb > 360:
                dir_rb = dir_rb - 360

            lb_distance = calc_channel_width(node_geom, lb_geom, dir_lb,
                                             out_dis, proj_nodes)
            rb_distance = calc_channel_width(node_geom, rb_geom, dir_rb,
                                             out_dis, proj_nodes)

            if lb_distance is None:
                warning(f"Warning: There is not a left bank polyline 90 degrees perpendicular to the aspect at Node ID {nodeID}.")
                lb_distance = 0.0

            if rb_distance is None:
                warning(f"Warning: There is not a right bank polyline 90 degrees perpendicular to the aspect at Node ID {nodeID}.")
                rb_distance = 0.0

            nodeDict[nodeID]["CHANWIDTH"] = (lb_distance + rb_distance) * con_to_m
            nodeDict[nodeID]["LEFT"] = lb_distance * con_to_m
            nodeDict[nodeID]["RIGHT"] = rb_distance * con_to_m

        message("Updating input nodes feature class")
        update_fc(nodeDict, nodes_fc, addFields)

        gc.collect()

        endTime = time.time()
        elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
        mspernode = timedelta(seconds=(endTime - startTime) / len(nodes)).microseconds
        message("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))

    except Exception:
        tbinfo = traceback.format_exc()
        pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
        raise
