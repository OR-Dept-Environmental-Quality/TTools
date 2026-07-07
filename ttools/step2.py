"""TTools Step 2: Measure Channel Widths"""

import sys
import gc
import time
import traceback
from datetime import timedelta
from math import ceil

from ttools.utils import to_meters_con, message, warning
from ttools.geo_package import (read_fc, get_crs, crs_equal, update_fc, fc_exists,
                            read_polyline_geometry, calc_channel_width)


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
        if not fc_exists(nodes_fc):
            raise ValueError("This output does not exist: \n" +
                     "{0}\n".format(nodes_fc))

        # Determine input spatial units
        proj_nodes = get_crs(nodes_fc)

        # Check to make sure the rb_fc/lb_fc and input points are
        # in the same projection.
        if not crs_equal(nodes_fc, rb_fc):
            raise ValueError("Input points and right bank feature class do not have " +
                     "the same projection. Please reproject your data.")

        if not crs_equal(nodes_fc, lb_fc):
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
        rb_geom = read_polyline_geometry(rb_fc)
        lb_geom = read_polyline_geometry(lb_fc)

        nodes = list(nodeDict.keys())
        nodes.sort()

        ten_percent = max(1, len(nodes) // 10)
        for n, nodeID in enumerate(nodes):
            if n > 0 and n % ten_percent == 0:
                message("Processed {0}% of nodes.".format(int(n / len(nodes) * 100)))
            node_x = float(nodeDict[nodeID]["POINT_X"])
            node_y = float(nodeDict[nodeID]["POINT_Y"])
            aspect = float(nodeDict[nodeID]["ASPECT"])

            # calculate right and left transect directions in degrees
            dir_lb = aspect - 90
            dir_rb = aspect + 90

            if dir_lb < 0:
                dir_lb = dir_lb + 360

            if dir_rb > 360:
                dir_rb = dir_rb - 360

            lb_distance = calc_channel_width(node_x, node_y, lb_geom,
                                             dir_lb, out_dis,
                                             proj_nodes)
            rb_distance = calc_channel_width(node_x, node_y, rb_geom,
                                             dir_rb, out_dis,
                                             proj_nodes)

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
