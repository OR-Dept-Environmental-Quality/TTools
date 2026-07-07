"""TTools Step 1: Create Stream Nodes"""

import sys
import os
import gc
import time
from datetime import timedelta
from math import ceil, atan2, degrees
from operator import itemgetter

import numpy as np

from ttools.utils import to_meters_con, from_meters_con, message
from ttools.geo_package import (write_fc, get_crs, crs_equal, sample_raster_at_point,
                            fc_exists, read_streamline_features,
                            position_along_line, transform_to_latlong)


def create_node_list(streamline_fc, sid_field, node_dx, checkDirection, z_raster):
    """Reads an input stream centerline file and returns the NODE ID,
    STREAM ID, and X/Y coordinates as a list"""
    nodeList = []
    nodeID = 0

    # Determine input projection and spatial units
    con_from_m = from_meters_con(streamline_fc)
    con_to_m = to_meters_con(streamline_fc)

    # Read the stream centerline feature class
    features = read_streamline_features(streamline_fc, sid_field)

    # Pull the stream IDs into a list
    sid_list = [f[0] for f in features]

    # Check for duplicate stream IDs
    dups = list(set([i for i in sid_list if sid_list.count(i) > 1]))
    if dups:
        raise ValueError("There are duplicate stream IDs in your input stream" +
                 "feature class." +
                 "\nHere are the duplicates:  \n" +
                 "{0}".format(dups))

    # Now create the nodes
    message("Creating Nodes")
    for streamID, geom, lineLength in features:
        numNodes = int(lineLength * con_to_m / node_dx)
        nodes = range(0, numNodes + 1)
        mid = range(0, numNodes)

        if checkDirection is True:
            flip = check_stream_direction(geom, z_raster, streamID)
        else:
            flip = 1

        # list of percentage of feature length to traverse
        positions = [n * node_dx * con_from_m / lineLength for n in nodes]
        segment_length = [node_dx] * numNodes + [lineLength * con_to_m % node_dx]
        mid_distance = node_dx * con_from_m / lineLength
        if mid_distance > 1:
            # this situation occurs when the stream < node_dx.
            # The azimuth is calculated for the entire stream line.
            mid_distance = 1

        i = 0
        for position in positions:
            node_x, node_y = position_along_line(
                geom, abs(flip - position))

            # Get the coordinates at the up/down midway point along
            # the line between nodes and calculate the stream azimuth
            if position == 0.0:
                mid_up_x, mid_up_y = position_along_line(
                    geom, abs(flip - (position + mid_distance)))
                mid_down_x, mid_down_y = node_x, node_y
            elif 0.0 < position + mid_distance < 1:
                mid_up_x, mid_up_y = position_along_line(
                    geom, abs(flip - (position + mid_distance)))
                mid_down_x, mid_down_y = position_along_line(
                    geom, abs(flip - (position - mid_distance)))
            else:
                mid_up_x, mid_up_y = node_x, node_y
                mid_down_x, mid_down_y = position_along_line(
                    geom, abs(flip - (position - mid_distance)))

            stream_azimuth = degrees(atan2((mid_down_x - mid_up_x),
                                           (mid_down_y - mid_up_y)))
            if stream_azimuth < 0:
                stream_azimuth = stream_azimuth + 360

            # list of "NODE_ID","STREAM_ID". "STREAM_KM", "LENGTH",
            # "POINT_X","POINT_Y", "ASPECT"
            nodeList.append([nodeID, streamID,
                             float(position * lineLength * con_to_m / 1000),
                             segment_length[i],
                             node_x, node_y, stream_azimuth])
            nodeID = nodeID + 1
            i = i + 1

    return(nodeList)


def create_nodes_fc(nodeList, nodes_fc, streamline_fc, sid_field, proj):
    """Create the output point feature class using
    the data from the nodes list"""
    message("Exporting Data")

    addFields = ["STREAM_ID", "STREAM_KM", "LENGTH",
                 "LONGITUDE", "LATITUDE", "ASPECT"]

    # Build nodeDict from the node list
    nodeDict = {}
    for row in nodeList:
        nodeID = row[0]
        nodeDict[nodeID] = {
            "STREAM_ID": row[1],
            "STREAM_KM": row[2],
            "LENGTH": row[3],
            "ASPECT": row[6],
            "POINT_X": row[4],
            "POINT_Y": row[5],
        }

    # Change X/Y from input spatial units to decimal degrees
    x_list = [nodeDict[nid]["POINT_X"] for nid in nodeDict]
    y_list = [nodeDict[nid]["POINT_Y"] for nid in nodeDict]
    lons, lats = transform_to_latlong(x_list, y_list, proj)

    for i, nodeID in enumerate(nodeDict):
        nodeDict[nodeID]["LONGITUDE"] = lons[i]
        nodeDict[nodeID]["LATITUDE"] = lats[i]

    write_fc(nodeDict, nodes_fc, addFields, proj)


def check_stream_direction(stream, z_raster, streamID):
    """Samples the elevation raster at both ends of the stream
    polyline to see which is the downstream end and returns flip = 1
    if the stream km need to be reversed"""

    down_x, down_y = position_along_line(stream, 0)
    up_x, up_y = position_along_line(stream, 1)

    z_down_val = sample_raster_at_point(z_raster, down_x, down_y)
    z_up_val = sample_raster_at_point(z_raster, up_x, up_y)

    z_down = float(z_down_val) if z_down_val is not None else -9999
    z_up = float(z_up_val) if z_up_val is not None else -9999

    if z_down <= z_up or z_down == -9999 or z_up == -9999:
        # do not reverse stream km
        flip = 0
    else:
        message("Reversing {0}".format(streamID))
        # reversed stream km
        flip = 1

    return flip


def step1(streamline_fc, sid_field, node_dx, cont_stream_km,
          nodes_fc, checkDirection=False, z_raster=None):
    """TTools Step 1: Create Stream Nodes

    This script will take an input polyline feature with unique stream IDs and generate evenly spaced points along each
    unique stream ID polyline at a user defined spacing measured from the downstream endpoint. The script can also check
    the digitized direction to determine the downstream end.

    Parameters:
        streamline_fc (str): Path to the stream centerline polyline feature class.
        sid_field (str): Name of the attribute field in streamline_fc holding the unique
            stream identifier such as the stream name or ID number.
        node_dx (float): Spacing between nodes in meters.
        cont_stream_km (bool): True/False flag to indicate that a continuous stream km
            should be used for all nodes regardless of the unique values in the sid_field.
        nodes_fc (str): Path and name of the output node feature class.
        checkDirection (bool): True/False flag to check if the stream was digitized in
            correct direction. If checkDirection = True, z_raster must be set.
        z_raster (str, optional): Path and name of the ground elevation raster.
            Ignored if checkDirection = False.

    Outputs:
        nodes_fc: New point feature class with the following fields:
            NODE_ID - Unique node ID.
            STREAM_ID - Field matching a unique stream identifier from the sid_field.
            STREAM_KM - Double measured from the downstream end of the stream for each STREAM ID.
            LONGITUDE - Decimal degrees X coordinate of the node using GCS_WGS_1984 datum.
            LATITUDE - Decimal degrees Y coordinate of the node using GCS_WGS_1984 datum.
            ASPECT - Stream aspect in the direction of flow.
    """

    # enable garbage collection
    gc.enable()

    try:
        # keeping track of time
        startTime = time.time()

        # Check if the output exists
        if fc_exists(nodes_fc):
            raise ValueError("This output already exists: \n" +
                     "{0}\n".format(nodes_fc) +
                     "Please rename your output.")

        # Get the spatial projection of the input stream lines
        proj = get_crs(streamline_fc)

        if checkDirection is True:
            if z_raster is None:
                raise ValueError("z_raster must be set when checkDirection is True.")

            # Check to make sure the elevation raster and input
            # streams are in the same projection.
            if not crs_equal(streamline_fc, z_raster):
                raise ValueError("Input stream line and elevation raster do not have "
                         "the same projection. Please reproject your data.")

        # Create the stream nodes and return them as a list
        nodeList = create_node_list(streamline_fc, sid_field, node_dx,
                                    checkDirection, z_raster)

        if cont_stream_km:
            # sort the list by stream ID and stream km
            nodeList = sorted(nodeList, key=itemgetter(1, 2))

            skm = 0.0
            for i in range(0, len(nodeList)):
                nodeList[i][2] = skm
                skm = skm + (node_dx * 0.001)

            # re sort the list by stream km (with downstream end at the top)
            nodeList = sorted(nodeList, key=itemgetter(2), reverse=True)

        else:
            # sort the list by stream ID and then stream km
            # (downstream end at the top)
            nodeList = sorted(nodeList, key=itemgetter(1, 2), reverse=True)

        # Create the output node feature class with the nodes list
        create_nodes_fc(nodeList, nodes_fc, streamline_fc, sid_field, proj)

        gc.collect()

        endTime = time.time()
        elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
        mspernode = timedelta(seconds=(endTime - startTime) / len(nodeList)).microseconds
        message("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))

    except Exception:
        import traceback
        tbinfo = traceback.format_exc()
        pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
        raise
