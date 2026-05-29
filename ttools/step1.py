"""TTools Step 1: Create Stream Nodes"""

import sys
import os
import gc
import time
import traceback
from datetime import timedelta
from math import ceil, atan2, degrees
from operator import itemgetter

import arcpy

from ttools.utils import to_meters_con, from_meters_con, message


def create_node_list(streamline_fc, sid_field, node_dx, checkDirection, z_raster):
    """Reads an input stream centerline file and returns the NODE ID,
    STREAM ID, and X/Y coordinates as a list"""
    nodeList = []
    incursorFields = ["SHAPE@", "SHAPE@LENGTH", sid_field]
    nodeID = 0

    # Determine input projection and spatial units
    proj = arcpy.Describe(streamline_fc).spatialReference
    con_from_m = from_meters_con(streamline_fc)
    con_to_m = to_meters_con(streamline_fc)

    # Pull the stream IDs into a list
    sid_list = []
    with arcpy.da.SearchCursor(streamline_fc, sid_field, "", proj) as Inrows:
        for row in Inrows:
            sid_list.append(row[0])

    # Check for duplicate stream IDs
    dups = list(set([i for i in sid_list if sid_list.count(i) > 1]))
    if dups:
        raise ValueError("There are duplicate stream IDs in your input stream" +
                 "feature class." +
                 "\nHere are the duplicates:  \n" +
                 "{0}".format(dups))

    # Now create the nodes. I'm pulling the fc data twice because on
    # speed tests it is faster compared to saving all the incursorFields
    # to a list and iterating over the list
    message("Creating Nodes")
    with arcpy.da.SearchCursor(streamline_fc, incursorFields, "", proj) as Inrows:
        for row in Inrows:
            lineLength = row[1]  # These units are in the units of projection
            numNodes = int(lineLength * con_to_m / node_dx)
            nodes = range(0, numNodes + 1)
            mid = range(0, numNodes)

            if checkDirection is True:
                flip = check_stream_direction(row[0], z_raster, row[2])
            else:
                flip = 1
            arcpy.SetProgressor("step", "Creating Nodes", 0, numNodes + 1, 1)
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
                node = row[0].positionAlongLine(abs(flip - position),
                                                True).centroid
                # Get the coordinates at the up/down midway point along
                # the line between nodes and calculate the stream azimuth
                if position == 0.0:
                    mid_up = row[0].positionAlongLine(
                        abs(flip - (position + mid_distance)), True).centroid
                    mid_down = node
                elif 0.0 < position + mid_distance < 1:
                    mid_up = row[0].positionAlongLine(
                        abs(flip - (position + mid_distance)), True).centroid
                    mid_down = row[0].positionAlongLine(
                        abs(flip - (position - mid_distance)), True).centroid
                else:
                    mid_up = node
                    mid_down = row[0].positionAlongLine(
                        abs(flip - (position - mid_distance)), True).centroid

                stream_azimuth = degrees(atan2((mid_down.X - mid_up.X),
                                               (mid_down.Y - mid_up.Y)))
                if stream_azimuth < 0:
                    stream_azimuth = stream_azimuth + 360

                # list of "NODE_ID","STREAM_ID". "STREAM_KM", "LENGTH",
                # "POINT_X","POINT_Y", "ASPECT", "SHAPE@X", "SHAPE@Y"
                nodeList.append([nodeID, row[2],
                                 float(position * lineLength * con_to_m / 1000),
                                 segment_length[i],
                                 node.X, node.Y, stream_azimuth, node.X, node.Y])
                nodeID = nodeID + 1
                i = i + 1

        arcpy.SetProgressorPosition()
    arcpy.ResetProgressor()
    return(nodeList)


def create_nodes_fc(nodeList, nodes_fc, streamline_fc, sid_field, proj):
    """Create the output point feature class using
    the data from the nodes list"""
    message("Exporting Data")

    # Determine Stream ID field properties
    sid_type = arcpy.ListFields(streamline_fc, sid_field)[0].type
    sid_precision = arcpy.ListFields(streamline_fc, sid_field)[0].precision
    sid_scale = arcpy.ListFields(streamline_fc, sid_field)[0].scale
    sid_length = arcpy.ListFields(streamline_fc, sid_field)[0].length

    # Create an empty output with the same projection as the input polyline
    cursorfields = ["NODE_ID",
                    "STREAM_ID",
                    "STREAM_KM",
                    "LENGTH",
                    "LONGITUDE",
                    "LATITUDE",
                    "ASPECT"]
    arcpy.CreateFeatureclass_management(os.path.dirname(nodes_fc),
                                        os.path.basename(nodes_fc),
                                        "POINT", "", "DISABLED", "DISABLED", proj)

    # Add attribute fields
    for f in cursorfields:
        if f == "STREAM_ID":
            arcpy.AddField_management(nodes_fc, f, sid_type, sid_precision,
                                      sid_scale, sid_length, "",
                                      "NULLABLE", "NON_REQUIRED")
        else:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(nodes_fc, cursorfields + ["SHAPE@X", "SHAPE@Y"]) as cursor:
        for row in nodeList:
            cursor.insertRow(row)

    # Change X/Y from input spatial units to decimal degrees
    # This selects the best geographic datum transformation for the
    # input projection and extent. ListTransformations returns
    # all the workable transformations with the best match first.
    # when the nodes and target share same datum (e.g. WGS 84)
    # the list is empty and no transformation needed.
    proj_dd = arcpy.SpatialReference(4326)  # GCS_WGS_1984
    extent = arcpy.Describe(nodes_fc).extent
    transforms = arcpy.ListTransformations(proj, proj_dd, extent)
    # Use the first transformation recommended by arcpy.
    transform_to_use = transforms[0] if transforms else ""

    with arcpy.da.UpdateCursor(nodes_fc, ["SHAPE@", "LONGITUDE", "LATITUDE"]) as cursor:
        for row in cursor:
            pt_dd = row[0].projectAs(proj_dd, transform_to_use)
            row[1] = pt_dd.centroid.X  # LONGITUDE
            row[2] = pt_dd.centroid.Y  # LATITUDE
            cursor.updateRow(row)


def check_stream_direction(stream, z_raster, streamID):
    """Samples the elevation raster at both ends of the stream
    polyline to see which is the downstream end and returns flip = 1
    if the stream km need to be reversed"""

    down = stream.positionAlongLine(0, True).centroid
    up = stream.positionAlongLine(1, True).centroid

    # when a single raster cell is sampled it is a little faster to
    # use arcpy compared to converting to an array and then sampling.
    # The code to use the array method is in an earlier version of the code.
    z_down = float(arcpy.GetCellValue_management(z_raster, str(down.X) + " " + str(down.Y), 1).getOutput(0))
    z_up = float(arcpy.GetCellValue_management(z_raster, str(up.X) + " " + str(up.Y), 1).getOutput(0))

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
        if arcpy.Exists(nodes_fc):
            raise ValueError("This output already exists: \n" +
                     "{0}\n".format(nodes_fc) +
                     "Please rename your output.")

        # Get the spatial projection of the input stream lines
        proj = arcpy.Describe(streamline_fc).spatialReference

        if checkDirection is True:
            if z_raster is None:
                raise ValueError("z_raster must be set when checkDirection is True.")

            proj_ele = arcpy.Describe(z_raster).spatialReference

            # Check to make sure the elevation raster and input
            # streams are in the same projection.
            if proj.name != proj_ele.name:
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
        tbinfo = traceback.format_exc()
        pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
        raise
