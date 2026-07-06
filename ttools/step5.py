"""TTools Step 5: Sample Landcover"""
import sys
import os
import gc
import time
import traceback
from datetime import timedelta
from math import radians, sin, cos, ceil
from collections import OrderedDict

import arcpy

from ttools.utils import (from_meters_con, from_z_units_to_meters_con,
                          coord_to_array, message, read_fc, update_fc,
                          get_raster_info, raster_to_array)


# -----------------------------------------------------------------------
# Shared functions (used by both orthogonal and star methods)
# -----------------------------------------------------------------------

def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the STREAM_ID,
    NODE_ID, and X/Y coordinates as a dictionary"""

    nodeDict = read_fc(nodes_fc)

    # Check to see if the last field exists if yes add it.
    # Grabs last field because often the first field, emergent, is zero
    check_field = addFields[len(addFields)-1]
    first_node = next(iter(nodeDict.values()))

    if not overwrite_data and check_field in first_node:
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


def update_lc_point_fc(lc_point_list, type, lc_point_fc,
                       nodes_fc, nodes_in_block, overwrite_data, proj):
    """Creates/updates the output landcover sample point feature
    class using the data from the landcover point list"""

    # Create an empty output with the same projection as the input polyline
    cursorfields = ["POINT_X", "POINT_Y"] + ["STREAM_ID", "NODE_ID", "SAMPLE_ID",
                                            "TRANS_DIR", "TRANSECT",
                                            "SAMPLE", "KEY"] + type

    # Check if the output exists and create if not
    if not arcpy.Exists(lc_point_fc):
        arcpy.CreateFeatureclass_management(os.path.dirname(lc_point_fc),
                                            os.path.basename(lc_point_fc),
                                            "POINT", "", "DISABLED", "DISABLED", proj)

        # Determine Stream ID field properties
        sid_type = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].type
        sid_precision = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].precision
        sid_scale = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].scale
        sid_length = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].length

        typeDict = {"POINT_X": "DOUBLE",
                    "POINT_Y": "DOUBLE",
                    "NODE_ID": "LONG",
                    "SAMPLE_ID": "LONG",
                    "TRANS_DIR": "DOUBLE",
                    "TRANSECT": "SHORT",
                    "SAMPLE": "SHORT",
                    "KEY": "TEXT",}

        for t in type:
            typeDict[t] = "DOUBLE"

        # Add attribute fields
        for f in cursorfields:
            if f == "STREAM_ID":
                arcpy.AddField_management(lc_point_fc, f, sid_type,
                                          sid_precision, sid_scale, sid_length,
                                          "", "NULLABLE", "NON_REQUIRED")
            else:
                arcpy.AddField_management(lc_point_fc, f, typeDict[f], "",
                                          "", "", "", "NULLABLE", "NON_REQUIRED")

    if not overwrite_data:
        # Build a query to retrieve existing rows from the nodes
        # that need updating
        whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_in_block))

        # delete those rows
        with arcpy.da.UpdateCursor(lc_point_fc, ["NODE_ID"], whereclause) as cursor:
            for row in cursor:
                cursor.deleteRow()

    with arcpy.da.InsertCursor(lc_point_fc, ["SHAPE@X", "SHAPE@Y"] +
                               cursorfields) as cursor:
        for row in lc_point_list:
            cursor.insertRow(row)


def setup_lcdata_headers_orthogonal(transsample_count):
    """Generates a list of the landcover data file column header names and data types"""

    type = ["ELE"]

    lcheaders = []
    otherheaders = []

    dirs = ["L", "R"]

    # Concatenate the type, dir, and sample and order in the correct way
    for d, dir in enumerate(dirs):
        for s, sample in enumerate(range(1, int(transsample_count) + 1)):
            if d==0 and s==0:
                lcheaders.append("LC_T0_S0") # add emergent
                lcheaders.append("LC_{0}_S{1}".format(dir, sample))
            else:
                lcheaders.append("LC_{0}_S{1}".format(dir, sample))

    for t in type:
        for d, dir in enumerate(dirs):
            for s, sample in enumerate(range(1, int(transsample_count) + 1)):
                if t !="ELE" and d==0 and s==0:
                    otherheaders.append(t+"T0_S0") # add emergent
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, sample))
                else:
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, sample))

    return lcheaders, otherheaders


def setup_lcdata_headers_star(transsample_count, trans_count):
    """Generates a list of the landcover data file
    column header names and data types"""

    type = ["ELE"]

    lcheaders = []
    otherheaders = []

    dirs = ["T{0}".format(x) for x in range(1, trans_count + 1)]

    zones = range(1,int(transsample_count)+1)

    # Concatenate the type, dir, and zone and order in the correct way

    for d, dir in enumerate(dirs):
        for z, zone in enumerate(zones):
            if d==0 and z==0:
                lcheaders.append("LC_T0_S0") # add emergent
                lcheaders.append("LC_{0}_S{1}".format(dir, zone))
            else:
                lcheaders.append("LC_{0}_S{1}".format(dir, zone))

    for t in type:
        for d, dir in enumerate(dirs):
            for z, zone in enumerate(zones):
                if t !="ELE" and d==0 and z==0:
                    otherheaders.append(t+"_T0_S0") # add emergent
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))
                else:
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))

    return lcheaders, otherheaders


def create_lc_point_list_orthogonal(nodeDict, nodes_in_block, transsample_count,
                                     transsample_distance, start_bank, con_from_m):
    """This builds a unique long form list of information for all the
    landcover samples in the block. This list is used to
    create/update the output feature class."""

    lc_point_list = []
    samplesPerNode = (2 * transsample_count) + 1

    if start_bank:
        # need to reduce sample number by 1 so first sample starts at bank
        bx = 1
    else:
        bx = 0

    for nodeID in nodes_in_block:
        origin_x = nodeDict[nodeID]["POINT_X"]
        origin_y = nodeDict[nodeID]["POINT_Y"]
        aspect = nodeDict[nodeID]["ASPECT"]
        streamID = nodeDict[nodeID]["STREAM_ID"]
        sampleID = nodeID * transsample_count

        # Determine the offset from the node to the start of the first transect.
        # Measured in meters here and converted to map units farther down
        if start_bank:
            offset_l = nodeDict[nodeID]["LEFT"]
            offset_r = nodeDict[nodeID]["RIGHT"]
        else:
            offset_l = 0
            offset_r = 0

        # calculate right and left transect directions in degrees
        dir_left = aspect - 90
        dir_right = aspect + 90

        if dir_left < 0:
             dir_left = dir_left + 360

        if dir_right > 360:
            dir_right = dir_right - 360

        dirs = [dir_left, dir_right]

        # This is the emergent/stream sample
        lc_point_list.append([origin_x, origin_y, origin_x, origin_y,
                              streamID, nodeID, sampleID,
                              0, 0, 0, "T0_S0"])

        for d, dir in enumerate(dirs):

            if d == 0:
                offset = offset_l
                prefix = "L"
            if d == 1:
                offset = offset_r
                prefix = "R"

            for sample in range(1, int(transsample_count + 1)):
                # Calculate the x and y coordinate of the
                # landcover sample location
                node_to_sample_dis = (offset + ((sample - bx) * transsample_distance)) * con_from_m
                pt_x = (node_to_sample_dis * sin(radians(dir))) + origin_x
                pt_y = (node_to_sample_dis * cos(radians(dir))) + origin_y

                key = '{0}_S{1}'.format(prefix, sample)
                sampleID = (nodeID * samplesPerNode) + (d * transsample_count) + sample

                # Add to the list
                lc_point_list.append([pt_x, pt_y, pt_x, pt_y,
                                      streamID, nodeID, sampleID,
                                      dir, d+1, sample, key])

    return lc_point_list


def create_lc_point_list_star(nodeDict, nodes_in_block, dirs, zones,
                               transsample_distance, zone_sample, con_from_m):
    """This builds a unique long form list of information for all the
    landcover samples in the block. This list is used to
    create/update the output feature class."""

    lc_point_list = []
    numDirs = len(dirs)
    numZones = len(zones)
    zonesPerNode = (numDirs * numZones) + 1

    if zone_sample:
        adjust = 0.5
    else:
        adjust = 0.0

    for nodeID in nodes_in_block:
        origin_x = nodeDict[nodeID]["POINT_X"]
        origin_y = nodeDict[nodeID]["POINT_Y"]
        streamID = nodeDict[nodeID]["STREAM_ID"]
        sampleID = nodeID * zonesPerNode

        # This is the emergent/stream sample
        lc_point_list.append([origin_x, origin_y, origin_x, origin_y,
                              streamID, nodeID, sampleID,
                              0, 0, 0, "T0_S0"])

        for d, dir in enumerate(dirs):
            for zone in zones:
                # Calculate the x and y coordinate of the
                # landcover sample location
                pt_x = ((zone - adjust) * transsample_distance * con_from_m *
                        sin(radians(dir))) + origin_x
                pt_y = ((zone - adjust) * transsample_distance * con_from_m *
                        cos(radians(dir))) + origin_y

                key = 'T{0}_S{1}'.format(d+1, zone)
                tran_num = '{:{}{}}'.format(d+1, 0, 3)
                samp_num = '{:{}{}}'.format(zone, 0, 2)
                sampleID = (nodeID * zonesPerNode) + (d * numZones) + zone

                # Add to the list
                lc_point_list.append([pt_x, pt_y, pt_x, pt_y,
                                      streamID, nodeID, sampleID,
                                      dir, d+1, zone, key])

    return lc_point_list


def create_block_list(nodeDict, nodes, block_size, buffer):
    """Returns two lists, one containing the coordinate extent
    for each block that will be iteratively extracted to an array
    and the other containing node IDs within each block extent."""

    message("Calculating block extents")
    x_coord_list = [nodeDict[nodeID]["POINT_X"] for nodeID in nodes]
    y_coord_list = [nodeDict[nodeID]["POINT_Y"] for nodeID in nodes]

    # calculate bounding box extent for samples
    x_min = min(x_coord_list)
    x_max = max(x_coord_list)
    y_min = min(y_coord_list)
    y_max = max(y_coord_list)

    x_width = int(x_max - x_min + 1)
    y_width = int(y_max - y_min + 1)

    block_extents = []
    block_nodes = []

    # Build data blocks
    for x in range(0, x_width, block_size):
        for y in range(0, y_width, block_size):

            # Lower left coordinate of block (in map units)
            block0_x_min = min([x_min + x, x_max])
            block0_y_min = min([y_min + y, y_max])
            # Upper right coordinate of block (in map units)
            block0_x_max = min([block0_x_min + block_size, x_max])
            block0_y_max = min([block0_y_min + block_size, y_max])

            block_x_max = block0_x_max
            block_x_min = block0_x_min
            block_y_max = block0_y_max
            block_y_min = block0_y_min

            nodes_in_block = []
            for nodeID in nodes:
                node_x = nodeDict[nodeID]["POINT_X"]
                node_y = nodeDict[nodeID]["POINT_Y"]
                if (block0_x_min <= node_x <= block0_x_max and
                    block0_y_min <= node_y <= block0_y_max):

                    nodes_in_block.append(nodeID)

            if nodes_in_block:

                # Minimize the size of the block0 by the true
                # extent of the nodes in the block
                node_x_min = min([nodeDict[nodeID]["POINT_X"] for nodeID in nodes_in_block])
                node_y_min = min([nodeDict[nodeID]["POINT_Y"] for nodeID in nodes_in_block])
                node_x_max = max([nodeDict[nodeID]["POINT_X"] for nodeID in nodes_in_block])
                node_y_max = max([nodeDict[nodeID]["POINT_Y"] for nodeID in nodes_in_block])

                if block0_x_min < node_x_min: block_x_min = node_x_min
                if block0_x_max > node_x_max: block_x_max = node_x_max
                if block0_y_min < node_y_min: block_y_min = node_y_min
                if block0_y_max > node_y_max: block_y_max = node_y_max

                # Add the block extent for processing and the buffer distance.
                # Because the node coordinates are used to determine if the node is within the block extent,
                # a buffer distance is added to ensure each block extent will include all samples around nodes
                # located at the edge of the block.

                # order 0 left,      1 bottom,    2 right,     3 top
                block_extents.append((block_x_min - buffer, block_y_min - buffer,
                                      block_x_max + buffer, block_y_max + buffer))
                block_nodes.append(nodes_in_block)

    return block_extents, block_nodes


def sample_raster(block, lc_point_list, raster, con):

    if con is not None:
        nodata_to_value = -9999 / con
    else:
        nodata_to_value = -9999

    raster_info = get_raster_info(raster)
    x_cellsize = raster_info["x_cellsize"]
    y_cellsize = raster_info["y_cellsize"]

    # localize the block extent values and add one cell distance to the size to ensure all cells that need to be
    # sampled are included in the array.
    block_x_min = block[0] - x_cellsize
    block_y_min = block[1] - y_cellsize
    block_x_max = block[2] + x_cellsize
    block_y_max = block[3] + y_cellsize

    # Get the coordinate extent of the input raster
    raster_x_min = raster_info["x_min"]
    raster_y_min = raster_info["y_min"]
    raster_x_max = raster_info["x_max"]
    raster_y_max = raster_info["y_max"]

    # Calculate the block x and y offset from the raster and adjust
    # the block coordinates so they are at the raster cell corners.
    block_x_min_corner = block_x_min - ((block_x_min - raster_x_min) % x_cellsize)
    block_y_min_corner = block_y_min - ((block_y_min - raster_y_min) % y_cellsize)
    block_x_max_corner = block_x_max + ((raster_x_max - block_x_max) % x_cellsize)
    block_y_max_corner = block_y_max + ((raster_y_max - block_y_max) % y_cellsize)

    # Construct the array. Note returned array is (row, col) so (y, x)
    try:
        raster_array, block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize = raster_to_array(
            raster, block_x_min_corner, block_y_min_corner,
            block_x_max_corner, block_y_max_corner, nodata_to_value)
    except:
        tbinfo = traceback.format_exc()
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"
        raise MemoryError(pymsg)

    # convert array values to meters if needed
    if con is not None:
        raster_array = raster_array * con

    lc_point_list_new = []
    if raster_array.max() > -9999:
        # There is at least one pixel of data
        n_rows, n_cols = raster_array.shape
        for point in lc_point_list:
            xy = coord_to_array(point[0], point[1], block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize)
            if 0 <= xy[1] < n_rows and 0 <= xy[0] < n_cols:
                point.append(raster_array[xy[1], xy[0]])
            else:
                # off raster sample
                point.append(-9999)
            lc_point_list_new.append(point)
    else:
        # No data, add -9999
        for point in lc_point_list:
            point.append(-9999)
            lc_point_list_new.append(point)
    return lc_point_list_new


# -----------------------------------------------------------------------
# Main step5 function
# -----------------------------------------------------------------------

def step5(nodes_fc, method, transsample_count, transsample_distance,
          lc_raster, lc_units, z_raster, z_units, lc_point_fc,
          block_size=10, overwrite_data=True,
          start_bank=True,
          trans_count=8, zone_sample=False, heatsource8=False):
    """TTools Step 5: Sample Landcover

    This module provides two sampling methods for landcover data:

    Orthogonal Method:
    The orthogonal sampling method will take an input point feature (from Step 1) and sample an input landcover
    raster along two transects that are orthogonal to the stream aspect (left and right looking downstream). The
    number of sample points along each transect and distance between samples is user defined. The transect can
    start at the stream node or at the right and left banks. The Orthogonal sampling method can be used to develop
    inputs for Heat Source 6, Washington Department of Ecology's Shade model, and the shade file for CE-QUAL-W2.

    Star Pattern Method:
    The star pattern sampling method will take an input point feature (from Step 1) and sample input landcover
    rasters along transects in any number of azimuth directions oriented outward from the stream node. The number
    of transects, sample points along each transect, and the distance between samples is user defined. The measured
    distance along each transect starts at the stream node. The star pattern sampling method can be used to develop
    landcover inputs for Heat Source 7 - 9.

    TTools steps 1 - 3 must be run before Step 5.

    Parameters:
        nodes_fc (str): Path to the TTools point feature class.
        method (str): Sampling method. Either "orthogonal" or "star".
        transsample_count (int): Number of samples per transect. The number DOES NOT include the
            sample at the stream node.
        transsample_distance (float): The distance between transect samples (meters).
        lc_raster (str): Path and name of the land cover code, height, or elevation raster.
        lc_units (str): z units of the lc_raster (aka units of height or elevation). Use "Feet",
            "Meters", or None if the lc_raster values are codes and do not represent elevation or
            height units.
        z_raster (str): Path and name of the ground elevation raster.
        z_units (str): z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is
            not in feet or meters the elevation values must be converted.
        lc_point_fc (str): Path and name of output sample point feature file.
        block_size (int): The x and y size in kilometers for each raster block pulled into an array.
            Start with 10 if you aren't sure and reduce if there is an error. To increase processing
            speed rasters are subdivided iteratively into smaller blocks and pulled into arrays for
            processing. Very large block sizes may use a lot of memory and result in an error.
        overwrite_data (bool): True/False flag if existing data in nodes_fc and lc_point_fc can be
            overwritten.
        start_bank (bool): method="orthogonal" only. True/False flag to indicate if the transect
            should start at the stream bank. If False the transect will start at the stream node.
            Normally set to True.
        trans_count (int): method="star" only. Number of transects per node. The degrees of
            separation between each transect is equal to 360 / trans_count. If trans_count = 8,
            the heading of the first transect is 45 degrees (Northeast), the second transect is
            90 degrees (East), and the eighth is 360 degrees (North). Transects are oriented
            clockwise from north. trans_count is ignored if heatsource8 = True.
        zone_sample (bool): method="star" only. True/False flag to indicate if the sample should
            represent a zone and be centered relative to the transsample_distance as measured from
            the stream node. If this is True the distance from the stream node to the sample
            location along each transect is equal to the transsample_distance * (sample number -
            0.5). This should be True if using heat source 7, heat source 8.0.1 - 8.0.5, or heat
            source 8.0.7 - 8.0.8 when the model is configured to use the zone method. This should
            be False if using heat source 9 or when using heat source 8.0.7 - 8.0.8 and the model
            is configured to use the point method.
        heatsource8 (bool): method="star" only. True/False flag to indicate if the star pattern
            with 7 transects should be used. Heat source version 7 and 8 use 7 transects around
            each stream node in the following directions: Northeast, East, Southeast, South,
            Southwest, West, Northwest.

    Outputs:
        nodes_fc: New fields are added into nodes_fc with the landcover and elevation values for
            each transect sample.
        lc_point_fc: New point feature class created with a point at each x/y sample and the
            sample raster values.
    """

    # enable garbage collection
    gc.enable()

    try:
        if method == "orthogonal":
            msg = "Step 5: Sample Landcover - Orthogonal Method"
        elif method == "star":
            msg = "Step 5: Sample Landcover - Star Pattern"
        else:
            raise ValueError("Invalid method. Use 'orthogonal' or 'star'.")
        message(msg)

        # keeping track of time
        startTime = time.time()

        # Check if the node fc exists
        if not arcpy.Exists(nodes_fc):
            raise ValueError("This output does not exist: \n" +
                     "{0}\n".format(nodes_fc))

        # Check if the lc point fc exists and delete if needed
        if arcpy.Exists(lc_point_fc) and overwrite_data:
            arcpy.Delete_management(lc_point_fc)

        # Determine input spatial units and set unit conversion factors
        proj = arcpy.Describe(nodes_fc).spatialReference
        proj_ele = arcpy.Describe(z_raster).spatialReference
        proj_lc = arcpy.Describe(lc_raster).spatialReference if lc_raster else None

        con_from_m = from_meters_con(nodes_fc)
        con_lc_to_m = from_z_units_to_meters_con(lc_units) if lc_units else None
        con_z_to_m = from_z_units_to_meters_con(z_units)

        # convert block size from km to meters to units of the node fc
        # in the future block size should be estimated based on available memory
        # memorysize = datatypeinbytes*nobands*block_size^2
        # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
        block_size = int(con_from_m * block_size * 1000)

        # Check to make sure the raster and input
        # points are in the same projection.
        if proj.name != proj_ele.name:
            raise ValueError("Input points and elevation raster do not have the " +
                     "same projection. Please reproject your data.")

        if proj_lc is not None and proj_lc.name != proj_ele.name:
            raise ValueError("The landcover and elevation rasters do not have the " +
                     "same projection. Please reproject your data.")

        # Setup the raster dictionary. It is ordered because
        # key list needs to correspond to the order of the attribute fields
        rasterDict = OrderedDict({"LC": lc_raster,
                                  "ELE": z_raster})

        # ---------------------------------------------------------------
        # Method-specific setup
        # ---------------------------------------------------------------
        if method == "orthogonal":
            lcheaders, otherheaders = setup_lcdata_headers_orthogonal(transsample_count)
            addFields = lcheaders + otherheaders

            # read the node data into the node dictionary
            nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)

        elif method == "star":
            # flag indicating the model should use the heat source 8 methods
            # (same as 8 directions but no north)
            if heatsource8:
                dirs = [45,90,135,180,225,270,315]
                trans_count = 7
            else:
                dirs = [x * 360.0 / trans_count for x in range(1, trans_count + 1)]

            zones = range(1, int(transsample_count + 1))

            lcheaders, otherheaders = setup_lcdata_headers_star(transsample_count, trans_count)
            addFields = lcheaders + otherheaders

            # read the node data into the node dictionary
            nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)

        # Get a list of the nodes, sort them
        nodes = list(nodeDict.keys())
        nodes.sort()

        # calculate the buffer distance (in raster spatial units) to add to
        # the base bounding box when extracting to an array. The buffer is
        # equal to the sample distance + 1 + max left/right channel width to make sure the block includes
        # all the landcover samples for each node.
        # Everything measured in meters here and converted to map units

        # Determine the max offset from the node to the start of the first transect sample.
        offset_l_max = [nodeDict[nodeID].get("LEFT", 0) for nodeID in nodeDict]
        offset_r_max = [nodeDict[nodeID].get("RIGHT", 0) for nodeID in nodeDict]
        offset_max = max(offset_l_max + offset_r_max)

        buffer = int(((transsample_count + 1) * transsample_distance + offset_max) * con_from_m)

        # Build the block list
        block_extents, block_nodes = create_block_list(nodeDict, nodes, block_size, buffer)

        # Iterate through each block, calculate sample coordinates,
        # convert raster to array, sample the raster
        total_samples = 0

        for p, block in enumerate(block_extents):
            nodes_in_block = block_nodes[p]
            nodes_in_block.sort()
            message("Processing block {0} of {1}".format(p + 1, len(block_extents)))

            # build the landcover sample list
            if method == "orthogonal":
                lc_point_list = create_lc_point_list_orthogonal(
                    nodeDict, nodes_in_block, transsample_count,
                    transsample_distance, start_bank, con_from_m)
            elif method == "star":
                lc_point_list = create_lc_point_list_star(
                    nodeDict, nodes_in_block, dirs, zones,
                    transsample_distance, zone_sample, con_from_m)

            for t, (type, raster) in enumerate(rasterDict.items()):
                if raster is None:
                    for i in range(0, len(lc_point_list)):
                        lc_point_list[i].append(-9999)
                else:
                    if raster == z_raster:
                        con = con_z_to_m
                    elif raster == lc_raster:
                        con = con_lc_to_m
                    else:
                        con = None

                    lc_point_list = sample_raster(block, lc_point_list, raster, con)

                # Update the node dict
                for row in lc_point_list:
                    key = "{0}_{1}".format(type, row[10])
                    nodeDict[row[5]][key] = row[11 + t]

            # Write the landcover data to the TTools point feature class
            update_fc(nodeDict, nodes_fc, addFields, nodes_in_block)

            # Build the output point feature class using the data
            update_lc_point_fc(lc_point_list, list(rasterDict.keys()),
                               lc_point_fc, nodes_fc, nodes_in_block,
                               overwrite_data, proj)

            total_samples = total_samples + len(lc_point_list)
            del lc_point_list
            gc.collect()

        endTime = time.time()

        elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
        mspersample = timedelta(seconds=(endTime - startTime) /
                                total_samples).microseconds
        message("Process Complete in {0} minutes. {1} microseconds per sample".format(elapsedmin, mspersample))

    except Exception:
        tbinfo = traceback.format_exc()
        pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
        raise


# -----------------------------------------------------------------------
# Wrapper functions
# -----------------------------------------------------------------------

def step5_star(nodes_fc, trans_count, transsample_count, transsample_distance,
               zone_sample, heatsource8,
               lc_raster, lc_units, z_raster, z_units, lc_point_fc,
               block_size=10, overwrite_data=True):
    """TTools Step 5: Sample Landcover - Star Pattern

    The star pattern sampling method will take an input point feature (from Step 1) and sample input landcover
    rasters along transects in any number of azimuth directions oriented outward from the stream node. The number
    of transects, sample points along each transect, and the distance between samples is user defined. The measured
    distance along each transect starts at the stream node. The star pattern sampling method can be used to develop
    landcover inputs for Heat Source 7 - 9.

    TTools steps 1 - 3 must be run before Step 5.

    Parameters:
        nodes_fc (str): Path to the TTools point feature class.
        trans_count (int): Number of transects per node. The degrees of
            separation between each transect is equal to 360 / trans_count. If trans_count = 8,
            the heading of the first transect is 45 degrees (Northeast), the second transect is
            90 degrees (East), and the eighth is 360 degrees (North). Transects are oriented
            clockwise from north. trans_count is ignored if heatsource8 = True.
        transsample_count (int): Number of samples per transect. The number DOES NOT include the
            sample at the stream node.
        transsample_distance (float): The distance between transect samples (meters).
        zone_sample (bool): True/False flag to indicate if the sample should
            represent a zone and be centered relative to the transsample_distance as measured from
            the stream node. If this is True the distance from the stream node to the sample
            location along each transect is equal to the transsample_distance * (sample number -
            0.5). This should be True if using heat source 7, heat source 8.0.1 - 8.0.5, or heat
            source 8.0.7 - 8.0.8 when the model is configured to use the zone method. This should
            be False if using heat source 9 or when using heat source 8.0.7 - 8.0.8 and the model
            is configured to use the point method.
        heatsource8 (bool): True/False flag to indicate if the star pattern
            with 7 transects should be used. Heat source version 7 and 8 use 7 transects around
            each stream node in the following directions: Northeast, East, Southeast, South,
            Southwest, West, Northwest.
        lc_raster (str): Path and name of the land cover code, height, or elevation raster.
        lc_units (str): z units of the lc_raster (aka units of height or elevation). Use "Feet",
            "Meters", or None if the lc_raster values are codes and do not represent elevation or
            height units.
        z_raster (str): Path and name of the ground elevation raster.
        z_units (str): z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is
            not in feet or meters the elevation values must be converted.
        lc_point_fc (str): Path and name of output sample point feature file.
        block_size (int): The x and y size in kilometers for each raster block pulled into an array.
            Start with 10 if you aren't sure and reduce if there is an error. To increase processing
            speed rasters are subdivided iteratively into smaller blocks and pulled into arrays for
            processing. Very large block sizes may use a lot of memory and result in an error.
        overwrite_data (bool): True/False flag if existing data in nodes_fc and lc_point_fc can be
            overwritten.

    Outputs:
        nodes_fc: New fields are added into nodes_fc with the landcover and elevation values for
            each transect sample.
        lc_point_fc: New point feature class created with a point at each x/y sample and the
            sample raster values.
    """
    step5(nodes_fc, method="star",
          transsample_count=transsample_count,
          transsample_distance=transsample_distance,
          lc_raster=lc_raster, lc_units=lc_units,
          z_raster=z_raster, z_units=z_units,
          lc_point_fc=lc_point_fc,
          block_size=block_size, overwrite_data=overwrite_data,
          trans_count=trans_count, zone_sample=zone_sample,
          heatsource8=heatsource8)


def step5_ortho(nodes_fc, start_bank, transsample_count, transsample_distance,
                lc_raster, lc_units, z_raster, z_units, lc_point_fc,
                block_size=10, overwrite_data=True):
    """TTools Step 5: Sample Landcover - Orthogonal Method

    The orthogonal sampling method will take an input point feature (from Step 1) and sample an input landcover
    raster along two transects that are orthogonal to the stream aspect (left and right looking downstream). The
    number of sample points along each transect and distance between samples is user defined. The transect can
    start at the stream node or at the right and left banks. The Orthogonal sampling method can be used to develop
    inputs for Heat Source 6, Washington Department of Ecology's Shade model, and the shade file for CE-QUAL-W2.

    TTools steps 1 - 3 must be run before Step 5.

    Parameters:
        nodes_fc (str): Path to the TTools point feature class.
        start_bank (bool): True/False flag to indicate if the transect
            should start at the stream bank. If False the transect will start at the stream node.
            Normally set to True.
        transsample_count (int): Number of samples per transect. The number DOES NOT include the
            sample at the stream node.
        transsample_distance (float): The distance between transect samples (meters).
        lc_raster (str): Path and name of the land cover code, height, or elevation raster.
        lc_units (str): z units of the lc_raster (aka units of height or elevation). Use "Feet",
            "Meters", or None if the lc_raster values are codes and do not represent elevation or
            height units.
        z_raster (str): Path and name of the ground elevation raster.
        z_units (str): z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is
            not in feet or meters the elevation values must be converted.
        lc_point_fc (str): Path and name of output sample point feature file.
        block_size (int): The x and y size in kilometers for each raster block pulled into an array.
            Start with 10 if you aren't sure and reduce if there is an error. To increase processing
            speed rasters are subdivided iteratively into smaller blocks and pulled into arrays for
            processing. Very large block sizes may use a lot of memory and result in an error.
        overwrite_data (bool): True/False flag if existing data in nodes_fc and lc_point_fc can be
            overwritten.

    Outputs:
        nodes_fc: New fields are added into nodes_fc with the landcover and elevation values for
            each transect sample.
        lc_point_fc: New point feature class created with a point at each x/y sample and the
            sample raster values.
    """
    step5(nodes_fc, method="orthogonal",
          transsample_count=transsample_count,
          transsample_distance=transsample_distance,
          lc_raster=lc_raster, lc_units=lc_units,
          z_raster=z_raster, z_units=z_units,
          lc_point_fc=lc_point_fc,
          block_size=block_size, overwrite_data=overwrite_data,
          start_bank=start_bank)
