"""TTools Step 4: Measure Topographic Shade Angles"""
import sys
import gc
import time
import traceback
from datetime import timedelta
from math import radians, sin, cos, ceil
from collections import defaultdict

import numpy as np

from ttools.utils import (to_meters_con, from_meters_con,
                          from_z_units_to_meters_con,
                          message)
from ttools.geo_package import (read_fc, update_fc, update_sample_fc,
                            get_crs, crs_equal, get_raster_info, raster_to_array,
                            fc_exists, delete_fc)


def nested_dict():
    """Build a nested dictionary"""
    return defaultdict(nested_dict)


def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the STREAM_ID,
    NODE_ID, and X/Y coordinates as a dictionary"""

    message("Reading nodes feature class")

    nodeDict = read_fc(nodes_fc)

    if overwrite_data:
        for nodeID in nodeDict:
            for field in addFields:
                nodeDict[nodeID].pop(field, None)

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


def update_topo_fc(topo_list, topo_fc, nodes_to_update, overwrite_data, proj):
    """Creates/updates the output topo point feature
    class using the data from the topo list"""

    fields_to_update = ["POINT_X", "POINT_Y", "STREAM_ID", "NODE_ID",
                        "AZIMUTH", "TOPOANGLE", "TOPO_ELE", "NODE_ELE",
                        "ELE_CHANGE", "TOPODIS", "SEARCHDIS", "NA_SAMPLES"]

    update_sample_fc(topo_list, topo_fc, fields_to_update,
                     nodes_to_update, overwrite_data, proj)


def build_search_array(searchDistance_min, searchDistance_max, cellsize, use_skippy):
    """Build a numpy array from the minimum to the max
    search distance by increments of the cellsize."""

    # use next cell over to avoid divide by zero errors
    if searchDistance_min <= 0: searchDistance_min = cellsize

    if use_skippy:
        # This is a modified version of the skippy
        # algorithm from Greg Pelletier. It is not being used but I'm
        # keeping it in here in case someone wants to turn it on.
        # It needs to be fixed so the ncells are adjusted based on distance
        searchDistance = searchDistance_min
        distanceList = [searchDistance_min]

        # ncells = max([int(ceil((block_x_max - block_x_min)/ x_cellsize)), 1])

        ncells = 1
        while not searchDistance > searchDistance_max:
            if ncells <= 10:
                searchDistance = searchDistance + (cellsize)
            if 10 < ncells <= 20:
                searchDistance = searchDistance + (cellsize * 3)
            if 20 < ncells <= 40:
                searchDistance = searchDistance + (cellsize * 6)
            if 40 < ncells <= 50:
                searchDistance = searchDistance + (cellsize * 12)
            if 50 < ncells <= 60:
                searchDistance = searchDistance + (cellsize * 25)
            if ncells > 60:
                searchDistance = searchDistance + (cellsize * 50)
            distanceList.append(searchDistance)
            ncells = ncells + 1
        distance_array = np.array(distanceList)
    else:
        if (searchDistance_max - searchDistance_min >= cellsize):
            distance_array = np.arange(searchDistance_min, searchDistance_max, cellsize)
        else:
            distance_array = np.array([searchDistance_min])
    return distance_array


def create_blocks(nodeDict, block_size, searchDistance_max, buffer, azimuths):
    """Returns two lists, one containing the coordinate extent
    for each block that will be iteratively extracted to an array
    and the other containing the start and stop distances for each
    topo line that is within the block extent."""

    message("Preparing blocks and topo sampling")

    # Create a dictionary to lookup which nodesIDs can be updated after
    # the block has been sampled.
    blockDict = nested_dict()

    # Get a list of the nodes, sort them
    nodes = list(nodeDict.keys())
    nodes.sort()

    topo_list = []
    x_coord_list = []
    y_coord_list = []

    for nodeID in nodes:
        node_x = nodeDict[nodeID]["POINT_X"]
        node_y = nodeDict[nodeID]["POINT_Y"]
        streamID = nodeDict[nodeID]["STREAM_ID"]
        z_node = nodeDict[nodeID]["Z_NODE"]

        for a in azimuths:
            # calculate x/y coordinates at max search distance
            end_x = ((searchDistance_max * sin(radians(a))) + node_x)
            end_y = ((searchDistance_max * cos(radians(a))) + node_y)

            x_coord_list.append(end_x)
            y_coord_list.append(end_y)

            topo_list.append([nodeID, streamID, a, z_node, node_x, node_y, end_x, end_y])

    # calculate bounding box extent for samples
    x_min = min(x_coord_list)
    x_max = max(x_coord_list)
    y_min = min(y_coord_list)
    y_max = max(y_coord_list)

    x_width = int(x_max - x_min + 1)
    y_width = int(y_max - y_min + 1)

    # Block Number
    b = 0

    # Build blocks
    for x in range(0, x_width, block_size):
        for y in range(0, y_width, block_size):

            # Lower left coordinate of block (in map units)
            block_x_min = min([x_min + x, x_max])
            block_y_min = min([y_min + y, y_max])
            # Upper right coordinate of block (in map units)
            block_x_max = min([block_x_min + block_size, x_max])
            block_y_max = min([block_y_min + block_size, y_max])

            topo_in_block = []

            block_segments = (((block_x_min, block_y_max), (block_x_min, block_y_min)),
                              ((block_x_min, block_y_min), (block_x_max, block_y_min)),
                              ((block_x_max, block_y_min), (block_x_max, block_y_max)),
                              ((block_x_max, block_y_max), (block_x_min, block_y_max)))

            # Now start iterating through the topo list to evaluate
            # if any part of the topo line is in the block extent
            for nodeID, streamID, a, z_node, node_x, node_y, end_x, end_y in topo_list:

                contains_node = False
                contains_end = False

                # check if the node is inside the block
                if (block_x_min <= node_x <= block_x_max and
                        block_y_min <= node_y <= block_y_max):
                    contains_node = True

                # check if the topo line end point is inside the block
                if (block_x_min <= end_x <= block_x_max and
                        block_y_min <= end_y <= block_y_max):
                    contains_end = True

                # check if the entire topo line segment is inside the block
                if contains_node and contains_end:
                    block_search_start = 0
                    block_search_end = searchDistance_max

                    topo_in_block.append([nodeID, streamID, a,
                                          z_node,
                                          node_x, node_y,
                                          end_x, end_y,
                                          block_search_start,
                                          block_search_end])

                # check if and where the topo segment cross the block
                else:
                    distance = []

                    # check each block segment for an intersection
                    for i, block_segment in enumerate(block_segments):

                        intersects, inter1_x, inter1_y, inter2_x, inter2_y = find_intersection(block_segment[0],
                                                                                               block_segment[1],
                                                                                               (node_x, node_y),
                                                                                               (end_x, end_y), True)

                        # if there is an intersection calculate the
                        # distance in fc units from node to
                        # that intersection
                        if intersects:
                            if a in [0, 180]:
                                # need to use the y direction for these
                                distance.append((inter1_y - node_y) / cos(radians(a)))
                            else:
                                distance.append((inter1_x - node_x) / sin(radians(a)))

                    if (len(distance) == 1 and
                            (0 < distance[0] < searchDistance_max) and
                            (contains_node or contains_end)):
                        # one intersection
                        if contains_node:
                            # This will be changed from zero to
                            # the cell size in build_search_array()
                            # Should probably just change it here
                            block_search_start = 0
                            block_search_end = distance[0]
                        elif contains_end:
                            # end of the topo line
                            block_search_start = distance[0]
                            block_search_end = searchDistance_max

                        # part of the topo line is in the block, add it
                        topo_in_block.append([nodeID, streamID, a,
                                              z_node,
                                              node_x, node_y,
                                              end_x, end_y,
                                              block_search_start,
                                              block_search_end])

                    elif len(distance) > 1:
                        # two intersections, crosses the block
                        # two intersections, end or start on block line
                        # three intersections, collinear over length of block
                        # three intersections, collinear w/ end or start on block line
                        block_search_start = min(i for i in distance if i is not None)
                        block_search_end = max(distance)

                        topo_in_block.append([nodeID, streamID, a,
                                              z_node,
                                              node_x, node_y,
                                              end_x, end_y,
                                              block_search_start,
                                              block_search_end])

                    del distance[:]

            if topo_in_block:
                # order 0 left,      1 bottom,    2 right,     3 top
                blockDict[b]["extent"] = (block_x_min - buffer, block_y_min - buffer,
                                      block_x_max + buffer, block_y_max + buffer)
                blockDict[b]["sample_extent"] = (block_x_min, block_y_min,
                                                 block_x_max, block_y_max)
                blockDict[b]["samples"] = topo_in_block
                blockDict[b]["nodes_to_update"] = []

            b = b + 1

    node_last_block = {}
    blockIDs = list(blockDict.keys())
    blockIDs.sort()
    for blockID in blockIDs:
        for sample in blockDict[blockID]["samples"]:
            node_last_block[sample[0]] = blockID

    for nodeID in node_last_block:
        blockDict[node_last_block[nodeID]]["nodes_to_update"].append(nodeID)

    return blockDict


def filter_distance_array(distance_array, node_x, node_y, a, sample_extent):
    """Return search distances and coordinates that fall within the block extent."""
    sin_a = sin(radians(a))
    cos_a = cos(radians(a))
    pt_x_array = distance_array * sin_a + node_x
    pt_y_array = distance_array * cos_a + node_y

    in_block = ((pt_x_array >= sample_extent[0]) &
                (pt_x_array <= sample_extent[2]) &
                (pt_y_array >= sample_extent[1]) &
                (pt_y_array <= sample_extent[3]))

    return distance_array[in_block], pt_x_array[in_block], pt_y_array[in_block]


def find_intersection(a, b, c, d, check_collinear=True):
    """Calculates 2D intersection coordinates of segments a-b and c-d.
    a,b,c,d are tuples in the form of (x,y). If checking for collinearity
    then segments that overlap or touch but do not cross will be
    considered an intersection. Returns the min and max points of
    intersection for overlap and the same points for intersections."""

    # a-b = block segment
    # c-d = topo line

    Dx_Cx = d[0] - c[0]
    Ay_Cy = a[1] - c[1]
    Dy_Cy = d[1] - c[1]
    Ax_Cx = a[0] - c[0]
    Bx_Ax = b[0] - a[0]
    By_Ay = b[1] - a[1]

    Cx_Bx = c[0] - b[0]
    Dx_Ax = d[0] - a[0]
    Cy_By = c[1] - b[1]
    Dy_Ay = d[1] - a[1]

    numerator_a = Dx_Cx * Ay_Cy - Dy_Cy * Ax_Cx
    numerator_b = Bx_Ax * Ay_Cy - By_Ay * Ax_Cx
    denominator = Dy_Cy * Bx_Ax - Dx_Cx * By_Ay

    if (check_collinear and numerator_a == 0 and numerator_b == 0 and denominator == 0):
        # Lines are collinear

        # if the signs are different the
        # segments have some overlap
        overlap_x = (Cx_Bx < 0) != (Dx_Ax < 0)
        overlap_y = (Cy_By < 0) != (Dy_Ay < 0)

        point_overlap = (a == d or b == c)

        # if any are True there is an intersection
        if (overlap_x or overlap_y or point_overlap):
            # There is overlap
            x = sorted((a[0], b[0], c[0], d[0]))
            y = sorted((a[1], b[1], c[1], d[1]))

            # min
            ixa = x[1]
            iya = y[1]

            # max
            ixb = x[2]
            iyb = y[2]

            return True, ixa, iya, ixb, iyb
        return False, None, None, None, None

    if denominator == 0:
        # Lines are parallel or
        # maybe collinear if check_collinear=False

        return False, None, None, None, None

    u_a = numerator_a / denominator
    u_b = numerator_b / denominator

    if (u_a >= 0) and (u_a <= 1) and (u_b >= 0) and (u_b <= 1):
        # segments intersect
        ixa = a[0] + Bx_Ax * u_a
        iya = a[1] + By_Ay * u_a

        ixb = c[0] + Dx_Cx * u_b
        iyb = c[1] + Dy_Cy * u_b

        return True, ixa, iyb, ixa, iyb
    return False, None, None, None, None


def get_topo_angles(block_extent, sample_extent, block_samples, z_raster, azimuthdisdict, searchDistance_max_m, con_z_to_m, con_to_m):
    """This gets the maximum topographic angle and other information for
    each topo line within the block. The data is saved to the nodeDict
    as a list."""

    if con_z_to_m is not None:
        nodata_to_value = -9999 / con_z_to_m
    else:
        nodata_to_value = -9999

    raster_info = get_raster_info(z_raster)
    x_cellsize = raster_info["x_cellsize"]
    y_cellsize = raster_info["y_cellsize"]

    # localize the block extent values and add one cell distance to the size to ensure all cells that need to be
    # sampled are included in the array.
    block_x_min = block_extent[0] - x_cellsize
    block_y_min = block_extent[1] - y_cellsize
    block_x_max = block_extent[2] + x_cellsize
    block_y_max = block_extent[3] + y_cellsize

    # Get the coordinates extent of the input raster
    # this could be in main so it doesn't have to run each time
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
        z_array, block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize = raster_to_array(
            z_raster, block_x_min_corner, block_y_min_corner,
            block_x_max_corner, block_y_max_corner, nodata_to_value)
    except:
        tbinfo = traceback.format_exc()
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"
        raise MemoryError(pymsg)

    # convert array values to meters if needed
    if con_z_to_m is not None:
        z_array = z_array * con_z_to_m

    topo_samples = []
    distance_array_dict = {}
    for sample in block_samples:
        a = sample[2]
        if a not in distance_array_dict:
            cellsize = azimuthdisdict[a]
            distance_array_dict[a] = build_search_array(0,
                                                        searchDistance_max_m,
                                                        cellsize,
                                                        use_skippy=False)

    # There is at least one pixel of data
    for (nodeID, streamID, a, z_node,
         node_x, node_y, end_x, end_y,
         block_search_start, block_search_end) in block_samples:

        # create list of distance movements along the topo line
        # that are within the block extent
        # this is in units of the fc
        distance_array, pt_x_array, pt_y_array = filter_distance_array(
            distance_array_dict[a], node_x, node_y, a, sample_extent)
        if len(distance_array) == 0:
            continue

        sin_a = sin(radians(a))
        cos_a = cos(radians(a))

        # Vectorized coord_to_array
        # dx is the x distance from the raster block's left edge to each topo sample point.
        # dy is the y distance from the raster block's top edge down to each topo sample point.
        dx = pt_x_array - block_x_min_corner
        dy = block_y_max_corner - pt_y_array
        col_x_array = ((dx - (dx % x_cellsize)) / x_cellsize).astype(int)
        row_y_array = ((dy - (dy % y_cellsize)) / y_cellsize).astype(int)

        # Clip to valid array indices
        col_x_array = np.clip(col_x_array, 0, z_array.shape[1] - 1)
        row_y_array = np.clip(row_y_array, 0, z_array.shape[0] - 1)

        # indices that are not off raster
        valid_index = ((0 <= row_y_array) & (row_y_array < z_array.shape[0]) &
                       (0 <= col_x_array) & (col_x_array < z_array.shape[1]))

        z_topo_array = np.full(len(row_y_array), -9999.0)
        z_topo_array[valid_index] = z_array[row_y_array[valid_index], col_x_array[valid_index]]

        # convert distances to meters
        distance_array_m = distance_array * con_to_m

        # Calculate the topo angles along the topo line
        angle_array = np.degrees(np.arctan((z_topo_array - z_node) / distance_array_m))
        # remove the off raster samples
        naindex = np.where(z_topo_array < -9998)
        for x in naindex[0]: angle_array[x] = -9999
        # Find the max topo angle
        topoAngle = angle_array.max()
        # array index at the max topo angle
        arryindex = np.where(angle_array == topoAngle)[0][0]
        z_topo = z_topo_array[arryindex]
        # elevation change between topo angle location and node elevation
        z_change = z_topo - z_node
        # distance from the node to topo angle location in units of fc
        topoAngleDistance = distance_array[arryindex]
        topoAngle_x = (topoAngleDistance * sin_a) + node_x
        topoAngle_y = (topoAngleDistance * cos_a) + node_y
        topoAngleDistance_m = topoAngleDistance * con_to_m
        off_rastersamples = (z_topo_array < -9998).sum()

        topo_samples.append([topoAngle_x, topoAngle_y,
                             topoAngle_x, topoAngle_y,
                             streamID, nodeID, a,
                             topoAngle, z_topo,
                             z_node, z_change,
                             topoAngleDistance_m,
                             searchDistance_max_m,
                             off_rastersamples])

    return topo_samples


def step4(nodes_fc, topo_directions, searchDistance_max_km, z_raster, z_units,
          topo_fc, block_size=10, overwrite_data=True):
    """TTools Step 4: Measure Topographic Shade Angles

    This script will take an input point feature (from Step 1) and calculate the maximum topographic
    shade angle for each stream node in different directions.

    TTools steps 1 - 3 must be run before Step 4.

    Parameters:
        nodes_fc (str): Path to the TTools point feature class.
        topo_directions (int): An integer value corresponding to the specific list of azimuth
            directions to sample (e.g. topo_directions = 1). Options listed below.
            1. Heat Source or Washington Department of Ecology's Shade model: [270, 180, 90]
            2. All cardinal and intercardinal directions: [45, 90, 135, 180, 225, 270, 315, 0/360]
            3. CE-QUAL-W2: [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280,
               300, 320, 340]
        searchDistance_max_km (float): The maximum distance in kilometers to search for the largest
            topographic shade angle.
        z_raster (str): Path and name of the ground elevation raster.
        z_units (str): z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is
            not in feet or meters the elevation values must be converted.
        topo_fc (str): Path and name of output topographic point feature. This feature identifies
            the location on the z_raster producing the largest topographic shade angle for each node
            and sample direction.
        block_size (int): The x and y size in kilometers for the z_raster blocks pulled into an
            array. Start with 10 if you aren't sure and reduce if there is an error. To increase
            processing the z_raster is subdivided iteratively into smaller blocks and pulled into
            arrays for processing. Very large block sizes may use a lot of memory and result in
            an error.
        overwrite_data (bool): True/False flag if existing data in nodes_fc and topo_fc can be
            overwritten.

    Outputs:
        nodes_fc: New fields are added into nodes_fc with the maximum topographic shade angles
            for each direction at each node.
        topo_fc: New point feature class created with a point at each x/y location that produces
            the largest topographic shade angle for each node and sample direction.
    """

    # enable garbage collection
    gc.enable()

    try:
        message("Step 4: Measure Topographic Angles")

        # keeping track of time
        startTime = time.time()

        # Check if the node fc exists
        if not fc_exists(nodes_fc):
            raise ValueError("This output does not exist: \n" +
                     "{0}\n".format(nodes_fc))

        # Check if the topo sample fc exists and delete if needed
        if fc_exists(topo_fc) and overwrite_data:
            delete_fc(topo_fc)

        # Determine input point spatial units
        proj = get_crs(nodes_fc)

        # Check to make sure the raster and input
        # points are in the same projection.
        if not crs_equal(nodes_fc, z_raster):
            raise ValueError("Input points and elevation raster do not have the " +
                     "same projection. Please reproject your data.")

        # Get the units conversion factors
        con_z_to_m = from_z_units_to_meters_con(z_units)
        con_to_m = to_meters_con(nodes_fc)
        con_from_m = from_meters_con(nodes_fc)
        searchDistance_max_m = searchDistance_max_km * 1000  # in meters

        # search distance in units of the fc
        searchDistance_max = int(con_from_m * searchDistance_max_m)

        # convert block size from km to meters to units of the node fc
        # in the future block size should be estimated based on available memory
        # memorysize = datatypeinbytes*nobands*block_size^2
        # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
        block_size = int(con_from_m * block_size * 1000)

        # Get the elevation raster cell size in units of the raster
        raster_info = get_raster_info(z_raster)
        x_cellsize = raster_info["x_cellsize"]
        y_cellsize = raster_info["y_cellsize"]

        if topo_directions == 1:
            azimuths = [270, 180, 90]

            azimuthdict = {90: "TOPO_E", 180: "TOPO_S", 270: "TOPO_W"}

        elif topo_directions == 2:  # All directions
            azimuths = [45, 90, 135, 180, 225, 270, 315, 0]

            azimuthdict = {45: "TOPO_NE", 90: "TOPO_E", 135: "TOPO_SE",
                           180: "TOPO_S", 225: "TOPO_SW", 270: "TOPO_W",
                           315: "TOPO_NW", 0: "TOPO_N"}

        elif topo_directions == 3:
            azimuths = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340]

            azimuthdict = {0: "TOPO01", 20: "TOPO02", 40: "TOPO03", 60: "TOPO04", 80: "TOPO05", 100: "TOPO06", 120: "TOPO07",
                           140: "TOPO08", 160: "TOPO09", 180: "TOPO10", 200: "TOPO11", 220: "TOPO12", 240: "TOPO13",
                           260: "TOPO14", 280: "TOPO15", 300: "TOPO16", 320: "TOPO17", 340: "TOPO18"}

        azimuthdisdict = dict.fromkeys(azimuths, None)

        for a in azimuthdisdict:
            if a <= 45:
                azimuthdisdict[a] = abs(y_cellsize / cos(radians(a)))
            elif 45 < a <= 135:
                azimuthdisdict[a] = abs(x_cellsize / sin(radians(a)))
            elif 135 < a <= 225:
                azimuthdisdict[a] = abs(y_cellsize / cos(radians(a)))
            elif 225 < a <= 315:
                azimuthdisdict[a] = abs(x_cellsize / sin(radians(a)))
            elif 315 < a <= 360:
                azimuthdisdict[a] = abs(y_cellsize / cos(radians(a)))

        # Build the topo field names
        addFields = [azimuthdict[a] for a in azimuths]

        # Read the feature class data into a dictionary
        nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)

        # Build the blockDict
        blockDict = create_blocks(nodeDict, block_size, searchDistance_max,
                                  x_cellsize, azimuths)

        # Iterate through each block
        all_topo_list = []
        all_nodes_updated = []
        blockIDs = list(blockDict.keys())
        blockIDs.sort()

        for p, blockID in enumerate(blockIDs):
            block_extent = blockDict[blockID]["extent"]
            sample_extent = blockDict[blockID]["sample_extent"]
            block_samples = blockDict[blockID]["samples"]
            block_samples.sort()

            message("Processing block {0} of {1}".format(p + 1, len(blockIDs)))

            # calculate coordinates along the
            # portion of the topo line in the block,
            # convert raster to array, sample the raster
            # calculate the topo angles and other info
            topo_samples = get_topo_angles(block_extent, sample_extent, block_samples,
                                           z_raster, azimuthdisdict,
                                           searchDistance_max, con_z_to_m,
                                           con_to_m)
            if topo_samples:
                # Update the nodeDict
                for sample in topo_samples:
                    nodeID = sample[5]
                    a = sample[6]
                    topoAngle = sample[7]

                    # Create a key to hold the topo list info for this block
                    topo_key = azimuthdict[a] + "_list"

                    if azimuthdict[a] in nodeDict[nodeID]:
                        if nodeDict[nodeID][azimuthdict[a]] < topoAngle:
                            nodeDict[nodeID][azimuthdict[a]] = topoAngle
                            nodeDict[nodeID][topo_key] = sample

                    else:
                        nodeDict[nodeID][azimuthdict[a]] = topoAngle
                        nodeDict[nodeID][topo_key] = sample

                del topo_samples

            # Check if any nodes can be updated in the node and topo fc
            if blockDict[blockID]["nodes_to_update"]:
                nodes_to_update = blockDict[blockID]["nodes_to_update"]

                # Write the topo data to the TTools point feature class
                update_fc(nodeDict, nodes_fc, addFields, nodes_to_update)

                # Build/add to the output topo feature class
                topo_list = []
                for nodeID in nodes_to_update:
                    for field in addFields:
                        topo_key = field + "_list"
                        topo_list.append(nodeDict[nodeID][topo_key])
                        # delete some data
                        nodeDict[nodeID].pop(field, None)
                        nodeDict[nodeID].pop(topo_key, None)

                all_topo_list.extend(topo_list)
                all_nodes_updated.extend(nodes_to_update)

                del topo_list
                gc.collect()

        # Write the topo feature class
        if all_topo_list:
            update_topo_fc(all_topo_list, topo_fc, all_nodes_updated,
                           overwrite_data, proj)

        endTime = time.time()
        elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
        mspernode = timedelta(seconds=(endTime - startTime) / len(nodeDict.keys())).microseconds
        message("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))

    except Exception:
        tbinfo = traceback.format_exc()
        pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
        raise
