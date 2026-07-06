"""TTools Step 3: Sample Elevation and Gradient"""
import sys
import gc
import time
import traceback
from datetime import timedelta
import itertools
from math import ceil
from collections import defaultdict

import arcpy

from ttools.utils import (from_meters_con, from_z_units_to_meters_con, coord_to_array,
                          message, read_fc, update_fc, get_raster_info,
                          raster_to_array)


def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the
    NODE_ID and X/Y coordinates as a dictionary"""

    nodeDict = read_fc(nodes_fc)

    if "GRADIENT" in addFields:
        # Check to see if the 1st field exists if yes add it to
        # the cursorfields to be retrieved.
        check_field = addFields[0]
    else:
        # Check to see if the last field exists if yes add it.
        # Grabs last field because often the first field, emergent, is zero
        check_field = addFields[len(addFields)-1]

    first_node = next(iter(nodeDict.values()))

    if not overwrite_data and check_field in first_node:
        filtered = {}
        for nodeID in nodeDict:
            val = nodeDict[nodeID].get(check_field)
            if "GRADIENT" in addFields:
                # 0 is a valid gradient value, don't treat as null
                if val is None or val < -9998:
                    filtered[nodeID] = nodeDict[nodeID]
            else:
                # if the data is null or zero (0 = default for shapefile),
                # it is retrieved and will be overwritten.
                if val is None or val == 0 or val < -9998:
                    filtered[nodeID] = nodeDict[nodeID]
        nodeDict = filtered

    if "GRADIENT" in addFields and len(nodeDict) == 0:
        raise ValueError("The gradient field checked in the input point feature class "+
                 "have existing data. There is nothing to process. Exiting")

    return nodeDict


def create_block_list(nodeDict, nodes, block_size, buffer):
    """Returns two lists, one containing the coordinate extent
    for each block that will be iteratively extracted to an array
    and the other containing the stream and node IDs within each
    block extent."""

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

                    nodes_in_block.append([nodeID, node_x, node_y])

            if nodes_in_block:

                # Minimize the size of the block0 by the true
                # extent of the nodes in the block
                node_x_min = min([nodeDict[nodeID[0]]["POINT_X"] for nodeID in nodes_in_block])
                node_y_min = min([nodeDict[nodeID[0]]["POINT_Y"] for nodeID in nodes_in_block])
                node_x_max = max([nodeDict[nodeID[0]]["POINT_X"] for nodeID in nodes_in_block])
                node_y_max = max([nodeDict[nodeID[0]]["POINT_Y"] for nodeID in nodes_in_block])

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


def calculate_gradient(zList, len_list, smooth_flag):

    skipupNodes = [0]
    gradientList = [0 for i in zList]

    for i in range(1, len(zList)):
        z = zList[i]
        zUp = zList[i - 1 - max(skipupNodes)]

        # Check if the gradient is <= 0, if yes keep going until
        # it is positive again. Then recalculate the gradient over
        # the longer distance
        if z > zUp and smooth_flag is True:
            skipupNodes.append(max(skipupNodes) + 1)
        else:
            dx_meters = sum(len_list[i:i + max(skipupNodes) + 1])
            #dx_meters = float(kmList[i] - kmList[i- max(skipdownNodes)]) * 1000
            gradient = (zUp - z) / dx_meters
            for Skip in skipupNodes:
                gradientList[i - Skip] = gradient
            skipupNodes = [0]

    return (gradientList)


def sample_raster(block, nodes_in_block, z_raster, cellcoords, con_z_to_m):

    if con_z_to_m is not None:
        nodata_to_value = -9999 / con_z_to_m
    else:
        nodata_to_value = -9999

    raster_info = get_raster_info(z_raster)
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
            z_raster, block_x_min_corner, block_y_min_corner,
            block_x_max_corner, block_y_max_corner, nodata_to_value)
    except:
        tbinfo = traceback.format_exc()
        pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"
        raise MemoryError(pymsg)

    # convert array values to meters if needed
    if con_z_to_m is not None:
        raster_array = raster_array * con_z_to_m

    z_list = []
    if raster_array.max() > -9999:
        # There is at least one pixel of data
        n_rows, n_cols = raster_array.shape
        for node in nodes_in_block:
            xy = coord_to_array(node[1], node[2], block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize)
            z_sampleList = []
            for coord in cellcoords:
                # Calculate the cell X/Y based on the base coordinate movement
                cell_x = xy[0] + coord[0]
                cell_y = xy[1] + coord[1]
                if 0 <= cell_y < n_rows and 0 <= cell_x < n_cols:
                    z_sampleList.append(raster_array[cell_y, cell_x])
                else:
                    # off raster sample
                    z_sampleList.append(-9999)

            # sample at node:
            if 0 <= xy[1] < n_rows and 0 <= xy[0] < n_cols:
                z_node = raster_array[xy[1], xy[0]]
            else:
                # off raster sample
                z_node = -9999

            # Remove no data values (-9999) unless they are all no data
            if not max(z_sampleList) < -9998:
                z_sampleList = [z for z in z_sampleList if z > -9999]
            # Get the lowest elevation
            node.append(min(z_sampleList))
            node.append(z_node)
            z_list.append(node)

    else:
        # No data, add -9999 for elevation and z_node
        for node in nodes_in_block:
            node.append(-9999)
            node.append(-9999)
            z_list.append(node)

    return z_list


def step3(nodes_fc, searchCells, smooth_flag, z_raster, z_units,
          block_size=10, overwrite_data=True):
    """TTools Step 3: Sample Elevation and Gradient

    This script will take an input point feature (from Step 1) and sample the input raster elevation to find the
    lowest elevation in a user defined search radius and calculate the gradient for each node in the downstream direction.

    TTools step 1 must be run before Step 3.

    Parameters:
        nodes_fc (str): Path to the TTools point feature class.
        searchCells (int): The number of cells to search around the node for the lowest elevation.
        smooth_flag (bool): True/False flag to indicate if smoothing should occur when the gradient
            is <= 0. When gradients are <= 0 the algorithm moves downstream until a positive gradient
            is found and then recalculates the gradient over the longer distance and applies the
            resulting value to all nodes over that distance.
        z_raster (str): Path and name of the ground elevation raster.
        z_units (str): z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is
            not in feet or meters the elevation values must be converted.
        block_size (int): The x and y size in kilometers for the z_raster blocks pulled into an
            array. Start with 10 if you aren't sure and reduce if there is an error. To increase
            processing the z_raster is subdivided iteratively into smaller blocks and pulled into
            arrays for processing. Very large block sizes may use a lot of memory and result in
            an error.
        overwrite_data (bool): True/False flag if existing data in nodes_fc can be overwritten.

    Outputs:
        nodes_fc: New fields ELEVATION and GRADIENT are added into nodes_fc.
    """

    # enable garbage collection
    gc.enable()

    try:
        message("Step 3: Sample Stream Elevations/Gradient")

        # keeping track of time
        startTime = time.time()

        # Check if the node fc exists
        if not arcpy.Exists(nodes_fc):
            raise ValueError("This output does not exist: \n" +
                     "{0}\n".format(nodes_fc))

        # Determine input point spatial units
        proj = arcpy.Describe(nodes_fc).spatialReference
        proj_ele = arcpy.Describe(z_raster).spatialReference

        # Check to make sure the raster and input
        # points are in the same projection.
        if proj.name != proj_ele.name:
            raise ValueError("Input points and elevation raster do not have the " +
                     "same projection. Please reproject your data.")

        # Get the units conversion factor
        con_z_to_m = from_z_units_to_meters_con(z_units)
        con_from_m = from_meters_con(nodes_fc)

        # convert block size from km to meters to units of the node fc
        # in the future block size should be estimated based on available memory
        # memorysize = datatypeinbytes*nobands*block_size^2
        # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
        block_size = int(con_from_m * block_size * 1000)

        # Get the elevation raster cell size
        raster_info = get_raster_info(z_raster)
        cellsize = raster_info["x_cellsize"]

        # calculate the buffer distance (in raster spatial units) to add to
        # the base bounding box when extracting to an array. The buffer is
        # equal to the cellsize * the searchcells to make sure the block includes
        # the surrounding cells at each corner
        buffer = int((searchCells + 1) * cellsize)

        # Make a list of the base x/y coordinate movements
        # from the node origin. These values will be
        # multiplied by the cell size.
        # searchCells = 0 samples at the node
        # searchCells = 1 cell width around node = 9 cells
        # searchCells = 2 cell widths around node = 25 cells ...
        cell_moves = [i for i in range(searchCells * -1, searchCells + 1, 1)]
        cellcoords = list(itertools.product(cell_moves, cell_moves))

        # read the data into a dictionary
        addFields = ["ELEVATION", "Z_NODE"]
        nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
        if len(nodeDict) != 0:

            # Get a list of the nodes, sort them
            nodes = list(nodeDict.keys())
            nodes.sort()
            n_nodes = len(nodes)

            # Build the block list
            block_extents, block_nodes = create_block_list(nodeDict, nodes, block_size, buffer)

            # Iterate through each block, calculate sample coordinates,
            # convert raster to array, sample the raster
            total_samples = 0

            for p, block in enumerate(block_extents):
                nodes_in_block = block_nodes[p]
                nodes_in_block.sort()

                message("Processing block {0} of {1}".format(p + 1, len(block_extents)))

                z_list = sample_raster(block, nodes_in_block, z_raster, cellcoords, con_z_to_m)

                # Update the node fc
                for row in z_list:
                    nodeDict[row[0]]["ELEVATION"] = row[3]
                    nodeDict[row[0]]["Z_NODE"] = row[4]

                nodes_to_query = [row[0] for row in z_list]

                # Write the elevation data to the TTools point feature class
                message("Updating input nodes feature class")
                update_fc(nodeDict, nodes_fc, addFields, nodes_to_query)

                total_samples = total_samples + len(z_list)
                del z_list
                gc.collect()

        else:
            message("The elevation field checked in the input point feature class " +
                    "have existing data. Advancing to gradient processing")
        del(nodeDict)

        # Start on gradients

        # read the data into a dictionary
        addFields = ["GRADIENT"]
        nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
        nodes_to_query = []

        # Group nodes by STREAM_ID for gradient calculation
        stream_groups = defaultdict(list)
        for nodeID in nodeDict:
            streamID = nodeDict[nodeID]["STREAM_ID"]
            stream_groups[streamID].append(nodeID)

        for n, streamID in enumerate(stream_groups):
            message("Calculating gradients stream {0} of {1}".format(n + 1, len(stream_groups)))

            # Sort nodes within stream by STREAM_KM descending (upstream to downstream)
            stream_nodes = sorted(stream_groups[streamID],
                                  key=lambda nodeID: nodeDict[nodeID]["STREAM_KM"],
                                  reverse=True)

            z_list = [nodeDict[nodeID]["ELEVATION"] for nodeID in stream_nodes]
            len_list = [nodeDict[nodeID]["LENGTH"] for nodeID in stream_nodes]

            # Calculate Gradient
            gradientList = calculate_gradient(z_list, len_list, smooth_flag)

            for i, nodeID in enumerate(stream_nodes):
                nodeDict[nodeID]["GRADIENT"] = gradientList[i]
                nodes_to_query.append(nodeID)

        message("Updating input nodes feature class")
        update_fc(nodeDict, nodes_fc, addFields, nodes_to_query)

        endTime = time.time()
        gc.collect()

        elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
        mspernode = timedelta(seconds=(endTime - startTime) / n_nodes).microseconds
        message("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))

    except Exception:
        tbinfo = traceback.format_exc()
        pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
        print(pymsg)
        raise
