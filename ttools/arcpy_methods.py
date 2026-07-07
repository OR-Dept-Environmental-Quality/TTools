"""
Arcpy methods for various TTools processing tasks. These methods are called
by geo_package.py when use_arcpy is True.
"""

import os
import numpy as np
import arcpy

# ---------------------------------------------------------------------------
# Vector
# ---------------------------------------------------------------------------

def fc_exists_arcpy(path):
    return arcpy.Exists(path)


def delete_fc_arcpy(path):
    arcpy.Delete_management(path)


def read_fc_arcpy(nodes_fc, addFields=None):
    """Reads the input point feature class and returns a dictionary
    with NODE_ID used as the dictionary key.

    Uses arcpy.da.SearchCursor with NODE_ID as the first field
    and SHAPE@X/SHAPE@Y for geometry coordinates.

    Parameters:
        nodes_fc (str): Path to the feature class.
        addFields (list of str, optional): Subset of attribute fields to read.
            If None all attribute fields are read (except NODE_ID which
            becomes the key). Geometry coordinates (POINT_X, POINT_Y) are
            always included.

    Returns:
        dict: {nodeID: {field: value, ..., "POINT_X": x, "POINT_Y": y}}
    """


    # Get all field names (excluding geometry and OID fields)
    all_fields = [f.name for f in arcpy.ListFields(nodes_fc)
                  if f.type not in ("Geometry", "OID")]

    if addFields:
        read_fields = [f for f in addFields if f in all_fields]
    else:
        read_fields = [f for f in all_fields if f != "NODE_ID"]

    # NODE_ID first, then requested fields, then geometry
    incursorFields = ["NODE_ID"] + read_fields + ["SHAPE@X", "SHAPE@Y"]

    proj = arcpy.Describe(nodes_fc).spatialReference

    nodeDict = {}
    with arcpy.da.SearchCursor(nodes_fc, incursorFields, "", proj) as Inrows:
        for row in Inrows:
            nodeID = row[0]
            nodeDict[nodeID] = {}
            for i, f in enumerate(read_fields):
                nodeDict[nodeID][f] = row[i + 1]
            nodeDict[nodeID]["POINT_X"] = row[-2]
            nodeDict[nodeID]["POINT_Y"] = row[-1]

    return nodeDict


def write_fc_arcpy(nodeDict, nodes_fc, addFields, proj):
    """Creates the output point node feature class using the data from
    the nodes dictionary.

    NODE_ID is always written as a field. Geometry is created
    from POINT_X and POINT_Y.

    Parameters:
        nodeDict (dict): {nodeID: {field: value, ..., "POINT_X": x,
            "POINT_Y": y}}
        nodes_fc (str): Output feature class path.
        addFields (list of str): Attribute field names to write (not
            including NODE_ID which is always written).
        proj (arcpy.SpatialReference): Coordinate reference system for
            the output.
    """

    out_path = os.path.dirname(nodes_fc)
    out_name = os.path.basename(nodes_fc)

    full_path = os.path.join(out_path, out_name)

    # Delete if exists
    if arcpy.Exists(full_path):
        arcpy.Delete_management(full_path)

    # Create the feature class
    arcpy.CreateFeatureclass_management(out_path, out_name, "POINT",
                                        "", "DISABLED", "DISABLED", proj)

    # Add NODE_ID field
    arcpy.AddField_management(full_path, "NODE_ID", "LONG", "", "", "",
                              "", "NULLABLE", "NON_REQUIRED")

    # Add attribute fields - infer type from data
    first_node = next(iter(nodeDict.values()))
    for f in addFields:
        val = first_node.get(f)
        if isinstance(val, str):
            arcpy.AddField_management(full_path, f, "TEXT", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")
        elif isinstance(val, (int, np.integer)):
            arcpy.AddField_management(full_path, f, "LONG", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")
        else:
            arcpy.AddField_management(full_path, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    # Insert rows
    cursorfields = ["NODE_ID"] + addFields + ["SHAPE@X", "SHAPE@Y"]
    with arcpy.da.InsertCursor(full_path, cursorfields) as cursor:
        for nodeID in nodeDict:
            row = [nodeID]
            for f in addFields:
                row.append(nodeDict[nodeID].get(f))
            row.append(nodeDict[nodeID]["POINT_X"])
            row.append(nodeDict[nodeID]["POINT_Y"])
            cursor.insertRow(row)


def update_fc_arcpy(nodeDict, nodes_fc, addFields, nodes_to_query=None):
    """Updates the input point feature class with data from the
    nodes dictionary. NODE_ID is used as the key field.

    Parameters:
        nodeDict (dict): {nodeID: {field: value, ...}}
        nodes_fc (str): Path to the existing feature class.
        addFields (list of str): Field names to update.
        nodes_to_query (list, optional): List of NODE_IDs to update.
            If None, updates all rows where NODE_ID is found in nodeDict.
    """


    # Add fields if they don't exist
    existingFields = [f.name for f in arcpy.ListFields(nodes_fc)]
    for f in addFields:
        if f not in existingFields:
            if f.startswith("LC"):
                arcpy.AddField_management(nodes_fc, f, "TEXT", "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")
            else:
                arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")

    # Build a query to retrieve just the nodes that need updating
    if nodes_to_query:
        whereclause = "{0} IN ({1})".format(
            "NODE_ID", ','.join(str(i) for i in nodes_to_query))
    else:
        whereclause = ""

    with arcpy.da.UpdateCursor(nodes_fc, ["NODE_ID"] + addFields,
                               whereclause) as cursor:
        for row in cursor:
            nodeID = row[0]
            if nodeID in nodeDict:
                for f, field in enumerate(addFields):
                    if field in nodeDict[nodeID]:
                        row[f + 1] = nodeDict[nodeID][field]
                cursor.updateRow(row)


# Field type mapping for sample feature classes (topo_fc, lc_point_fc)
_SAMPLE_FIELD_TYPES = {
    "POINT_X": "DOUBLE",
    "POINT_Y": "DOUBLE",
    "NODE_ID": "LONG",
    "SAMPLE_ID": "LONG",
    "TRANS_DIR": "DOUBLE",
    "TRANSECT": "SHORT",
    "SAMPLE": "SHORT",
    "KEY": "TEXT",
    "AZIMUTH": "DOUBLE",
    "TOPOANGLE": "DOUBLE",
    "TOPO_ELE": "DOUBLE",
    "NODE_ELE": "DOUBLE",
    "ELE_CHANGE": "DOUBLE",
    "TOPODIS": "DOUBLE",
    "SEARCHDIS": "DOUBLE",
    "NA_SAMPLES": "SHORT",
}


def _get_sample_field_type(field_name, sample_value=None):
    """Get the arcpy field type for a sample feature class field.

    Uses a fixed mapping for known fields. STREAM_ID type is
    inferred from the data. LC fields are TEXT. All other
    raster fields (ELE, LAI, k, OH, CAN) are DOUBLE.
    """
    if field_name in _SAMPLE_FIELD_TYPES:
        return _SAMPLE_FIELD_TYPES[field_name]
    if field_name == "STREAM_ID":
        if isinstance(sample_value, str):
            return "TEXT"
        return "LONG"
    # LC fields are TEXT
    if field_name.startswith("LC"):
        return "TEXT"
    # Other raster fields (ELE, LAI, k, OH, CAN) are DOUBLE
    return "DOUBLE"


def update_sample_fc_arcpy(point_list, sample_fc, fields_to_update,
                           nodes_to_update, overwrite_data, proj):
    """Creates/updates the output sample point feature class using
    the data from the point list.

    Creates the feature class if it doesn't exist, adds fields,
    optionally deletes existing rows for the specified nodes,
    and inserts new rows.

    Parameters:
        point_list (list of list): Each row: [x, y, field1_val, field2_val, ...].
            row[0:2] are geometry coordinates (SHAPE@X, SHAPE@Y).
            row[2:] map to fields_to_update in order.
        sample_fc (str): Path to the output sample feature class.
        fields_to_update (list of str): Field names for the attribute data
            (maps to row[2:]).
        nodes_to_update (list): List of NODE_IDs being updated.
        overwrite_data (bool): If True, overwrite all existing data. If False,
            delete existing rows for nodes_to_update before inserting.
        proj (arcpy.SpatialReference): Coordinate reference system for
            the output.
    """

    out_path = os.path.dirname(sample_fc)
    out_name = os.path.basename(sample_fc)

    full_path = os.path.join(out_path, out_name)

    # Create the feature class if it doesn't exist
    if not arcpy.Exists(full_path):
        arcpy.CreateFeatureclass_management(out_path, out_name, "POINT",
                                            "", "DISABLED", "DISABLED", proj)

        # Get a sample row to infer STREAM_ID type
        sample_row = point_list[0]

        # Add attribute fields
        for i, f in enumerate(fields_to_update):
            field_type = _get_sample_field_type(
                f, sample_value=sample_row[i + 2])
            arcpy.AddField_management(full_path, f, field_type, "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    if not overwrite_data:
        # Delete existing rows for nodes being updated
        whereclause = "{0} IN ({1})".format(
            "NODE_ID", ','.join(str(i) for i in nodes_to_update))

        with arcpy.da.UpdateCursor(full_path, ["NODE_ID"],
                                   whereclause) as cursor:
            for row in cursor:
                cursor.deleteRow()

    # Insert new rows
    with arcpy.da.InsertCursor(full_path, ["SHAPE@X", "SHAPE@Y"] +
                               fields_to_update) as cursor:
        for row in point_list:
            cursor.insertRow(row)


def get_crs_arcpy(path):
    """Get the CRS from a vector or raster.

    Parameters:
        path (str): Path to the file.

    Returns:
        arcpy.SpatialReference: The coordinate reference system.
    """


    desc = arcpy.Describe(path)
    return desc.spatialReference


def crs_equal_arcpy(fc1, fc2):
    """Compare CRS names for two feature class vector or raster paths using arcpy."""
    crs_a = get_crs_arcpy(fc1)
    crs_b = get_crs_arcpy(fc2)
    return crs_a.name == crs_b.name


def get_field_info_arcpy(path, field_name):
    """Get the data type for a feature class attribute field using arcpy.

    Parameters:
        path (str): Path to the feature class.
        field_name (str): Name of the field.

    Returns:
        str: The arcpy field type (e.g. "Double", "String", "Integer").
    """


    fields = arcpy.ListFields(path, field_name)
    if fields:
        return fields[0].type
    raise ValueError(f"Field '{field_name}' not found in {path}")


# ---------------------------------------------------------------------------
# Raster
# ---------------------------------------------------------------------------

def get_raster_info_arcpy(raster_path):
    """Get raster metadata using arcpy.

    Parameters:
        raster_path (str): Path to the raster file.

    Returns:
        dict: Keys: 'x_cellsize', 'y_cellsize', 'x_min', 'y_min',
            'x_max', 'y_max', 'crs', 'nodata'.
    """


    desc = arcpy.Describe(raster_path)
    sr = desc.spatialReference
    extent = desc.extent

    raster = arcpy.Raster(raster_path)

    return {
        "x_cellsize": raster.meanCellWidth,
        "y_cellsize": raster.meanCellHeight,
        "x_min": extent.XMin,
        "y_min": extent.YMin,
        "x_max": extent.XMax,
        "y_max": extent.YMax,
        "crs": sr,
        "nodata": raster.noDataValue,
    }


def raster_to_array_arcpy(raster_path, x_min, y_min, x_max, y_max,
                            nodata_value=-9999):
    """Read a raster into a numpy array for a specified block extent.

    Uses arcpy.RasterToNumPyArray with a lower-left origin point
    and column/row counts derived from the bounding box.

    Parameters:
        raster_path (str): Path to the raster file.
        x_min, y_min, x_max, y_max (float): Bounding box of the block
            in map coordinates.
        nodata_value (float): Value to use for nodata pixels in the
            output array.

    Returns:
        tuple: (array, x_min_actual, y_max_actual, x_cellsize, y_cellsize)
    """


    raster = arcpy.Raster(raster_path)
    x_cellsize = raster.meanCellWidth
    y_cellsize = raster.meanCellHeight

    # Clamp to raster extent
    desc = arcpy.Describe(raster_path)
    extent = desc.extent
    x_min = max(x_min, extent.XMin)
    y_min = max(y_min, extent.YMin)
    x_max = min(x_max, extent.XMax)
    y_max = min(y_max, extent.YMax)

    # Calculate ncols and nrows
    ncols = int((x_max - x_min) / x_cellsize)
    nrows = int((y_max - y_min) / y_cellsize)

    if ncols <= 0 or nrows <= 0:
        array = np.full((max(1, nrows), max(1, ncols)),
                        nodata_value, dtype=np.float64)
        return (array, x_min, y_max, x_cellsize, y_cellsize)

    # arcpy.RasterToNumPyArray takes the lower-left corner
    lower_left = arcpy.Point(x_min, y_min)
    array = arcpy.RasterToNumPyArray(raster_path, lower_left,
                                     ncols, nrows, nodata_value)

    # y_max_actual is the top of the block
    y_max_actual = y_min + nrows * y_cellsize

    return (array, x_min, y_max_actual, x_cellsize, y_cellsize)


def sample_raster_at_point_arcpy(raster_path, x, y):
    """Sample a raster value at a single x/y coordinate using arcpy.

    Uses arcpy.GetCellValue_management to sample the raster.

    Parameters:
        raster_path (str): Path to the raster file.
        x, y (float): Coordinates in the raster's CRS.

    Returns:
        float: The raster value at (x, y), or None if nodata.
    """


    result = arcpy.GetCellValue_management(
        raster_path, "{0} {1}".format(x, y), "1")
    val_str = result.getOutput(0)

    if val_str == "NoData" or val_str == "":
        return None
    return float(val_str)


# ---------------------------------------------------------------------------
# Spatial processing
# ---------------------------------------------------------------------------

def read_streamline_features_arcpy(streamline_fc, sid_field):
    """Read streamline features and return stream IDs, geometries,
    and lengths.

    Parameters:
        streamline_fc (str): Path to the stream centerline polyline
            feature class.
        sid_field (str): Name of the attribute field holding the unique
            stream identifier.

    Returns:
        list of tuple: Each tuple is (streamID, geometry, length) where
            geometry is an arcpy Polyline and length is in the units
            of projection.
    """
    proj = arcpy.Describe(streamline_fc).spatialReference
    incursorFields = ["SHAPE@", "SHAPE@LENGTH", sid_field]
    features = []
    with arcpy.da.SearchCursor(streamline_fc, incursorFields, "", proj) as Inrows:
        for row in Inrows:
            features.append((row[2], row[0], row[1]))
    return features


def position_along_line_arcpy(geom, fraction):
    """Return the (x, y) coordinates at a fractional position along a line.

    Parameters:
        geom (arcpy.Polyline): The line geometry.
        fraction (float): Fractional position along the line (0.0 to 1.0).

    Returns:
        tuple: (x, y) coordinates.
    """
    pt = geom.positionAlongLine(fraction, True).centroid
    return (pt.X, pt.Y)


def transform_to_latlong_arcpy(x_list, y_list, proj):
    """Transform coordinates from a projected CRS to decimal degrees (EPSG:4326).

    Parameters:
        x_list (list of float): X coordinates in the projected CRS.
        y_list (list of float): Y coordinates in the projected CRS.
        proj (arcpy.SpatialReference): The source coordinate reference system.

    Returns:
        tuple: (longs, lats) lists of decimal degree coordinates.
    """
    
    # This selects the best geographic datum transformation for the
    # input projection and extent. ListTransformations returns
    # all the workable transformations with the best match first.
    # when the nodes and target share same datum (e.g. WGS 84)
    # the list is empty and no transformation needed.
    proj_dd = arcpy.SpatialReference(4326)  # GCS_WGS_1984
    extent = arcpy.Extent(min(x_list), min(y_list),
                          max(x_list), max(y_list),
                          spatial_reference=proj)
    transforms = arcpy.ListTransformations(proj, proj_dd, extent)
    transform_name = transforms[0] if transforms else ""

    longs = []
    lats = []
    for x, y in zip(x_list, y_list):
        pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), proj)
        pt_dd = pt_geom.projectAs(proj_dd, transform_name)
        longs.append(pt_dd.centroid.X)
        lats.append(pt_dd.centroid.Y)
    return (longs, lats)


def read_polyline_geometry_arcpy(polyline_fc):
    """Reads an input polyline into an arcpy polyline geometry object.

    Parameters:
        polyline_fc (str): Path to the polyline feature class.

    Returns:
        arcpy.Polyline: A single polyline geometry representing all features.
    """
    proj_polyline = arcpy.Describe(polyline_fc).spatialReference
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


def calc_channel_width_arcpy(node_x, node_y, bank_geom, aspect, line_dis,
                             proj):
    """Calculate the distance from a node to a bank along a
    perpendicular transect.

    Parameters:
        node_x (float): Node X coordinate.
        node_y (float): Node Y coordinate.
        bank_geom (arcpy.Polyline): Bank geometry object.
        aspect (float): Transect direction in degrees.
        line_dis (float): Maximum search distance in spatial units.
        proj (arcpy.SpatialReference): Coordinate reference system.

    Returns:
        float: Distance from node to bank in spatial units.
    """
    node_geom = arcpy.PointGeometry(arcpy.Point(node_x, node_y), proj)
    pt1 = node_geom.pointFromAngleAndDistance(aspect, line_dis, "PLANAR")
    line = arcpy.Polyline(arcpy.Array([node_geom.centroid, pt1.centroid]), proj)
    pt2 = line.intersect(bank_geom, 1)

    if pt2.centroid:
        to_bank_distance = node_geom.distanceTo(pt2)
    else:
        # No intersection = no bank 90 deg from aspect
        # Find the minimum distance to bank, regardless of the angle
        near_distance = bank_geom.queryPointAndDistance(node_geom.centroid)
        to_bank_distance = near_distance[2]

    return (to_bank_distance)
