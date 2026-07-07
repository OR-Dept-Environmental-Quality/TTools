"""
Geo package router

This routes function calls from the various TTools steps and utils to either
the open source packages (geopandas/rasterio/fiona/shapely/pyproj) or ESRI's arcpy package
methods. Default is open source. Using arcpy requires a valid ESRI ArcPro license.

Call set_geo_package(use_arcpy=True) to use arcpy.
set_geo_package is typically done once at import by the Arc toolbox.
"""

_use_arcpy = False


def set_geo_package(use_arcpy=False):
    """Sets the core geo package/s to use for TTools.
    When use_arcpy is True, all functions use ESRI's arcpy package for
    geoprocessing. Using arcpy requires a valid ESRI ArcPro license.

    When use_arcpy is False (default), the open source packages are
    used (geopandas/rasterio/fiona/shapely/pyproj)."""

    global _use_arcpy
    _use_arcpy = use_arcpy


# ---------------------------------------------------------------------------
# Vector Read/Write
# ---------------------------------------------------------------------------

def fc_exists(path):
    """Check whether a feature class path exists.

    Checks .shp, .gpkg/layer, or .gdb/layer paths by checking that the
    container file exists and the layer is present.

    Parameters:
        path (str): Path to the feature class.

    Returns:
        bool: True if the feature class exists.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import fc_exists_arcpy
        return fc_exists_arcpy(path)
    else:
        from ttools.open_methods import fc_exists_open
        return fc_exists_open(path)


def delete_fc(path):
    """Delete a feature class."""
    if _use_arcpy:
        from ttools.arcpy_methods import delete_fc_arcpy
        return delete_fc_arcpy(path)
    else:
        from ttools.open_methods import delete_fc_open
        return delete_fc_open(path)


def read_fc(nodes_fc, addFields=None):
    """Reads the input point feature class and returns a dictionary
    with NODE_ID used as the dictionary key.

    Parameters:
        nodes_fc (str): Path to the feature class.
        addFields (list of str, optional): Subset of attribute fields to read.
            If None all attribute fields are read. Geometry coordinates
            (POINT_X, POINT_Y) are always included.

    Returns:
        dict: {nodeID: {field: value, ..., "POINT_X": x, "POINT_Y": y}}
    """
    if _use_arcpy:
        from ttools.arcpy_methods import read_fc_arcpy
        return read_fc_arcpy(nodes_fc, addFields)
    else:
        from ttools.open_methods import read_fc_open
        return read_fc_open(nodes_fc, addFields)


def write_fc(nodeDict, nodes_fc, addFields, proj):
    """Create the output point feature class using the data from
    the nodes dictionary.

    NODE_ID is always written as a field. Geometry is created
    from POINT_X and POINT_Y values in the nodeDict.

    Parameters:
        nodeDict (dict): {nodeID: {field: value, ..., "POINT_X": x, "POINT_Y": y}}
        nodes_fc (str): Output feature class path.
        addFields (list of str): Attribute field names to write (not including
            NODE_ID which is always written).
        proj (CRS): Coordinate reference system for the output.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import write_fc_arcpy
        return write_fc_arcpy(nodeDict, nodes_fc, addFields, proj)
    else:
        from ttools.open_methods import write_fc_open
        return write_fc_open(nodeDict, nodes_fc, addFields, proj)


def update_fc(nodeDict, nodes_fc, addFields, nodes_to_query=None):
    """Updates the input point feature class with data from the
    nodes dictionary.

    NODE_ID is used as the key field.

    Parameters:
        nodeDict (dict): {nodeID: {field: value, ...}}
        nodes_fc (str): Path to the existing feature class.
        addFields (list of str): Field names to update.
        nodes_to_query (list, optional): List of NODE_IDs to update. If None,
            updates all rows where NODE_ID is found in nodeDict.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import update_fc_arcpy
        return update_fc_arcpy(nodeDict, nodes_fc, addFields, nodes_to_query)
    else:
        from ttools.open_methods import update_fc_open
        return update_fc_open(nodeDict, nodes_fc, addFields, nodes_to_query)


def update_sample_fc(point_list, sample_fc, fields_to_update,
                     nodes_to_update, overwrite_data, proj):
    """Creates/updates the output sample point feature class using
    the data from the point list.

    Used for topo_fc (Step 4) and lc_point_fc (Step 5).
    Each row in point_list has geometry in [0:2] and attribute
    values in [2:] mapping to fields_to_update.

    Parameters:
        point_list (list of list): Each row: [x, y, field1_val, field2_val, ...].
        sample_fc (str): Path to the output sample feature class.
        fields_to_update (list of str): Field names for the attribute data.
        nodes_to_update (list): List of NODE_IDs being updated.
        overwrite_data (bool): If True, overwrite all existing data.
        proj (CRS): Coordinate reference system for the output.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import update_sample_fc_arcpy
        return update_sample_fc_arcpy(point_list, sample_fc, fields_to_update,
                                      nodes_to_update, overwrite_data, proj)
    else:
        from ttools.open_methods import update_sample_fc_open
        return update_sample_fc_open(point_list, sample_fc, fields_to_update,
                                     nodes_to_update, overwrite_data, proj)


# ---------------------------------------------------------------------------
# Spatial processing
# ---------------------------------------------------------------------------

def read_streamline_features(streamline_fc, sid_field):
    """Read streamline features and return stream IDs, geometries,
    and lengths.

    Parameters:
        streamline_fc (str): Path to the stream centerline polyline
            feature class.
        sid_field (str): Name of the attribute field holding the unique
            stream identifier.

    Returns:
        list of tuple: Each tuple is (streamID, geometry, length).
    """
    if _use_arcpy:
        from ttools.arcpy_methods import read_streamline_features_arcpy
        return read_streamline_features_arcpy(streamline_fc, sid_field)
    else:
        from ttools.open_methods import read_streamline_features_open
        return read_streamline_features_open(streamline_fc, sid_field)


def position_along_line(geom, fraction):
    """Return the (x, y) coordinates at a fractional position along a line.

    Parameters:
        geom: The line geometry (Shapely or arcpy).
        fraction (float): Fractional position along the line (0.0 to 1.0).

    Returns:
        tuple: (x, y) coordinates.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import position_along_line_arcpy
        return position_along_line_arcpy(geom, fraction)
    else:
        from ttools.open_methods import position_along_line_open
        return position_along_line_open(geom, fraction)


def transform_to_latlong(x_list, y_list, proj):
    """Transform coordinates from a projected CRS to decimal degrees.

    Parameters:
        x_list (list of float): X coordinates in the projected CRS.
        y_list (list of float): Y coordinates in the projected CRS.
        proj (CRS): The source coordinate reference system.

    Returns:
        tuple: (lons, lats) lists of decimal degree coordinates.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import transform_to_latlong_arcpy
        return transform_to_latlong_arcpy(x_list, y_list, proj)
    else:
        from ttools.open_methods import transform_to_latlong_open
        return transform_to_latlong_open(x_list, y_list, proj)


def read_polyline_geometry(polyline_fc):
    """Reads an input polyline into a single geometry object.

    Parameters:
        polyline_fc (str): Path to the polyline feature class.

    Returns:
        geometry object: A single geometry representing all polyline features.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import read_polyline_geometry_arcpy
        return read_polyline_geometry_arcpy(polyline_fc)
    else:
        from ttools.open_methods import read_polyline_geometry_open
        return read_polyline_geometry_open(polyline_fc)


def calc_channel_width(node_x, node_y, bank_geom, aspect, line_dis, proj):
    """Calculate the distance from a node to a bank along a
    perpendicular transect.

    Parameters:
        node_x (float): Node X coordinate.
        node_y (float): Node Y coordinate.
        bank_geom: Bank geometry object.
        aspect (float): Transect direction in degrees.
        line_dis (float): Maximum search distance in spatial units.
        proj (CRS): Coordinate reference system.

    Returns:
        float: Distance from node to bank in spatial units.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import calc_channel_width_arcpy
        return calc_channel_width_arcpy(node_x, node_y, bank_geom, aspect,
                                        line_dis, proj)
    else:
        from ttools.open_methods import calc_channel_width_open
        return calc_channel_width_open(node_x, node_y, bank_geom, aspect,
                                       line_dis, proj)


# ---------------------------------------------------------------------------
# CRS and field info
# ---------------------------------------------------------------------------

def get_crs(path):
    """Get the CRS from a vector or raster file.

    Parameters:
        path (str): Path to the file.

    Returns:
        CRS: The coordinate reference system.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import get_crs_arcpy
        return get_crs_arcpy(path)
    else:
        from ttools.open_methods import get_crs_open
        return get_crs_open(path)


def crs_equal(fc1, fc2):
    """Compare CRS names for two feature class vector or raster paths using selected geo package."""
    if _use_arcpy:
        from ttools.arcpy_methods import crs_equal_arcpy
        return crs_equal_arcpy(fc1, fc2)
    else:
        from ttools.open_methods import crs_equal_open
        return crs_equal_open(fc1, fc2)


def get_field_info(path, field_name):
    """Get the data type for a feature class attribute field.

    Parameters:
        path (str): Path to the feature class.
        field_name (str): Name of the field.

    Returns:
        str: The dtype of the field.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import get_field_info_arcpy
        return get_field_info_arcpy(path, field_name)
    else:
        from ttools.open_methods import get_field_info_open
        return get_field_info_open(path, field_name)


# ---------------------------------------------------------------------------
# Raster Processing
# ---------------------------------------------------------------------------

def get_raster_info(raster_path):
    """Get raster metadata.

    Parameters:
        raster_path (str): Path to the raster file.

    Returns:
        dict: Keys: 'x_cellsize', 'y_cellsize', 'x_min', 'y_min',
            'x_max', 'y_max', 'crs', 'nodata'.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import get_raster_info_arcpy
        return get_raster_info_arcpy(raster_path)
    else:
        from ttools.open_methods import get_raster_info_open
        return get_raster_info_open(raster_path)


def raster_to_array(raster_path, x_min, y_min, x_max, y_max,
                      nodata_value=-9999):
    """Read a raster into a numpy array for a specified block extent.

    Parameters:
        raster_path (str): Path to the raster file.
        x_min, y_min, x_max, y_max (float): Bounding box of the block
            in map coordinates.
        nodata_value (float): Value to use for nodata pixels in the
            output array.

    Returns:
        tuple: (array, x_min_actual, y_max_actual, x_cellsize, y_cellsize)
    """
    if _use_arcpy:
        from ttools.arcpy_methods import raster_to_array_arcpy
        return raster_to_array_arcpy(raster_path, x_min, y_min, x_max, y_max,
                                       nodata_value)
    else:
        from ttools.open_methods import raster_to_array_open
        return raster_to_array_open(raster_path, x_min, y_min, x_max, y_max,
                                      nodata_value)


def sample_raster_at_point(raster_path, x, y):
    """Sample a raster value at a single x/y coordinate.

    Parameters:
        raster_path (str): Path to the raster file.
        x, y (float): Coordinates in the raster's CRS.

    Returns:
        float: The raster value at (x, y), or None if nodata.
    """
    if _use_arcpy:
        from ttools.arcpy_methods import sample_raster_at_point_arcpy
        return sample_raster_at_point_arcpy(raster_path, x, y)
    else:
        from ttools.open_methods import sample_raster_at_point_open
        return sample_raster_at_point_open(raster_path, x, y)
