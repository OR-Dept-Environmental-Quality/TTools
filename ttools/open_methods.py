"""
Open source geo methods for various TTools processing tasks using
pandas, geopandas, shapely, fiona, numpy, rasterio, and pyproj. Underlying
many of these methods is GDAL which is packaged with fiona. These methods are called
by geo_package.py when use_arcpy is False (the default).
"""

import os
from math import radians, sin, cos

import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio
import rasterio.windows
from shapely.geometry import Point, LineString
from shapely.ops import nearest_points, unary_union
from pyproj import CRS, Transformer


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _get_driver(path):
    """Determine the GDAL/Fiona driver from the file path extension.

    Parameters:
        path (str): File path.

    Returns:
        str: Driver name for Fiona/GeoPandas.
    """
    ext = os.path.splitext(path)[-1].lower()

    if ext == ".shp":
        return "ESRI Shapefile"
    elif ext == ".gpkg":
        return "GPKG"
    elif ext == ".gdb" or ".gdb" in path:
        return "OpenFileGDB"
    else:
        raise ValueError(
            f"Unsupported file format: {ext}. "
            "Use .shp, .gpkg, or .gdb"
        )


def _get_fc_path(path):
    """Get the path to the feature class which may be inside a .gdb or .gpkg.

    Handles paths like:
        C:/data/my.gdb/layer_name
        C:/data/my.gpkg/layer_name

    Parameters:
        path (str): File path, possibly including a layer name after
            .gdb or .gpkg.

    Returns:
        tuple: (file_path, layer_name) where layer_name may be None.
    """
    # Clean up the path separators
    path_clean = path.replace("\\", "/")

    # Check for .gdb or .gpkg
    for ext in [".gdb", ".gpkg"]:
        ext_idx = path_clean.lower().find(ext)
        if ext_idx != -1:
            container_path = path_clean[:ext_idx + len(ext)]
            remainder = path_clean[ext_idx + len(ext):].strip("/")
            layer_name = remainder if remainder else None
            return container_path, layer_name

    return path, None


def _open_raster(raster_path):
    """Open a raster path, including rasters in .gdb or .gpkg."""
    file_path, layer_name = _get_fc_path(raster_path)

    if layer_name and file_path.lower().endswith(".gdb"):
        open_path = 'OpenFileGDB:"{0}":{1}'.format(file_path, layer_name)
        return rasterio.open(open_path)

    if layer_name and file_path.lower().endswith(".gpkg"):
        return rasterio.open(file_path, TABLE=layer_name)

    return rasterio.open(raster_path)


# ---------------------------------------------------------------------------
# Vector
# ---------------------------------------------------------------------------

def fc_exists_open(path):
    """Checks whether a feature class exists.

    Parameters:
        path (str): Path to the feature class.

    Returns:
        bool: True if the feature class exists.
    """
    import fiona

    file_path, layer_name = _get_fc_path(path)
    if not os.path.exists(file_path):
        return False
    if layer_name is not None:
        return layer_name in fiona.listlayers(file_path)
    return True


def delete_fc_open(path):
    """Delete a feature class using fiona."""
    import fiona

    file_path, layer_name = _get_fc_path(path)
    driver = _get_driver(file_path)
    if layer_name is not None:
        fiona.remove(file_path, driver=driver, layer=layer_name)
    else:
        fiona.remove(file_path, driver=driver)


def read_fc_open(nodes_fc, addFields=None):
    """Reads the input point feature class and returns a dictionary
    with NODE_ID used as the dictionary key.

    Parameters:
        nodes_fc (str): Path to the nodes feature class.
        addFields (list of str, optional): Subset of attribute fields to read.
            If None all attribute fields are read (except NODE_ID which
            becomes the key). Geometry coordinates (POINT_X, POINT_Y) are
            always included.

    Returns:
        dict: {nodeID: {field: value, ..., "POINT_X": x, "POINT_Y": y}}
    """
    file_path, layer_name = _get_fc_path(nodes_fc)

    if layer_name:
        gdf = gpd.read_file(file_path, layer=layer_name)
    else:
        gdf = gpd.read_file(nodes_fc)

    nodeDict = {}

    if addFields:
        read_fields = [f for f in addFields if f in gdf.columns]
    else:
        read_fields = [c for c in gdf.columns
                       if c not in ("geometry", "NODE_ID")]

    for idx, row in gdf.iterrows():
        nodeID = row["NODE_ID"]
        nodeDict[nodeID] = {}
        for f in read_fields:
            nodeDict[nodeID][f] = row[f]
        if row.geometry is not None:
            nodeDict[nodeID]["POINT_X"] = row.geometry.x
            nodeDict[nodeID]["POINT_Y"] = row.geometry.y

    return nodeDict


def write_fc_open(nodeDict, nodes_fc, addFields, proj):
    """Create the output point feature class using the data from
    the nodes dictionary.

    NODE_ID is always written as a field. Geometry is created
    from POINT_X and POINT_Y.

    Parameters:
        nodeDict (dict): {nodeID: {field: value, ..., "POINT_X": x,
            "POINT_Y": y}}
        nodes_fc (str): Output feature class path.
        addFields (list of str): Attribute field names to write (not
            including NODE_ID which is always written).
        proj (CRS): Coordinate reference system for the output.
    """
    data = []
    for nodeID in nodeDict:
        row_data = {"NODE_ID": nodeID}
        for f in addFields:
            row_data[f] = nodeDict[nodeID].get(f)
        row_data["geometry"] = Point(nodeDict[nodeID]["POINT_X"],
                                     nodeDict[nodeID]["POINT_Y"])
        data.append(row_data)

    gdf = gpd.GeoDataFrame(data, crs=proj)

    file_path, layer_name = _get_fc_path(nodes_fc)
    driver = _get_driver(file_path)

    if driver == "OpenFileGDB" and layer_name is None:
        raise ValueError(
            "A layer name is required when writing to a .gdb. "
            "Either append it to the path (e.g. 'output.gdb/my_layer') "
            "or pass layer_name='my_layer'."
        )

    kwargs = {"driver": driver}
    if layer_name:
        kwargs["layer"] = layer_name

    gdf["NODE_ID"] = gdf["NODE_ID"].astype("int32")
    gdf.to_file(file_path, **kwargs)


def update_fc_open(nodeDict, nodes_fc, addFields, nodes_to_query=None):
    """Updates the input point feature class with data from the
    nodes dictionary.

    Reads the existing data, updates the specified fields from
    nodeDict, and writes back. NODE_ID is used as the key field.

    Parameters:
        nodeDict (dict): {nodeID: {field: value, ...}}
        nodes_fc (str): Path to the existing feature class.
        addFields (list of str): Field names to update.
        nodes_to_query (list, optional): List of NODE_IDs to update.
            If None, updates all rows where NODE_ID is found in nodeDict.
    """
    file_path, layer_name = _get_fc_path(nodes_fc)

    if layer_name:
        gdf = gpd.read_file(file_path, layer=layer_name)
    else:
        gdf = gpd.read_file(nodes_fc)

    # Add fields if they don't exist
    missing_fields = [f for f in addFields if f not in gdf.columns]
    new_fields = {}
    for f in missing_fields:
        if f.startswith("LC"):
            new_fields[f] = pd.Series("", index=gdf.index)
        else:
            new_fields[f] = pd.Series(np.nan, index=gdf.index)

    if new_fields:
        gdf = gpd.GeoDataFrame(
            pd.concat([gdf, pd.DataFrame(new_fields)], axis=1),
            geometry=gdf.geometry.name, crs=gdf.crs)

    string_fields = set()
    for field in addFields:
        if (hasattr(gdf[field], 'dtype') and
                str(gdf[field].dtype) in ("object", "str", "string")):
            string_fields.add(field)

    # Build lookup set
    if nodes_to_query:
        query_set = set(nodes_to_query)
    else:
        query_set = set(nodeDict.keys())

    rows = []
    present_rows = []
    index = []
    for nodeID in query_set:
        if nodeID in nodeDict:
            row = []
            present_row = []
            for field in addFields:
                if field in nodeDict[nodeID]:
                    val = nodeDict[nodeID][field]
                    # Convert to string if the column has string dtype
                    if field in string_fields:
                        val = str(val) if val is not None else ""
                    present_row.append(True)
                else:
                    val = np.nan
                    present_row.append(False)
                row.append(val)
            rows.append(row)
            present_rows.append(present_row)
            index.append(nodeID)

    if rows:
        update_df = pd.DataFrame(rows, index=index, columns=addFields)
        present_df = pd.DataFrame(present_rows, index=index,
                                  columns=addFields)
        target = gdf["NODE_ID"].isin(index)
        ordered = update_df.reindex(gdf.loc[target, "NODE_ID"])
        present = present_df.reindex(gdf.loc[target, "NODE_ID"])
        current = gdf.loc[target, addFields]
        values = np.where(present.to_numpy(), ordered.to_numpy(),
                          current.to_numpy())
        gdf.loc[target, addFields] = values

    # Write back
    driver = _get_driver(file_path)
    kwargs = {"driver": driver}
    if layer_name:
        kwargs["layer"] = layer_name

    gdf["NODE_ID"] = gdf["NODE_ID"].astype("int32")
    gdf.to_file(file_path, **kwargs)


def update_sample_fc_open(point_list, sample_fc, fields_to_update,
                          nodes_to_update, overwrite_data, proj):
    """Creates/updates the output sample point feature class using
    the data from the point list.

    Creates the feature class if it doesn't exist, optionally
    removes existing rows for the specified nodes, and adds
    new rows.

    Parameters:
        point_list (list of list): Each row: [x, y, field1_val, field2_val, ...].
            row[0:2] are geometry coordinates. row[2:] map to
            fields_to_update in order.
        sample_fc (str): Path to the output sample feature class.
        fields_to_update (list of str): Field names for the attribute data
            (maps to row[2:]).
        nodes_to_update (list): List of NODE_IDs being updated.
        overwrite_data (bool): If True, overwrite all existing data. If False,
            remove existing rows for nodes_to_update before adding new ones.
        proj (CRS): Coordinate reference system for the output.
    """
    # Build GeoDataFrame from point_list
    data = []
    for row in point_list:
        d = {"geometry": Point(row[0], row[1])}
        for i, f in enumerate(fields_to_update):
            d[f] = row[i + 2]
        data.append(d)

    new_gdf = gpd.GeoDataFrame(data, crs=proj)

    file_path, layer_name = _get_fc_path(sample_fc)

    # If not overwriting and the file exists, merge with existing data
    if not overwrite_data and fc_exists_open(sample_fc):
        if layer_name:
            existing_gdf = gpd.read_file(file_path, layer=layer_name)
        else:
            existing_gdf = gpd.read_file(sample_fc)

        # Remove rows for nodes being updated
        update_set = set(nodes_to_update)
        existing_gdf = existing_gdf[~existing_gdf["NODE_ID"].isin(update_set)]

        new_gdf = gpd.GeoDataFrame(
            pd.concat([existing_gdf, new_gdf], ignore_index=True),
            crs=proj
        )

    driver = _get_driver(file_path)
    kwargs = {"driver": driver}
    if layer_name:
        kwargs["layer"] = layer_name

    typeDict = {
        "NODE_ID": "int32",
        "SAMPLE_ID": "int32",
        "TRANSECT": "int16",
        "SAMPLE": "int16",
        "NA_SAMPLES": "int16",
        "AZIMUTH": "float64",
        "SEARCHDIS": "float64",
    }

    # Use type in typeDict, or string for LC, and float64 for all other fields.
    for field in fields_to_update:
        if field in typeDict:
            new_gdf[field] = new_gdf[field].astype(typeDict[field])
        elif field == "STREAM_ID":
            if str(new_gdf[field].dtype) in ("object", "str", "string"):
                new_gdf[field] = new_gdf[field].astype(str)
            else:
                new_gdf[field] = new_gdf[field].astype("int32")
        elif field.startswith("LC") or field == "KEY":
            new_gdf[field] = new_gdf[field].astype(str)
        else:
            new_gdf[field] = new_gdf[field].astype("float64")

    new_gdf.to_file(file_path, **kwargs)


def get_crs_open(path):
    """Get the CRS from a vector or raster file.

    Parameters:
        path (str): Path to the file.

    Returns:
        CRS: The coordinate reference system.
    """
    gdb_path, layer_name = _get_fc_path(path)

    # Try as vector first
    try:
        if layer_name:
            gdf = gpd.read_file(gdb_path, layer=layer_name, rows=0)
        else:
            gdf = gpd.read_file(path, rows=0)
        return CRS(gdf.crs)
    except Exception:
        pass

    # Try as raster
    try:
        with _open_raster(path) as src:
            return CRS(src.crs)
    except Exception:
        pass

    raise ValueError(f"Could not determine CRS from: {path}")


def crs_equal_open(fc1, fc2):
    """Compare CRS names for two feature class vector or raster paths using open source packages."""
    crs_a = get_crs_open(fc1)
    crs_b = get_crs_open(fc2)
    return crs_a.equals(crs_b)


def get_field_info_open(path, field_name):
    """Get the data type for a feature class attribute field.

    Parameters:
        path (str): Path to the feature class.
        field_name (str): Name of the field.

    Returns:
        str: The numpy/pandas dtype of the field.
    """
    gdf = read_fc_open(path)
    if field_name in gdf.columns:
        return str(gdf[field_name].dtype)
    raise ValueError(f"Field '{field_name}' not found in {path}")


# ---------------------------------------------------------------------------
# Raster
# ---------------------------------------------------------------------------

def get_raster_info_open(raster_path):
    """Get raster extent, cell size, and crs, and no data value.

    Parameters:
        raster_path (str): Path to the raster file.

    Returns:
        dict: Keys: 'x_cellsize', 'y_cellsize', 'x_min', 'y_min',
            'x_max', 'y_max', 'crs', 'nodata'.
    """
    with _open_raster(raster_path) as src:
        transform = src.transform
        bounds = src.bounds

        return {
            "x_cellsize": abs(transform.a),
            "y_cellsize": abs(transform.e),
            "x_min": bounds.left,
            "y_min": bounds.bottom,
            "x_max": bounds.right,
            "y_max": bounds.top,
            "crs": CRS(src.crs),
            "nodata": src.nodata,
        }


def raster_to_array_open(raster_path, x_min, y_min, x_max, y_max,
                      nodata_value=-9999):
    """Read a raster into a numpy array for a specified block extent.

    Uses rasterio.windows.

    Parameters:
        raster_path (str): Path to the raster file.
        x_min, y_min, x_max, y_max (float): Bounding box of the block
            in map coordinates.
        nodata_value (float): Value to use for nodata pixels in the
            output array.

    Returns:
        tuple: (array, x_min_actual, y_max_actual, x_cellsize, y_cellsize)
            where array is a 2D numpy array (row, col) and the min/max
            values are the actual cell-aligned coordinates of the block.
    """
    with _open_raster(raster_path) as src:
        x_cellsize = abs(src.transform.a)
        y_cellsize = abs(src.transform.e)

        # Clamp to raster extent
        r_bounds = src.bounds
        x_min = max(x_min, r_bounds.left)
        y_min = max(y_min, r_bounds.bottom)
        x_max = min(x_max, r_bounds.right)
        y_max = min(y_max, r_bounds.top)

        # Convert to row/col indices for precise windowing.
        # Use round() for start indices because callers (e.g. step3)
        # pre-align coordinates to cell corners — floating-point
        # arithmetic can land just below the true integer (e.g.
        # 1918.9999999999 instead of 1919), and int() truncation
        # would shift the window by one cell.
        col_start = round((x_min - r_bounds.left) / x_cellsize)
        col_end = int(np.ceil((x_max - r_bounds.left) / x_cellsize))
        row_start = round((r_bounds.top - y_max) / y_cellsize)
        row_end = int(np.ceil((r_bounds.top - y_min) / y_cellsize))

        # Clamp to valid raster dimensions
        col_start = max(0, col_start)
        row_start = max(0, row_start)
        col_end = min(col_end, src.width)
        row_end = min(row_end, src.height)

        # Handle empty blocks (block entirely outside raster)
        n_cols = max(0, col_end - col_start)
        n_rows = max(0, row_end - row_start)

        # Compute the actual aligned coordinates of the block
        x_min_aligned = r_bounds.left + col_start * x_cellsize
        y_max_aligned = r_bounds.top - row_start * y_cellsize

        if n_cols == 0 or n_rows == 0:
            array = np.full((max(1, n_rows), max(1, n_cols)),
                            nodata_value, dtype=np.float64)
            return (array, x_min_aligned, y_max_aligned,
                    x_cellsize, y_cellsize)

        window = rasterio.windows.Window(
            col_start, row_start, n_cols, n_rows)

        array = src.read(1, window=window)

        # Replace nodata
        if src.nodata is not None:
            array = array.astype(np.float64)
            array[array == src.nodata] = nodata_value

        return (array, x_min_aligned, y_max_aligned,
                x_cellsize, y_cellsize)


def sample_raster_at_point_open(raster_path, x, y):
    """Sample a raster value at a single x/y coordinate.

    Parameters:
        raster_path (str): Path to the raster file.
        x, y (float): Coordinates in the raster's CRS.

    Returns:
        float: The raster value at (x, y), or None if nodata.
    """
    with _open_raster(raster_path) as src:
        row, col = src.index(x, y)
        val = src.read(1, window=((row, row + 1), (col, col + 1)))
        result = float(val[0, 0])

        if src.nodata is not None and result == src.nodata:
            return None
        return result


# ---------------------------------------------------------------------------
# Spatial processing
# ---------------------------------------------------------------------------

def read_streamline_features_open(streamline_fc, sid_field):
    """Read streamline features and return stream IDs, geometries,
    and lengths.

    Parameters:
        streamline_fc (str): Path to the stream centerline polyline
            feature class.
        sid_field (str): Name of the attribute field holding the unique
            stream identifier.

    Returns:
        list of tuple: Each tuple is (streamID, geometry, length) where
            geometry is a Shapely LineString and length is in the units
            of projection.
    """
    file_path, layer_name = _get_fc_path(streamline_fc)
    if layer_name:
        gdf = gpd.read_file(file_path, layer=layer_name)
    else:
        gdf = gpd.read_file(streamline_fc)
    features = []
    for idx, row in gdf.iterrows():
        features.append((row[sid_field], row.geometry, row.geometry.length))
    return features


def position_along_line_open(geom, fraction):
    """Return the (x, y) coordinates at a fractional position along a line.

    Parameters:
        geom (LineString): The line geometry.
        fraction (float): Fractional position along the line (0.0 to 1.0).

    Returns:
        tuple: (x, y) coordinates.
    """
    pt = geom.interpolate(fraction, normalized=True)
    return (pt.x, pt.y)


def transform_to_latlong_open(x_list, y_list, proj):
    """Transform coordinates from a projected CRS to decimal degrees (EPSG:4326).

    Parameters:
        x_list (list of float): X coordinates in the projected CRS.
        y_list (list of float): Y coordinates in the projected CRS.
        proj (CRS): The source coordinate reference system.

    Returns:
        tuple: (longs, lats) lists of decimal degree coordinates.
    """
    proj_dd = CRS("EPSG:4326")  # GCS_WGS_1984
    transformer = Transformer.from_crs(proj, proj_dd, always_xy=True)
    longs, lats = transformer.transform(x_list, y_list)
    return (longs, lats)


def read_polyline_geometry_open(polyline_fc):
    """Reads an input polyline into a single Shapely geometry object.

    Parameters:
        polyline_fc (str): Path to the polyline feature class.

    Returns:
        geometry: A single geometry representing all polyline features.
    """
    file_path, layer_name = _get_fc_path(polyline_fc)
    if layer_name:
        gdf = gpd.read_file(file_path, layer=layer_name)
    else:
        gdf = gpd.read_file(polyline_fc)
    # Union all geometries into a single geometry object
    poly_geom = unary_union(gdf.geometry)
    return poly_geom


def calc_channel_width_open(node_x, node_y, bank_geom, aspect, line_dis,
                            proj):
    """Calculate the distance from a node to a stream bank along a
    perpendicular transect.

    Parameters:
        node_x (float): Node X coordinate.
        node_y (float): Node Y coordinate.
        bank_geom: Bank geometry object.
        aspect (float): Transect direction in degrees.
        line_dis (float): Maximum search distance in spatial units.
        proj: Coordinate reference system (unused in open-source backend).

    Returns:
        float: Distance from node to bank in spatial units.
    """
    # Create a point at the given angle and distance from the node
    pt1_x = node_x + line_dis * sin(radians(aspect))
    pt1_y = node_y + line_dis * cos(radians(aspect))

    line = LineString([(node_x, node_y), (pt1_x, pt1_y)])
    node_pt = Point(node_x, node_y)

    pt2 = line.intersection(bank_geom)

    if not pt2.is_empty:
        # If intersection returns multiple points, use the nearest one
        if pt2.geom_type == 'MultiPoint':
            to_bank_distance = min(node_pt.distance(p) for p in pt2.geoms)
        elif pt2.geom_type == 'Point':
            to_bank_distance = node_pt.distance(pt2)
        else:
            # Intersection is a line segment or other geometry
            to_bank_distance = node_pt.distance(pt2)
    else:
        # No intersection = no bank 90 deg from aspect
        # Find the minimum distance to bank, regardless of the angle
        nearest_pt = nearest_points(node_pt, bank_geom)[1]
        to_bank_distance = node_pt.distance(nearest_pt)

    return (to_bank_distance)
