"""
Shared utility functions for TTools.

These functions were duplicated across multiple TTools scripts
and are consolidated here.
"""

from math import ceil
import numpy as np
import arcpy


def message(msg):
    """Print a message."""
    try:
        arcpy.AddMessage(msg)
    except Exception:
        print(msg)


def warning(msg):
    """Print a warning message."""
    try:
        arcpy.AddWarning(msg)
    except Exception:
        print(msg)


# -----------------------------------------------------------------------
# ArcPy feature class and raster helpers
# -----------------------------------------------------------------------

def read_fc(nodes_fc, addFields=None):
    """Reads the input point feature class and returns a dictionary
    keyed by NODE_ID."""

    # Get all field names except geometry and OID fields
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


def update_fc(nodeDict, nodes_fc, addFields, nodes_to_query=None):
    """Updates the input point feature class with data from the
    nodes dictionary."""

    # Add fields if they do not exist
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


def get_raster_info(raster_path):
    """Get raster metadata using arcpy."""
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


def raster_to_array(raster_path, block_x_min_corner, block_y_min_corner,
                    block_x_max_corner, block_y_max_corner,
                    nodata_value=-9999):
    """Construct the raster array from the lower left block corner coordinates."""
    raster = arcpy.Raster(raster_path)
    x_cellsize = raster.meanCellWidth
    y_cellsize = raster.meanCellHeight

    # calculate the number of cols/rows from the lower left
    ncols = int((block_x_max_corner - block_x_min_corner) / x_cellsize)
    nrows = int((block_y_max_corner - block_y_min_corner) / y_cellsize)

    lower_left = arcpy.Point(block_x_min_corner, block_y_min_corner)
    array = arcpy.RasterToNumPyArray(raster_path, lower_left,
                                     ncols, nrows, nodata_value)

    return (array, block_x_min_corner, block_y_max_corner,
            x_cellsize, y_cellsize)


# -----------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------

def to_meters_con(inFeature):
    """Returns the conversion factor to get from the
    input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except Exception:
        raise ValueError(
            "{0} has a coordinate system ".format(inFeature) +
            "that is not projected or not recognized. Use a projected "
            "coordinate system preferably in linear units of feet or meters."
        )
    return con_to_m


def from_meters_con(inFeature):
    """Returns the conversion factor to get from meters to the
    spatial units of the input feature class"""
    try:
        con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except Exception:
        raise ValueError(
            "{0} has a coordinate system ".format(inFeature) +
            "that is not projected or not recognized. Use a projected "
            "coordinate system preferably in linear units of feet or meters."
        )
    return con_from_m


def from_z_units_to_meters_con(zUnits):
    """Returns the conversion factor to get from the input z
    units to meters"""
    try:
        con_z_to_m = float(zUnits)
    except (ValueError, TypeError):
        if zUnits == "Meters":
            con_z_to_m = 1.0
        elif zUnits == "Feet":
            con_z_to_m = 0.3048
        else:
            con_z_to_m = None
    return con_z_to_m


def coord_to_array(easting, northing, block_x_min, block_y_max,
                   x_cellsize, y_cellsize):
    """converts x/y coordinates to col and row of the array"""
    col_x = int(((easting - block_x_min) - ((easting - block_x_min) % x_cellsize)) / x_cellsize)
    row_y = int(((block_y_max - northing) - ((block_y_max - northing) % y_cellsize)) / y_cellsize)
    return [col_x, row_y]
