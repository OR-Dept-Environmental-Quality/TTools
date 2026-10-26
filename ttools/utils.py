"""
Shared utility functions for TTools.
"""

from ttools.geo_package import to_meters_con


# -----------------------------------------------------------------------
# Messages
# -----------------------------------------------------------------------

_message_method = None
_warning_method = None


def set_message_method(method):
    """Set which method handles messages (e.g. arcpy.AddMessage)."""
    global _message_method
    _message_method = method


def set_warning_method(method):
    """Set which method prints warning messages (e.g. arcpy.AddWarning)."""
    global _warning_method
    _warning_method = method


def message(msg):
    """Print a message."""
    print(msg)
    if _message_method:
        _message_method(msg)


def warning(msg):
    """Print a warning message."""
    print(msg)
    if _warning_method:
        _warning_method(msg)


# -----------------------------------------------------------------------
# Shared utility functions
# -----------------------------------------------------------------------

def from_meters_con(inFeature):
    """Returns the conversion factor to get from meters to the
    spatial units of the input feature class"""
    return 1.0 / to_meters_con(inFeature)


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


def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize, raster_x_min, raster_y_max):
    """converts x/y coordinates to col and row of the array"""
    dx = easting - raster_x_min
    dy = raster_y_max - northing
    col_x = round((dx - (dx % x_cellsize)) / x_cellsize) - round((block_x_min - raster_x_min) / x_cellsize)
    row_y = round((dy - (dy % y_cellsize)) / y_cellsize) - round((raster_y_max - block_y_max) / y_cellsize)
    return [col_x, row_y]
