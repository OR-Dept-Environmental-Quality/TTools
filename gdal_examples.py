# Some example code showing how to use gdal to do various tasks.

# Install gdal for windows from here:
# http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal
# this will install into primary python install

# For mac
# http://www.kyngchaos.com/software/frameworks
# this will install into primary python install

from osgeo import ogr
from osgeo import gdal

# Read ESRI Filegeodatabase
driver = ogr.GetDriverByName("OpenFileGDB")
db = driver.Open(r"D:\Projects\TTools_9\Example_data.gdb", 0)
fc = db.GetLayer("out_nodes")
sr = fc.GetSpatialRef()
sr.ExportToProj4()

# Read a Shapefile, and get field names
source = ogr.Open(r"D:\Projects\TTools_9\out_nodes.shp", 1)
layer = source.GetLayer()
layer_defn = layer.GetLayerDefn()
field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
source = None # Close the Shapefile

# Add a new field
source = ogr.Open(r"D:\Projects\TTools_9\out_nodes.shp", 1)
layer = source.GetLayer()
new_field = ogr.FieldDefn("NewFieldName", ogr.OFTInteger)
layer.CreateField(new_field)
source = None # Close the Shapefile

def RasterToArray(raster):
    """Converts a rasters to a numpy array, returns the geotransform data and number or rows and cols"""

    data = gdal.Open(raster)
    band = data.GetRasterBand(1)
    arry = band.ReadAsArray()
    gt = data.GetGeoTransform()
    proj = data.GetProjection
    rows = data.RasterYSize
    cols = data.RasterYSize
    
    return arry, gt, proj, rows, cols

def CoordToArray(easting, northing, gt, rows):
    '''converts x/y coordinates to col and row of the array'''
    xy = []
    xy.append((easting - gt[0])/ gt[1])  # col
    xy.append((northing - (gt[3] - (rows * gt[1]))) / gt[1])  # row
    
    return xy


