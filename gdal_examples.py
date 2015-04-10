# Some example code showing how to use gdal to do various tasks.

# Install gdal for windows from here:
# http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal
# this will install into primary python install

# For mac
# http://www.kyngchaos.com/software/frameworks
# this will install into primary python install

from __future__ import print_function
from osgeo import ogr
from osgeo import gdal

#raster = r"/Users/rmichie/Desktop/Test/44123C4319_Rasters/veght44123c4319/veght44123c4319.img"
raster = r"D:\Projects\TTools_9\Raster\LiDAR\be\be_lidar"
gdbpath = r"D:\Projects\TTools_9\Example_data.gdb"
fc = "out_nodes"
shpfile =  r"D:\Projects\TTools_9\out_nodes.shp"
newfield = "Test"

# raster value = 149.96002197265625
easting = 551997.077
northing = 945129.128

# raster value = 3.480010986328125
easting = 551979.017
northing = 945125.731

def ReadFileGDB(gdbpath, gdb_fc):
    # Read ESRI Filegeodatabase
    driver = ogr.GetDriverByName("OpenFileGDB")
    db = driver.Open(gdbpath, 0)
    fc = db.GetLayer(gdb_fc)
    sr = fc.GetSpatialRef()
    sr.ExportToProj4()
    return fc, sr

def ReadSHP(shpfile):
    # Read a Shapefile, and get field names
    source = ogr.Open(shpfile, 1)
    layer = source.GetLayer()
    layer_defn = layer.GetLayerDefn()
    field_names = [layer_defn.GetFieldDefn(i).GetName() for i in range(layer_defn.GetFieldCount())]
    source = None # Close the Shapefile
    return layer, field_names


def AddFields(shpfile, AddFieldsList):
    # Add new fields to a shapefile
    source = ogr.Open(shpfile, 1)
    layer = source.GetLayer()
    for field in AddFieldsList:
        new_field = ogr.FieldDefn(field, ogr.OFTInteger)
        layer.CreateField(new_field)
    source = None # Close the Shapefile


def ReadRaster(raster):
    # Read a Raster
    data = gdal.Open(raster)
    return data

def RasterToArray(raster):
    """Converts a rasters to a numpy array, returns the geotransform data and number or rows and cols"""

    data = gdal.Open(raster)
    band = data.GetRasterBand(1)
    arry = band.ReadAsArray(xoff=0, yoff=0, buf_obj=None)
    gt = data.GetGeoTransform()
    proj = data.GetProjection()
    rows = data.RasterYSize
    cols = data.RasterXSize
    return arry, gt, proj, rows, cols

def ReadRasterDirectly(raster, xy):
    '''Reads the rasters directly, returns data in binary and must be unpacked'''
    import struct
    data = gdal.Open(raster)
    band = data.GetRasterBand(1)
    val = band.ReadRaster(xoff= int(xy[0]), yoff= int(xy[1]), xsize=1, ysize=1, buf_type=gdal.GDT_Float32)
    uval = struct.unpack('f' , val) # unpack as float
    return(uval[0])

def CoordToArray(easting, northing, gt):
    '''converts x/y coordinates to col and row of the array'''
    xy = []
    xy.append((easting - gt[0])/ gt[1])  # col
    xy.append((northing - gt[3]) / gt[5])  # row
    return xy


data = ReadRaster(raster)
arry, gt, proj, rows, cols = RasterToArray(raster)
xy = CoordToArray(easting, northing, gt)
print(arry[xy[1], xy[0]])
val = ReadRasterDirectly(raster, xy)
print(val)
print("done")
