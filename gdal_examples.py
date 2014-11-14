# Some example code showing how to use gdal to do various tasks.

# Install gdal from here:
# http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal
# this will install into primary python install

from osgeo import ogr
from osgeo import gdal

# Read ESRI Filegeodatabase
driver = ogr.GetDriverByName("OpenFileGDB")
db = driver.Open(r"D:\Projects\TTools_9\Example_data.gdb", 0)
fc = db.GetLayer("out_nodes")
sr = fc.GetSpatialRef()
sr.ExportToProj4()

