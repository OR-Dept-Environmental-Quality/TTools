#######################################################################################
# TTools
# Convert to csv - v 0.1
# Ryan Michie

# This script will take an input point feature and output a csv file in the landcover
# data format for heatsource 9.

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from collections import defaultdict
import csv

# parameter values
# 0: input TTools point file (inPoint)
# 1: output csv file to write to name/path outcsv_final)

def tree(): 
    """Build a nested dictionary"""
    return defaultdict(tree)

def read_pointfile(pointfile, readfields):
    """Reads an input point file and returns the NODE ID and X/Y coordinates as a nested dictionary"""
    pnt_dict = tree()
    Incursorfields = ["STREAM_ID", "NODE_ID", "SHAPE@X","SHAPE@Y"] + readfields
    # Determine input point spatial units
    proj = arcpy.Describe(inPoint).spatialReference
    with arcpy.da.SearchCursor(pointfile,Incursorfields,"",proj) as Inrows:
	for row in Inrows:
	    pnt_dict[row[0]]["POINT_X"] = row[1]
	    pnt_dict[row[0]]["POINT_Y"] = row[2] 
	    for f in xrange(0,len(readfields)):
		pnt_dict[row[0]][readfields[f]] = row[3+f]
    return(pnt_dict)

def read_csv(csvfile):
  """Reads an input csv file and returns the header row as a list and the data as a nested dictionary"""
  csv_dict = tree()
  with open(csvfile, "rb") as f:
    reader = csv.reader(f)
    header = reader.next()
    if header[0] != "NODE_ID":
      sys.exit("csv file does not have NODE_ID as first column")
    for row in reader:
      for key in xrange(0,len(header)):
	csv_dict[row[0]][header[key]] = row[key]		
  return(header,csv_dict)

try:

	#inPoint = arcpy.GetParameterAsText(0)
	#outcsv_final = arcpy.GetParameterAsText(1)
	
	# Start Fill in Data
	inPoint = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
	outcsv_final = "D:\Projects\TTools_9\out_nodes.csv"
	# End Fill in Data
	removelist = [u"OBJECTID",u"Shape",u"NODE_ID"]
	
	# Get all the column headers in the point file and remove the ones in removelist
	header = [field.name for field in arcpy.Describe(inPoint).fields]
	header_clean = [h for h in header if h not in removelist]
	
	NODES = read_pointfile(inPoint, header_clean)

	####################################################################################################### 
	# Output the csv file
	
	# list of Topo headers
	
	csv_headers = ["NODE_ID","STREAM_KM", "LONGITUDE", "LATITUDE", "TOPO_270", "TOPO_180", "TOPO_90"]

	# Get the stream km dictionary keys and sort them
	NODE_keys = NODES.keys()
	NODE_keys.sort()    
	
	# make a wide format list by node ID from the nested dictionary 
	csv_out = [[NODES[nodeID][h] for h in csv_header] for nodeID in NODE_keys]

	# Add the header row
	csv_out.insert(0,csv_header)
	
	# Export to csv
	with open(outcsv_final, "wb") as f:
		writer = csv.writer(f)
		writer.writerows(csv_out3)	



	
	
# For arctool errors
except arcpy.ExecuteError:
	msgs = arcpy.GetMessages(2)
	#arcpy.AddError(msgs)
	print(msgs)
	
# For other errors
except:
	import traceback, sys
	tb = sys.exc_info()[2]
	tbinfo = traceback.format_tb(tb)[0]
	    
	pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
	msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
	
	#arcpy.AddError(pymsg)
	#arcpy.AddError(msgs)
	
	print(pymsg)
	print(msgs)	