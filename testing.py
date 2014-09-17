# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
import itertools
from arcpy import env
from math import sqrt, pow, ceil
from collections import defaultdict

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

pointfile = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"

env.overwriteOutput = True
def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

try:

    #"""Reads the input point file and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]
    AddFields = ["TEST","TOPO_S","TOPO_E"]
    #AddFields = ["TOPO_W","TOPO_S","TOPO_E"]
    OverwriteData = False
    
    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(pointfile):
	ExistingFields.append(f.name)     
    
    # Check to see if the 1st field exists if yes add it.
    if OverwriteData == False and (AddFields[0] in ExistingFields) == True:
	Incursorfields.append(AddFields[0])
    else:
	OverwriteData = True
    
    # Determine input point spatial units
    proj = arcpy.Describe(pointfile).spatialReference
    
    with arcpy.da.SearchCursor(pointfile,Incursorfields,"",proj) as Inrows:
	if OverwriteData == True:
	    for row in Inrows:
		pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
		pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
		pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
	else:
	    for row in Inrows:
		# Is the data null or zero, if yes grab it.
		if row[5] == None or row[5] == 0:
		    pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
		    pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
		    pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]		    

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