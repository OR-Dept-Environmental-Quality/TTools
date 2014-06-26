#######################################################################################
# TTools
# Step 1: Create Stream Nodes  version 0.5
# Ryan Michie

# This script will take an input polyline feature with unique stream IDs and generate 
# evenly spaced points along each unique stream ID line at a user defined spacing 
# measured from the downstream endpoint.

# The outputs include a point feature class and matching csv file with 
# a unique "NODE_ID" field,
# a "STREAM_ID" field matching a unique polyline ID field identifed by the user,
# a "LENGTH" field with the distance in meters from the downstream end of the stream, and
# "STREAM_KM", X, and Y.

# Future Updates
# make many of these tasks into methods.

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# parameter values
# 0: input stream centerline polyline (inLine)
# 1: unique StreamID field (StreamIDfield)
# 2: spacing between nodes (NodeDistance)
# 3: output point file name/path (outpoint_final)
# 4: output csv file name/path (outcsv_final)

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from arcpy import env
from math import radians, sin, cos
from collections import defaultdict
from operator import itemgetter
import csv

env.overwriteOutput = True

#enable garbage collection
gc.enable()

try:

	#inLine = arcpy.GetParameterAsText(0)
	#StreamIDfield = arcpy.GetParameterAsText(1)
	#NodeDistance = arcpy.GetParameterAsText(2)
	#outpoint_final = arcpy.GetParameterAsText(3)
	#outcsv_final = arcpy.GetParameterAsText(4)
	
	# Start Fill in Data
	inLine = r"D:\Projects\RestorationExample\Shapefiles\McFee_Streamline_Dissolve.shp"
	SIDname = "GNIS_ID"
	NodeDistance = 50
	outpoint_final = r"D:\Projects\RestorationExample\out_nodes.shp"
	outcsv_final = r"D:\Projects\RestorationExample\out_nodes.csv"
	# End Fill in Data
	Incursorfields = ["SHAPE@",SIDname]
	NODES_shp = []
	
	#keeping track of time
	startTime= time.time()	
	
	# Determine input spatial units and the Stream ID field properties
	proj = arcpy.Describe(inLine).SpatialReference
	SIDtype = arcpy.ListFields(inLine,SIDname)[0].type
	SIDprecision = arcpy.ListFields(inLine,SIDname)[0].precision
	SIDscale = arcpy.ListFields(inLine,SIDname)[0].scale
	SIDlength = arcpy.ListFields(inLine,SIDname)[0].length
	
	NID = 0
	with arcpy.da.SearchCursor(inLine,Incursorfields) as Inrows:
		#arcpy.SetProgressor("step", "Creating Nodes", 0, Inrows, 1)
		print("Creating Nodes")
		for row in Inrows:
			LineLength = row[0].getLength("PRESERVE_SHAPE")
			numNodes = int(LineLength / NodeDistance)
			nodes = range(0,numNodes+1)
			positions = [node * NodeDistance / LineLength for node in nodes] # list of Lengths in meters
			for position in positions:
				XY = row[0].positionAlongLine(position,True).centroid
				# list of "NODE_ID",SIDname,"LENGTH","STREAM_KM","POINT_X","POINT_Y","SHAPE@X","SHAPE@Y"
				NODES_shp.append((NID, row[1], int(position * LineLength), float(position * LineLength /1000), XY.X, XY.Y, XY.X, XY.Y ))
				NID = NID + 1
		#arcpy.SetProgressorPosition()
		print("Adding X/Y")
	
	del(Inrows, row, positions, position, NID, XY, LineLength, numNodes, nodes)		
	gc.collect()
		
	####################################################################################################### 
	# build the output= point feature class using the data from the NODES_shp
	
	#arcpy.AddMessage("Exporting Data")
	print("Exporting Data")	
	
	#Create an empty output with the same projection as the input polyline
	cursorfields = ["NODE_ID",SIDname,"LENGTH","STREAM_KM","POINT_X","POINT_Y","SHAPE@X","SHAPE@Y"]
	arcpy.CreateFeatureclass_management(os.path.dirname(outpoint_final),os.path.basename(outpoint_final), "POINT","","DISABLED","DISABLED",proj)
	
	# Add attribute fields
	arcpy.AddField_management(outpoint_final, "NODE_ID", "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, SIDname, SIDtype, SIDprecision, SIDscale, SIDlength, "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "LENGTH", "LONG", "","", "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "STREAM_KM", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "POINT_X", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "POINT_Y", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	
	cursor = arcpy.da.InsertCursor(outpoint_final, cursorfields)                  
	
	for row in NODES_shp:
		cursor.insertRow(row)
	del(cursor,row)
	
	####################################################################################################### 
	# Output the csv file
	
	# Remove the last pair of X and Y
	NODES_csv = [NODES_shp[row][0:6] for row in range(0,len(NODES_shp))]	
	
	# Add the header row
	LC_Header = ["NODE_ID",SIDname,"LENGTH","Stream_KM","POINT_X","POINT_Y"]
	NODES_csv.insert(0,LC_Header)
	
	# Export to csv
	with open(outcsv_final, "wb") as f:
		writer = csv.writer(f)
		writer.writerows(NODES_csv)	

	####################################################################################################### 
	
	endTime = time.time()
	elapsedmin= (endTime - startTime) / 60	
	print("Process Complete in %s minutes" % (elapsedmin))
	#arcpy.AddMessage("Process Complete at %s, %s minutes" % (endTime, elapsedmin))

	
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