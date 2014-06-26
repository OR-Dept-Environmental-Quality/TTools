#######################################################################################
# TTools
# Step 1: Create Stream Nodes  0.1
# Ryan Michie

# This script will take an input polyline feature and generate evenly spaced points along 
# the line at a user defined spacing distance measured from the downstream endpoint.

# The outputs include a point feature class and matching csv file with 
# a unique "NODE_ID" field,
# a "STREAM_ID" field matching a unique polyline ID field identifed by the user,
# a "LENGTH" value mesured in meters from the downstream end of the stream, a Stream KM value, 
# and any other attribute field in the input polyline feature.

# Future Updates
# Topo is not included in the csv output. Might have to read in the topo csv output and concatenate it to landcover output.
# Add the sampled raster values back into the input point feature class like previous versions of TTools.

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# parameter values
# 0: input stream centerline polyline (inLine)
# 1: spacing between nodes (NodeDistance)
# 2: output point file name/path (outpoint_final)
# 3: output csv file name/path (outcsv_final)

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from arcpy import env
from math import radians, sin, cos
from collections import defaultdict
from operator import itemgetter
import csv


# Check out the ArcGIS Spatial Analyst extension license
#arcpy.CheckOutExtension("Spatial")
#from arcpy.sa import *
from arcpy.management import *

env.overwriteOutput = True

#enable garbage collection
gc.enable()

try:

	#inLine = arcpy.GetParameterAsText(0)
	#NodeDistance = arcpy.GetParameterAsText(1)
	#outpoint_final = arcpy.GetParameterAsText(2)
	#outcsv_final = arcpy.GetParameterAsText(3)
	
	# Start Fill in Data
	inLine = r"D:\Projects\RestorationExample\Shapefiles\McFee_Streamline_Dissolve.shp"
	NodeDistance = 50
	outpoint_final = r"D:\Projects\RestorationExample\out_samplepoint.shp"
	outcsv_final = r"D:\Projects\RestorationExample\out_sampledata.csv"
	# End Fill in Data
	
	pts = []
	with arcpy.da.SearchCursor(inLine,"SHAPE@") as rows:
		for row in rows:
			theLength_m = row[0].length
			numNodes = int(theLength_m / NodeDistance)
			nodes = range(0,numNodes+1)
			position_m = [node * NodeDistance for node in nodes]
			for position in position_m:
				pts.append(row[0].positionAlongLine(position))
	#del pts	

	#Add X and Y fields to inpoints
	#arcpy.AddMessage("Adding X/Y")
	print("Adding X/Y")
	arcpy.AddXY_management(inPoint)
	
	length = []
	origin_x = []
	origin_y = []	

	InRows = arcpy.SearchCursor(inPoint,"","","LENGTH; POINT_X; POINT_Y","")	
	
	# Determine input spatial units and set conversion factor to get from meters to the input spatial units
	proj = arcpy.Describe(inPoint).SpatialReference
	unitCode = arcpy.Describe(inPoint).SpatialReference.linearUnitCode
	if unitCode == 9001: #International meter
		units_con = 1 
	if unitCode == 9002: #International foot
		units_con = 3.280839895013123359580052493
	if unitCode == 9003: #US Survey foot
		units_con = 3937/1200
	if unitCode == 9005: #Clarke's foot
		units_con = 3.280869330266635653352551371	
	
	if CanopyDataType == "Codes":        
		type = ['LC','ELE']
	if CanopyDataType == "LAI":  #Use LAI methods
		type = ['HEIGHT','LAI','k','ELE']
		emergentlabel ='LAI_EMERGENT'
	if CanopyDataType == "CanopyCover":  #Use Canopy Cover methods
		type = ['HEIGHT','CCV','ELE']
		emergentlabel = 'CCV_EMERGENT'			
	
	# Set the converstion factor to make sure the iput elevation z units are in meters
	if EleUnits == "Meters": #International meter
		ele_con = 1 
	if EleUnits == "Feet": #International foot
		ele_con = 3.280839895013123359580052493
	if EleUnits == "Other": #Some other units
		sys.exit("Please modify your raster elevation units so they are either in meters or feet")	
	
			      
	if NumDirections == 999:  #999 is a flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
		azimuths = [45,90,135,180,225,270,315]
	else:        
		azimuths = [x * 360.0 / NumDirections for x in range(1,NumDirections+ 1)]
	
	if StreamSample == "TRUE":
		zone = range(0,int(NumZones))
	else:
		zone = range(1,int(NumZones+1))
		
	for row in InRows:
	# Get the raw values from the input points
	# Should put an X/Y field checker here and add/calculate those fields if not present
		length.append(row.getValue("LENGTH"))
		origin_x.append(row.getValue("POINT_X"))
		origin_y.append(row.getValue("POINT_Y"))
		
	del(InRows,row)	
	#Calculate the unique number of dictionary values
	#N = len(length) * len(azimuths) * len(zone) * len(type)
	
	# build the NODES nested dictionary and save a generic i value
	def tree(): return defaultdict(tree)
	
	i=0
	NODES = tree()	
	type2 = type + ['SAMPLE_X','SAMPLE_Y']
	for l in length:
		for a in azimuths:
			for z in zone:
				for t in type2:
					NODES[l][a][z][t] = i
					i = i +1


	del(i,x,l,t,a,z)	
	
	#keeping track of time
	startTime= time.time()
	
	#arcpy.AddMessage("Extracting values")
	print("Extracting raster values")	

	for l in range(0,len(length)):
		print("Processing Node %s of %s" % (l+1, len(length)))
		#arcpy.AddMessage("Process Node %s of %s" % (l+1, len(length)))	
		for a in azimuths:
			for z in zone:
				for t in type2:		
								
					# Calculate the x and y sample location
					_X_ = (z * TransDistance * units_con * sin(radians(a))) + origin_x[l]
					_Y_ = (z * TransDistance * units_con * cos(radians(a))) + origin_y[l]
					xypoint = str(_X_) + " " + str(_Y_) # xy string requirement for arcpy.GetCellValue
					
					# Sample the point value from the appropriate raster and add to NODES dictionary
					if t == "SAMPLE_X":
						NODES[length[l]][a][z]['SAMPLE_X'] = _X_
					if t == "SAMPLE_Y":
						NODES[length[l]][a][z]['SAMPLE_Y'] = _Y_
					if t == "ELE":
						thevalue = arcpy.GetCellValue_management(EleRaster, xypoint)
						NODES[length[l]][a][z][t] = float(thevalue.getOutput(0)) * ele_con
					if t in ["LC","HEIGHT"]:
						thevalue = arcpy.GetCellValue_management(LCRaster, xypoint)
						NODES[length[l]][a][z][t] = float(thevalue.getOutput(0))
					if t in ["LAI","CCV"]:
						thevalue = arcpy.GetCellValue_management(CanopyRaster, xypoint)
						NODES[length[l]][a][z][t] = float(thevalue.getOutput(0))
					if t == "k":
						thevalue =arcpy.GetCellValue_management(kRaster, xypoint)
						NODES[length[l]][a][z][t] = float(thevalue.getOutput(0))
			
	del(l,a,z,t,thevalue,_X_,_Y_)			
	gc.collect()
		
	
	####################################################################################################### 
	# build the output= point feature class using the data from the NODES dictionary
	#arcpy.AddMessage("Exporting Data")
	print("Exporting Data")	
	
	#Create an empty output with the same projection as the input points
	cursorfields = ["LENGTH","AZIMUTH","ZONE"] + type2 + ["SHAPE@X","SHAPE@Y"]
	arcpy.CreateFeatureclass_management(os.path.dirname(outpoint_final),os.path.basename(outpoint_final), "POINT","","DISABLED","DISABLED",proj)
	
	# Add attribute fields
	arcpy.AddField_management(outpoint_final, "LENGTH", "DOUBLE", 16, 3, "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "AZIMUTH", "DOUBLE", 16, 2, "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "ZONE", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")	
	for t in type2:
		arcpy.AddField_management(outpoint_final, t, "DOUBLE", 16, "", "", "", "NULLABLE", "NON_REQUIRED")	
	
	# put the dictionary in list form needed to build the output point feature class
	type3 = type2 + ["SAMPLE_X","SAMPLE_Y"]
	laz = [[l,a,z] for l in length for a in azimuths for z in zone]
	NODES_shp = [laz[row] + [NODES[laz[row][0]][laz[row][1]][laz[row][2]][t] for t in type3] for row in range(0,len(laz))]	
	
	cursor = arcpy.da.InsertCursor(outpoint_final, cursorfields)                  
	
	for row in NODES_shp:
		cursor.insertRow(row)
	del(cursor,l,a,z,laz,row)
	
	####################################################################################################### 
	# Output the csv file
	
	# Get the stream km dictionary keys and sort them
	NODE_keys = NODES.keys()
	NODE_keys.sort()
	
	NODES_csv = [[NODES[l][a][z][t] for t in type for a in azimuths for z in zone] for l in length]
	
	# add in the stream km at the beginning of the list
	for l in range(0,len(NODE_keys)):
		NODES_csv[l].insert(0,NODE_keys[l])	
	
	# Add the header row
	LC_Header = ["Stream_KM"]
	
	for t in range(0,len(type)):
		for a in range(0,len(azimuths)):
			for z in range(0,len(zone)):
				LC_Header.append(type[t]+'_A'+str(a+1)+'_Z'+str(zone[z]))		
	
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