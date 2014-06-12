#######################################################################################
# TTools
# Step 5: Sample Landcover - Point Method v 0.93
#
# This script will take a point input (from Step 1) and sample a landcover raster
# in a user specificed number of cardianal directions with point samples spaced at a user defined distance
# moving away from the stream.

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.0 SP3 or higher to run.
#
# Ryan Michie
#######################################################################################

# parameter values
# 0: input TTools point file (inPoint)
# 1: input number of directions to sample (NumDirections)
# 2: input number of vegetation (transverse) samples in each direction (NumZones)
# 3: input The distance between transverse samples (TransDistance)
# 4: input canopy data type. 1."Codes", 2."CanopyCover", or 3."LAI" (CanopyData)
# 5: input landcover code or height raster (LCRaster)
# 6: input (optional) canopy cover or LAI raster (CanopyRaster)
# 7: input (optional) k coeffcient raster (kRaster)
# 8: input elevation raster (EleRaster)
# 9: output sample point file name/path (outpoint_final)
# 10: output csv file name/path (outcsv_final)

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
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
from arcpy.management import *

env.overwriteOutput = True

#enable garbage collection
gc.enable()


try:

	#inPoint = arcpy.GetParameterAsText(0)
	#NumDirections = long(arcpy.GetParameterAsText(1))
	#NumZones = long(arcpy.GetParameterAsText(2))
        #TransDistance = long(arcpy.GetParameterAsText(3))
	#CanopyData = = arcpy.GetParameterAsText(4) # One of these: 1."Codes", 2."CanopyCover", or 3."LAI"
	#LCRaster = arcpy.GetParameterAsText(5) # This is either landcover height or codes
	#CanopyRaster = arcpy.GetParameterAsText(6) # OPTIONAL This is either canopy cover or LAI raster
	#kRaster = arcpy.GetParameterAsText(7) # OPTIONAL The k value raster for LAI
	#EleRaster = arcpy.GetParameterAsText(8)
	#outpoint_final = arcpy.GetParameterAsText(9)
	#outcsv_final = arcpy.GetParameterAsText(10)
	
	# Start Fill in Data
	inPoint = "D:/Projects/RestorationExample/Shapefiles/V8_Star/McFee_TTools756_post_star_HARN.shp"
	NumDirections = 8
	NumZones = 6 # includes stream sample
	TransDistance = 8
	LCRaster = "D:/Projects/RestorationExample/Raster/LiDAR/veg_ht_int/veght_lidar" # This is either landcover height or codes
	CanopyDataType = "Codes"
	CanopyRaster = "" # OPTIONAL This is either canopy cover or a LAI raster
	kRaster = "" # OPTIONAL This is the k value for LAI
	EleRaster = "D:/Projects/RestorationExample/Raster/LiDAR/be/be_lidar" 
	outpoint_final = "D:/Projects/RestorationExample/out_samplepoint.shp"
	outcsv_final = "D:/Projects/RestorationExample/out_sampledata.csv"
	# End Fill in Data

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
			      
	if NumDirections == 999:  #999 is a flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
		azimuths = [45,90,135,180,225,270,315]
	else:        
		azimuths = [x * 360.0 / NumDirections for x in range(1,NumDirections+ 1)]
			
	zone = range(0,int(NumZones))	
	
	i=0		
	for row in InRows:
	# Get the raw values from the input points
	# Should put an X/Y field checker here and add/calculate those fields if not present
		length.append(row.getValue("LENGTH"))
		origin_x.append(row.getValue("POINT_X"))
		origin_y.append(row.getValue("POINT_Y"))
		i= i + 1		
	del(InRows,row)	
	#Calculate the unique number of dictionary values
	N = len(length) * len(azimuths) * len(zone) * len(type)
	
	# Build the Dictionarys
	dictkeys = ["LENGTH","SAMPLE_X","SAMPLE_Y","VARIABLE","AZIMUTH","ZONE","VALUE"]
	DATA = defaultdict(list)
	for d in dictkeys:
		for i in range(0,N):
			DATA[d].append(i)	
	
	def tree(): return defaultdict(tree)
	
	i=0
	NODES = tree()
	#for l in length:
		#for t in type:
			#for a in azimuths:
				#for z in zone:
					#NODES[l][t][a][z] = i
					#i = i +1

	# New test	
	type2 = type + ['SAMPLE_X','SAMPLE_Y']
	for l in length:
		for a in azimuths:
			for z in zone:
				for t in type2:
					NODES[l][a][z][t] = i
					i = i +1


	del(i,x,d,l,t,a,z)	
	
	#keeping track of time
	startTime= time.time()
	end = len(length)
	
	#arcpy.AddMessage("Extracting values")
	print("Extracting raster values")	
	i = 0
	for l in range(0,len(length)):
		print("Processing Node %s of %s" % (l+1, end))
		#arcpy.AddMessage("Process Node %s of %s" % (l+1, len(length)))	
		for a in range(0,len(azimuths)):
			for z in range(0,len(zone)):
				for t in range(0,len(type)):		
				
					DATA['LENGTH'][i] = length[l]					
					# Determine Sample location,
					_X_ = (zone[z] * TransDistance * units_con * sin(radians(azimuths[a]))) + origin_x[l]
					_Y_ = (zone[z] * TransDistance * units_con * cos(radians(azimuths[a]))) + origin_y[l]
					
					DATA['SAMPLE_X'][i] = _X_
					DATA['SAMPLE_Y'][i] = _Y_
					NODES[length[l]][azimuths[a]][z]['SAMPLE_X'] = _X_
					NODES[length[l]][azimuths[a]][z]['SAMPLE_Y'] = _Y_
					xypoint = str(_X_) + " " + str(_Y_) # string requiremetns for GetCellValue
					
					# Sample the point value from the raster
					if type[t] == "ELE":
						thevalue = arcpy.GetCellValue_management(EleRaster, xypoint)
					if type[t] in ["LC","HEIGHT"]:
						thevalue = arcpy.GetCellValue_management(LCRaster, xypoint)
					if type[t] in ["LAI","CCV"]:
						thevalue = arcpy.GetCellValue_management(CanopyRaster, xypoint)
					if type[t] == "k":
						thevalue =arcpy.GetCellValue_management(kRaster, xypoint)[_dict_]
					DATA['VALUE'][i] = float(thevalue.getOutput(0))
					#NODES[length[l]][type[t]][azimuths[a]][z] = float(thevalue.getOutput(0))
					NODES[length[l]][azimuths[a]][z][type[t]] = float(thevalue.getOutput(0))
					# other data
					DATA['VARIABLE'][i] = type[t]
					DATA['AZIMUTH'][i] = azimuths[a]
					DATA['ZONE'][i] = zone[z]
	
					i = i +1
			
	del(l,a,z,t,i,thevalue)			
	gc.collect()
		
	
	####################################################################################################### 
	# build the point feature class uisng the data from the dictionary
	#arcpy.AddMessage("Exporting Data")
	print("Exporting Data")	
	
	#Create an empty output with the same projection as the input points
	cursorfields = ["LENGTH","SAMPLE_X","SAMPLE_Y","VARIABLE","AZIMUTH","ZONE","VALUE","SHAPE@X","SHAPE@Y"]
	shapekeys = ["LENGTH","SAMPLE_X","SAMPLE_Y","VARIABLE","AZIMUTH","ZONE","VALUE","SAMPLE_X","SAMPLE_Y"]

	arcpy.CreateFeatureclass_management(os.path.dirname(outpoint_final),os.path.basename(outpoint_final), "POINT","","DISABLED","DISABLED",proj)
	
	#Add Fields
	arcpy.AddField_management(outpoint_final, "LENGTH", "DOUBLE", 16, 3, "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "SAMPLE_X", "DOUBLE", 16, "", "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "SAMPLE_Y", "DOUBLE", 16, "", "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "VARIABLE", "TEXT", "", "", "30", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "AZIMUTH", "DOUBLE", 16, 2, "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "ZONE", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoint_final, "VALUE", "DOUBLE", 16, 3, "", "", "NULLABLE", "NON_REQUIRED")	
	
	# put the dictionary in a list form needed to build the output point feature class
	data = [[DATA[k][row] for k in shapekeys] for row in range(0,N)]
	
	cursor = arcpy.da.InsertCursor(outpoint_final, cursorfields)                  
		
	for row in data:
		cursor.insertRow(row)
	del cursor
	
	####################################################################################################### 
	# output csv and point files
	
	# Get the stream km dictionary keys and sort them
	NODE_keys = NODES.keys()
	NODE_keys.sort()
	
	#NODES_csv = [[NODES[l][t][a][z] for t in type for a in azimuths for z in zone] for l in length]
	# new test
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
	
	
	#csvout = sorted(data_csv, key=itemgetter(0,1,2,3))	
	####################################################################################################### 
	endTime = time.time()
	elapsedmin= (endTime - startTime) / 60	
	# output csv and point files
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