#######################################################################################
# TTools
# Step 5: Sample Landcover - Point Method v 0.90
#
# This script will take a point input (from Step 1) and sample a landcover raster
# in a user specificed number of cardianal directions with point samples spaced at a user defined distance
# moving away from the stream.

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.0 SP3 or higher to run.
# The script will run with ArcGIS 10.0 SP1 or SP 2 but crashes after about 15-30 loops due to an ArcGIS bug
# See http://support.esri.com/en/bugs/nimbus/TklNMDYzODE0
#
# Ryan Michie
#######################################################################################

# parameter values
# 0: input TTools point file (inPoint)
# 1: input number of directions to sample (NumDirections)
# 2: input number of vegetation (transverse) samples in each direction (NumZones)
# 3: input The distance between transverse samples (TransDistance)
# 4: input canopy data type. 1.Codes, 2.Canopy Cover, or 3.LAI (CanopyData)
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


# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

env.overwriteOutput = True

#enable garbage collection
gc.enable()

try:

	#inPoint = arcpy.GetParameterAsText(0)
	#NumDirections = long(arcpy.GetParameterAsText(1))
	#NumZones = long(arcpy.GetParameterAsText(2))
        #TransDistance = long(arcpy.GetParameterAsText(3))
	#CanopyData = = arcpy.GetParameterAsText(4) # One of these: 1.Codes, 2.Canopy Cover, or 3.LAI
	#LCRaster = arcpy.GetParameterAsText(5) # This is either landcover height or codes
	#CanopyRaster = arcpy.GetParameterAsText(6) # OPTIONAL This is either canopy cover or LAI raster
	#kRaster = arcpy.GetParameterAsText(7) # OPTIONAL The k value raster for LAI
	#EleRaster = arcpy.GetParameterAsText(8)
	#outpoint_final = arcpy.GetParameterAsText(9)
	#outcsv_final = arcpy.GetParameterAsText(10)
	
	# Start Fill in Data
	inPoint = "D:/Projects/RestorationExample/Shapefiles/V8_Star/McFee_TTools756_post_star_HARN.shp"
	NumDirections = 8
	NumZones = 5
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
	
	length = {}
	origin_x = {}
	origin_y = {}	

	InRows = arcpy.SearchCursor(inPoint,"","","LENGTH; POINT_X; POINT_Y","")
		
	i=0		
	for row in InRows:
	# Get the raw values from the input points
	# Should put an X/Y field checker here and add/calculate those fields if not present
		length[i] = [row.getValue("LENGTH")]
		origin_x[i] = row.getValue("POINT_X")
		origin_y[i] = row.getValue("POINT_Y")
		i= i + 1		
	del InRows
	
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
		dir = ['NE','E','SE','S','SW','W','NW']
		dirnum = [45,90,135,180,225,270,315]
	else:        
		dir = ['DIR' + str(x) for x in range(1,NumDirections+ 1)]
		dirnum = [x for x in range(1,NumDirections+ 1)]
			
	zone = range(1,int(NumZones)+1)	
	
	#Calculate the unique number of dictionary values
	en = i * (len(dirnum)+1) * (len(zone)+1) * (len(type)+1)
	
	# TODO need to fix dictionary so it is easier to extract values based on i for each key.
	# TODO needs to be is initilized with number of unique values (en) .
	dictkeys = ["LENGTH","VARIABLE","DIRECTION","ZONE","SAMPLE_X","SAMPLE_Y","VALUE"]
	DATA = [dict.fromkeys(dictkeys)]
	
	Angle_Incr = 360.0 / NumDirections	
	
	# for maunal starts
	start = 0 #use the first wedge node FID
	end = i # use the last wedge node FID + 1
	#end = i
	n = end - start
	
	#keeping track of time
	startTime= time.time()
	
	i = 0
	for l in range(start,end):
		for d in range(0,len(dirnum)):
			if NumDirections == 999:
				Angle = dirnum[d]
			else:
				Angle = dirnum[d] * Angle_Incr
			 
			for z in range(0,len(zone)):
				for t in range(0,len(type)):
					# Determine Sample location, create point object
					
					SAMPLE_X[i] = (zone[z] * TransDistance * units_con * sin(radians(Angle))) + origin_x[l]
					SAMPLE_Y[i] = (zone[z] * TransDistance * units_con * cos(radians(Angle))) + origin_y[l]
					xypoint = str(SAMPLE_X[i]) + " " + str(SAMPLE_Y[i]) # TODO i don't think this works. check weird arc requiremetns for GetCellValue
					
					#pointList = [arcpy.Point(SAMPLE_X[i], SAMPLE_Y[i])]					
					
					# Sample the point
					if type[t] == "ELE":
						VALUE[i] = arcpy.GetCellValue_management(EleRaster, xypoint)					
					if type[t] in ["LC","HEIGHT"]:
						VALUE[i] = arcpy.GetCellValue_management(LCRaster, xypoint)					
					if type[t] in ["LAI","CCV"]:
						VALUE[i] = arcpy.GetCellValue_management(CanopyRaster, xypoint)					
					if type[t] == "k":
						VALUE[i] =arcpy.GetCellValue_management(kRaster, xypoint)					
					# Extract all the other info
					COLNAME[i] = type[t]+'_'+dir[d]+'_'+str(zone[z])
					VARIABLE[i] = type[t]
					DIRECTION[i] = dirnum[d]
					ZONE[i] = zone[z]
					LDZ[i] = (Length[l] * 10000) + ((dirnum[d] +1) * 100) + (zone[z] + 1)
					ID[i] = i
					i = i +1
		
	del(l,d,z,t,i)			
	gc.collect()
	
	# build the point file for each sample using sample x/y
	# rehape dictionary into wide format for heatsource import
	# output csv and point files.

	
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