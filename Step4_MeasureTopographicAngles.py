#######################################################################################
# TTools
# Step 4: Measure Topographic Angles - v 0.1
# Ryan Michie

# This script will take an input point feature (from Step 1) and determine the maximum
# topographic elevation and the the slope angle from the point in different directions.

# The outputs include the slope angle saved into the point feature class and csv file
# created in step 1, a new point feature class for each x/y location 
# of the maximum elevation, and a csv file formatted for input into heat source 9.

# Future Updates
#

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# parameter values
# 0: input TTools point file (inPoint)
# 1: input the directions to sample (Directions) 1. [W,S,E], 2. [NE,E,SE,S,SW,W,NW]
# 2: input the maximum km distance to search (MaxSearchDistance_km)
# 3: input elevation raster (EleRaster)
# 4: input elvation raster units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 5: output sample point file name/path (outpoint_final)
# 6: output csv file to write to name/path outcsv_final) - OPTIONAL



# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from arcpy import env
from math import radians, sin, cos, atan
from collections import defaultdict
from operator import itemgetter
import csv


# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
#from arcpy.sa import *
#from arcpy.management import *

env.overwriteOutput = True

#enable garbage collection
gc.enable()

def GetLinearUnitConversion(inFeature):
  """Get the conversion factor to get from meters to the input spatial units"""
  unitCode = arcpy.Describe(inFeature).SpatialReference.linearUnitCode
  if unitCode == 9001: #International meter
	  units_con = 1 
  if unitCode == 9002: #International foot
	  units_con = 3.280839895013123359580052493
  if unitCode == 9003: #US Survey foot
	  units_con = 3937/1200
  if unitCode == 9005: #Clarke's foot
	  units_con = 3.280869330266635653352551371 
  if unitCode not in [9001,9002,9003,9005]:
          system.exit("Unreconized spatial reference. Use projection with units of feet or meters.") 
  return units_con

try:

	#inPoint = arcpy.GetParameterAsText(0)
	#Directions = long(arcpy.GetParameterAsText(1))
	#MaxSearchDistance_km = long(arcpy.GetParameterAsText(2))
	#EleRaster = arcpy.GetParameterAsText(3)
	#EleUnits = arcpy.GetParameterAsText(4)
	#outpoint_final = arcpy.GetParameterAsText(11)
	#outcsv_final = arcpy.GetParameterAsText(12)
	
	# Start Fill in Data
	inPoint = r"D:\Projects\RestorationExample\out_nodes.shp"
	Directions = 1
	MaxSearchDistance_km = 2
	EleRaster = r"D:\Projects\RestorationExample\Raster\LiDAR\be\be_lidar"
	EleUnits = "Feet"
	outpoint_final = r"D:\Projects\RestorationExample\topo_samplepoint.shp"
	outcsv_final = r"D:\Projects\RestorationExample\out_nodes.csv"
	# End Fill in Data

	#Add X and Y fields to inpoints
	#arcpy.AddMessage("Adding X/Y")
	print("Adding X/Y")
	arcpy.AddXY_management(inPoint)

	# Determine input point spatial units and set conversion factor to get from the raster units to the input spatial units
	proj = arcpy.Describe(inPoint).SpatialReference	
	rasterproj = arcpy.Describe(EleRaster).SpatialReference
	mto_node_con = GetLinearUnitConversion(inPoint)
	mto_Ele_con = GetLinearUnitConversion(EleRaster)
	units_con=  mto_node_con / mto_Ele_con
	MaxSearchDistance = MaxSearchDistance_km * mto_node_con * 1000
	
	if EleUnits = "Feet":
	       eleZ_to_m = 3.280839895013123359580052493
	else:
	       eleZ_to_m = 1
	
	# Get the elevation raster cell size
	CellSizeResult = arcpy.GetRasterProperties_management(EleRaster, "CELLSIZEX")
	CellSize = float(CellSizeResult.getOutput(0))	      
	if Directions == 2: # All directions
		azimuths = [45,90,135,180,225,270,315]
	else:        
		azimuths = [270,180,90]
	
	
	#keeping track of time
	startTime= time.time()	
		
	Incursorfields = ["SHAPE@XY","NODE_ID"]
	NODES = []
	
	with arcpy.da.SearchCursor(inPoint,Incursorfields) as Inrows:
		#arcpy.SetProgressor("step", "Creating Nodes", 0, Inrows, 1)
		print("Sampling Elevations")
		for row in Inrows:
			# Get the raw values from the input points
			NODES.append(row)
		del(Inrows,row)
	
	#arcpy.AddMessage("Extracting values")
	print("Extracting raster values")	

	n = 0
	for node in NODES:
		print("Processing Node %s of %s" % (n+1, len(NODES)))
		#arcpy.AddMessage("Process Node %s of %s" % (n+1, len(NODES)))
		node_x = node[0][0]
		node_y = node[0][1]
		nodexy = str(node_x) + " " + str(node_y) # xy string requirement for arcpy.GetCellValue	
		thevalue = arcpy.GetCellValue_management(EleRaster, nodexy)
		nodeZ = float(thevalue.getOutput(0)) * ele_con		
		SearchDistance = 0
		MaxShadeAngle = 0
		for a in azimuths:
			while not SearchDistance > MaxSearchDistance:
			  # This is the Skippy algorithm from Greg Pelletier
			  if i <= 10:
				  SearchDistance = SearchDistance + (CellSize * units_con)
			  if i <= 20:
				  SearchDistance = SearchDistance + (CellSize * units_con * 3)
			  if i <= 40:
				  SearchDistance = SearchDistance + (CellSize * units_con * 6)
			  if i <= 50:
				  SearchDistance = SearchDistance + (CellSize * units_con * 12)
			  if i <= 60:
				  SearchDistance = SearchDistance + (CellSize * units_con * 25)
			  else:
				  SearchDistance = SearchDistance + (CellSize * units_con * 50)
			  i = i + 1			  
			  
			  # Calculate the x and y sample location.
			  sample_x = (SearchDistance * sin(radians(a))) + node_x
			  sample_y = (SearchDistance * cos(radians(a))) + node_y
			  samplexy = str(sample_x) + " " + str(sample_y) # xy string requirement for arcpy.GetCellValue
					
			  # Sample the elevation value from the elevation raster
			  thevalue = arcpy.GetCellValue_management(EleRaster, samplexy)
			  sampleZ= float(thevalue.getOutput(0)) * ele_con
			    
			  #calculate the topographic shade angle (in degrees)
			  ShadeAngle = atn((sampleZ - nodeZ) / SearchDistance)
			  if ShadeAngle > MaxShadeAngle:
			        MaxShadeAngle = ShadeAngle
				MaxShadeAngle_X = sample_x
				MaxShadeAngle_Y = sample_y
			
			Node
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