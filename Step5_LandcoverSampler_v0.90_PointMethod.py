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
# 1: input Number of directions to sample (WedgeZones)
# 2: input Number of vegetation (transverse) samples in each direction (VegZones)
# 3: Input vegetation height Raster (VegRaster)
# 4: input The distance between transverse samples (TransDistance)
# 5: output polygon file name/path (outpoly_final)
# 6: output table file name/path (outtable_final)

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
	#VegZones = long(arcpy.GetParameterAsText(2))
        #TransDistance = long(arcpy.GetParameterAsText(3))
	#CanopyData = = arcpy.GetParameterAsText(4) # One of these: 1.Codes, 2.Canopy Cover, or 3.LAI
	#LCRaster = arcpy.GetParameterAsText(5) # This is either landcover height or codes
	#CanopyRaster = arcpy.GetParameterAsText(6) # OPTIONAL This is either canopy cover or LAI raster
	#kRaster = arcpy.GetParameterAsText(7) # OPTIONAL The k value raster for LAI
	#EleRaster = arcpy.GetParameterAsText(8)

	#outpoly_final = arcpy.GetParameterAsText(6)
	#outtable_veg_final = arcpy.GetParameterAsText(7) #os.path.dirname(outpoly_final) + "/out_table_veg_final.dbf"
	
	# Start Fill in Data
	inPoint = "D:/Projects/RestorationExample/Shapefiles/V8_Star/McFee_TTools756_post_star.shp"
	NumDirections = 8
	NumZones = 5
	TransDistance = 8
	LCRaster = "D:/Projects/RestorationExample/Raster/LiDAR/veg_ht_int" # This is either landcover height or codes
	CanopyDataType = "CanopyCover"
	CanopyRaster = "" # OPTIONAL This is either canopy cover or a LAI raster
	kRaster = "" # OPTIONAL This is the k value for LAI
	EleRaster = "D:/Projects/RestorationExample/Raster/LiDAR/be" 
	outtable_final = "D:/ArcGis/LIDAR/yachats/output/yachats_wedges_CCC.shp"
	# End Fill in Data

	#Add X and Y fields to inpoints
	#arcpy.AddMessage("Adding X/Y")
	print("Adding X/Y")
	arcpy.AddXY_management(inPoint)

	Angle_Incr = 360.0 / NumDirections / 2

	pointArray = arcpy.Array()
	pntObj = arcpy.Point()
	Length = {}
	origin_x = {}
	origin_y = {}
	COLNAME = []
	VARIABLE = []
	DIRECTION = []
	ZONE = []
	SAMPLE_X = []
	SAMPLE_Y = []
	LDZ = {}
	VALUE = []
	ID = []
	
	i=0
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
	if CanopyDataType == "Canopy Cover":  #Use Canopy Cover methods
		type = ['HEIGHT','CCV','ELE']
		emergentlabel = 'CCV_EMERGENT'			
			      
	if NumDirections == 999:  #999 is a flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
		dir = ['NE','E','SE','S','SW','W','NW']
	else:        
		dir = ['DIR' + str(x) for x in range(1,NumDirections+ 1)]
			
	zone = range(1,int(NumZones)+1)
			
	for row in InRows:
	# Get the raw values from the input points
	# Should put an X/Y field checker here and add/calculate those fields if not present
		Length[i] = row.getValue("LENGTH")
		origin_x[i] = row.getValue("POINT_X")
		origin_y[i] = row.getValue("POINT_Y")
		i= i + 1		
	del InRows

	# for maunal starts
	start = 0 #use the first wedge node FID
	end = i # use the last wedge node FID + 1
	#end = i
	n = end - start
	
	#keeping track of time
	startTime= time.time()
	
	i = 0
	for l in range(start,end):
		for d in range(0,len(dir)):
			Angle = d * Angle_Incr
			for z in range(0,len(zone)):
				for t in range(0,len(type)):
					# Determine Sample location, create point object
					Dis = (z + 0) * TransDistance * units_con
					SAMPLE_X[i] = (Dis * sin(radians(Angle))) + origin_x[i]
					SAMPLE_Y[i] = (Dis * cos(radians(Angle))) + origin_y[i]
					pointList = arcpy.Point(SAMPLE_X[i], SAMPLE_Y[i])					
					
					# Sample the point
					if type[t] == "ELE":
						VALUE[i] = ExtractByPoints(EleRaster, pointList,"INSIDE")					
					if type[t] in ["LC","HEIGHT"]:
						VALUE[i] = ExtractByPoints(LCRaster, pointList,"INSIDE")					
					if type[t] in ["LAI","CCV"]:
						VALUE[i] = ExtractByPoints(CanopyRaster, pointList,"INSIDE")					
					if type[t] = "k":
						VALUE[i] = ExtractByPoints(kRaster, pointList,"INSIDE")					
					# Extract all the other info
					COLNAME[i] = type[t]+'_'+dir[d]+'_'+str(zone[z])
					VARIABLE[i] = type[t]
					DIRECTION[i] = dir[d]
					ZONE[i] = zone[z]
					LDZ[i] = (Length[i] * 10000) + ((dir[i] +1) * 100) + (zone[i] + 1)
					ID[i] = i
					i = i +1
	
	del(l,d,z,t,i)			
	gc.collect()
	

	
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