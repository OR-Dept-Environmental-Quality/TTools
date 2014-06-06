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
	LCRaster = "D:/Projects/RestorationExample/Raster/LiDAR/veg_ht_int" # This would either be Landcover height or codes
	CanopyDataType = "CanopyCover"
	#CanopyRaster = "D:/ArcGis/LIDAR/yachats/VegHt/meters/veght_yach_m_i.img # OPTIONAL This would either be canopy cover or LAI
	#kRaster = # OPTIONAL This would be the k value for LAI
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
	VALUE = []
	
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

	if CanopyDataTyupe == "Codes":        
		type = ['LC']
	if CanopyDataType == "LAI":  #Use LAI methods
		type = ['LC','LAI','k']
		emergentlabel ='LAI_EMERGENT'
	if CanopyDataType == "Canopy Cover":  #Use Canopy Cover methods
		type = ['LC','CCV']
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
	
	for i in range(start,end):
		for d in range(0,len(dir)):
			Angle = d * Angle_Incr
			for z in range(0,len(zone)):
				Dis = (z + 0) * TransDistance * units_con
				pntObj.X = (Dis * sin(radians(Angle))) + origin_x[i]
				pntObj.Y = (Dis * cos(radians(Angle))) + origin_y[i]				
		

		#Create an empty output polygon with the same projection as the input points
		outpoly_temp = os.path.dirname(outpoly_final) + "/out_temp_wedge_%s.shp" % i
		arcpy.CreateFeatureclass_management(os.path.dirname(outpoly_temp),os.path.basename(outpoly_temp), "POINT","","DISABLED","DISABLED",proj)
		arcpy.AddField_management(outpoly_temp, "LWV", "TEXT","","", 30, "", "NULLABLE", "NON_REQUIRED")
		arcpy.AddField_management(outpoly_temp, "LENGTH", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
		arcpy.AddField_management(outpoly_temp, "DIR", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
		arcpy.AddField_management(outpoly_temp, "ZONE", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
				
		# read in iterate through each wedge and vegzone to create a polygon
		OutRows = arcpy.InsertCursor(outpoly_temp)
		for wzone in xrange(WedgeZones):
			AngleStart = ((wzone + 1) * Angle_Incr) - (Angle_Incr / 2)
			Angle = AngleStart
			for vzone in xrange(VegZones):
				BotDis = (vzone + 0) * TransDistance * units_con
				TopDis = (vzone + 1) * TransDistance * units_con
				#iterate through each vertex X/Y to complete the polygon
				if vzone == 0: # First veg zone that is a triangle shape
					for vertex in range(7):
						if vertex in [0]: #bottom Start Angle
							pntObj.X = (BotDis * sin(radians(AngleStart))) + origin_x[i]
							pntObj.Y = (BotDis * cos(radians(AngleStart))) + origin_y[i]
						if vertex in [1]: #top Start Angle
							pntObj.X = (TopDis * sin(radians(AngleStart))) + origin_x[i]
							pntObj.Y = (TopDis * cos(radians(AngleStart))) + origin_y[i]
						if vertex in [2,3,4,5]: #top 2-5
							Angle = Angle + Angle_Incr_sub
							pntObj.X = (TopDis * sin(radians(Angle))) + origin_x[i]
							pntObj.Y = (TopDis * cos(radians(Angle))) + origin_y[i]
						if vertex in [6]: #bottom End Angle
							pntObj.X = (BotDis * cos(radians(Angle))) + origin_x[i]
							pntObj.Y = (BotDis * cos(radians(Angle))) + origin_y[i]
							Angle = Angle - (Angle_Incr_sub * 4)
						polyArray.add(pntObj)
				if vzone != 0:
					for vertex in range(11):
						if vertex in [0]: #bottom Start Angle
							pntObj.X = (BotDis * sin(radians(AngleStart))) + origin_x[i]
							pntObj.Y = (BotDis * cos(radians(AngleStart))) + origin_y[i]
						if vertex in [1]: #top Start Angle
							pntObj.X = (TopDis * sin(radians(AngleStart))) + origin_x[i]
							pntObj.Y = (TopDis * cos(radians(AngleStart))) + origin_y[i]
						if vertex in [2,3,4,5]: #top 2-5
							Angle = Angle + Angle_Incr_sub
							pntObj.X = (TopDis * sin(radians(Angle))) + origin_x[i]
							pntObj.Y = (TopDis * cos(radians(Angle))) + origin_y[i]
						if vertex in [6]: #bottom End Angle
							pntObj.X = (BotDis * sin(radians(Angle))) + origin_x[i]
							pntObj.Y = (BotDis * cos(radians(Angle))) + origin_y[i]
						if vertex in [7,8,9,10]: #bottom 7-10
							Angle =  Angle - Angle_Incr_sub
							pntObj.X = (BotDis * sin(radians(Angle))) + origin_x[i]
							pntObj.Y = (BotDis * cos(radians(Angle))) + origin_y[i]
						polyArray.add(pntObj)
				del vertex
				gc.collect()
				feat = OutRows.newRow()
				feat.LWV = "000000" + str((Length[i] * 10000) + ((wzone +1) * 100) + vzone + 1)
				feat.LENGTH = Length[i]
				feat.WZONE = wzone + 1
				feat.VZONE = vzone + 1
				feat.shape = polyArray
				OutRows.insertRow(feat)
				polyArray.removeAll()
			del vzone
			gc.collect()
		del wzone
		del OutRows
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
