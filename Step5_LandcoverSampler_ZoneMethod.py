#######################################################################################
# TTools
# Step 5: Sample Landcover - Zone Method v 0.96
#
# This script will take a point input (from Step 1) and generate a series of pie shapped
# polygon zones extending outward from the point. The number of polygon zones and the 
# distance the polygones extend from the point is user defined.
# The script will then sample the landcover raster grid and run spatial statistics on each zone.
# The output is the zone polygons and a table file of the zonal summary statistics.
#
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
	#WedgeZones = long(arcpy.GetParameterAsText(1))
	#VegZones = long(arcpy.GetParameterAsText(2))
        #TransDistance = long(arcpy.GetParameterAsText(3))
	#VegRaster = arcpy.GetParameterAsText(4)
	#EleRaster = arcpy.GetParameterAsText(5)

	#outpoly_final = arcpy.GetParameterAsText(6)
	#outtable_veg_final = arcpy.GetParameterAsText(7) #os.path.dirname(outpoly_final) + "/out_table_veg_final.dbf"
        #outtable_ele_final = arcpy.GetParametersAsText(8) #os.path.dirname(outpoly_final) + "/out_table_ele_final.dbf"
	
	# Start Fill in Data
	inPoint = "D:/ArcGis/LIDAR/yachats/Yachats_TTools_756_HARN.shp"
	WedgeZones = 8
	VegZones = 5
	TransDistance = 8
	VegRaster = "D:/ArcGis/LIDAR/yachats/VegHt/meters/veght_yach_m_i.img"
	EleRaster = "D:/ArcGis/LIDAR/yachats/BareEarth/meters/be_yach_m"
	outpoly_final = "D:/ArcGis/LIDAR/yachats/output/yachats_wedges_CCC.shp"
	outtable_veg_final = "D:/ArcGis/LIDAR/yachats/output/yachats_outtable_CCC_veght_m.dbf"
	outtable_ele_final = "D:/ArcGis/LIDAR/yachats/output/yachats_outtable_CCC_be_m.dbf"
	# End Fill in Data

	#Add X and Y fields to inpoints
	#arcpy.AddMessage("Adding X/Y")
	print("Adding X/Y")
	arcpy.AddXY_management(inPoint)

	Angle_Incr = 360.0 / WedgeZones
	Angle_Incr_sub = Angle_Incr / 4.0

	polyArray = arcpy.Array()
	pntObj = arcpy.Point()
	length = []
	origin_x = []
	origin_y = []
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

	#Create an empty output with the same projection as the input points
	#arcpy.AddMessage("Creating empty output")
	print("Creating empty output")
	arcpy.CreateFeatureclass_management(os.path.dirname(outpoly_final),os.path.basename(outpoly_final), "POLYGON","","DISABLED","DISABLED",proj)
	arcpy.AddField_management(outpoly_final, "LWV", "TEXT","","", 30, "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoly_final, "LENGTH", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoly_final, "WZONE", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
	arcpy.AddField_management(outpoly_final, "VZONE", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
		
	# Get the raw values from the input points
	# Should put an X/Y field checker here and add/calculate those fields if not present
	for row in InRows:
			length.append(row.getValue("LENGTH"))
			origin_x.append(row.getValue("POINT_X"))
			origin_y.append(row.getValue("POINT_Y"))
			i= i + 1		
	del(InRows,row)
     
	# for maunal starts
	start = 0 #use the first wedge node FID
	end = i # use the last wedge node FID + 1
	#end = i
	n = end - start
	for i in range(start,end):

		#keeping track of time
		startTime= time.time()

		#Create an empty output polygon with the same projection as the input points
		outpoly_temp = os.path.dirname(outpoly_final) + "/out_temp_wedge_%s.shp" % i
		arcpy.CreateFeatureclass_management(os.path.dirname(outpoly_temp),os.path.basename(outpoly_temp), "POLYGON","","DISABLED","DISABLED",proj)
		arcpy.AddField_management(outpoly_temp, "LWV", "TEXT","","", 30, "", "NULLABLE", "NON_REQUIRED")
		arcpy.AddField_management(outpoly_temp, "LENGTH", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
		arcpy.AddField_management(outpoly_temp, "WZONE", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
		arcpy.AddField_management(outpoly_temp, "VZONE", "DOUBLE", 16, 0, "", "", "NULLABLE", "NON_REQUIRED")
		
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

		#Name output tables and Execute ZonalStatisticsAsTable 
		if i > start:
			outtable_veg = os.path.dirname(outpoly_final) + "out_temp_veg_table_%s.dbf" % i
			outtable_ele = os.path.dirname(outpoly_final) + "out_temp_ele_table_%s.dbf" % i
		else:
			outtable_veg = outtable_veg_final
			outtable_ele = outtable_ele_final
			
		arcpy.sa.ZonalStatisticsAsTable(outpoly_temp, "LWV", VegRaster, outtable_veg, "DATA", "ALL")
		arcpy.sa.ZonalStatisticsAsTable(outpoly_temp, "LWV", EleRaster, outtable_ele, "DATA", "ALL")

		#Append temp table to final table and delete
		if i > start:
			arcpy.Append_management(outtable_veg, outtable_veg_final, "NO_TEST")
			arcpy.Delete_management(outtable_veg)
			arcpy.Append_management(outtable_ele, outtable_ele_final, "NO_TEST")
			arcpy.Delete_management(outtable_ele)

		#Append temp shapefile to final and delete
		arcpy.Append_management(outpoly_temp, outpoly_final, "NO_TEST")
		arcpy.Delete_management(outpoly_temp)

		## Display progress
		endTime = time.time()
		elapsedsec= endTime - startTime
		a = Length[i]
		z = i + 1
		minremain = int(((n-z)* elapsedsec/60) + 1)
		#arcpy.AddMessage("Completed record %s of %s in - %s minutes remaining" % (z, n, minremain))
		print("Completed record %s of %s - %s minutes remaining" % (z, n, minremain))


	del Length
	del origin_x
	del origin_y
	gc.collect()

	#arcpy.AddMessage("Process complete")
	print("Process complete")

#This bit of code captures errors and sends them to the window for viewing

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
