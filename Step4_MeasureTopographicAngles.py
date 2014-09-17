#######################################################################################
# TTools
# Step 4: Measure Topographic Angles - v 0.93
# Ryan Michie

# Measure_Topographic_Angles will take an input point feature (from Step 1) and calculate the maximum
# topographic elevation and the the slope angle from each node in different directions.

# INPUTS
# 0: Input TTools point feature class (inPoint)
# 1: input the directions to sample (Directions) 1. [W,S,E], 2. [NE,E,SE,S,SW,W,NW]
# 2: input the maximum km distance to search (MaxSearchDistance_km)
# 3: input elevation raster (EleRaster)
# 4: input elevation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 5: output sample point file name/path (outpoint_final)

# OUTPUTS
# point feature class (edit inPoint) - Added fields with topographic shade angles for each direction at each node
# point feature class (new) - point for each x/y location of the maximum elevation.

# Future Updates
#

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from arcpy import env
from math import degrees, radians, sin, cos, atan, ceil
from collections import defaultdict
#from operator import itemgetter

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#inPoint = parameters[0].valueAsText
#Directions = parameters[1].valueAsText # Needs to be a long
#MaxSearchDistance_km = parameters[2].valueAsText
#EleRaster = parameters[3].valueAsText
#EleUnits = parameters[4].valueAsText
#outpoint_final = parameters[5].valueAsText

# Start Fill in Data
inPoint = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
Directions = 1
MaxSearchDistance_km = .005
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = 1
outpoint_final = r"D:\Projects\TTools_9\Example_data.gdb\topo_samplepoint"
# End Fill in Data

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def ReadPointFileOLD(pointfile):
    """Reads the input point file and returns the NODE_ID and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["NODE_ID", "SHAPE@X","SHAPE@Y"]
    # Determine input point spatial units
    proj = arcpy.Describe(inPoint).spatialReference
    with arcpy.da.SearchCursor(pointfile,Incursorfields,"",proj) as Inrows:
	for row in Inrows:
	    pnt_dict[row[0]]["POINT_X"] = row[1]
	    pnt_dict[row[0]]["POINT_Y"] = row[2] 
    return(pnt_dict)

def ReadPointFile(pointfile, OverwriteData, AddFields):
    """Reads the input point file and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y",]
    
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
		# if the data is null or zero (0 = default for shapefile), it is retreived and will be overwritten.
		if row[5] == None or row[5] == 0:
		    pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
		    pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
		    pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4] 

def CreateTopoPointFile(pointList, pointfile, proj):
    """Creates the output topo point feature class using the data from the nodes list"""
    arcpy.AddMessage("Exporting Data")
    print("Exporting Data")
    
    #Create an empty output with the same projection as the input polyline
    cursorfields = ["NODE_ID","AZIMUTH","TOPOANGLE","TOPO_ELE","ELE_CHANGE","DISTANCE_m","SEARCH_m","NA_SAMPLES","POINT_X","POINT_Y"]
    arcpy.CreateFeatureclass_management(os.path.dirname(pointfile),os.path.basename(pointfile), "POINT","","DISABLED","DISABLED",proj)

    # Add attribute fields # TODO add dictionary of field types so they aren't all double
    for f in cursorfields:
	if f == "AZIMUTH":
	    arcpy.AddField_management(pointfile, f, "TEXT", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	else:
	    arcpy.AddField_management(pointfile, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	    
    with arcpy.da.InsertCursor(pointfile, cursorfields + ["SHAPE@X","SHAPE@Y"]) as cursor:
	for row in pointList:
	    cursor.insertRow(row)

def UpdatePointFile(pointDict, pointfile, AddFields): 
    """Updates the input point feature class with data from the nodes dictionary"""
    
    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(pointfile):
	ExistingFields.append(f.name)     
    
    # Check to see if the field exists and add it if not
    for f in AddFields:
	if (f in ExistingFields) == False:
	    arcpy.AddField_management(pointfile, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")      
    
    with arcpy.da.UpdateCursor(pointfile,["NODE_ID"] + AddFields) as cursor:
	for row in cursor:
	    for f in xrange(0,len(AddFields)):
		node =row[0]
		row[f+1] = pointDict[node][AddFields[f]]
		cursor.updateRow(row)

def UpdatePointFile2(pointDict, pointfile, AddFields): 
    """Updates the input point feature class with data from the nodes dictionary"""
    # Add attribute fields # TODO add a check to se if the field already exists. if yes ask to overwrite.
    for f in AddFields:
	arcpy.AddField_management(pointfile, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	
    whereclause = """%s = %s""" % arcpy.AddFieldDelimiters(pointfile, "NODE_ID"), nodeID
    
    with arcpy.da.UpdateCursor(pointfile,["NODE_ID"] + AddFields, whereclause) as cursor:
	for row in pointDict:
	    for f in xrange(0,len(AddFields)):
		node =row[0]
		row[f+1] = pointDict[node][AddFields[f]]
		cursor.updateRow(row)

def ToMetersUnitConversion(inFeature):
    """Returns the conversion factor to get from the input spatial units to meters"""
    unitCode = arcpy.Describe(inFeature).SpatialReference.linearUnitCode
    if unitCode == 9001: #International meter
	units_con = 1 
    if unitCode == 9002: #International foot
	units_con = 0.3048
    if unitCode == 9003: #US Survey foot
        units_con = 1200/3937
    if unitCode == 9005: #Clarke's foot
	units_con =  0.3047972654 
    if unitCode not in [9001,9002,9003,9005]:
	arcpy.AddError("{0} has an unrecognized spatial reference. Use projection with units of feet or meters.".format(inFeature))
	system.exit("Unrecognized spatial reference. Use projection with units of feet or meters.")
    return units_con

#enable garbage collection
gc.enable()
  
try:
    #keeping track of time
    startTime= time.time()    
    
    # Determine input point spatial units
    proj = arcpy.Describe(inPoint).spatialReference
    proj_ele = arcpy.Describe(EleRaster).spatialReference
	
    # Check to make sure the raster and input points are in the same projection.
    if proj.name != proj_ele.name:
	arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(inPoint,EleRaster))
	sys.exit("Input points and elevation raster do not have the same projection. Please reproject your data.")
	# reproject the inpoints
	#inPoint_rp = os.path.dirname(inPoint) + "\\rp_" + os.path.basename(inPoint)
	#arcpy.Project_management(inPoint,inPoint_rp,proj_ele)
	#inPoint = inPoint_rp	

    nodexy_to_m = ToMetersUnitConversion(inPoint)
    Elexy_to_m = ToMetersUnitConversion(EleRaster)
    units_con=  nodexy_to_m / Elexy_to_m
    MaxSearchDistance = MaxSearchDistance_km * 1/nodexy_to_m * 1000
	
    # Determine the elevation Z units conversion into meters
    if EleUnits == 1: # Feet
	eleZ_to_m = 0.3048
    if EleUnits == 2: # Meters
	eleZ_to_m = 1
    if EleUnits == 3: # Other
	sys.exit("Please modify your raster elevtion units to feet or meters.") 	
	
    # Get the elevation raster cell size
    CellSizeResult = arcpy.GetRasterProperties_management(EleRaster, "CELLSIZEX")
    CellSize = float(CellSizeResult.getOutput(0))	      
    if Directions == 2: # All directions
	azimuths = [45,90,135,180,225,270,315]
    else:        
	azimuths = [270,180,90]	
    
    azimuthdict = {45:"NE",90:"E",135:"SE",180:"S",225:"SW",270:"W",315:"NW",365:"N"}
	
    # read the data into a nested dictionary
    NODES = ReadPointFile(inPoint)
    TOPO_fc = []
    n = 0
    elapsedNodemin = 0.2
    for nodeID in NODES:
	startNodeTime= time.time()
	print("Processing Node %s of %s. %s minutes remaining" % (n+1, len(NODES), elapsedNodemin * (len(NODES) - n)))
	# Set the progressor
	arcpy.SetProgressor("step", "",0, len(NODES), 1)
	arcpy.SetProgressorLabel("Processing Node %s of %s. %s minutes remaining" % (n+1, len(NODES), elapsedNodemin * (len(NODES) - n)))
	#arcpy.AddMessage("Processing Node %s of %s" % (nodeID+1, len(NODES)))
	node_x = float(NODES[nodeID]["POINT_X"])
	node_y = float(NODES[nodeID]["POINT_Y"])
    
	nodexy = str(node_x) + " " + str(node_y) # xy string requirement for arcpy.GetCellValue	
	thevalue = arcpy.GetCellValue_management(EleRaster, nodexy)
	if str(thevalue.getOutput(0)) == "NoData":
	    sys.exit("There is no elevation to sample at sample node %s Please check your elevation raster and the spatial coordinates of your input data.") % (SreamNode_ID)
	else:
	    sampleZ= float(thevalue.getOutput(0)) *  eleZ_to_m		
	    nodeZ = float(thevalue.getOutput(0)) *  eleZ_to_m		
	for a in xrange(0,len(azimuths)):
	    SearchDistance = 0
	    MaxShadeAngle = 0
	    offRasterSamples = 0
	    i = 0
	    MaxZChange = 0
	    SampleZFinal = nodeZ
	    MaxShadeAngle = 0
	    FinalSearchDistance = 0
	    MaxShadeAngle_X = node_x
	    MaxShadeAngle_Y = node_y		  
	    while not SearchDistance > MaxSearchDistance:
		# This is the skippy algorithm from Greg Pelletier
		if i <= 10:
		    SearchDistance = SearchDistance + (CellSize * units_con)
		if 10 < i <= 20:
		    SearchDistance = SearchDistance + (CellSize * units_con * 3)
		if 20 < i <= 40:
		    SearchDistance = SearchDistance + (CellSize * units_con * 6)
		if 40 < i <= 50:
		    SearchDistance = SearchDistance + (CellSize * units_con * 12)
		if 50 < i <= 60:
		    SearchDistance = SearchDistance + (CellSize * units_con * 25)
		if i > 60:
		    SearchDistance = SearchDistance + (CellSize * units_con * 50)
		i = i + 1			  
			  
		# Calculate the x and y sample location.
		sample_x = (SearchDistance * sin(radians(azimuths[a]))) + node_x
		sample_y = (SearchDistance * cos(radians(azimuths[a]))) + node_y
		samplexy = str(sample_x) + " " + str(sample_y) # xy string requirement for arcpy.GetCellValue
					
		# Sample the elevation value from the elevation raster
		thevalue = arcpy.GetCellValue_management(EleRaster, samplexy)
		if str(thevalue.getOutput(0)) == "NoData":
		    offRasterSamples = offRasterSamples + 1
		else:
		    sampleZ= float(thevalue.getOutput(0)) *  eleZ_to_m
		
		#calculate the topographic shade angle (in degrees)
		ShadeAngle = degrees(atan((sampleZ - nodeZ) / SearchDistance * nodexy_to_m))
		if ShadeAngle > MaxShadeAngle:
		    MaxZChange = sampleZ - nodeZ
		    SampleZFinal = sampleZ
		    MaxShadeAngle = ShadeAngle
		    FinalSearchDistance = SearchDistance * nodexy_to_m
		    MaxShadeAngle_X = sample_x
		    MaxShadeAngle_Y = sample_y
	    
	    TOPO_fc.append([nodeID, azimuthdict[azimuths[a]], MaxShadeAngle, SampleZFinal, MaxZChange,FinalSearchDistance, MaxSearchDistance * nodexy_to_m, offRasterSamples, MaxShadeAngle_X, MaxShadeAngle_Y, MaxShadeAngle_X,MaxShadeAngle_Y])
	    NODES[nodeID]["TOPO_"+ str(azimuthdict[azimuths[a]])]= MaxShadeAngle
	endNodeTime = time.time()
	elapsedNodemin = ceil(((endNodeTime - startNodeTime) / 60)* 10)/10	
	arcpy.SetProgressorPosition()
	n = n + 1
    arcpy.ResetProgressor()
    del(i,nodeID,a,ShadeAngle,MaxShadeAngle,MaxZChange, MaxShadeAngle_X,MaxShadeAngle_Y,thevalue,SampleZFinal,offRasterSamples, FinalSearchDistance)			
    gc.collect()
	
    ####################################################################################################### 
    # Build the output point feature class using the data from the NODES list
    CreateTopoPointFile(TOPO_fc,outpoint_final,proj)

    ####################################################################################################### 
    # Write the topo angles to the TTools point feature class
    AddFields = ["TOPO_"+ azimuthdict[a] for a in azimuths]
    UpdatePointFile(NODES, inPoint, AddFields)		
    
    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10	
    print("Process Complete in %s minutes" % (elapsedmin))
    arcpy.AddMessage("Process Complete in %s minutes" % (elapsedmin))

	
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