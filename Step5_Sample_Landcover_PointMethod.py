#######################################################################################
# TTools
# Step 5: Sample Landcover - Star Pattern, Point Method v 0.98
# Ryan Michie

# Sample_Landcover_PointMethod will take an input point feature (from Step 1) and sample input landcover rasters
# in a user specificed number of cardianal directions with point samples spaced at a user defined distance
# moving away from the stream.

# INPUTS
# 0: input TTools point file (inPoint)
# 1: input number of directions to sample (NumDirections)
# 2: input number of transverse vegetation samples in each azimuth direction (NumZones)
# 3: include stream sample (True/False)
# 4: input The distance between transverse samples (TransDistance)
# 5: input canopy data type. 1."Codes", 2."CanopyCover", or 3."LAI" (CanopyData)
# 6: input landcover code or height raster (LCRaster)
# 7: input (optional) canopy cover or LAI raster (CanopyRaster)
# 8: input (optional) k coeffcient raster (kRaster)
# 9: input elevation raster (EleRaster)
# 10: input elvation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 11: output sample point file name/path (outpoint_final)

# OUTPUTS
# point feature class (edit inPoint) - added fields with Landcover and elevation data for each azimuth direction at each node
# point feature class (new) - point at each x/y sample and the sample raster values

# Future Updates

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from arcpy import env
from math import radians, sin, cos
from collections import defaultdict
from operator import itemgetter

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
from arcpy.management import *

env.overwriteOutput = True

# Parameter fields for python toolbox
#inPoint = parameters[0].valueAsText
#NumDirections = parameters[1].valueAsText # LONG
#NumZones = parameters[2].valueAsText # LONG
#StreamSample = parameters[3].valueAsText True/False
#TransDistance = parameters[4].valueAsText # LONG
#CanopyData = = parameters[5].valueAsText One of these: 1."Codes", 2."CanopyCover", or 3."LAI"
#LCRaster = parameters[6].valueAsText # This is either landcover height or codes
#CanopyRaster = parameters[7].valueAsText # OPTIONAL This is either canopy cover or LAI raster
#kRaster = parameters[8].valueAsText # OPTIONAL The k value raster for LAI
#EleRaster = parameters[9].valueAsText
#EleUnits = parameters[10].valueAsText
#outpoint_final = aparameters[11].valueAsText
    
# Start Fill in Data
inPoint = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
NumDirections = 8
NumZones = 4 
StreamSample = 'True' # include stream sample in number of zones? (True/False)
TransDistance = 8
LCRaster = r"D:\Projects\TTools_9\Example_data.gdb\veght_lidar_ft" # This is either landcover height or codes
CanopyDataType = "Codes"
CanopyRaster = "" # OPTIONAL This is either canopy cover or a LAI raster
kRaster = "" # OPTIONAL This is the k value for LAI
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = "Feet"
outpoint_final = r"D:\Projects\TTools_9\Example_data.gdb\LC_samplepoint"
# End Fill in Data

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def ReadPointFile(pointfile):
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

def CreateLCPointFile(pointList, LCFields, pointfile, proj):
    """Creates the output landcover sample point feature class using the data from the nodes list"""
    arcpy.AddMessage("Exporting Data")
    print("Exporting Data")
    
    arcpy.CreateFeatureclass_management(os.path.dirname(pointfile),os.path.basename(pointfile), "POINT","","DISABLED","DISABLED",proj)
    
    AddFields = ["NODE_ID","AZIMUTH","ZONE"] + LCFields + ["POINT_X","POINT_Y"]

    # Add attribute fields # TODO add dictionary of field types so they aren't all double
    for f in AddFields:
	arcpy.AddField_management(pointfile, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	    
    with arcpy.da.InsertCursor(pointfile, AddFields + ["SHAPE@X","SHAPE@Y"]) as cursor:
	for row in pointList:
	    cursor.insertRow(row)

def UpdatePointFile(pointDict, pointfile, AddFields): 
    """Updates the input point feature class with data from the nodes dictionary"""
    # Add attribute fields # TODO add a check to se if the field already exists. if yes ask to overwrite.
    for f in AddFields:
	arcpy.AddField_management(pointfile, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")    
    
    with arcpy.da.UpdateCursor(pointfile,["NODE_ID"] + AddFields) as cursor:
	for row in cursor:
	    for f in range(0,len(AddFields)):
		node =row[0]
		row[f+1] = pointDict[node][AddFields[f]]
		cursor.updateRow(row)

def FromMetersUnitConversion(inFeature):
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
	system.exit("Unrecognized spatial reference. Use projection with units of feet or meters.")
    return units_con

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()
    
    NODES = ReadPointFile(inPoint)	
	
    # Determine input spatial units and set conversion factor to get from meters to the input spatial units
    proj = arcpy.Describe(inPoint).SpatialReference
    units_con = FromMetersUnitConversion(inPoint)
	
    # Determine which type of rasters we are using. This helps build the dictionary keys and ouput headers
    if CanopyDataType == "Codes":        
	type = ['LC','ELE']
    if CanopyDataType == "LAI":  #Use LAI methods
	type = ['LC','LAI','k','ELE']
    if CanopyDataType == "CanopyCover":  #Use Canopy Cover methods
	type = ['LC','CCV','ELE']
    
    # Set the converstion factor to get from the input elevation z units to meters
    if EleUnits == "Meters": #International meter
	ele_con = 1 
    if EleUnits == "Feet": #International foot
	ele_con = 0.3048
    if EleUnits == "Other": #Some other units
	sys.exit("Please modify your raster elevation z units so they are either in meters or feet")	
	
    if NumDirections == 999:  #999 is a flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
	azimuths = [45,90,135,180,225,270,315]
    else:        
	azimuths = [x * 360.0 / NumDirections for x in range(1,NumDirections+ 1)]
	
    zone = range(1,int(NumZones+1))
    # TODO This is a future function that may replace the emergent methods.
    # this method there is a regular landcover sample for each azimuth direction at the stream node.
    #if StreamSample == "TRUE":
	#zone = range(0,int(NumZones))
    #else:
	#zone = range(1,int(NumZones+1))
	
    #arcpy.AddMessage("Extracting values")
    print("Extracting raster values")	

    for n in NODES:
	print("Processing Node %s of %s" % (n+1, len(NODES)))
	# Set the progressor
	arcpy.SetProgressor("step","",0, len(NODES), 1)
	arcpy.SetProgressorLabel("Process Node %s of %s" % (n+1, len(NODES)))	
	#arcpy.AddMessage("Process Node %s of %s" % (n+1, len(NODES)))
	
	#Add input metadata
	NODES[n]["NUM_DIR"] = NumDirections
	NODES[n]["NUM_ZONES"] = NumZones
	NODES[n]["SAMPLE_DIS"] = TransDistance
	NODES[n]["LC_RASTER"] = LCRaster
	NODES[n]["ELE_RASTER"] = EleRaster
	NODES[n]["CAN_RASTER"] = CanopyRaster
	NODES[n]["k_RASTER"] = kRaster
	
	# xy string requirement for arcpy.GetCellValue
	origin_x = NODES[n]["POINT_X"]
	origin_y = NODES[n]["POINT_Y"]
	xypoint = str(origin_x) + " " + str(origin_y) 
	
	# Sample the stream/emergent
	if StreamSample == "True":
	    for t in type:
		LCkey = t+"_EMERGENT"
		# Sample the point value from the appropriate raster and add to NODES dictionary
		if t == "ELE":
		    thevalue = arcpy.GetCellValue_management(EleRaster, xypoint)
		    NODES[n][LCkey] = float(thevalue.getOutput(0)) * ele_con
		    
		if t in ["LC"]:
		    thevalue = arcpy.GetCellValue_management(LCRaster, xypoint)
		    NODES[n][LCkey] = float(thevalue.getOutput(0))
		
		if t in ["LAI","CCV"]:
		    thevalue = arcpy.GetCellValue_management(CanopyRaster, xypoint)
		    NODES[n][LCkey] = float(thevalue.getOutput(0))
		
		if t == "k":
		    thevalue =arcpy.GetCellValue_management(kRaster, xypoint)
		    NODES[n][LCkey] = float(thevalue.getOutput(0))	    
	
	# Sample moving away from the stream
	for a in azimuths:
	    for z in zone:
		# Calculate the x and y coordinate of the landcover sample location
		_X_ = (z * TransDistance * units_con * sin(radians(a))) + origin_x
		_Y_ = (z * TransDistance * units_con * cos(radians(a))) + origin_y
		xypoint = str(_X_) + " " + str(_Y_) # xy string requirement for arcpy.GetCellValue
		
		# Add the coordinates into the NODES dictionary
		NODES[n]["POINT_X_"+str(int(a))+'_Z'+str(z)] = _X_
		NODES[n]["POINT_Y_"+str(int(a))+'_Z'+str(z)] = _Y_			
		
		for t in type:
		    LCkey = t+'_'+str(int(a))+'_Z'+str(z)
		    # Sample the point value from the appropriate raster and add to NODES dictionary
		    if t == "ELE":
			thevalue = arcpy.GetCellValue_management(EleRaster, xypoint)
			NODES[n][LCkey] = float(thevalue.getOutput(0)) * ele_con
		    if t in ["LC"]:
			thevalue = arcpy.GetCellValue_management(LCRaster, xypoint)
			NODES[n][LCkey] = float(thevalue.getOutput(0))
		    if t in ["LAI","CCV"]:
			thevalue = arcpy.GetCellValue_management(CanopyRaster, xypoint)
			NODES[n][LCkey] = float(thevalue.getOutput(0))
		    if t == "k":
			thevalue =arcpy.GetCellValue_management(kRaster, xypoint)
			NODES[n][LCkey] = float(thevalue.getOutput(0))
	arcpy.SetProgressorPosition()
    arcpy.ResetProgressor()
    del(n,a,z,t,thevalue,origin_x,origin_y,xypoint,_X_,_Y_)			
    gc.collect()
    
    ####################################################################################################### 
    # Build the output point feature class using the data from the NODES dictionary
	    
    # these are the fields for the output list. the extra XY at the end is for arcpy's "SHAPE@X","SHAPE@Y"
    type2 = type + ["POINT_X","POINT_Y","POINT_X","POINT_Y"]
    
    # This builds a unique long form list of all node, azimuth, and zone values (naz). 
    # This will be used in the next step to build NODES dictionary keys in order to iterate through each sample 
    # and get the data. Azimuth 0 and Zone 0 are the emergent samples. Those are stored in a seperate list (naz_e).
    naz= [[n,a,z] for n in NODES for a in azimuths for z in zone]
    
    
    # Build the NODE dictonary keys, get the values, and save it as a list. 
    # This list is used to create the output point feature class
    LCpointlist = [naz[row] + [NODES[naz[row][0]][t+'_'+str(int(naz[row][1]))+'_Z'+str(naz[row][2])] for t in type2] for row in range(0,len(naz))]
    # add the emergent samples
    
    if StreamSample == "True":
	type_emergent = [t+"_EMERGENT" for t in type] + ["POINT_X","POINT_Y","POINT_X","POINT_Y"]
	naz_e= [[n,a,z] for n in NODES for a in [0] for z in [0]]
	LCpointlist = LCpointlist + [naz_e[row] + [NODES[naz_e[row][0]][t] for t in type_emergent] for row in range(0,len(naz_e))]
    
    CreateLCPointFile(LCpointlist, type, outpoint_final, proj)
    ####################################################################################################### 
    
    # Create the landcover headers by concatenating the type, azimuth direction, and zone
    
    if StreamSample == "True":
	LCHeaders = ["LC_EMERGENT"]
    else:
	LCHeaders = []
    for t in type:
	for a in range(0,len(azimuths)):
	    for z in range(0,len(zone)):
		if StreamSample == "True" and t in ["CCV","LAI","k"] and a==0 and z==0:
		    LCHeaders.append(t+'_EMERGENT') # add this when we are at the beginning of the LAI/Canopy cover series
		    LCHeaders.append(t+'_'+str(int(azimuths[a]))+'_Z'+str(zone[z]))
		else:
		    LCHeaders.append(t+'_'+str(int(azimuths[a]))+'_Z'+str(zone[z]))     
    
    LCHeaders = LCHeaders + ["NUM_DIR","NUM_ZONES","SAMPLE_DIS"]
    
    # Write the landcover data to the TTools point feature class
    UpdatePointFile(NODES, inPoint, LCHeaders)

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
