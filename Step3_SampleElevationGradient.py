#######################################################################################
# TTools
# Step 3: Sample Stream Elevations/ Gradient - v 0.1
# Ryan Michie

# Sample_ElevationsGradient will take an input point feature (from Step 1) and sample the input raster elevation
# to find the lowest elevation in a user defined search radius and calculate the gradient for each node in the downstream direction.

# INPUTS
# 0: Input TTools point feature class (inPoint)
# 1: input the number of samples around the node (searchCells) 1. [0],  2. [9], 3. [25]
# 2: input flag for smoothing if gradient is zero or negative (smooth_flag) 1. True, 2. False
# 2: input elevation raster (EleRaster)
# 3: input elevation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"

# OUTPUTS
# point feature class (edit inPoint) - Added fields with ELEVATION and GRADIENT for node

# Future Updates

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
import itertools
from arcpy import env
from math import sqrt, pow, ceil
from collections import defaultdict

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#inPoint = parameters[0].valueAsText
#searchCells = parameters[1].valueAsText # Needs to be a int
#smooth_flag = parameters[2].valueAsText # Needs to be a int
#EleRaster = parameters[3].valueAsText
#EleUnits = parameters[4].valueAsText

# Start Fill in Data
inPoint = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
searchCells = 9
smooth_flag = "True"
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = 1
# End Fill in Data

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def DistanceBetweenPoints(Xa,Ya,Xb,Yb):
    """Determines the distance between two points using the pythagorean theorem and returns distance between them in map units"""
    from math import sqrt, pow
    dist = math.sqrt(math.pow((Xa - Xb),2) + math.pow((Ya - Yb),2))
    return dist

def ReadPointFile(pointfile):
    """Reads the input point file and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y",]
    # Determine input point spatial units
    proj = arcpy.Describe(inPoint).spatialReference
    with arcpy.da.SearchCursor(pointfile,Incursorfields,"",proj) as Inrows:
	for row in Inrows:
	    pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
	    pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
	    pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
    return(pnt_dict)

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

def GetElevation(samplexy, EleRaster, LowElev, eleZ_to_m):
    """"""
    thevalue = arcpy.GetCellValue_management(EleRaster, samplexy)
    if str(thevalue.getOutput(0)) == "NoData":
	offRasterSamples = offRasterSamples + 1
    else:
	SampleElev= float(thevalue.getOutput(0)) *  eleZ_to_m
		
    #See if the sample elevation is the lowest for this stream node
    if SampleElev < LowElev:
	LowElev = SampleElev
    return(LowElev)

def CalculateGradient(LowElev,LowElevUp):
    return (grad)

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
    
    with arcpy.da.UpdateCursor(pointfile,["STREAM_ID","NODE_ID"] + AddFields) as cursor:
	for row in cursor:
	    for f in xrange(0,len(AddFields)):
		streamID = row[0]
		nodeID =row[1]
		row[f+2] = pointDict[streamID][nodeID][AddFields[f]]
		cursor.updateRow(row)

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
    
    # Make a list of the base x/y coordinate movments. These values will be multipled by the cell size.
    if searchCells == 0: 
	cellcoords = list(itertools.product([0],[0]))
	
    if searchCells == 9:
	cellcoords = list(itertools.product([-1,0,1],[-1,0,1]))

    if searchCells == 25:
	cellcoords = list(itertools.product([-2,-1,0,1,2],[-2,-1,0,1,2]))
        
    # read the data into a nested dictionary
    NODES = ReadPointFile(inPoint)
    n = 0
    for streamID in NODES:
	SkipDownNodes = [1]
	
	for nodeID in NODES[streamID]:
	    node_x = float(NODES[streamID][nodeID]["POINT_X"])
	    node_y = float(NODES[streamID][nodeID]["POINT_Y"])
	    LowElev = 99999
	    
	    for coord in cellcoords:
		# Calculate the cell X/Y based on the base coordinate movement
		cell_x = node_x + (coord[0] * CellSize)
		cell_y = node_y + (coord[1] * CellSize)
		
		samplexy = str(cell_x) + " " + str(cell_y) # xy string requirement for arcpy.GetCellValue
		
		#GetElevation(samplexy, EleRaster, LowElev, eleZ_to_m)
		# Sample the elevation value from the elevation raster
		thevalue = arcpy.GetCellValue_management(EleRaster, samplexy)
		if str(thevalue.getOutput(0)) == "NoData":
		    offRasterSamples = offRasterSamples + 1
		else:
		    SampleElev= float(thevalue.getOutput(0)) *  eleZ_to_m
		
		#See if the sample elevation is the lowest for this stream node
		if SampleElev < LowElev:
		    LowElev = SampleElev
	    
	    # Save the elevation back into the dictionary
	    NODES[streamID][nodeID]["ELEVATION"] = LowElev
	    
	    # CalculateGradient
	    if nodeID > 0:
		# Calculate gradient between the current node and downstream node
		Ele = LowElev
		EleDown = NODES[streamID][nodeID - 1]["ELEVATION"]
		dx_meters = float(NODES[streamID][nodeID]["STREAM_KM"] - NODES[streamID][nodeID - 1]["STREAM_KM"]) * 1000
		#GradientDown = (Ele - EleDown) / (dx_meters)
		# Check if the gradient is <= 0, if yes keep going until
		# it is positive again. Then recalculate the gradient over the longer distance		
		if Ele <= EleDown and smooth_flag == "True":
		    EleDownSkip = EleDown
		    SkipDownNodes.append(max(SkipDownNodes) + 1)
		    NODES[streamID][nodeID -1]["GRADIENT"] = GradientDown
		else:
		    EleDownSkip = NODES[streamID][nodeID - max(SkipDownNodes)]["ELEVATION"]
		    dx_meters = float(NODES[streamID][nodeID]["STREAM_KM"] - NODES[streamID][nodeID - max(SkipDownNodes)]["STREAM_KM"]) * 1000
		    GradientDown = (Ele - EleDownSkip) / dx_meters
		    for Skip in SkipDownNodes:
			NODES[streamID][nodeID - Skip]["GRADIENT"] = GradientDown
		    SkipDownNodes = [1]

    # Write the Elevation and Gradient to the TTools point feature class
    AddFields = ["ELEVATION", "GRADIENT"]
    UpdatePointFile(NODES, inPoint, AddFields)	    

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