#######################################################################################
# TTools for ArcGIS 10.1 and 10.2
# python toolbox - v 9.0.0 beta 1
# Ryan Michie

# CHANGE THIS TO Toolbox.pyt

# This is the master TTools code for the python toolbox used in ArcGIS 10.1 and 10.2. The code for each 
# tool class comes from the individual python files.

# Step 1: Create Stream Nodes  version 0.92
# Step 2: Measure Channel Widths TBA
# Step 3: Sample Stream Elevations/ Gradient TBA
# Step 4: Measure Topographic Angles - v 0.92
# Step 5: Sample Landcover - Star Pattern, Point Method v 0.98
# Convert to csv - v 0.91

# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
import arcpy
from arcpy import env
from __future__ import division
import sys, os, string, gc, shutil, time
from math import radians, sin, cos, atan
from collections import defaultdict
from operator import itemgetter
import csv

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
from arcpy.management import *

env.overwriteOutput = True

class TTools(object):
    def __init__(self):
        """TTools is a series of ArcGIS arcsripts used to sample geospatial data and assemble high-resolution inputs for the Heat Source model or other water quality analysis."""
        self.label = "TTools"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Step1_Create_Stream_Nodes, Step4_Measure_Topographic_Angles, Step5_Sample_Landcover_PointMethod, Output_To_csv]

class Step1_Create_Stream_Nodes(object):
    def __init__(self):
        """This script will take an input polyline feature with unique stream IDs and generate evenly spaced points along each unique stream ID line at a user defined spacing measured from the downstream endpoint"""
        self.label = "Step1_Create_Stream_Nodes"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        inLine =arcpy.Parameter(
            name="inLine",
            displayName="Input Stream Centerlines/s",
            direction="Input",
            datatype="GPFeatureLayer",
            parameterType="Required")
        inLine.filter.list = ["LINE"])
        
        SIDname =arcpy.Parameter(
            name="SIDname",
            displayName="Stream Identifier Field Name",
            direction="Input",
            datatype="Field",
            parameterType="Required")      
        
        NodeDistance =arcpy.Parameter(
            name="NodeDistance",
            displayName="Desired distance between nodes (meters)",
            direction="Input",
            datatype="Double",
            parameterType="Required")
        
        outpoint_final =arcpy.Parameter(
            name="outpoint_final",
            displayName="Output point features",
            direction="Output", 
            datatype="DEFeatureClass",
            parameterType="Required")
        
        parameters = [inLine, SIDname, NodeDistance, outpoint_final]
        
        return parameters

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        if parameters[0].altered:
            parameters[1].value = arcpy.ValidateFieldName(parameters[1].value, parameters[0].value)        
        
        if parameters[1].altered:
            parameters[2].value = arcpy.ValidateFieldName(parameters[2].value, parameters[1].value)        
        
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        inLine = parameters[0].valueAsText
        StreamIDfield = parameters[1].valueAsText
        NodeDistance = parameters[2].valueAsText
        outpoint_final = parameters[3].valueAsText        
        
	def CreateNodes(inLine):
	    """Reads an input stream centerline file and returns the NODE ID, STREAM ID, and X/Y coordinates as a list"""
	    NODES = []
	    Incursorfields = ["SHAPE@",SIDname]
	    NID = 0
	    # Determine input point spatial units
	    proj = arcpy.Describe(inLine).spatialReference
	    with arcpy.da.SearchCursor(inLine,Incursorfields,"",proj) as Inrows:
		print("Creating Nodes")	
		for row in Inrows:
		    LineLength = row[0].getLength("PRESERVE_SHAPE")
		    numNodes = int(LineLength / NodeDistance)
		    nodes = range(0,numNodes+1)
		    arcpy.SetProgressor("step", "Creating Nodes", 0, numNodes+1, 1)
		    positions = [n * NodeDistance / LineLength for n in nodes] # list of Lengths in meters
		    for position in positions:
			node = row[0].positionAlongLine(position,True).centroid
			# list of "NODE_ID",STREAM_ID,"STREAM_KM","POINT_X","POINT_Y","SHAPE@X","SHAPE@Y"
			NODES.append((NID, row[1], float(position * LineLength /1000), node.X, node.Y, node.X, node.Y ))
			NID = NID + 1
		    arcpy.SetProgressorPosition()
	    arcpy.ResetProgressor()
	    return(NODES)
	
	def CreatePointFile(pointList,pointfile, SIDname, proj):
	    """Create the output point feature class using the data from the nodes list"""
	    arcpy.AddMessage("Exporting Data")
	    
	    # Determine Stream ID field properties
	    SIDtype = arcpy.ListFields(inLine,SIDname)[0].type
	    SIDprecision = arcpy.ListFields(inLine,SIDname)[0].precision
	    SIDscale = arcpy.ListFields(inLine,SIDname)[0].scale
	    SIDlength = arcpy.ListFields(inLine,SIDname)[0].length    
	    
	    #Create an empty output with the same projection as the input polyline
	    cursorfields = ["NODE_ID","STREAM_ID","STREAM_KM","LONGITUDE","LATITUDE","SHAPE@X","SHAPE@Y"]
	    arcpy.CreateFeatureclass_management(os.path.dirname(pointfile),os.path.basename(pointfile), "POINT","","DISABLED","DISABLED",proj)
	
	    # Add attribute fields
	    arcpy.AddField_management(pointfile, "NODE_ID", "LONG", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	    arcpy.AddField_management(pointfile, "STREAM_ID", SIDtype, SIDprecision, SIDscale, SIDlength, "", "NULLABLE", "NON_REQUIRED")
	    arcpy.AddField_management(pointfile, "STREAM_KM", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	    arcpy.AddField_management(pointfile, "LONGITUDE", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
	    arcpy.AddField_management(pointfile, "LATITUDE", "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")
		    
	    with arcpy.da.InsertCursor(pointfile, cursorfields) as cursor:
		for row in pointList:
		    cursor.insertRow(row)
	    
	    #Change X/Y from input spatial units to decimal degrees
	    proj_dd = arcpy.SpatialReference(4269) #GCS_North_American_1983 
	    with arcpy.da.UpdateCursor(pointfile,["SHAPE@X","SHAPE@Y","LONGITUDE","LATITUDE"],"",proj_dd) as cursor:
		for row in cursor:
		    row[2] = row[0] # LONGITUDE
		    row[3] = row[1] # LATITUDE
		    cursor.updateRow(row)
	
	#enable garbage collection
	gc.enable()
	
	#keeping track of time
	startTime= time.time()   
		
	# Create the stream nodes and return them as a list
	NODES = CreateNodes(inLine)
	    
	#sort the list by stream ID and then stream km
	NODES = sorted(NODES, key=itemgetter(1,2), reverse=True)
		
	# Get the spatial projecton of the input stream lines
	proj = arcpy.Describe(inLine).SpatialReference
		
	# Create the output point feature class with the nodes list
	CreatePointFile(NODES,outpoint_final, SIDname, proj)
		
	gc.collect()
		
	endTime = time.time()
	elapsedmin= (endTime - startTime) / 60	
	arcpy.AddMessage("Process Complete at %s, %s minutes" % (endTime, elapsedmin))        
        
        return

class Step4_Measure_Topographic_Angles(object):
    def __init__(self):
        """Measure_Topographic_Angles will take an input point feature (from Step 1) and calculate the maximum topographic elevation and the the slope angle from each node in different directions."""
        self.label = "Step4_Measure_Topographic_Angles"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return

class Step5_Sample_Landcover_PointMethod(object):
    def __init__(self):
        """Sample_Landcover_PointMethod will take an input point feature (from Step 1) and sample input landcover rasters in a user specificed number of cardianal directions with point samples spaced at a user defined distance moving away from the stream."""

        self.label = "Step5_Sample_Landcover"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return
        
class Output_To_csv(object):
    def __init__(self):
        """Output_To_csv will take the point feature created from using steps 1-5 and output a landcover data csv file formatted for heat source 9"""

        self.label = "Output_To_csv"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return
        
      