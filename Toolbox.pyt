#######################################################################################
# TTools for ArcGIS 10.1 and 10.2
# python toolbox - v 9.0.0 beta 1
# Ryan Michie

# This is the master TTools code for the python toolbox used in ArcGIS 10.1 and 10.2. The code for each 
# tool class comes from the individual python files.

# Step 1: Create Stream Nodes  version 0.92
# Step 2: Measure Channel Widths TBA
# Step 3: Sample Stream Elevations/ Gradient TBA
# Step 4: Measure Topographic Angles - v 0.92
# Step 5: Sample Landcover - Star Pattern, Point Method v 0.98
# Output To csv - v 0.91

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

class Toolbox:
    def __init__(self):
        """TTools is a series of ArcGIS arcsripts used to sample geospatial data and assemble high-resolution inputs for the Heat Source model or other water quality analysis."""
        self.label = "TTools"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Step1_Create_Stream_Nodes]
        #self.tools = [Step1_Create_Stream_Nodes, Step4_Measure_Topographic_Angles, Step5_Sample_Landcover_PointMethod, Output_To_csv]

class Step1_Create_Stream_Nodes:
    def __init__(self):
        """This script will take an input polyline feature with unique stream IDs and generate evenly spaced points along each unique stream ID line at a user defined spacing measured from the downstream endpoint"""
        self.label = "Step1_Create_Stream_Nodes"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self, parameters):
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

        node_dx =arcpy.Parameter(
            name="node_dx",
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

        parameters = [inLine, SIDname, node_dx, outpoint_final]

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
        SIDname = parameters[1].valueAsText
        node_dx = parameters[2].valueAsText
        outpoint_final = parameters[3].valueAsText
        
        #enable garbage collection
        gc.enable()

        #keeping track of time
        startTime= time.time()   
    
        # Create the stream nodes and return them as a list
        NodeList = CreateNodeList(inLine)
    
        #sort the list by stream ID and then stream km
        NodeList = sorted(NodeList, key=itemgetter(1,2), reverse=True)
    
        # Get the spatial projecton of the input stream lines
        proj = arcpy.Describe(inLine).SpatialReference
    
        # Create the output point feature class with the nodes list
        CreateNodesFC(NodeList, outpoint_final, SIDname, proj)
    
        gc.collect()
    
        endTime = time.time()
        elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
        mspernode = timedelta(seconds=(endTime - startTime) / len(NodeList)).microseconds
        print("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))
        #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))       

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

        NodesFC =arcpy.Parameter(
            name="NodesFC",
            displayName="Input TTools point feature class",
            direction="Input",
            datatype="GPFeatureLayer",
            parameterType="Required")	

        multiplecsv = arcpy.Parameter(
            name="multiplecsv",
            displayName="Create seperate csv files for each STREAM_ID",
            direction="Input",
            datatype="GPBoolean",
            parameterType="Required")

        outcsv_dir = arcpy.Parameter(
            name="outcsv_dir",
            displayName="Path directory where output csv file will be saved",
            direction="Input",
            datatype="GPString",
            parameterType="Required")

        outcsv_file = arcpy.Parameter(
            name="outcsv_file",
            displayName="Name of the csv file",
            direction="Input",
            datatype="GPString",
            parameterType="Required")	

        parameters = [NodesFC, multiplecsv, outcsv_dir, outcsv_file]

        return parameters

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

        NodesFC = parameters[0].valueAsText
        multiplecsv = parameters[1].valueAsText
        outcsv_dir = parameters[2].valueAsText
        outcsv_file = parameters[3].valueAsText	

        def ReadNodesFC(NodesFC, readfields):
            """Reads an input point file and returns the NODE ID and X/Y coordinates as a nested dictionary"""
            pnt_dict = NestedDictTree()
            Incursorfields = ["STREAM_ID","NODE_ID"] + readfields
            # Determine input point spatial units
            proj = arcpy.Describe(NodesFC).spatialReference
            with arcpy.da.SearchCursor(NodesFC,Incursorfields,"",proj) as Inrows:
                for row in Inrows:
                    for f in xrange(0,len(readfields)):
                        pnt_dict[row[0]][row[1]][readfields[f]] = row[2+f]
            return(pnt_dict)

        def write_csv(csvlist, csvfile):
            """write the input list to csv"""
            with open(outcsv_final, "wb") as f:
                writer = csv.writer(f)
                writer.writerows(csv_out)    
        
        def NestedDictTree(): 
            """Build a nested dictionary"""
            return defaultdict(NestedDictTree)
        
        try:
            #keeping track of time
            startTime= time.time()
            arcpy.AddMessage("Export to csv") 

            removelist = [u"OBJECTID",u"Id",u"Shape",u"ELEVATION",u"GRADIENT",u"NUM_DIR",u"NUM_ZONES",u"SAMPLE_DIS"]

            # Get all the column headers in the point file and remove the ones in removelist
            header = [field.name for field in arcpy.Describe(NodesFC).fields]
            header_clean = [h for h in header if h not in removelist]

            NODES = read_pointfile(NodesFC, header_clean)

            # make a wide format list by node ID from the nested dictionary
            if multiplecsv == "True":
                for streamID in NODES:
                    csv_out = [[NODES[streamID][nodeID][h] for h in header_clean] for nodeID in NODES[streamID]]
                    outcsv_final = outcsv_dir + "\\" + outcsv_file.replace(".csv", "") + "_" + str(streamID) + ".csv"

                    #sort the list by stream km
                    csv_out = sorted(csv_out, key=itemgetter(1), reverse=True)

                    # Add the header row
                    csv_out.insert(0,header_clean)	    

                    # write it
                    write_csv(csv_out,outcsv_final)

            else:
                csv_out = [[NODES[streamID][nodeID][h] for h in header_clean] for streamID in NODES for nodeID in NODES[streamID]]
                outcsv_final = outcsv_dir+ "\\" + outcsv_file

                #sort the list by stream ID and then stream km
                csv_out = sorted(csv_out, key=itemgetter(1,2), reverse=True)	

                # Add the header row
                csv_out.insert(0,header_clean)

                # write it
                write_csv(csv_out,outcsv_final)

            gc.collect()

            endTime = time.time()
            elapsedmin= (endTime - startTime) / 60	
            arcpy.AddMessage("Process Complete at %s, %s minutes" % (endTime, elapsedmin))    

        # For arctool errors
        except arcpy.ExecuteError:
            msgs = arcpy.GetMessages(2)
            arcpy.AddError(msgs)

        # For other errors
        except:
            import traceback, sys
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]

            pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)	
        return