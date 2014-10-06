#######################################################################################
# TTools
# Step 5: Sample Landcover - Star Pattern, Point Method v 0.98
# Ryan Michie

# Sample_Landcover_PointMethod will take an input point feature (from Step 1) and sample input landcover rasters
# in a user specificed number of cardianal directions with point samples spaced at a user defined distance
# moving away from the stream.

# INPUTS
# 0: Input TTools point feature class (PointFile)
# 1: input number of transects per node (trans_count)
# 2: input number of samples per transect (transsample_count)
# 3: include stream sample in transect count (True/False)
# 4: input The distance between transect samples (transsample_distance)
# 5: input landcover data type. 1."Codes" 2. "CanopyCover", or 3."LAI" (CanopyData)
# 6: use heatsource 8 methods 1. True 2. False
# 7: input landcover code or height raster (LCRaster)
# 8: input (optional) canopy cover or LAI raster (CanopyRaster)
# 9: input (optional) k coeffcient raster (kRaster)
# 10: input elevation raster (EleRaster)
# 11: input elvation raster z units (EleUnits) 1. "Feet", 2. "Meters" 3. "Other"
# 12: output sample point file name/path (outpoint_final)
# 13: input flag if existing data can be over written (OverwriteData) 1. True, 2. False

# OUTPUTS
# point feature class (edit PointFile) - added fields with Landcover and elevation data for each azimuth direction at each node
# point feature class (new) - point at each x/y sample and the sample raster values

# Future Updates

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys
import os
import string 
import gc
import shutil
import time
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
#PointFile = parameters[0].valueAsText
#trans_count = parameters[1].valueAsText # LONG
#transsample_count = parameters[2].valueAsText # LONG
#StreamSample = parameters[3].valueAsText True/False
#transsample_distance = parameters[4].valueAsText # LONG
#CanopyDataType = parameters[5].valueAsText One of these: 1."CanopyCover", or 2."LAI"
#LCRaster = parameters[6].valueAsText # This is either landcover height or codes
#CanopyRaster = parameters[7].valueAsText # OPTIONAL This is either canopy cover or LAI raster
#kRaster = parameters[8].valueAsText # OPTIONAL The k value raster for LAI
#EleRaster = parameters[9].valueAsText
#EleUnits = parameters[10].valueAsText
#outpoint_final = aparameters[11].valueAsText

# Start Fill in Data
PointFile = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
trans_count = 8 
transsample_count = 4 # does not include a sample at the stream node
StreamSample = True # include a sample at the stream node (emergent sample)? (True/False)
transsample_distance = 8
LCRaster = r"D:\Projects\TTools_9\Example_data.gdb\veght_lidar_ft" # This is either landcover height or codes
CanopyDataType = "Codes"
heatsource8 = False
CanopyRaster = "" # OPTIONAL This is either canopy cover or a LAI raster
kRaster = "" # OPTIONAL This is the k value for LAI
EleRaster = r"D:\Projects\TTools_9\Example_data.gdb\be_lidar_ft"
EleUnits = "Feet"
outpoint_final = r"D:\Projects\TTools_9\Example_data.gdb\LC_samplepoint"
OverwriteData = False
# End Fill in Data

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

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

    with arcpy.da.SearchCursor(pointfile, Incursorfields,"",proj) as Inrows:
        if OverwriteData == True:
            for row in Inrows:
                pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
                pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[5] == None or row[5] == 0:
                    pnt_dict[row[0]][row[1]]["STREAM_KM"] = row[2] 
                    pnt_dict[row[0]][row[1]]["POINT_X"] = row[3]
                    pnt_dict[row[0]][row[1]]["POINT_Y"] = row[4]
    return pnt_dict

def CreateLCPointFile(pointList, LCFields, pointfile, proj):
    """Creates the output landcover sample point feature class using the data from the nodes list"""
    arcpy.AddMessage("Exporting Data")
    print("Exporting Data")

    arcpy.CreateFeatureclass_management(os.path.dirname(pointfile),os.path.basename(pointfile), "POINT","","DISABLED","DISABLED",proj)

    #AddFields = ["NODE_ID","AZIMUTH","ZONE"] + LCFields + ["POINT_X","POINT_Y"]
    AddFields = ["STREAM_ID","NODE_ID","AZIMUTH","ZONE"] + LCFields + ["POINT_X","POINT_Y"]

    # Add attribute fields # TODO add dictionary of field types so they aren't all double
    for f in AddFields:
        arcpy.AddField_management(pointfile, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(pointfile, AddFields + ["SHAPE@X","SHAPE@Y"]) as cursor:
        for row in pointList:
            cursor.insertRow(row)

def SetupLCDataHeaders(transsample_count, trans_count, CanopyDataType, StreamSample, heatsource8):
    """Generates a list of the landcover data file column header names and data types"""
    
    if CanopyDataType == "Codes":        
        type = ["LC","ELE"]
        
    if CanopyDataType == "LAI":  #Use LAI methods
        type = ["HT","LAI", "k", "OH", "ELE"]
        
    if CanopyDataType == "CanopyCover":  #Use Canopy Cover methods
        type = ["HT","CAN", "OH", "ELE"]
    
    lcdataheaders =[] 
    if heatsource8 == True: # a flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
        dir = ['NE','E','SE','S','SW','W','NW']
    else:        
        dir = ['D' + str(x) for x in range(1, trans_count + 1)]

    zone = range(1,int(transsample_count)+1)
    
    # Concatenate the type, dir, and zone and order in the correct way
    for t in type:
        for d in range(0,len(dir)):
            for z in range(0,len(zone)):
                if StreamSample == "True" and t !="ELE" and d==0 and z==0:
                    lcdataheaders.append(t+"_EMERGENT") # add emergent
                    lcdataheaders.append(t+"_"+dir[d]+"_"+str(zone[z]))
                else:
                    lcdataheaders.append(t+"_"+dir[d]+"_"+str(zone[z]))
    
    return lcdataheaders, type

def SampleLandcover(NodeDict, dir, zone):
    """This method through the process of sampling the landcover rasters"""
    #arcpy.AddMessage("Extracting values")
    print("Extracting raster values")
    n = 0

    for streamID in NodeDict:
        for nodeID in NodeDict[streamID]:
            print("Processing stream %s of %s, Node %s of %s" % (n+1, len(NodeDict), nodeID+1, len(NodeDict[streamID])))
            # Set the progressor
            arcpy.SetProgressor("step","",0, len(NodeDict), 1)
            arcpy.SetProgressorLabel("Process Node %s of %s" % (n+1, len(NodeDict)))	
            #arcpy.AddMessage("Process Node %s of %s" % (n+1, len(NodeDict)))

            #Add input metadata
            #NodeDict[streamID][nodeID]["NUM_DIR"] = trans_count
            #NodeDict[streamID][nodeID]["NUM_ZONES"] = transsample_count
            #NodeDict[streamID][nodeID]["SAMPLE_DIS"] = transsample_distance
            #NodeDict[streamID][nodeID]["LC_RASTER"] = LCRaster
            #NodeDict[streamID][nodeID]["ELE_RASTER"] = EleRaster
            #NodeDict[streamID][nodeID]["CAN_RASTER"] = CanopyRaster
            #NodeDict[streamID][nodeID]["k_RASTER"] = kRaster

            # xy string requirement for arcpy.GetCellValue
            origin_x = NodeDict[streamID][nodeID]["POINT_X"]
            origin_y = NodeDict[streamID][nodeID]["POINT_Y"]
            xypoint = str(origin_x) + " " + str(origin_y)            

            # Sample the stream/emergent
            if StreamSample == "True":
                for t in type:
                    LCkey = t+"_EMERGENT"
                    # Sample the point value from the appropriate raster and add to Nodes dictionary
                    if t == "ELE":
                        thevalue = arcpy.GetCellValue_management(EleRaster, xypoint)
                        NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0)) * ele_con

                    if t in ["LC", "HT"]:
                        thevalue = arcpy.GetCellValue_management(LCRaster, xypoint)
                        NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0))

                    if t in ["LAI","CCV"]:
                        thevalue = arcpy.GetCellValue_management(CanopyRaster, xypoint)
                        NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0))

                    if t == "k":
                        thevalue =arcpy.GetCellValue_management(kRaster, xypoint)
                        NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0))	    

            for d in range(0,len(dir)):
                for z in zone:
                    # Calculate the x and y coordinate of the landcover sample location
                    _X_ = (z * transsample_distance * units_con * sin(radians(dir[d]))) + origin_x
                    _Y_ = (z * transsample_distance * units_con * cos(radians(dir[d]))) + origin_y
                    xypoint = str(_X_) + " " + str(_Y_) # xy string requirement for arcpy.GetCellValue

                    # Add the coordinates into the NodeDict dictionary
                    NodeDict[streamID][nodeID]["POINT_X_D"+str(d+1)+'_'+str(z)] = _X_
                    NodeDict[streamID][nodeID]["POINT_Y_D"+str(d+1)+'_'+str(z)] = _Y_			

                    for t in type:
                        LCkey = t+'_D'+str(d+1)+'_'+str(z)
                        # Sample the point value from the appropriate raster and add to Nodes dictionary
                        if t == "ELE":
                            thevalue = arcpy.GetCellValue_management(EleRaster, xypoint)
                            NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0)) * ele_con
                        if t in ["LC"]:
                            thevalue = arcpy.GetCellValue_management(LCRaster, xypoint)
                            NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0))
                            #NodeDict[n][LCkey] = float(thevalue.getOutput(0))
                        if t in ["LAI","CCV"]:
                            thevalue = arcpy.GetCellValue_management(CanopyRaster, xypoint)
                            NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0))
                        if t == "k":
                            thevalue =arcpy.GetCellValue_management(kRaster, xypoint)
                            NodeDict[streamID][nodeID][LCkey] = float(thevalue.getOutput(0))
        n = n + 1
        arcpy.SetProgressorPosition()
    return NodeDict

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
    
    # Determine input spatial units and set conversion factor to get from meters to the input spatial units
    proj = arcpy.Describe(PointFile).SpatialReference
    units_con = FromMetersUnitConversion(PointFile)

    # Set the converstion factor to get from the input elevation z units to meters
    if EleUnits == "Meters": #International meter
        ele_con = 1 
    if EleUnits == "Feet": #International foot
        ele_con = 0.3048
    if EleUnits == "Other": #Some other units
        sys.exit("Please modify your raster elevation z units so they are either in meters or feet")	

    if heatsource8 == True: # flag indicating the model should use the heat source 8 methods (same as 8 directions but no north)
        dir = [45,90,135,180,225,270,315]
    else:        
        dir = [x * 360.0 / trans_count for x in range(1,trans_count+ 1)]

    zone = range(1,int(transsample_count+1))
    
    # TODO This is a future function that may replace the emergent methods.
    # If True there is a regular landcover sample at the stream node
    # for each azimuth direction vs a single emergent sample at the stream node.
    #if StreamSample == "TRUE":
        #zone = range(0,int(transsample_count))
    #else:
        #zone = range(1,int(transsample_count+1))    
    
    AddFields, type = SetupLCDataHeaders(transsample_count, trans_count, CanopyDataType, StreamSample, heatsource8)
    NodeDict = ReadPointFile(PointFile, OverwriteData, AddFields)
    NodeDict = SampleLandcover(NodeDict, dir, zone)

    arcpy.ResetProgressor()
    #del(n,a,z,t,thevalue,origin_x,origin_y,xypoint,_X_,_Y_)			
    gc.collect()

    ####################################################################################################### 
    # Build the output point feature class using the data from the Nodes dictionary

    # these are the fields for the output list. the extra XY at the end is for arcpy's "SHAPE@X","SHAPE@Y"
    type2 = type + ["POINT_X","POINT_Y","POINT_X","POINT_Y"]

    # This builds a unique long form list of all node, azimuth directions, and zone values (ndz). 
    # This will be used in the next step to build node dictionary keys in order to iterate through each sample 
    # and get the data. Azimuth direction 0 and Zone 0 are the emergent samples. Those are stored in a seperate list (ndz_e).
    #ndz= [[n,d,z] for n in NodeDict for d in dir for z in zone ]
    ndz= [[s,n,d,z] for s in NodeDict for n in NodeDict[s] for d in dir for z in zone ]

    # Build the node dictonary keys, get the values, and save it as a list. 
    # This list is used to create the output point feature class
    #LCpointlist = [ndz[row] + [NodeDict[ndz[row][0]][t+'_'+str(int(ndz[row][1]))+'_Z'+str(ndz[row][2])] for t in type2] for row in range(0,len(ndz))]
    LCpointlist = [ndz[row] + [NodeDict[ndz[row][0]][ndz[row][1]][t+'_'+str(int(ndz[row][2]))+'_Z'+str(ndz[row][3])] for t in type2] for row in range(0,len(ndz))]
    # add the emergent samples

    if StreamSample == True:
        type_emergent = [t+"_EMERGENT" for t in type] + ["POINT_X","POINT_Y","POINT_X","POINT_Y"]
        #ndz_e= [[n,d,z] for n in NodeDict for d in [0] for z in [0]]
        ndz_e= [[s,n,d,z] for s in NodeDict for n in NodeDict[s] for d in [0] for z in [0]]
        #LCpointlist = LCpointlist + [ndz_e[row] + [NodeDict[ndz_e[row][0]][t] for t in type_emergent] for row in range(0,len(ndz_e))]
        LCpointlist = LCpointlist + [ndz_e[row] + [NodeDict[ndz[row][0]][ndz_e[row][1]][t] for t in type_emergent] for row in range(0,len(ndz_e))]

    CreateLCPointFile(LCpointlist, type, outpoint_final, proj)
    ####################################################################################################### 

    # Create the landcover headers to be added to the TTools point feature class
    AddFields = AddFields + ["NUM_DIR","NUM_ZONES","SAMPLE_DIS"]

    # Write the landcover data to the TTools point feature class
    UpdatePointFile(NodeDict, PointFile, AddFields)

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