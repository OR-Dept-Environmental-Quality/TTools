# Import system modules
from __future__ import division, print_function
import sys
import os
import string
import gc
import shutil
import time
from datetime import timedelta
import arcpy
import itertools
from arcpy import env
from math import degrees, atan2, sqrt, pow, ceil, pi, radians
from collections import defaultdict
#######################################################################################
# TTools
# Step 2: Measure Channel Widths / Aspects - v 0.1
# Ryan Michie

# INPUTS
# 0: Input TTools point feature class (NodesFC)
# 1: input Right Bank feature class(rb)
# 2: input Left Bank feature class(lb)
# 3: input flag if existing data can be over written (OverwriteData) 1. True, 2. False

# OUTPUTS
# point feature class (edit NodesFC) - Added fields with ASPECT and CHANWIDTH for each node

# Future Updates

# This version is for manual starts from within python.
# This script requires Python 2.6, comptypes package, and ArcGIS 10.1 or higher to run.


# if using ArcObjects requires installation of the comtypes package
# https://pypi.python.org/pypi/comtypes
# 1) Download zip
# 2) open terminal and CD to directory
# 3) > python install setup.py

# Once comtypes is installed, the following modifications
# need to be made for compatibility with ArcGIS 10.1 or 10.2:
# 1) find the comtypes folder in site packages
# 1) Delete automation.pyc, automation.pyo, safearray.pyc, safearray.pyo
# 2) Edit automation.py
# 3) Add the following entry to the _ctype_to_vartype dictionary (line 794):
#    POINTER(BSTR): VT_BYREF|VT_BSTR,


#######################################################################################

#Check out the ArcGIS Spatial Analyst extension license
#arcpy.CheckOutExtension("Spatial")

env.overwriteOutput = True

# Parameter fields for python toolbox
#NodesFC = parameters[0].valueAsText
#rb = parameters[1].valueAsText
#lb = parameters[2].valueAsText
#OverwriteData = parameters[3].valueAsText

# Start Fill in Data
NodesFC = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
rb = r"D:\Projects\TTools_9\Example_data.gdb\McFee_RB"
lb = r"D:\Projects\TTools_9\Example_data.gdb\McFee_LB"
OverwriteData = True
# End Fill in Data

def Downloadcomtypes():
    import urllib
    
    print("Downloading comtypes")
    url = "https://pypi.python.org/packages/source/c/comtypes/comtypes-1.1.1.zip"
    file_name = url.split('/')[-1]
    out_filepath = "C:/" + file_name
    urllib.urlretrieve (url, r"C:\comtypes-1.1.1.zip")

def GetLibPath():
    """Return location of ArcGIS type libraries as string.
    Based on work from Matt Wilkie and others
    https://bitbucket.org/maphew/canvec/src/eaf2678de06f/Canvec/Scripts/parco.py"""
    import _winreg

    # this looks for a file which should always be present in the ESRI  com library folder
    # first in the regkey (and directory) expected for v9x, then for v10
    # this method of testing is ugly, reeks of code smell, but works
    # after 2 hours of trying to make it happen elegantly I give up, for now. (mhw)

    keyESRI = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE,"SOFTWARE\\ESRI\\ArcGIS")
    esriComDir = _winreg.QueryValueEx(keyESRI,"InstallDir") [0] + "com\\"

    if not os.path.isfile(esriComDir+'\\esriSystemUI.olb'):
        keyESRI = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE,"SOFTWARE\\ESRI\\Desktop10.1")
        esriComDir = _winreg.QueryValueEx(keyESRI,"InstallDir") [0] + "com\\"
    
    if not os.path.isfile(esriComDir+'\\esriSystemUI.olb'):
            keyESRI = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE,"SOFTWARE\\ESRI\\Desktop10.2")
            esriComDir = _winreg.QueryValueEx(keyESRI,"InstallDir") [0] + "com\\"    

    return esriComDir
    
def GetModule(sModuleName):
    """Import ArcGIS module"""
    from comtypes.client import GetModule
    sLibPath = GetLibPath()
    GetModule(sLibPath + sModuleName)

def GetStandaloneModules():
    """Import commonly used ArcGIS libraries for standalone scripts"""
    GetModule("esriSystem.olb")
    GetModule("esriGeometry.olb")
    GetModule("esriCarto.olb")
    GetModule("esriDisplay.olb")
    GetModule("esriGeoDatabase.olb")
    GetModule("esriDataSourcesGDB.olb")
    GetModule("esriDataSourcesFile.olb")
    GetModule("esriOutput.olb")

def DistanceBetweenPoints(Xa,Ya,Xb,Yb):
    """Determines the distance between two points using the pythagorean theorem and returns distance between them in map units"""
    from math import sqrt, pow
    dist = math.sqrt(math.pow((Xa - Xb),2) + math.pow((Ya - Yb),2))
    return dist

def NestedDictTree(): 
    """Build a nested dictionary"""
    return defaultdict(NestedDictTree)

def ReadNodesFC(NodesFC, OverwriteData, AddFields):
    """Reads the input point feature class and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y",]

    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(NodesFC):
        ExistingFields.append(f.name)     

    # Check to see if the 1st field exists if yes add it.
    if OverwriteData == False and (AddFields[0] in ExistingFields) == True:
        Incursorfields.append(AddFields[0])
    else:
        OverwriteData = True

    # Determine input point spatial units
    proj = arcpy.Describe(NodesFC).spatialReference

    with arcpy.da.SearchCursor(NodesFC,Incursorfields,"",proj) as Inrows:
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
    if len(pnt_dict) == 0:
        sys.exit("The fields checked in the input point feature class have existing data. There is nothing to process. Exiting")
            
    return(pnt_dict)

def ToMetersUnitConversion(inFeature):
    """Returns the conversion factor to get from the input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that is not projected or not recognized. Use a projected coordinate system preferably in linear units of feet or meters.".format(inFeature))
        sys.exit("Coordinate system is not projected or not recognized. Use a projected coordinate system, preferably in linear units of feet or meters.")   
    return con_to_m

def CalcChannelWidthAdvancedLic(nodexy, rb, lb):
    rb_near_result = arcpy.analysis.GenerateNearTable(
        nodexy, rb, r'in_memory\neartable', '', 'LOCATION','NO_ANGLE', 'CLOSEST')
    
    with arcpy.da.SearchCursor(rb_near_result, ['NEAR_DIST']) as rows:
            row = rows.next()
            rb_distance = row[0]    
    
    lb_near_result = arcpy.analysis.GenerateNearTable(
        nodexy, lb, r'in_memory\neartable', '', 'LOCATION','NO_ANGLE', 'CLOSEST')
    
    with arcpy.da.SearchCursor(lb_near_result, ['NEAR_DIST']) as rows:
                row = rows.next()
                lb_distance = row[0]    
    return(rb_distance + lb_distance)
    

def UpdateNodesFC(pointDict, NodesFC, AddFields): 
    """Updates the input point feature class with data from the nodes dictionary"""
    print("Updating input point feature class")

    # Get a list of existing fields
    ExistingFields = []
    for f in arcpy.ListFields(NodesFC):
        ExistingFields.append(f.name)     

    # Check to see if the field exists and add it if not
    for f in AddFields:
        if (f in ExistingFields) == False:
            arcpy.AddField_management(NodesFC, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")   

    with arcpy.da.UpdateCursor(NodesFC,["STREAM_ID","NODE_ID"] + AddFields) as cursor:
        for row in cursor:
            for f in xrange(0,len(AddFields)):
                streamID = row[0]
                nodeID =row[1]
                row[f+2] = pointDict[streamID][nodeID][AddFields[f]]
                cursor.updateRow(row)
                
def CalcAspect(dx, dy):
    
    RadAngle = atan2(dY,dX)
    aspect = (RadAngle * 180) / pi
    aspect2 = degrees(atan2(dY,dX))

def DistanceBetweenPoints(Xa,Ya,Xb,Yb):
    """Determines the distance between two points using the pythagorean theorem and returns distance between them in map units"""
    from math import sqrt, pow
    dist = math.sqrt(math.pow((Xa - Xb),2) + math.pow((Ya - Yb),2))
    return dist

#enable garbage collection
gc.enable()

try:
    print("Step 2: Measure Channel Width and Aspect") 
    
    #Downloadcomtypes()
    #Installcomtypes()
    #from comtypes.client import GetModule
    #GetModule(GetLibPath() + "esriSystem.olb")
    #GetModule(GetLibPath() + "esriGeometry.olb")
    #from comtypes.gen.esriGeometry import IPoint
    
    #keeping track of time
    startTime= time.time()    

    # Determine input spatial units
    proj = arcpy.Describe(NodesFC).spatialReference
    proj_rb = arcpy.Describe(rb).spatialReference
    proj_lb = arcpy.Describe(lb).spatialReference
    
    # Check to make sure the rb/lb and input points are in the same projection.
    if proj.name != proj_rb.name:
        arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(NodesFC, rb))
        sys.exit("Input points and right bank feature class do not have the same projection. Please reproject your data.")
    
    if proj.name != proj_lb.name:
            arcpy.AddError("{0} and {1} do not have the same projection. Please reproject your data.".format(NodesFC, lb))
            sys.exit("Input points and left bank feature class do not have the same projection. Please reproject your data.")     
    
    nodexy_to_m = ToMetersUnitConversion(NodesFC)
    
    AddFields = ["CHANWIDTH"]

    # Read the feature class data into a nested dictionary
    NodeDict = ReadNodesFC(NodesFC, OverwriteData, AddFields)
    
    
    n = 1
    for streamID in NodeDict:
        print("Processing stream %s of %s" % (n, len(NodeDict)))
        i =0 
        for nodeID in NodeDict[streamID]:
            
            node_x = float(NodeDict[streamID][nodeID]["POINT_X"])
            node_y = float(NodeDict[streamID][nodeID]["POINT_Y"])
            node_geom = arcpy.PointGeometry(arcpy.Point(node_x, node_y), proj)
            
            #NodeDict[streamID][nodeID]["ASPECT"] = CalcAspect(NodeDict[streamID], streamID)
            NodeDict[streamID][nodeID]["CHANWIDTH"] = CalcChannelWidth(node_geom, rb_geom, lb_geom) * nodexy_to_m
        n = n + 1         
        
    UpdateNodesFC(NodeDict, NodesFC, AddFields)

    gc.collect()

    endTime = time.time()
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / n).microseconds
    print("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))
    #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))

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