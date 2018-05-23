#######################################################################################
# TTools
# Step 5: Sample Landcover - Star Pattern, Zone Method v 0.97 
#
# WARNING: This script is still under development.
#
# This script will take a point input (from Step 1) and generate a series of pie shapped
# polygon zones extending outward from the point. The number of polygon zones and the 
# distance the polygons extend from the point is user defined.
# The script will then sample the landcover raster and run spatial statistics on each zone.
# The output is the zone polygons and a table file of the zonal summary statistics.

#
# Ryan Michie
#######################################################################################

# parameter values
# 0: input TTools point file (nodes_fc)
# 1: input Number of directions to sample (trans_count)
# 2: input Number of vegetation (transverse) samples in each direction (transsample_count)
# 3: Input vegetation height Raster (lc_raster)
# 4: input The distance between transverse samples (transsample_distance)
# 5: output polygon file name/path (zones_fc)
# 6: output table file name/path (outtable_final)

# Import system modules
from __future__ import division
from __future__ import print_function
import sys
import os
import gc
import time
import traceback
from datetime import timedelta
from math import radians, sin, cos, ceil
from collections import defaultdict, OrderedDict
import arcpy
from arcpy import env

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *

env.overwriteOutput = True

#enable garbage collection
gc.enable()

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
trans_count = 8 
transsample_count = 5 # does not include a sample at the stream node
transsample_distance = 8
heatsource8 = False
sampleID_for_code = True  # Use the sampleID instead of height for the landcover code
lc_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_vght_m_mosaic"
lc_units = "Meters"
canopy_data_type = "#" # OPTIONAL This is either 1. "CanopyCover", or 2."LAI"
canopy_raster = "#" # OPTIONAL This is either canopy cover or a LAI raster
k_raster = "#" # OPTIONAL This is the k value for LAI
oh_raster = "#" # OPTIONAL This is the overhang raster
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
z_units = "Meters"
zones_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\LC_zones_two"
overwrite_data = True
# End Fill in Data
# ----------------------------------------------------------------------

env.workspace = os.path.dirname(nodes_fc)

def nested_dict(): 
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the STREAM_ID,
    NODE_ID, and X/Y coordinates as a nested dictionary"""
    nodeDict = nested_dict()
    incursorFields = ["NODE_ID", "STREAM_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]

    # Get a list of existing fields
    existingFields = [f.name for f in arcpy.ListFields(nodes_fc)] 

    # Check to see if the last field exists if yes add it. 
    # Grabs last field becuase often the first field, emergent, is zero
    if overwrite_data is False and (addFields[len(addFields)-1] in existingFields) is True:
        incursorFields.append(addFields[len(addFields)-1])
    else:
        overwrite_data = True

    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference

    with arcpy.da.SearchCursor(nodes_fc, incursorFields,"",proj) as Inrows:
        if overwrite_data:
            for row in Inrows:
                # NodeID should always be int.
                nodeID = int(row[0])
                nodeDict[nodeID]["STREAM_ID"] = row[1] 
                nodeDict[nodeID]["STREAM_KM"] = row[2] 
                nodeDict[nodeID]["POINT_X"] = row[3]
                nodeDict[nodeID]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[5] is None or row[5] == 0 or row[5] < -9998:
                    # NodeID should always be int.
                    nodeID = int(row[0])                    
                    nodeDict[nodeID]["STREAM_ID"] = row[1] 
                    nodeDict[nodeID]["STREAM_KM"] = row[2] 
                    nodeDict[nodeID]["POINT_X"] = row[3]
                    nodeDict[nodeID]["POINT_Y"] = row[4]
    
    if len(nodeDict) == 0:
        sys.exit("The fields checked in the input point feature class " +
                 "have existing data. There is nothing to process. Exiting")
              
    return nodeDict


def sample_raster(zones_fc, node, raster, con):
                            
    out_table = "temp_zonal_stats"
    
    arcpy.MakeFeatureLayer_management(in_features=zones_fc,
                                      out_layer="temp_zone_lyr",
                                      where_clause="""NODE_ID = {0}""".format(node))    
    

    arcpy.sa.ZonalStatisticsAsTable("temp_zone_lyr", "SAMPLE_ID", raster, out_table, "DATA", "MEAN_STD")
    
    with arcpy.da.SearchCursor(out_table, ["SAMPLE_ID", "MEAN", "STD"]) as cursor:
        data = [row for row in cursor]
                
    arcpy.Delete_management(out_table)
                
    return data  

def setup_lcdata_headers(transsample_count, trans_count,
                          canopy_data_type, stream_sample):
    """Generates a list of the landcover data file
    column header names and data types"""

    type = ["ELE"]

    #Use LAI methods   
    if canopy_data_type == "LAI":
        type = type + ["LAI","k","OH"]

    #Use Canopy Cover methods    
    if canopy_data_type == "CanopyCover":  
        type = type + ["CAN","OH"]

    lcheaders = []
    otherheaders = []
    
    dirs = ["T{0}".format(x) for x in range(1, trans_count + 1)]

    zones = range(1,int(transsample_count)+1)

    # Concatenate the type, dir, and zone and order in the correct way

    for d, dir in enumerate(dirs):
        for z, zone in enumerate(zones):
            if stream_sample and d==0 and z==0:
                lcheaders.append("LC_T0_S0") # add emergent
                lcheaders.append("LC_{0}_S{1}".format(dir, zone))
            else:
                lcheaders.append("LC_{0}_S{1}".format(dir, zone))

    for t in type:
        for d, dir in enumerate(dirs):
            for z, zone in enumerate(zones):
                if stream_sample and t !="ELE" and d==0 and z==0:
                    otherheaders.append(t+"_T0_S0") # add emergent
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))
                else:
                    otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))    

    #type = ["LC"] + type

    return lcheaders, otherheaders

def update_zones_fc(zonesDict, type, zones_fc):
    """Updates the zone feature class with data from the dictionary"""
    
    # Build a query to retreive just the samples that needs updating
    whereclause = """{0} IN ({1})""".format("SAMPLE_ID", ','.join(str(i) for i in zonesDict.keys()))

    with arcpy.da.UpdateCursor(zones_fc,["SAMPLE_ID"] + type, whereclause) as cursor:  
        for row in cursor:
            for f, field in enumerate(type):
                sampleID = row[0]
                row[f+1] = zonesDict[sampleID][f]
                cursor.updateRow(row)

def from_meters_con(inFeature):
    """Returns the conversion factor to get from meters to the
    spatial units of the input feature class"""
    try:
        con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system that".format(inFeature)+
        "is not projected or not recognized. Use a projected"+
        "coordinate system preferably"+
        "in linear units of feet or meters.")
        sys.exit("Coordinate system is not projected or not recognized. "+
                 "Use a projected coordinate system, preferably in linear "+
                 "units of feet or meters.")   
    return con_from_m

def from_z_units_to_meters_con(zUnits):
    """Returns the converstion factor to get from
    the input z units to meters"""
        
    try:
        con_z_to_m = float(zUnits)
    except:
        if zUnits == "Meters":
            con_z_to_m = 1.0 
        elif zUnits == "Feet":
            con_z_to_m = 0.3048
        else: con_z_to_m = None # The conversion factor will not be used
    
    return con_z_to_m

def make_zones_fc(nodeDict, zones_fc, nodes, dirs, zones, type,
                  transsample_distance, heatsource8, proj):
    """This builds the zones feature class and returns a dictionary of
    of the samples IDs as the key and a list of the node ID and LC key
    for the values"""
    
    print("Making zones fc. May take awhile for big datasets")
           
    # Check if the output exists and create if not
    if not arcpy.Exists(zones_fc):
        arcpy.CreateFeatureclass_management(os.path.dirname(zones_fc),
                                            os.path.basename(zones_fc),
                                            "POLYGON","","DISABLED","DISABLED", proj)
        
        # Determine Stream ID field properties
        sid_type = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].type
        sid_precision = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].precision
        sid_scale = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].scale
        sid_length = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].length
        
        typeDict = {"NODE_ID": "LONG",
                    "SAMPLE_ID": "LONG",
                    "TRANS_AZIMUTH": "DOUBLE",
                    "TRANSECT": "SHORT",
                    "ZONE": "SHORT",
                    "KEY": "TEXT",}
        
        # Add MEAN and STD fields
        stat_fields = ["{0}_{1}".format(t, s) for t in type for s in ["MEAN", "STD"]]
        
        for t in stat_fields:
            typeDict[t] = "DOUBLE"

        # Add attribute fields
        # STREAM_ID
        arcpy.AddField_management(zones_fc, "STREAM_ID", sid_type,
                                  sid_precision, sid_scale, sid_length,
                                  "", "NULLABLE", "NON_REQUIRED")
        
        addFields = ["NODE_ID", "SAMPLE_ID", "TRANS_AZIMUTH",
                        "TRANSECT", "ZONE", "KEY"] + stat_fields

        for f in addFields:
            arcpy.AddField_management(zones_fc, f, typeDict[f], "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")
    else:        
        # Check to see if all the stat attribute fields are there
        field_names = [f.name for f in arcpy.ListFields(zones_fc)]
        
        # Add MEAN and STD fields
        stat_fields = ["{0}_{1}".format(t, s) for t in type for s in ["MEAN", "STD"]]
        
        for f in stat_fields:
            if f not in field_names:
                arcpy.AddField_management(zones_fc, f, typeDict[f], "", "", "",
                                              "", "NULLABLE", "NON_REQUIRED")

    if not overwrite_data:
        # Build a query to retreive existing rows
        whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes))      

        # delete those rows
        with arcpy.da.UpdateCursor(zones_fc,["NODE_ID"], whereclause) as cursor:
            for row in cursor:
                cursor.deleteRow()
                
    numDirs = len(dirs)
    numZones = len(zones)
    zonesPerNode = (numDirs * numZones) + 1
    
    if heatsource8:
        angleIncr = 45.0
    else:
        angleIncr = 360.0 / numDirs
        
    angleIncr_sub = angleIncr / 3.0

    polyArray = arcpy.Array()
    pntObj = arcpy.Point()
    
    sampleDict = {}
    
    with arcpy.da.InsertCursor(zones_fc, ["SHAPE@", "STREAM_ID"] +
                               addFields) as cursor:    
        for nodeID in nodes:
            #print("making zone fc: {0:.0f}% complete".format((nodeID+1)/len(nodes) *100))
            origin_x = nodeDict[nodeID]["POINT_X"]
            origin_y = nodeDict[nodeID]["POINT_Y"]
            streamID = nodeDict[nodeID]["STREAM_ID"]
            
            sampleID = (nodeID * zonesPerNode)
            sampleDict[sampleID] = [nodeID, "T0_S0"]
        
            for d, dir in enumerate(dirs):
                angleStart = dir - (angleIncr / 2)
                angle = angleStart
                
                for z, zone in enumerate(zones):
                    # each sample zone
                    
                    key = 'T{0}_S{1}'.format(d+1, zone)
                    tran_num = '{:{}{}}'.format(d+1, 0, 3)
                    samp_num = '{:{}{}}'.format(zone, 0, 2)
                    sampleID = (nodeID * zonesPerNode) + (d * numZones) + zone
                    
                    sampleDict[sampleID] = [nodeID, key]
                
                    botDis = (z + 0) * transsample_distance * con_from_m
                    topDis = (z + 1) * transsample_distance * con_from_m
                    
                    #iterate through each vertex X/Y to complete the polygon
                    if zone == 1: # First veg zone that is a triangle shape
                        for vertex in range(6):
                            if vertex in [0]: #bottom start angle
                                pntObj.X = (botDis * sin(radians(angleStart))) + origin_x
                                pntObj.Y = (botDis * cos(radians(angleStart))) + origin_y
                            if vertex in [1]: #top start angle
                                pntObj.X = (topDis * sin(radians(angleStart))) + origin_x
                                pntObj.Y = (topDis * cos(radians(angleStart))) + origin_y
                            if vertex in [2,3,4]: # top 2-5
                                angle = angle + angleIncr_sub
                                pntObj.X = (topDis * sin(radians(angle))) + origin_x
                                pntObj.Y = (topDis * cos(radians(angle))) + origin_y
                            if vertex in [5]: #bottom end angle
                                pntObj.X = (botDis * cos(radians(angle))) + origin_x
                                pntObj.Y = (botDis * cos(radians(angle))) + origin_y
                                angle = angle - (angleIncr_sub * 3)
                            polyArray.add(pntObj)
                    if zone > 1:
                        for vertex in range(9):
                            if vertex in [0]: #bottom start angle
                                pntObj.X = (botDis * sin(radians(angleStart))) + origin_x
                                pntObj.Y = (botDis * cos(radians(angleStart))) + origin_y
                            if vertex in [1]: #top start angle
                                pntObj.X = (topDis * sin(radians(angleStart))) + origin_x
                                pntObj.Y = (topDis * cos(radians(angleStart))) + origin_y
                            if vertex in [2,3,4]: # top 2-5
                                angle = angle + angleIncr_sub
                                pntObj.X = (topDis * sin(radians(angle))) + origin_x
                                pntObj.Y = (topDis * cos(radians(angle))) + origin_y
                            if vertex in [5]: #bottom end angle
                                pntObj.X = (botDis * sin(radians(angle))) + origin_x
                                pntObj.Y = (botDis * cos(radians(angle))) + origin_y
                            if vertex in [6,7,8]: # bottom 7-10
                                angle =  angle - angleIncr_sub
                                pntObj.X = (botDis * sin(radians(angle))) + origin_x
                                pntObj.Y = (botDis * cos(radians(angle))) + origin_y
                            polyArray.add(pntObj)
                    del vertex        
                    
                    this_zone = [arcpy.Polygon(polyArray), streamID, nodeID,
                                 sampleID, dir, d+1, zone, key]
                    
                    this_zone = this_zone + [-9999 for t in stat_fields]
                    cursor.insertRow(this_zone) 
                    
                    polyArray.removeAll()
    return sampleDict
                    
def update_nodes_fc(nodeDict, nodes_fc, addFields, node_to_query):
    """Updates the input point feature class with data from the
    nodes dictionary"""
    
    # Build a query to retreive just the nodes that needs updating
    whereclause = """{0} IN ({1})""".format("NODE_ID", node_to_query)

    with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:  
        for row in cursor:
            nodeID = int(row[0])
            for f, field in enumerate(addFields):
                row[f+1] = nodeDict[nodeID][field]
                cursor.updateRow(row)
                
try:
    print("Step 5: Sample Landcover - Star Pattern, Zone Method")
    
    #keeping track of time
    startTime= time.time()
    
    # Check if the node fc exists
    if not arcpy.Exists(nodes_fc):
        arcpy.AddError("\nThis output does not exist: \n" +
                       "{0}\n".format(nodes_fc))
        sys.exit("This output does not exist: \n" +
                 "{0}\n".format(nodes_fc))
    
    # Check if the zone fc exists and delete if needed
    if arcpy.Exists(zones_fc) and overwrite_data:
        arcpy.Delete_management(zones_fc)
    
    
    # Determine input spatial units and set unit conversion factors
    proj = arcpy.Describe(nodes_fc).SpatialReference
    proj_ele = arcpy.Describe(z_raster).spatialReference
    proj_lc = arcpy.Describe(lc_raster).spatialReference
    
    con_from_m = from_meters_con(nodes_fc)
    con_lc_to_m = from_z_units_to_meters_con(lc_units)
    con_z_to_m = from_z_units_to_meters_con(z_units)
        
    # Check to make sure the raster and input 
    # points are in the same projection.
    if proj.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not ".format(nodes_fc,z_raster)+
                       "have the same projection."+
                       "Please reproject your data.")
        sys.exit("Input points and elevation raster do not have the "+
                 "same projection. Please reproject your data.")    
    
    if proj_lc.name != proj_ele.name:
        arcpy.AddError("{0} and {1} do not ".format(proj_lc,z_raster)+
                       "have the same projection."+
                       "Please reproject your data.")
        sys.exit("The landcover and elevation rasters do not have the "+
                 "same projection. Please reproject your data.")    
       
    # Setup the raster dictionary. It is ordered because
    # key list needs to correspond to the order of the attribute fields
    rasterDict = OrderedDict({"LC": lc_raster,
                              "ELE": z_raster})
    
    if canopy_data_type == "LAI":  #Use LAI methods
        if canopy_raster is not "#":
            rasterDict["LAI"] = canopy_raster
            
        if k_raster is not "#":
            rasterDict["k"] = k_raster

        if oh_raster is not "#":
            rasterDict["OH"] = oh_raster
        
    if canopy_data_type == "CanopyCover":  #Use Canopy Cover methods
        if canopy_raster is not "#":
            rasterDict["CAN"] = canopy_raster

        if oh_raster is not "#":
            rasterDict["OH"] = oh_raster
    
    # There is no stream sample for zones
    stream_sample = False
    
    # flag indicating the model should use the heat source 8 methods 
    # same as 8 directions but no north
    if heatsource8:
        dirs = [45,90,135,180,225,270,315]
        trans_count = 7
    else:        
        dirs = [x * 360.0 / trans_count for x in range(1,trans_count+ 1)]

    zones = range(1,int(transsample_count+1))
    
    lcheaders, otherheaders = setup_lcdata_headers(transsample_count, trans_count,
                                            canopy_data_type, stream_sample)
    
    addFields = lcheaders + otherheaders
    
    # Get a list of existing fields
    # Check to see if all the type attribute fields are there
    existingFields = [f.name for f in arcpy.ListFields(nodes_fc)]    
        
    # Check to see if the field exists and add it if not
    for f in lcheaders:
        if not (f in existingFields):
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")    

    # Check to see if the field exists and add it if not
    for f in otherheaders:
        if not (f in existingFields):
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")   
    
    # read the node data into the node dictionary
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
    
    # Get a list of the nodes, sort them
    nodes = nodeDict.keys()
    nodes.sort()
    
    # build the zone list
    sampleDict = make_zones_fc(nodeDict, zones_fc, nodes, dirs, zones,
                               rasterDict.keys(), transsample_distance,
                               heatsource8, proj)
    
    stat_fields = ["{0}_{1}".format(t, s) for t in rasterDict.keys() for s in ["MEAN", "STD"] ]
    
    for n, nodeID in enumerate(nodes):
        print("Processing node {0} of {1}".format(n + 1, len(nodes)))
        
        sampleDict2 = defaultdict(list)
        
        for type, raster in rasterDict.iteritems():
            if raster == z_raster:
                con = con_z_to_m
            elif raster == lc_raster:
                con = con_lc_to_m
            else:
                con = None  

            data_list = sample_raster(zones_fc, nodeID, raster, con) 
            
            # Update the node fc
            if sampleID_for_code:
                for row in data_list:
                    key = "{0}_{1}".format(type, sampleDict[row[0]][1])
                    nodeDict[nodeID][key] = row[0]
                    sampleDict2[row[0]].append(row[1])
                    sampleDict2[row[0]].append(row[2])
            
            else:
                for row in data_list:
                    key = "{0}_{1}".format(type, sampleDict[row[0]][1])
                    nodeDict[nodeID][key] = row[1]
                    sampleDict2[row[0]].append(row[1])
                    sampleDict2[row[0]].append(row[2])
            
        # update zones fc with new sampled values       
        update_zones_fc(sampleDict2, stat_fields, zones_fc)

        # Write the landcover data to the TTools point feature class 
        update_nodes_fc(nodeDict, nodes_fc, addFields, nodeID)      
    
    
    endTime = time.time()
    
    gc.collect()

    total_samples = trans_count * transsample_count * len(nodes)
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspersample = timedelta(seconds=(endTime - startTime) /
                            total_samples).microseconds
    print("Process Complete in {0} minutes. {1} microseconds per sample".format(elapsedmin, mspersample))
    #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per sample" % (elapsedmin, mspersample))

# For arctool errors
except arcpy.ExecuteError:
    msgs = arcpy.GetMessages(2)
    #arcpy.AddError(msgs)
    print(msgs)

# For other errors
except:
    tbinfo = traceback.format_exc()

    pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" +str(sys.exc_info()[1])
    msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

    #arcpy.AddError(pymsg)
    #arcpy.AddError(msgs)

    print(pymsg)
    print(msgs)
