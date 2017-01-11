#######################################################################################
# TTools
# Step 5: Sample Landcover - Star Pattern, Zone Method v 0.97 
#
#
# WARNING: DO NOT USE THIS SCRIPT. Still under development.
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

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes_major"
trans_count = 8 
transsample_count = 5 # does not include a sample at the stream node
transsample_distance = 8
heatsource8 = False
lc_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_vght_m_mosaic"
lc_units = "Meters"
canopy_data_type = "#" # OPTIONAL This is either 1. "CanopyCover", or 2."LAI"
canopy_raster = "#" # OPTIONAL This is either canopy cover or a LAI raster
k_raster = "#" # OPTIONAL This is the k value for LAI
oh_raster = "#" # OPTIONAL This is the overhang raster
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
z_units = "Meters"
zones_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\LC_samplepoint_major"
overwrite_data = True
# End Fill in Data
# ----------------------------------------------------------------------


def nested_dict(): 
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc(nodes_fc, overwrite_data, addFields):
    """Reads the input point feature class and returns the STREAM_ID,
    NODE_ID, and X/Y coordinates as a nested dictionary"""
    nodeDict = nested_dict()
    incursorFields = ["NODE_ID", "STREAM_ID", "STREAM_KM", "SHAPE@X","SHAPE@Y"]

    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)

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
                nodeDict[row[0]]["STREAM_ID"] = row[1] 
                nodeDict[row[0]]["STREAM_KM"] = row[2] 
                nodeDict[row[0]]["POINT_X"] = row[3]
                nodeDict[row[0]]["POINT_Y"] = row[4]
        else:
            for row in Inrows:
                # Is the data null or zero, if yes grab it.
                if row[5] is None or row[5] == 0 or row[5] < -9998:
                    nodeDict[row[0]]["STREAM_ID"] = row[1] 
                    nodeDict[row[0]]["STREAM_KM"] = row[2] 
                    nodeDict[row[0]]["POINT_X"] = row[3]
                    nodeDict[row[0]]["POINT_Y"] = row[4]
    
    if len(nodeDict) == 0:
        sys.exit("The fields checked in the input point feature class " +
                 "have existing data. There is nothing to process. Exiting")
              
    return nodeDict


def sample_raster(zone_list, raster, con):
                
    zone_list_new = []
    
    if con is None:
        con = 1
    
    for row in zone_list:
        arcpy.sa.ZonalStatisticsAsTable(row[0], "SAMPLE_ID", raster, out_table, "DATA", "MEAN")
        
        with arcpy.da.SearchCursor(out_table, ["MEAN"]) as cursor:
            for row in cursor:            
                zone_list_new.append(row[0] * con)
                
    arcpy.Delete_management(out_table)
                
    return zone_list_new  



def setup_LC_data_headers(transsample_count, trans_count,
                          canopy_data_type, stream_sample, heatsource8):
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
    # a flag indicating the model should use the heat source 8 methods 
    # (same as 8 directions but no north)
    if heatsource8:
        dirs = ["T{0}".format(x) for x in range(1, 8)]
    else:        
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

    type = ["LC"] + type

    return lcheaders, otherheaders, type

def update_zones_fc(zone_list, type, zones_fc, nodes_fc,
                       nodes, overwrite_data, proj):
    
    """Creates/updates the output zones feature
    class using the data from the zone list"""

    cursorfields = ["STREAM_ID","NODE_ID", "SAMPLE_ID", 
                    "TRANS_AZIMUTH","TRANSECT", "SAMPLE", "LC_KEY"] + type    

    # Check if the output exists and create if not
    if not arcpy.Exists(zones_fc):
        arcpy.CreateFeatureclass_management(os.path.dirname(zones_fc),
                                            os.path.basename(zones_fc),
                                            "POINT","","DISABLED","DISABLED",proj)

        # Determine Stream ID field properties
        sid_type = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].type
        sid_precision = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].precision
        sid_scale = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].scale
        sid_length = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].length    

        # Add attribute fields # TODO add dictionary of field types 
        # so they aren't all double
        for f in cursorfields:
            if f == "STREAM_ID":
                arcpy.AddField_management(zones_fc, f, sid_type,
                                          sid_precision, sid_scale, sid_length,
                                          "", "NULLABLE", "NON_REQUIRED")

            elif f in ["LC_KEY","SAMPLE_ID"]:
                arcpy.AddField_management(zones_fc, f, "TEXT", "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")                

            else:
                arcpy.AddField_management(zones_fc, f, "DOUBLE", "", "", "",
                                          "", "NULLABLE", "NON_REQUIRED")

    if not overwrite_data:
        # Build a query to retreive existing rows from the nodes 
        # that need updating
        whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_in_block))        

        # delete those rows
        with arcpy.da.UpdateCursor(zones_fc,["NODE_ID"], whereclause) as cursor:  
            for row in cursor:
                cursor.deleteRow()     

    with arcpy.da.InsertCursor(zones_fc, ["SHAPE@"] +
                               cursorfields) as cursor:
        for row in zone_list:
            cursor.insertRow(row)


def make_zones_list(nodeDict, dirs, zones, trans_count, transsample_distance):
    """This builds a unique long form list of information for all the
    landcover zones. This list is used to
    create/update the output feature class."""
    
    # Get a list of the nodes, sort them
    nodes = nodeDict.keys()
    nodes.sort()    
        
    zone_list = []
    
    angleIncr = 360.0 / trans_count
    angleIncr_sub = angleIncr / 4.0

    polyArray = arcpy.Array()
    pntObj = arcpy.Point()
    
    for nodeID in nodes:
        origin_x = nodeDict[nodeID]["POINT_X"]
        origin_y = nodeDict[nodeID]["POINT_Y"]
        streamID = nodeDict[nodeID]["STREAM_ID"]
    
        for d, dir in enumerate(dirs):
            angleStart = dir - (angleIncr / 2)
            angle = angleStart
            
            for z, zone in enumerate(zones):
                # each sample zone
                
                lc_key = 'LC_T{0}_S{1}'.format(d+1, z)
                tran_num = '{:{}{}}'.format(d+1, 0, 3)
                samp_num = '{:{}{}}'.format(zone, 0, 2)
                sampleID = '{0}{1}{2}'.format(nodeID, tran_num, samp_num)
            
            
                botDis = (z + 0) * transsample_distance * con_from_m
                topDis = (z + 1) * transsample_distance * con_from_m
                
                #iterate through each vertex X/Y to complete the polygon
                if zone == 1: # First veg zone that is a triangle shape
                    for vertex in range(7):
                        if vertex in [0]: #bottom start angle
                            pntObj.X = (botDis * sin(radians(angleStart))) + origin_x
                            pntObj.Y = (botDis * cos(radians(angleStart))) + origin_y
                        if vertex in [1]: #top start angle
                            pntObj.X = (topDis * sin(radians(angleStart))) + origin_x
                            pntObj.Y = (topDis * cos(radians(angleStart))) + origin_y
                        if vertex in [2,3,4,5]: # top 2-5
                            angle = angle + angleIncr_sub
                            pntObj.X = (topDis * sin(radians(angle))) + origin_x
                            pntObj.Y = (topDis * cos(radians(angle))) + origin_y
                        if vertex in [6]: #bottom end angle
                            pntObj.X = (botDis * cos(radians(angle))) + origin_x
                            pntObj.Y = (botDis * cos(radians(angle))) + origin_y
                            angle = angle - (angleIncr_sub * 4)
                        polyArray.add(pntObj)
                if zone > 1:
                    for vertex in range(11):
                        if vertex in [0]: #bottom start angle
                            pntObj.X = (botDis * sin(radians(angleStart))) + origin_x
                            pntObj.Y = (botDis * cos(radians(angleStart))) + origin_y
                        if vertex in [1]: #top start angle
                            pntObj.X = (topDis * sin(radians(angleStart))) + origin_x
                            pntObj.Y = (topDis * cos(radians(angleStart))) + origin_y
                        if vertex in [2,3,4,5]: # top 2-5
                            angle = angle + angleIncr_sub
                            pntObj.X = (topDis * sin(radians(angle))) + origin_x
                            pntObj.Y = (topDis * cos(radians(angle))) + origin_y
                        if vertex in [6]: #bottom end angle
                            pntObj.X = (botDis * sin(radians(angle))) + origin_x
                            pntObj.Y = (botDis * cos(radians(angle))) + origin_y
                        if vertex in [7,8,9,10]: # bottom 7-10
                            angle =  angle - angleIncr_sub
                            pntObj.X = (botDis * sin(radians(angle))) + origin_x
                            pntObj.Y = (botDis * cos(radians(angle))) + origin_y
                        polyArray.add(pntObj)
                del vertex        
                
                # Add to the list 
                zone_list.append([arcpy.Polygon(polyArray),
                                  streamID, nodeID, sampleID,
                                  dir, d+1, z, lc_key])
                
                polyArray.removeAll()
         
    return zone_list    
      

try:
    print("Step 5: Sample Landcover - Star Pattern, Point Method")
    
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
       
    # Setup the raster list
    typeraster = [lc_raster, z_raster]
    if canopy_data_type == "LAI":  #Use LAI methods
        if canopy_raster is not "#": typeraster.append(canopy_raster)
        else: typeraster.append(None)
        if k_raster is not "#": typeraster.append(k_raster)
        else: typeraster.append(None)
        if oh_raster is not "#": typeraster.append(oh_raster)
        else: typeraster.append(None)
        
    if canopy_data_type == "CanopyCover":  #Use Canopy Cover methods
        if canopy_raster is not "#": typeraster.append(canopy_raster)
        else: typeraster.append(None)
        if oh_raster is not "#": typeraster.append(oh_raster)
        else: typeraster.append(None)   
    
    stream_sample = True
    
    # flag indicating the model should use the heat source 8 methods 
    # (same as 8 directions but no north)
    if heatsource8:
        dirs = [45,90,135,180,225,270,315]
    else:        
        dirs = [x * 360.0 / trans_count for x in range(1,trans_count+ 1)]

    zones = range(1,int(transsample_count+1))
    
    lcheaders, otherheaders, type = setup_LC_data_headers(transsample_count, trans_count,
                                            canopy_data_type, stream_sample,
                                            heatsource8)
    
    addFields = lcheaders + otherheaders
    
    # Get a list of existing fields
    existingFields = []
    for f in arcpy.ListFields(nodes_fc):
        existingFields.append(f.name)
        
    # Check to see if the field exists and add it if not
    for f in lcheaders:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "TEXT", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")    

    # Check to see if the field exists and add it if not
    for f in otherheaders:
        if (f in existingFields) is False:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")    
    
    # read the node data into the node dictionary
    nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)
    
    # Get a list of the nodes, sort them
    nodes = nodeDict.keys()
    nodes.sort()
    
    
    # build the zone list
    zone_list = make_zone_list(nodeDict, dirs, zones, trans_count, transsample_distance)

    for raster in typeraster:
        if raster is None:
            for i in range(0, len(zone_list)):
                zone_list[i].append(-9999)
        else: 
            if raster == z_raster:
                con = con_z_to_m
            elif raster == lc_raster:
                con = con_lc_to_m
            else:
                con = None  

            zone_list = sample_raster(nodes, zone_list, raster, con)

   
    # Update the node fc
    for row in zone_list:
        for t, item, in enumerate(type):
            lc_key = '{0}_T{1}_S{2}'.format(item, row[8], row[9])
            nodeDict[row[5]][lc_key] = row[11 + t]
    
    # Write the landcover data to the TTools point feature class 
    update_nodes_fc(nodeDict, nodes_fc, addFields, nodes)    

    # Build the output point feature class using the data         
    update_zone_fc(zone_list, type, zone_fc,
                       nodes_fc, nodes, overwrite_data, proj)

    total_samples = total_samples + len(lc_point_list)
    del zone_list
    gc.collect()

    endTime = time.time()

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