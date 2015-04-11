########################################################################
# TTools
# Step 1: Create Stream Nodes  version 0.94
# Ryan Michie

# This script will take an input polyline feature with unique 
# stream IDs and generate evenly spaced points along each 
# unique stream ID polyline at a user defined spacing measured from 
# the downstream endpoint. The script can also check the digitized
# direction to determine the downstream end.

# INPUTS
# 0: Stream centerline polyline (streamline_fc)
# 1: Unique StreamID field (sid_field)
# 2: Spacing between nodes in meters (node_dx)
# 3: OPTIONAL True/False flag to check if the stream was digitized in 
#    correct direction (checkDirection)
# 4: OPTIONAL Elevation Raster used in the check stream 
#    direction procedure (z_raster)
# 5: Path/Name of output node feature class (nodes_fc)

# OUTPUTS
# point feature class

# The output point feature class has the following fields: 
# 0: NODE_ID - unique node ID
# 1: STREAM_ID"- field matching a unique polyline ID field identifed 
#     by the user,
# 2: STREAM_KM - double measured from the downstream end of the stream 
#     for each STREAM ID
# 3: LONGITUDE - decimal degrees X coordinate of the node using GCS_WGS_1984 datum.
# 4: LATITUDE - decimal degrees Y coordinate of the node using GCS_WGS_1984 datum.
# 5. STREAM_AZMTH - stream azimuth in the direction of flow"

# Future Updates
# eliminate arcpy and use gdal for reading/writing feature class data

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

########################################################################

# Import system modules
from __future__ import division, print_function
import sys
import os
import gc
import time
from datetime import timedelta
from math import ceil, atan2, degrees
from operator import itemgetter
import arcpy
from arcpy import env

env.overwriteOutput = True

# Parameter fields for python toolbox
#streamline_fc = parameters[0].valueAsText
#sid_field = parameters[1].valueAsText
#node_dx = parameters[2].valueAsText
#checkDirection = parameters[3].valueAsText
#z_raster = = parameters[4].valueAsText
#nodes_fc = parameters[5].valueAsText

# ----------------------------------------------------------------------
# Start Fill in Data
streamline_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_streams"
sid_field = "NAME"
node_dx = 50
checkDirection = True
z_raster = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_be_m_mosaic"
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
# End Fill in Data
# ----------------------------------------------------------------------

def create_node_list(streamline_fc, checkDirection, z_raster):
    """Reads an input stream centerline file and returns the NODE ID,
    STREAM ID, and X/Y coordinates as a list"""
    nodeList = []
    incursorFields = ["SHAPE@","SHAPE@LENGTH", sid_field]
    nodeID = 0
    # Determine input projection and spatial units
    proj = arcpy.Describe(streamline_fc).spatialReference
    con_from_m = from_meters_con(streamline_fc)
    con_to_m = to_meters_con(streamline_fc)
    
    # Pull the stream IDs into a list
    sid_list = []
    with arcpy.da.SearchCursor(streamline_fc, sid_field,"",proj) as Inrows:
        for row in Inrows:
            sid_list.append(row[0])
    
    # Check for duplicate stream IDs
    dups = list(set([i for i in sid_list if sid_list.count(i) > 1]))
    if dups:
        sys.exit("There are duplicate stream IDs in your input stream"+
                 "feature class."+
                 "\nHere are the duplicates:  \n"+
                 "{0}".format(dups))        
    
    # Now create the nodes. I'm pulling the fc data twice because on 
    # speed tests it is faster compared to saving all the incursorFields 
    # to a list and iterating over the list        
    print("Creating Nodes")
    with arcpy.da.SearchCursor(streamline_fc, incursorFields,"",proj) as Inrows:
        for row in Inrows:    
            lineLength = row[1]  # These units are in the units of projection
            numNodes = int(lineLength * con_to_m / node_dx)
            nodes = range(0,numNodes+1)
            mid = range(0,numNodes)
            
            if checkDirection is True:
                flip = check_stream_direction(row[0], z_raster, row[2])
            else:
                flip = 1
            arcpy.SetProgressor("step", "Creating Nodes", 0, numNodes+1, 1)
            # list of percentage of feature length to traverse
            positions = [n * node_dx * con_from_m / lineLength for n in nodes]
            mid_distance = node_dx * con_from_m / lineLength
            for position in positions:
                node = row[0].positionAlongLine(abs(flip - position),
                                                True).centroid
                
                # Get the coordinates at the up/down midway point along
                # the line between nodes and calculate the stream azimuth
                if position == 0.0:
                    mid_up = row[0].positionAlongLine(
                        abs(flip - (position + mid_distance)),True).centroid
                    mid_down = node
                elif 0.0 < position + mid_distance < 1:
                    mid_up = row[0].positionAlongLine(
                        abs(flip - (position + mid_distance)),True).centroid
                    mid_down = row[0].positionAlongLine(
                        abs(flip - (position - mid_distance)),True).centroid
                else:
                    mid_up = node
                    mid_down = row[0].positionAlongLine(
                        abs(flip - (position - mid_distance)),True).centroid
            
                stream_azimuth = degrees(atan2((mid_down.X - mid_up.X),
                                               (mid_down.Y - mid_up.Y)))
                if stream_azimuth < 0:
                    stream_azimuth = stream_azimuth + 360
                
                # list of "NODE_ID",STREAM_ID,"STREAM_KM",
                # "POINT_X","POINT_Y", STREAM_AZMTH, "SHAPE@X","SHAPE@Y"
                nodeList.append((nodeID, row[2],
                                 float(position * lineLength * con_to_m /1000),
                                 node.X, node.Y, stream_azimuth, node.X, node.Y))
                nodeID = nodeID + 1
        arcpy.SetProgressorPosition()
    arcpy.ResetProgressor()
    return(nodeList)

def create_nodes_fc(nodeList, nodes_fc, sid_field, proj):
    """Create the output point feature class using
    the data from the nodes list"""
    #arcpy.AddMessage("Exporting Data")
    print("Exporting Data")

    # Determine Stream ID field properties
    sid_type = arcpy.ListFields(streamline_fc,sid_field)[0].type
    sid_precision = arcpy.ListFields(streamline_fc,sid_field)[0].precision
    sid_scale = arcpy.ListFields(streamline_fc,sid_field)[0].scale
    sid_length = arcpy.ListFields(streamline_fc,sid_field)[0].length    

    #Create an empty output with the same projection as the input polyline
    cursorfields = ["NODE_ID",
                    "STREAM_ID",
                    "STREAM_KM",
                    "LONGITUDE",
                    "LATITUDE",
                    "STREAM_AZMTH"]
    arcpy.CreateFeatureclass_management(os.path.dirname(nodes_fc),
                                        os.path.basename(nodes_fc),
                                        "POINT","","DISABLED","DISABLED",proj)

    # Add attribute fields
    for f in cursorfields:
        if f == "STREAM_ID":
            arcpy.AddField_management(nodes_fc, f, sid_type, sid_precision,
                                      sid_scale, sid_length, "",
                                      "NULLABLE", "NON_REQUIRED")
        else:
            arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                      "", "NULLABLE", "NON_REQUIRED")

    with arcpy.da.InsertCursor(nodes_fc, cursorfields + ["SHAPE@X","SHAPE@Y"]) as cursor:
        for row in nodeList:
            cursor.insertRow(row)

    #Change X/Y from input spatial units to decimal degrees
    proj_dd = arcpy.SpatialReference(4326) # GCS_WGS_1984
    with arcpy.da.UpdateCursor(nodes_fc,["SHAPE@X","SHAPE@Y","LONGITUDE",
                                         "LATITUDE"],"",proj_dd) as cursor:
        for row in cursor:
            row[2] = row[0] # LONGITUDE
            row[3] = row[1] # LATITUDE
            cursor.updateRow(row)

def check_stream_direction(stream, z_raster, streamID):
    """Samples the elevation raster at both ends of the stream
    polyline to see which is the downstream end and returns flip = 1
    if the stream km need to be reversed"""
    
    down = stream.positionAlongLine(0,True).centroid
    up = stream.positionAlongLine(1,True).centroid
    
    # when a single raster cell is sampled it is a little faster to 
    # use arcpy compared to converting to an array and then sampling. 
    # I left the code just in case though
    z_down = float(arcpy.GetCellValue_management (z_raster, str(down.X)+ " "+ str(down.Y),1).getOutput(0))
    z_up = float(arcpy.GetCellValue_management (z_raster, str(up.X) + " "+ str(up.Y),1).getOutput(0))
    #z_down = arcpy.RasterToNumPyArray(z_raster, arcpy.Point(down.X, down.Y), 1, 1, -9999)[0][0]
    #z_up = arcpy.RasterToNumPyArray(z_raster, arcpy.Point(up.X, up.Y), 1, 1, -9999)[0][0]
    
    if z_down <= z_up or z_down == -9999 or z_up == -9999:
        # do not reverse stream km
        flip = 0
    else:
        print("Reversing " + streamID)
        # reversed stream km
        flip = 1
        
    return flip
    
def to_meters_con(inFeature):
    """Returns the conversion factor to get from the
    input spatial units to meters"""
    try:
        con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system ".format(inFeature)+
                       "that is not projected or not recognized. "+
                       "Use a projected coordinate system "
                       "preferably in linear units of feet or meters.")
        sys.exit("Coordinate system is not projected or not recognized. "+
                 "Use a projected coordinate system, preferably in linear "+
                 "units of feet or meters.")
    return con_to_m
    
def from_meters_con(inFeature):
    """Returns the conversion factor to get from meters to the
    spatial units of the input feature class"""
    try:
        con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
    except:
        arcpy.AddError("{0} has a coordinate system ".format(inFeature)+
                       "that is not projected or not recognized. "+
                       "Use a projected coordinate system "
                       "preferably in linear units of feet or meters.")
        sys.exit("Coordinate system is not projected or not recognized. "+
                 "Use a projected coordinate system, preferably in linear "+
                 "units of feet or meters.")
    return con_from_m

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()
    
    # Check if the output exists
    #if arcpy.Exists(nodes_fc):
        #arcpy.AddError("\nThis output already exists: \n" +
                       #"{0}\n".format(nodes_fc) + 
                       #"Please rename your output.")
        #sys.exit("This output already exists: \n" +
                 #"{0}\n".format(nodes_fc) + 
                 #"Please rename your output.")    
    
    # Get the spatial projecton of the input stream lines
    proj = arcpy.Describe(streamline_fc).SpatialReference    
    
    if checkDirection is True:
        proj_ele = arcpy.Describe(z_raster).spatialReference
    
        # Check to make sure the  elevatiohn raster and input 
        # streams are in the same projection.
        if proj.name != proj_ele.name:
            arcpy.AddError("{0} and {1} do not ".format(nodes_fc,z_raster)+
                           "have the same projection."+
                           "Please reproject your data.")
            sys.exit("Input stream line and elevation raster do not have "
                     "the same projection. Please reproject your data.")
    
    # Create the stream nodes and return them as a list
    nodeList = create_node_list(streamline_fc, checkDirection, z_raster)

    #sort the list by stream ID and then stream km
    nodeList = sorted(nodeList, key=itemgetter(1,2), reverse=True)
    
# Create the output node feature class with the nodes list
    create_nodes_fc(nodeList, nodes_fc, sid_field, proj)

    gc.collect()

    endTime = time.time()
    elapsedmin = ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / len(nodeList)).microseconds
    print("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))
    #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))	

# For arctool errors
except arcpy.ExecuteError:
    msgs = arcpy.GetMessages(2)
    #arcpy.AddError(msgs)
    print(msgs)

# For other errors
except:
    import traceback
    tb = sys.exc_info()[2]
    tbinfo = traceback.format_tb(tb)[0]

    pymsg = "PYTHON ERRORS:\nTraceback info:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
    msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

    #arcpy.AddError(pymsg)
    #arcpy.AddError(msgs)

    print(pymsg)
    print(msgs)