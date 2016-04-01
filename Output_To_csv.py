########################################################################
# TTools
# Output_To_csv - v 0.93
# Ryan Michie

# Output_To_csv will take the node point feature created from using 
# steps 1-5 and output landcover and morphology data csv files formatted
# for heat source 9.

# INPUTS
# 0: Input TTools point feature class (nodes_fc)
# 1: create seperate csv files for each STREAM_ID (multiplecsv) 
#    True/False. if True will append stream ID to csv name.
# 2: path directory where the output csv file will be saved (outputdir)
# 3: name of the landcover data csv file (lcdatafile)
# 4: name of the morphology data csv file (morphfile)

# OUTPUTS
# landcover and morphology data csv files formatted for input 
# into heatsource 9

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

# UPDATES
# These cases result in the incorrect outputs and need to fixed:

# 1. Not all the columns are in the node fc -> missing columns in 
# the csv. Instead the script should output the missing col 
# filled with None.

# 2. The columns in node_fc are out of order b/c the steps were not run 
# in sequential order or rerun with different arguments -> columns in 
# the csv are in the wrong order.

########################################################################

# Import system modules
from __future__ import division, print_function
import sys
import os
import gc
import time
import arcpy
import traceback
from os.path import join
from datetime import timedelta
from math import ceil
from collections import defaultdict
from operator import itemgetter
import csv

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
multiplecsv = False
topo_directions = 1
trans_count = 8 
transsample_count = 5 # does not include a sample at the stream node
transsample_distance = 8
canopy_data_type = "#" # OPTIONAL This is either 1. "CanopyCover", or 2."LAI"
heatsource8 = False

outputdir = r"D:\Projects\TTools_9"
lcdatafile = "lcdata.csv"
morphfile = "morphdata.csv"
# End Fill in Data
# ----------------------------------------------------------------------

# Parameter fields for python toolbox
#nodes_fc = parameters[0].valueAsText
#multiplecsv = parameters[1].valueAsText
#outputdir = parameters[2].valueAsText
#lcdatafile = parameters[3].valueAsText
#morphfile = parameters[4].valueAsText

def nested_dict(): 
    """Build a nested dictionary"""
    return defaultdict(nested_dict)

def read_nodes_fc(nodes_fc, readfields):
    """Reads an input point file and returns the NODE ID and X/Y
    coordinates as a nested dictionary"""
    nodeDict = nested_dict()
    incursorFields = ["STREAM_ID","NODE_ID"] + readfields
    # Determine input point spatial units
    proj = arcpy.Describe(nodes_fc).spatialReference
    with arcpy.da.SearchCursor(nodes_fc,incursorFields,"",proj) as Inrows:
        for row in Inrows:
            for f in xrange(0,len(readfields)):
                nodeDict[row[0]][row[1]][readfields[f]] = row[2+f]
    return(nodeDict)

def write_csv(outputdir, filename, colnames, outlist):
    """write the output list to csv"""
    
    # insert column header names
    outlist.insert(0, colnames)
    
    with open(join(outputdir, filename), "wb") as file_object:
        writer = csv.writer(file_object,  dialect= "excel")
        writer.writerows(outlist)
        
def setup_LC_data_headers(transsample_count, trans_count,
                          canopy_data_type, heatsource8, 
                          topo_directions):
    """Generates a list of the landcover data file
    column header names and data types"""
    
    type = ["LC","ELE"]
    
    lcdataheaders = ["STREAM_KM","LONGITUDE","LATITUDE"]
    
    if topo_directions == 2: # All directions
        azimuths = [45,90,135,180,225,270,315,365]
        last_azimuth = 45
    else:        
        azimuths = [270,180,90]
        last_azimuth = 90
    
        azimuthdict = {45:"TOPO_NE",90:"TOPO_E",135:"TOPO_SE",
                       180:"TOPO_S",225:"TOPO_SW",270:"TOPO_W",
                       315:"TOPO_NW",365:"TOPO_N"}
    
    # Add/build the topo field names
    lcdataheaders = lcdataheaders + [azimuthdict[a] for a in azimuths]

    #Use LAI methods   
    if canopy_data_type == "LAI":
        type = type + ["LAI","k","OH"]

    #Use Canopy Cover methods    
    if canopy_data_type == "CanopyCover":  
        type = type + ["CAN","OH"]

    # a flag indicating the model should use the heat source 8 methods 
    # (same as 8 directions but no north)
    if heatsource8:
        dirs = ["T{0}".format(x) for x in range(1, 8)]
    else:        
        dirs = ["T{0}".format(x) for x in range(1, trans_count + 1)]

    zones = range(1,int(transsample_count)+1)
    stream_sample = True

    # Concatenate the type, dir, and zone and order in the correct way
    for t in type:
        for d, dir in enumerate(dirs):
            for z, zone in enumerate(zones):
                if stream_sample is True and t !="ELE" and d==0 and z==0:
                    #lcdataheaders.append(t+"_EMERGENT") # add emergent
                    lcdataheaders.append(t+"_T0_S0") # add emergent
                    lcdataheaders.append("{0}_{1}_S{2}".format(t, dir, zone))
                else:
                    lcdataheaders.append("{0}_{1}_S{2}".format(t, dir, zone))

    return lcdataheaders

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()
    #arcpy.AddMessage("Export to csv")
    print("Export to csv")

    # Get all the column headers in the point file
    nodes_fc_headers = [field.name for field in arcpy.Describe(nodes_fc).fields]
    
    id_headers = ["STREAM_ID", "NODE_ID"]    
    morphheaders = ["STREAM_KM","ELEVATION","GRADIENT","BOTTOM_WIDTH",
                    "CHANNEL_ANGLE_Z","MANNINGS_n","SED_THERMAL_CONDUCTIVITY",
                    "SED_THERMAL_DIFFUSIVITY","SED_HYPORHEIC_THICKNESSS",
                    "%HYPORHEIC_EXCHANGE","POROSITY"]
        
    lcdataheaders = setup_LC_data_headers(transsample_count, trans_count, 
                                            canopy_data_type, heatsource8,
                                            topo_directions)
    
    
    file_header_list = [lcdataheaders, morphheaders]
    file_name_list = [lcdatafile, morphfile]    
    
    # This checks if the landcover sample colnames are in the node_fc.
    headers_to_query = []
    for header in nodes_fc_headers:
        if any(txt in header for txt in lcdataheaders):
            headers_to_query.append(header)
    

    nodeDict = read_nodes_fc(nodes_fc, headers_to_query)
    
    # Output the csv file 
    # make a wide format list by node ID from the nested dictionary
    n_nodes = 0
    if multiplecsv is True:
        streamIDs = nodeDict.keys()
        streamIDs.sort()        
        for i, streamID in enumerate(streamIDs):
            for ouput_header_list in file_header_list:
                outlist = []
                nodeIDs = nodeDict[streamID].keys()
                nodeIDs.sort()                
                n_nodes = n_nodes + len(nodeIDs)
                for nodeID in nodeIDs:
                    row_list = [streamID, nodeID]
                    for header in ouput_header_list:
                        if header in nodeDict[streamID][nodeID].keys():
                            val = nodeDict[streamID][nodeID][header]
                            row_list.append(val)
                        else:
                            # use None when there is no data in the Nodes 
                            # feature class
                            row_list.append(None)
                    outlist.append(row_list)
                       
                #sort the list by stream km
                outlist = sorted(outlist, key=itemgetter(1), reverse=True)
            
                # name the output
                filename = file_name_list[i].replace(".csv", "") + "_" + str(streamID) + ".csv"
                
                colnames = id_headers + ouput_header_list
    
                # write it
                write_csv(outputdir, filename, colnames, outlist)

    else:
                
        i = 0
        for ouput_header_list in file_header_list:
            outlist = []
            streamIDs = nodeDict.keys()
            streamIDs.sort()             
            for streamID in streamIDs:
                nodeIDs = nodeDict[streamID].keys()
                nodeIDs.sort()                 
                n_nodes = n_nodes + len(nodeIDs)
                for nodeID in nodeIDs:
                    row_list = [streamID, nodeID]
                    for header in ouput_header_list:
                        if header in nodeDict[streamID][nodeID].keys():
                            val = nodeDict[streamID][nodeID][header]
                            row_list.append(val)
                        else:
                            # use None when there is no data in the 
                            # Nodes feature class
                            row_list.append(None)
                    outlist.append(row_list)

            #sort the list by stream ID and then stream km
            outlist = sorted(outlist, key=itemgetter(1,2), reverse=True)
    
            # name the output
            filename = file_name_list[i]
            
            colnames = id_headers + ouput_header_list
    
            # write it
            write_csv(outputdir, filename, colnames, outlist)
            i = i + 1

    endTime = time.time()
    gc.collect()
    
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / n_nodes / 2).microseconds
    print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))    

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