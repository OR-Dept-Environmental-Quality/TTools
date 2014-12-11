########################################################################
# TTools
# Output_To_csv - v 0.92
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

########################################################################

# Import system modules
from __future__ import division, print_function
import sys
import os
import gc
import time
import arcpy
from os.path import join
from datetime import timedelta
from math import ceil
from collections import defaultdict
from operator import itemgetter
import csv

# Parameter fields for python toolbox
#nodes_fc = parameters[0].valueAsText
#multiplecsv = parameters[1].valueAsText
#outputdir = parameters[2].valueAsText
#lcdatafile = parameters[3].valueAsText
#morphfile = parameters[4].valueAsText

# ----------------------------------------------------------------------
# Start Fill in Data
nodes_fc = r"D:\Projects\TTools_9\JohnsonCreek.gdb\jc_stream_nodes"
multiplecsv = False
outputdir = r"D:\Projects\TTools_9"
lcdatafile = "lcdata.csv"
morphfile = "morphdata.csv"
# End Fill in Data
# ----------------------------------------------------------------------

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

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()
    arcpy.AddMessage("Export to csv") 

    # Get all the column headers in the point file
    nodes_fc_headers = [field.name for field in arcpy.Describe(nodes_fc).fields]
    
    id_headers = ["STREAM_ID", "NODE_ID"]    
    lcdataheaders =["STREAM_KM","LONGITUDE","LATITUDE","TOPO_W","TOPO_S","TOPO_E"]
    morphheaders = ["STREAM_KM","ELEVATION","GRADIENT","BOTTOM_WIDTH",
                    "CHANNEL_ANGLE_Z","MANNINGS_n","SED_THERMAL_CONDUCTIVITY",
                    "SED_THERMAL_DIFFUSIVITY","SED_HYPORHEIC_THICKNESSS",
                    "%HYPORHEIC_EXCHANGE","POROSITY"]
    
    lcsample_headers = ["LC_T","HT_T","ELE_T","LAI_T","k_T","CAN_T","OH_T"]
    
    # This gets the landcover sample colnames by looking for 
    # them iteratively    
    for header in nodes_fc_headers:
        if any(txt in header for txt in lcsample_headers):
            lcdataheaders.append(header)    
    file_header_list = [lcdataheaders, morphheaders]
    file_name_list = [lcdatafile, morphfile]

    nodeDict = read_nodes_fc(nodes_fc, nodes_fc_headers)

    # Output the csv file 
    # make a wide format list by node ID from the nested dictionary
    n_nodes = 0
    if multiplecsv is True:
        for streamID in nodeDict:
            i = 0
            for ouput_header_list in file_header_list:
                outlist = []
                colnames = id_headers + ouput_header_list
                n_nodes = n_nodes + len(nodeDict[streamID])
                for nodeID in nodeDict[streamID]:
                    row_list = []
                    for header in colnames:
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
    
                # write it
                write_csv(outputdir, filename, colnames, outlist)
                i = i + 1

    else:
        i = 0
        for ouput_header_list in file_header_list:
            outlist = []
            for streamID in nodeDict:
                colnames = id_headers + ouput_header_list
                n_nodes = n_nodes + len(nodeDict[streamID])
                for nodeID in nodeDict[streamID]:
                    row_list = []
                    for header in colnames:
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
    
            # write it
            write_csv(outputdir, filename, colnames, outlist)
            i = i + 1

    endTime = time.time()
    gc.collect()
    
    elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
    mspernode = timedelta(seconds=(endTime - startTime) / n_nodes / 2).microseconds
    print("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))    

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