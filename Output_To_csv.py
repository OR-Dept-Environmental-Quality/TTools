#######################################################################################
# TTools
# Output_To_csv - v 0.91
# Ryan Michie

# Output_To_csv will take the point feature created from using steps 1-5 and output a landcover data csv file formatted 
# for heat source 9.

# INPUTS
# 0: Input TTools point feature class (inPoint)
# 1: create seperate csv files for each STREAM_ID (multiplecsv) True/False. if True will appends stream ID to csv name.
# 2: path directory where output csv file will be saved (outcsv_dir)
# 3: name of the csv file (outcsv_file)

# OUTPUTS
# csv file, with the same columns and order as the input point file (minus a few feature class specfic columns)

# This version is for manual starts from within python.
# This script requires Python 2.6 and ArcGIS 10.1 or higher to run.

#######################################################################################

# Import system modules
from __future__ import division, print_function
import sys, os, string, gc, shutil, time
import arcpy
from collections import defaultdict
from operator import itemgetter
import csv

# Parameter fields for python toolbox
#inPoint = parameters[0].valueAsText
#multiplecsv = parameters[1].valueAsText
#outcsv_dir = parameters[2].valueAsText
#outcsv_file = parameters[3].valueAsText

# Start Fill in Data
inPoint = r"D:\Projects\TTools_9\Example_data.gdb\out_nodes"
multiplecsv = "True"
outcsv_dir = r"D:\Projects\TTools_9"
outcsv_file = "out_nodes.csv"
# End Fill in Data

def read_pointfile(pointfile, readfields):
    """Reads an input point file and returns the NODE ID and X/Y coordinates as a nested dictionary"""
    pnt_dict = NestedDictTree()
    Incursorfields = ["STREAM_ID","NODE_ID"] + readfields
    # Determine input point spatial units
    proj = arcpy.Describe(inPoint).spatialReference
    with arcpy.da.SearchCursor(pointfile,Incursorfields,"",proj) as Inrows:
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

def read_csv(csvfile):
    """Reads an input csv file and returns the header row as a list and the data as a nested dictionary"""
    csv_dict = NestedDictTree()
    with open(csvfile, "rb") as f:
        reader = csv.reader(f)
        header = reader.next()
        if header[0] != "NODE_ID":
            sys.exit("csv file does not have NODE_ID as first column")
        for row in reader:
            for key in xrange(0,len(header)):
                csv_dict[row[0]][header[key]] = row[key]
    return(header,csv_dict)

#enable garbage collection
gc.enable()

try:
    #keeping track of time
    startTime= time.time()
    arcpy.AddMessage("Export to csv") 

    removelist = [u"OBJECTID",u"Id",u"Shape",u"ELEVATION",u"GRADIENT",u"NUM_DIR",u"NUM_ZONES",u"SAMPLE_DIS"]

    # Get all the column headers in the point file and remove the ones in removelist
    header = [field.name for field in arcpy.Describe(inPoint).fields]
    header_clean = [h for h in header if h not in removelist]

    NODES = read_pointfile(inPoint, header_clean)

    ####################################################################################################### 
    # Output the csv file

    # Get the stream km dictionary keys and sort them
    #NODE_keys = NODES.keys()
    #NODE_keys.sort()    

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
    print("Process Complete in %s minutes" % (elapsedmin))
    arcpy.AddMessage("Process Complete at %s, %s minutes" % (endTime, elapsedmin))    

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