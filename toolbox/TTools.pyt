#!/usr/bin/env python3
"""
TTools
Version: [v9.2]

Primary Author: Ryan Michie
Organization: [Oregon Department of Environmental Quality]

Adapted from older versions of TTools authored by: Brian Kasper & Matt Boyd
Boyd, M., and Kasper, B. 2003. Analytical methods for dynamic open channel heat and mass transfer: Methodology for heat source model Version 7.0.

Most Python 2 to Python 3 updates by:
Tommy Franzen
Organization: [The Freshwater Trust]

Arc Toolbox TTools.pyt developed by Tetra Tech, with funding from the U.S. Environmental Protection Agency, June 2025

"""

import os
import datetime

def log_run(tool_name, params):
    """
    Every time a particular tool is executed log_run will:
    Append a line to a .log file.
    The log filename includes the current date.
    The log line includes the timestamp and parameter values.
    """

    # Note, any None, empty string, or "#" is turned into "<null>". Blank or no entries cause A crash

    # Directory to store logs (in the same folder as this script)
    log_dir = os.path.join(os.path.dirname(__file__), "logs")
    os.makedirs(log_dir, exist_ok=True)

    # Log filename with date
    today_str = datetime.datetime.now().strftime("%Y-%m-%d")
    log_file = os.path.join(log_dir, f"ttools_log_{today_str}.log")

    # Timestamp
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Extract parameter values safely
    param_values = []
    for p in params:
        try:
            if hasattr(p, "valueAsText"):
                val = p.valueAsText
            else:
                val = str(p)

            # Replace None, blank, or "#" with a placeholder so it accepts as string
            if val in (None, "", "#"):
                val = "<null>"

            param_values.append(str(val))
        except Exception as e:
            param_values.append(f"<error: {e}>")

    log_entry = f"{timestamp} | TTOOLS: {tool_name} | " + " | ".join(param_values) + "\n"

    # Append to log file
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(log_entry)
#-----------
import arcpy
import sys
import os
import gc
import time
import traceback
from datetime import timedelta
from math import ceil, atan2, degrees
from operator import itemgetter
from collections import Counter

class Toolbox:
    def __init__(self):
        self.label = "StreamTools"
        self.alias = "streamtools"
        self.tools = [CreateStreamNodes, MeasureChannelWidth, SampleElevationGradient, MeasureTopographicAngles, SampleLandcoverStartPattern, SampleLandcoverOrthogonalMethod]  # Other tools will be added here later


class CreateStreamNodes:
    # This constructor method is called when ArcGIS creates an instance of the tool.
    # It defines the tool's label and description, which appear in the Geoprocessing pane.
    def __init__(self):
        self.label = "Step 1: Create Stream Nodes"
        self.description = "Step 1: Generate evenly spaced nodes along stream centerlines."

    # This method defines the input and output parameters for the tool.
    # ArcGIS Pro uses these to build the UI in the Geoprocessing pane.
    #For a complete list of available datatypes e.g. DEFeatureClass, GPRasterlayer etc see
    #  <https://pro.arcgis.com/en/pro-app/latest/arcpy/geoprocessing_and_python/defining-parameter-data-types-in-a-python-toolbox.htm>

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            displayName="Input Stream Feature Class",
            name="streamline_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Stream ID Field",
            name="sid_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        param1.parameterDependencies = [param0.name]

        param2 = arcpy.Parameter(
            displayName="Node Spacing (meters)",
            name="node_dx",
            datatype="Double",
            parameterType="Required",
            direction="Input")

        param3 = arcpy.Parameter(
            displayName="Continuous Stream KM",
            name="cont_stream_km",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        param3.value = False

        param4 = arcpy.Parameter(
            displayName="Check Stream Direction?",
            name="checkDirection",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        param4.value = False

        param5 = arcpy.Parameter(
            displayName="Elevation Raster (Required if Check Direction is True)",
            name="z_raster",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
        
        param6 = arcpy.Parameter(
            displayName="Output Nodes Feature Class",
            name="nodes_fc",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        return [param0, param1, param2, param3, param4, param5, param6]

    def isLicensed(self):
        # This method tells ArcGIS whether this tool requires a special license.
        # Returning True means the tool is always available to run (no license restrictions).
        return True

    # This method can dynamically update parameter values or visibility
    # before the tool is run. Not used here, so it simply returns.
    def updateParameters(self, parameters):
        return

    # This method customizes validation messages shown in the tool dialog.
    # Useful for warning users about bad combinations or missing inputs.
    def updateMessages(self, parameters):
        return

    # This is the method that ArcGIS calls when the user clicks 'Run'.
    # It extracts parameter values and passes them to the tool's main logic.

    def execute(self, parameters, messages):
    # Log input parameters to file
        log_run(self.label, parameters)
        """
        TTools Step 1: Create Stream Nodes

        This script will take an input polyline feature with unique stream IDs and generate evenly spaced points along each
        unique stream ID polyline at a user defined spacing measured from the downstream endpoint. The script can also check
        the digitized direction to determine the downstream end.

        REQUIREMENTS
        ESRI ArcPro
        Python 3.7+

        INPUT VARIABLES
        0: streamline_fc:
        Path to the stream centerline polyline feature class.

        1: sid_field:
        Name of the attribute field in streamline_fc holding the unique stream identifier such as the stream name or ID number.

        2: node_dx
        Spacing between nodes in meters.

        3: cont_stream_km
        Boolean (True/False) flag to indicate that a continuous stream km should be used for all nodes regardless of the
        unique the values in the sid_field

        4: checkDirection
        Boolean (True/False) flag to check if the stream was digitized in correct direction. If checkDirection = True
        z_raster must be set.

        5: z_raster:
        Path and name of the ground elevation raster. Ignored if checkDirection = False.

        6: nodes_fc
        Path and name of the output to node feature class.

        OUTPUTS
        0. nodes_fc:
        New point feature class with the following fields:
        NODE_ID - Unique node ID.
        STREAM_ID - Field matching a unique stream identifier from the sid_field.
        STREAM_KM - Double measured from the downstream end of the stream for each STREAM ID
        LONGITUDE - Decimal degrees X coordinate of the node using GCS_WGS_1984 datum.
        LATITUDE - Decimal degrees Y coordinate of the node using GCS_WGS_1984 datum.
        ASPECT - Stream aspect in the direction of flow.

        """

        # Import system modules
        import sys
        import os
        import gc
        import time
        import traceback
        from datetime import timedelta
        from math import ceil, atan2, degrees
        from operator import itemgetter
        import arcpy
        from arcpy import env

        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------

        # Future Updates
        # eliminate arcpy and use gdal for reading/writing feature class data

        env.overwriteOutput = True

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
            arcpy.AddMessage("Creating Nodes")
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
                    segment_length = [node_dx] * numNodes + [lineLength * con_to_m % node_dx]
                    mid_distance = node_dx * con_from_m / lineLength
                    if mid_distance > 1:
                        # this situation occurs when the stream < node_dx.
                        # The azimuth is calculated for the entire stream line.
                        mid_distance = 1

                    i = 0
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

                        # list of "NODE_ID","STREAM_ID". "STREAM_KM", "LENGTH",
                        # "POINT_X","POINT_Y", "ASPECT", "SHAPE@X", "SHAPE@Y"
                        nodeList.append([nodeID, row[2],
                                         float(position * lineLength * con_to_m /1000),
                                         segment_length[i],
                                         node.X, node.Y, stream_azimuth, node.X, node.Y])
                        nodeID = nodeID + 1
                        i = i + 1

                arcpy.SetProgressorPosition()
            arcpy.ResetProgressor()
            return(nodeList)

        def create_nodes_fc(nodeList, nodes_fc, sid_field, proj):
            """Create the output point feature class using
            the data from the nodes list"""
            #arcpy.AddMessage("Exporting Data")
            arcpy.AddMessage("Exporting Data")
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
                            "LENGTH",
                            "LONGITUDE",
                            "LATITUDE",
                            "ASPECT"]
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
                arcpy.AddMessage("Reversing {0}".format(streamID))
                print("Reversing {0}".format(streamID))
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

# Variable assignment based on user input - MSF
        streamline_fc = parameters[0].valueAsText
        sid_field = parameters[1].valueAsText
        node_dx = float(parameters[2].value)
        cont_stream_km = bool(parameters[3].value)
        checkDirection = bool(parameters[4].value)
        z_raster = parameters[5].valueAsText if parameters[5].value else None
        nodes_fc = parameters[6].valueAsText

        try:
            #keeping track of time
            startTime= time.time()
            # Check if the output exists
            if arcpy.Exists(nodes_fc):
                arcpy.AddError("\nThis output already exists: \n" +
                               "{0}\n".format(nodes_fc) +
                               "Please rename your output.")
                sys.exit("This output already exists: \n" +
                         "{0}\n".format(nodes_fc) +
                         "Please rename your output.")

            # Get the spatial projection of the input stream lines
            proj = arcpy.Describe(streamline_fc).SpatialReference

            if checkDirection is True:
                proj_ele = arcpy.Describe(z_raster).spatialReference

                # Check to make sure the  elevation raster and input
                # streams are in the same projection.
                if proj.name != proj_ele.name:
                    arcpy.AddError("{0} and {1} do not ".format(nodes_fc,z_raster)+
                                   "have the same projection."+
                                   "Please reproject your data.")
                    sys.exit("Input stream line and elevation raster do not have "
                             "the same projection. Please reproject your data.")

            # Create the stream nodes and return them as a list
            nodeList = create_node_list(streamline_fc, checkDirection, z_raster)

            if cont_stream_km:
                #sort the list by stream ID and stream km
                nodeList = sorted(nodeList, key=itemgetter(1, 2))

                skm = 0.0
                for i in range(0, len(nodeList)):
                    nodeList[i][2] = skm
                    skm = skm + (node_dx * 0.001)

                # re sort the list by stream km (with downstream end at the top)
                nodeList = sorted(nodeList, key=itemgetter(2), reverse=True)

            else:
                #sort the list by stream ID and then stream km (downstream end at the top)
                nodeList = sorted(nodeList, key=itemgetter(1,2), reverse=True)

            # Create the output node feature class with the nodes list
            create_nodes_fc(nodeList, nodes_fc, sid_field, proj)

            gc.collect()

            endTime = time.time()
            elapsedmin = ceil(((endTime - startTime) / 60)* 10)/10
            mspernode = timedelta(seconds=(endTime - startTime) / len(nodeList)).microseconds
            
            arcpy.AddMessage("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))
            print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))
            #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))

        # For arctool errors
        except arcpy.ExecuteError:
            msgs = arcpy.GetMessages(2)
            arcpy.AddError(msgs)
            print(msgs)

        # For other errors
        except:
            tbinfo = traceback.format_exc()

            pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" +str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)

            print(pymsg)
            print(msgs)

class MeasureChannelWidth:
    def __init__(self):
        self.label = "Step 2: Measure Channel Widths"
        self.description = "Step 2: Measure the distance from each stream node to left and right banks, and compute total channel width."

    def getParameterInfo(self):
        param0 = arcpy.Parameter(
            displayName="Stream Nodes Feature Class",
            name="nodes_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Right Bank Feature Class",
            name="rb_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        param2 = arcpy.Parameter(
            displayName="Left Bank Feature Class",
            name="lb_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
            
        param3 = arcpy.Parameter(
            displayName="Overwrite Existing Width Fields",
            name="overwrite_data",
            datatype="Boolean",
            parameterType="Required",
            direction="Input")
        param3.value = True    


        return [param0, param1, param2, param3]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
    # Log input parameters to file
        log_run(self.label, parameters)
        #!/usr/bin/python

        """
        TTools Step 2: Measure Channel Widths

        This script will measure the channel width and distance to the right and left banks 90 degrees perpendicular
        to the stream aspect at each stream node. If there isn't a channel bank polyline 90 degrees perpendicular to the stream
        aspect, the script will return the distance to the nearest edge. Output distances are in meters.

        REQUIREMENTS
        TTools steps 1 must be run before Step 2.
        ESRI ArcPro
        Python 3.7+

        INPUT VARIABLES
        0: nodes_fc:
        Path to the TTools point feature class.

        1: rb_fc
        The right bank feature class. The right bank is on the right looking downstream.

        2: lb_fc
        The left bank feature class. The left bank is on the left looking downstream.

        3: overwrite_data:
        True/False flag if existing data in nodes_fc can be overwritten.

        OUTPUTS
        0. nodes_fc:
        New fields listed below are added into nodes_fc.
        CHANWIDTH: distance in meters between left and right banks.
        LEFT: distance in meters from the stream node to the closest edge of the left bank feature 90 degrees
        perpendicular to teh stream aspect.
        RIGHT: distance in meters from the stream ode to the closest edge of the right bank feature 90 degrees
        perpendicular to teh stream aspect.
        """
        # Import system modules
        import sys
        import gc
        import time
        import traceback
        from datetime import timedelta
        import arcpy
        from arcpy import env
        from math import ceil
        from collections import defaultdict

        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------

        # Future Updates
        # eliminate arcpy and use gdal for reading/writing feature class data

        def nested_dict():
            """Build a nested dictionary"""
            return defaultdict(nested_dict)

        def read_nodes_fc(nodes_fc, overwrite_data, addFields):
            """Reads the input point feature class and returns the STREAM_ID, NODE_ID, and X/Y coordinates as a nested dictionary"""
            nodeDict = nested_dict()
            incursorFields = ["STREAM_ID", "NODE_ID", "STREAM_KM", "ASPECT", "SHAPE@X", "SHAPE@Y"]

            # Get a list of existing fields
            existingFields = []
            for f in arcpy.ListFields(nodes_fc):
                existingFields.append(f.name)

            # Check to see if the 1st field exists if yes add it.
            if overwrite_data is False and (addFields[0] in existingFields) is True:
                incursorFields.append(addFields[0])
            else:
                overwrite_data = True

            # Determine input point spatial units
            proj_nodes = arcpy.Describe(nodes_fc).spatialReference

            with arcpy.da.SearchCursor(nodes_fc, incursorFields, "", proj_nodes) as Inrows:
                if overwrite_data is True:
                    for row in Inrows:
                        nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2]
                        nodeDict[row[0]][row[1]]["ASPECT"] = row[3]
                        nodeDict[row[0]][row[1]]["POINT_X"] = row[4]
                        nodeDict[row[0]][row[1]]["POINT_Y"] = row[5]
                else:
                    for row in Inrows:
                        # if the data is null or zero (0 = default for shapefile),
                        # it is retrieved and will be overwritten.
                        if row[6] is None or row[6] == 0 or row[6] < -9998:
                            nodeDict[row[0]][row[1]]["STREAM_KM"] = row[2]
                            nodeDict[row[0]][row[1]]["ASPECT"] = row[3]
                            nodeDict[row[0]][row[1]]["POINT_X"] = row[4]
                            nodeDict[row[0]][row[1]]["POINT_Y"] = row[5]
            if len(nodeDict) == 0:
                sys.exit("The fields checked in the input point feature class " +
                         "have existing data. There is nothing to process. Exiting")

            return (nodeDict)

        def to_meters_con(inFeature):
            """Returns the conversion factor to get from the
            input spatial units to meters"""
            try:
                con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
            except:
                arcpy.AddError("{0} has a coordinate system that ".format(inFeature) +
                               "is not projected or not recognized. Use a " +
                               "projected coordinate system preferably in linear " +
                               "units of feet or meters.")
                sys.exit("Coordinate system is not projected or not recognized. " +
                         "Use a projected coordinate system, preferably in " +
                         "linear units of feet or meters.")
            return con_to_m

        def read_polyline_geometry(polyline_fc, proj_polyline):
            """Reads an input polyline into an arcpy polyline geometry object"""
            poly_list = []
            # Get the x and y of each vertex in the polyline and save
            # it as a list.
            for row in arcpy.da.SearchCursor(polyline_fc, ["SHAPE@"]):
                for part in row[0]:
                    for pnt in part:
                        poly_list.append(arcpy.Point(pnt.X, pnt.Y))
            poly_array = arcpy.Array(poly_list)
            # put it into a geometry object.
            poly_geom = arcpy.Polyline(poly_array, proj_polyline)
            del row
            return (poly_geom)

        def calc_channel_width(node_geom, bank_geom, aspect, line_dis, proj_nodes):
            pt1 = node_geom.pointFromAngleAndDistance(aspect, line_dis, "PLANAR")
            line = arcpy.Polyline(arcpy.Array([node_geom.centroid, pt1.centroid]), proj_nodes)
            pt2 = line.intersect(bank_geom, 1)
    
            if pt2.centroid:
                to_bank_distance = node_geom.distanceTo(pt2)
            else:
                # No intersection = no bank 90 deg from aspect
                # Find the minimum distance to bank, regardless of the angle
                near_distance = bank_geom.queryPointAndDistance(node_geom.centroid)
                to_bank_distance = near_distance[2]

            return (to_bank_distance)

        def update_nodes_fc(nodeDict, nodes_fc, addFields):
            """Updates the input point feature class with
            data from the nodes dictionary"""
            arcpy.AddMessage("Updating input point feature class")
            print("Updating input point feature class")

            # Get a list of existing fields
            existingFields = []
            for f in arcpy.ListFields(nodes_fc):
                existingFields.append(f.name)

            # Check to see if the field exists and add it if not
            for f in addFields:
                if (f in existingFields) is False:
                    arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                              "", "NULLABLE", "NON_REQUIRED")

            with arcpy.da.UpdateCursor(nodes_fc, ["STREAM_ID", "NODE_ID"] +
                                                 addFields) as cursor:
                for row in cursor:
                    for f, field in enumerate(addFields):
                        streamID = row[0]
                        nodeID = row[1]
                        row[f + 2] = nodeDict[streamID][nodeID][field]
                        cursor.updateRow(row)

        # enable garbage collection
        gc.enable()

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        rb_fc = parameters[1].valueAsText
        lb_fc = parameters[2].valueAsText
        overwrite_data = bool(parameters[3].value)

        try:
            arcpy.AddMessage("Step 2: Measure Channel Width")
            print("Step 2: Measure Channel Width")

            # keeping track of time
            startTime = time.time()

            # Check if the output exists
            if not arcpy.Exists(nodes_fc):
                arcpy.AddError("\nThis output does not exist: \n" +
                               "{0}\n".format(nodes_fc))
                sys.exit("This output does not exist: \n" +
                         "{0}\n".format(nodes_fc))

            if overwrite_data is True:
                env.overwriteOutput = True
            else:
                env.overwriteOutput = False

            # Determine input spatial units
            proj_nodes = arcpy.Describe(nodes_fc).spatialReference
            proj_rb = arcpy.Describe(rb_fc).spatialReference
            proj_lb = arcpy.Describe(lb_fc).spatialReference

            # Check to make sure the rb_fc/lb_fc and input points are
            # in the same projection.
            if proj_nodes.name != proj_rb.name:
                arcpy.AddError("{0} and {1} do not have ".format(nodes_fc, rb_fc) +
                               "the same projection. Please reproject your data.")
                sys.exit("Input points and right bank feature class do not have " +
                         "the same projection. Please reproject your data.")

            if proj_nodes.name != proj_lb.name:
                arcpy.AddError("{0} and {1} do not have ".format(nodes_fc, lb_fc) +
                               "the same projection. Please reproject your data.")
                sys.exit("Input points and left bank feature class do not have " +
                         "the same projection. Please reproject your data.")

            con_to_m = to_meters_con(nodes_fc)

            # distance in the spatial units of nodes_fc to look for the right or left bank
            # 5000 meters should be long enough
            out_dis = 5000 * 1 / con_to_m

            addFields = ["CHANWIDTH", "LEFT", "RIGHT"]

            # Read the feature class data into a nested dictionary
            nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)

            # Read each of the bank polylines into geometry objects
            rb_geom = read_polyline_geometry(rb_fc, proj_rb)
            lb_geom = read_polyline_geometry(lb_fc, proj_lb)

            for n, streamID in enumerate(nodeDict):
                arcpy.AddMessage("Processing stream {0} of {1}".format(n + 1, len(nodeDict)))
                print("Processing stream {0} of {1}".format(n + 1, len(nodeDict)))

                for nodeID in nodeDict[streamID]:
                    node_x = float(nodeDict[streamID][nodeID]["POINT_X"])
                    node_y = float(nodeDict[streamID][nodeID]["POINT_Y"])
                    aspect = float(nodeDict[streamID][nodeID]["ASPECT"])
                    node_geom = arcpy.PointGeometry(arcpy.Point(node_x, node_y), proj_nodes)

                    # calculate right and left transect directions in degrees
                    dir_lb = aspect - 90
                    dir_rb = aspect + 90

                    if dir_lb < 0:
                        dir_lb = dir_lb + 360

                    if dir_rb > 360:
                        dir_rb = dir_rb - 360

                    lb_distance = calc_channel_width(node_geom, lb_geom, dir_lb, out_dis, proj_nodes)
                    rb_distance = calc_channel_width(node_geom, rb_geom, dir_rb, out_dis, proj_nodes)

                    if lb_distance is None:
                        arcpy.AddMessage(f"Warning: There is not a left bank polyline 90 degrees perpendicular to the aspect at Node ID {nodeID}.")                      
                        print(f"Warning: There is not a left bank polyline 90 degrees perpendicular to the aspect at Node ID {nodeID}.")
                        lb_distance = 0.0

                    if rb_distance is None:
                        arcpy.AddMessage(f"Warning: There is not a right bank polyline 90 degrees perpendicular to the aspect at Node ID {nodeID}.")
                        print(f"Warning: There is not a right bank polyline 90 degrees perpendicular to the aspect at Node ID {nodeID}.")
                        rb_distance = 0.0

                    nodeDict[streamID][nodeID]["CHANWIDTH"] = (lb_distance + rb_distance) * con_to_m
                    nodeDict[streamID][nodeID]["LEFT"] = lb_distance * con_to_m
                    nodeDict[streamID][nodeID]["RIGHT"] = rb_distance * con_to_m
            
            update_nodes_fc(nodeDict, nodes_fc, addFields)
    
            gc.collect()

            endTime = time.time()
            elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
            mspernode = timedelta(seconds=(endTime - startTime) / (n + 1)).microseconds
            arcpy.AddMessage("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))
            print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))
            # arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))

        # For arctool errors
        except arcpy.ExecuteError:
            msgs = arcpy.GetMessages(2)
            arcpy.AddError(msgs)
            print(msgs)

        # For other errors
        except:
            tbinfo = traceback.format_exc()

            pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)

            print(pymsg)
            print(msgs)


class SampleElevationGradient:
    def __init__(self):
        self.label = "Step 3: Sample Elevation Gradient"
        self.description = "Step 3: Sample elevation from raster and calculate downstream gradient for each node."

    def getParameterInfo(self):
        param0 = arcpy.Parameter(
            displayName="Stream Nodes Feature Class",
            name="nodes_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
            
        param1 = arcpy.Parameter(
            displayName="Search Radius (cells)",
            name="searchCells",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param1.value = 2

        param2 = arcpy.Parameter(
            displayName="Smooth Flat or Negative Gradients",
            name="smooth_flag",
            datatype="Boolean",
            parameterType="Required",
            direction="Input")
        param2.value = True

        param3 = arcpy.Parameter(
            displayName="Elevation Raster",
            name="z_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        param4 = arcpy.Parameter(
            displayName="Elevation Units",
            name="z_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param4.filter.type = "ValueList"
        param4.filter.list = ["Meters", "Feet"] 
        
        param5 = arcpy.Parameter(
            displayName="Raster Block Size (km)",
            name="block_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param5.value = 5

        param6 = arcpy.Parameter(
            displayName="Overwrite Existing Elevation/Gradient Fields",
            name="overwrite_data",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param6.value = True

        return [param0, param1, param2, param3, param4, param5, param6]
        
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return
    
    def updateMessages(self, parameters):
        return
    
    def execute(self, parameters, messages):
    
    # Log input parameters to file
        log_run(self.label, parameters)
        """
        TTools Step 3: Sample Elevation and Gradient
    
        This script will take an input point feature (from Step 1) and sample the input raster elevation to find the
        lowest elevation in a user defined search radius and calculate the gradient for each node in the downstream direction.
    
        REQUIREMENTS
        TTools steps 1 must be run before Step 3.
        ESRI ArcGIS Pro
        Python 3.7+
    
        INPUT VARIABLES
        0: nodes_fc:
        Path to the TTools point feature class.
    
        1: searchCells
        The number of cells to search around the node for the lowest elevation.
    
        2: smooth_flag
        Boolean (True/False) flag to indicate if smoothing should occur when the gradient is <= 0. When gradients are <= 0
        the algorithm moves downstream until it a positive gradient is found and then recalculate the gradient over the longer
        distance and applies the resulting value to all nodes over that distance.
    
        3: z_raster:
        Path and name of the ground elevation raster.
    
        4: z_units:
        z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is not in feet or meters the elevation
        values must be converted.
    
        6: block_size:
        The x and y size in kilometers for the z_raster blocks pulled into an array. Start with 5 if you aren't sure and reduce
        if there is an error. To increase processing the z_raster is subdivided iteratively into smaller blocks and pulled
        into arrays for processing. Very large block sizes may use a lot of memory and result in an error.
    
        7: overwrite_data:
        True/False flag if existing data in nodes_fc can be overwritten.
    
        OUTPUTS
        0. nodes_fc:
        New fields ELEVATION and GRADIENT are added into nodes_fc.
    
        """
        # Import system modules
        import sys
        import gc
        import time
        import traceback
        from datetime import timedelta
        import arcpy
        import itertools
        from arcpy import env
        from math import ceil
        from collections import defaultdict
    
        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------
    
        # Future Updates
        # eliminate arcpy and use gdal for reading/writing feature class data
    
        def nested_dict(): 
            """Build a nested dictionary"""
            return defaultdict(nested_dict)
    
        def read_nodes_fc1(nodes_fc, overwrite_data, addFields):
            """Reads the input point feature class and returns the
            NODE_ID and X/Y coordinates as a nested dictionary"""
            nodeDict = nested_dict()
            incursorFields = ["NODE_ID", "SHAPE@X","SHAPE@Y"]
    
            # Get a list of existing fields
            existingFields = []
            for f in arcpy.ListFields(nodes_fc):
                existingFields.append(f.name)
    
            # Check to see if the last field exists if yes add it. 
            # Grabs last field because often the first field, emergent, is zero
            if overwrite_data is False and (addFields[len(addFields)-1] in existingFields) is True:
                incursorFields.append(addFields[len(addFields)-1])
            else:
                overwrite_data = True
        
            # Check to see if all the new fields exist and add them if not
            for f in addFields:
                if (f in existingFields) is False:
                    arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")    
    
            # Determine input point spatial units
            proj = arcpy.Describe(nodes_fc).spatialReference
    
            with arcpy.da.SearchCursor(nodes_fc, incursorFields,"",proj) as Inrows:
                if overwrite_data is True:
                    for row in Inrows:
                        nodeID =row[0]
                        nodeDict[nodeID]["POINT_X"] = row[1]
                        nodeDict[nodeID]["POINT_Y"] = row[2]
                else:
                    for row in Inrows:
                        # Is the data null or zero, if yes grab it.
                        if row[3] is None or row[3] == 0 or row[3] < -9998:
                            nodeID =row[0]
                            nodeDict[nodeID]["POINT_X"] = row[1]
                            nodeDict[nodeID]["POINT_Y"] = row[2]
              
            return nodeDict
    
        def read_nodes_fc2(nodes_fc, overwrite_data, addFields):
            """Reads the input point file, adds new fields, and returns the
            STREAM_ID, STREAM_KM, NODE_ID, LENGTH, ELEVATION, and X/Y coordinates
            as a nested dictionary"""
    
            nodeDict = nested_dict()
            incursorFields = ["STREAM_ID", "STREAM_KM", "NODE_ID", "LENGTH", "ELEVATION", "SHAPE@X","SHAPE@Y"]
    
            # Get a list of existing fields
            existingFields = []
            for f in arcpy.ListFields(nodes_fc):
                existingFields.append(f.name)    
    
            # Check to see if the 1st field exists if yes add it to 
            # the cursorfields to be retrieved.
            if overwrite_data is False and (addFields[0] in existingFields) is True:
                incursorFields.append(addFields[0])
            else:
                overwrite_data = True
    
            # Check to see if all the new fields exist and add them if not
            for f in addFields:
                if (f in existingFields) is False:
                    arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "", "", "NULLABLE", "NON_REQUIRED")    
    
            # Determine input point spatial units
            proj = arcpy.Describe(nodes_fc).spatialReference
    
            with arcpy.da.SearchCursor(nodes_fc,incursorFields,"",proj) as Inrows:
                if overwrite_data is True:
                    for row in Inrows:
                        streamID = row[0]
                        stream_km = row[1]
                        nodeDict[streamID][stream_km]["NODE_ID"] = row[2]
                        nodeDict[streamID][stream_km]["LENGTH"] = row[3]
                        nodeDict[streamID][stream_km]["ELEVATION"] = row[4]
                        nodeDict[streamID][stream_km]["POINT_X"] = row[5]
                        nodeDict[streamID][stream_km]["POINT_Y"] = row[6]
                else:
                    for row in Inrows:
                        # if the data is null or zero (0 = default for shapefile),
                        # it is retrieved and will be overwritten.
                        if row[7] is None or row[7] < -9998:
                            streamID = row[0]
                            stream_km = row[1]                    
                            nodeDict[streamID][stream_km]["NODE_ID"] = row[2] 
                            nodeDict[streamID][stream_km]["LENGTH"] = row[3]
                            nodeDict[streamID][stream_km]["ELEVATION"] = row[4]
                            nodeDict[streamID][stream_km]["POINT_X"] = row[5]
                            nodeDict[streamID][stream_km]["POINT_Y"] = row[6]
            if len(nodeDict) == 0:
                sys.exit("The gradient field checked in the input point feature class "+
                         "have existing data. There is nothing to process. Exiting")    
            return nodeDict
    
        def calculate_gradient(zList, len_list, smooth_flag):
    
            skipupNodes = [0]
            gradientList = [0 for i in zList]
    
            for i in range(1,len(zList)):
                z = zList[i]
                zUp = zList[i - 1 - max(skipupNodes)]
    
                # Check if the gradient is <= 0, if yes keep going until
                # it is positive again. Then recalculate the gradient over 
                # the longer distance		
                if z > zUp and smooth_flag is True:
                    skipupNodes.append(max(skipupNodes) + 1)
                else:
                    dx_meters = sum(len_list[i:i+max(skipupNodes)+1])
                    #dx_meters = float(kmList[i] - kmList[i- max(skipdownNodes)]) * 1000
                    gradient = (zUp - z) / dx_meters
                    for Skip in skipupNodes:
                        gradientList[i-Skip] = gradient
                    skipupNodes = [0]
            
            return (gradientList)
    
        def update_nodes_fc1(nodeDict, nodes_fc, addFields, nodes_to_query):
            """Updates the input point feature class with data from
            the nodes dictionary with node_id as the primary key"""
            #print("Updating input point feature class")
    
            # Build a query to retrieve just the nodes that needs updating
            whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))
    
            with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:
                for row in cursor:
                    for f, field in enumerate(addFields):
                        nodeID =row[0]
                        row[f+1] = nodeDict[nodeID][field]
                        cursor.updateRow(row)
    
        def update_nodes_fc2(nodeDict, nodes_fc, addFields, nodes_to_query):
            """Updates the input point feature class with data from
            the nodes dictionary with stream IDs as the primary key"""
            arcpy.AddMessage("Updating input point feature class")
            print("Updating input point feature class")
    
            # Build a query to retreive just the nodes that needs updating
            whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))
    
            with arcpy.da.UpdateCursor(nodes_fc,["STREAM_ID","NODE_ID","STREAM_KM"] + addFields, whereclause) as cursor:  
                for row in cursor:
                    for f, field in enumerate(addFields):
                        streamID = row[0]
                        stream_km =row[2]
                        row[f+3] = nodeDict[streamID][stream_km][field]
                        cursor.updateRow(row)

    def create_block_list(nodeDict, nodes, block_size, buffer):
        """Returns two lists, one containing the coordinate extent
        for each block that will be iteratively extracted to an array
        and the other containing the stream and node IDs within each
        block extent."""

        print("Calculating block extents")
        x_coord_list = [nodeDict[nodeID]["POINT_X"] for nodeID in nodes]
        y_coord_list = [nodeDict[nodeID]["POINT_Y"] for nodeID in nodes]

        # calculate bounding box extent for samples
        x_min = min(x_coord_list)
        x_max = max(x_coord_list)
        y_min = min(y_coord_list)
        y_max = max(y_coord_list)

        x_width = int(x_max - x_min + 1)
        y_width = int(y_max - y_min + 1)

        block_extents = []
        block_nodes = []

        # Build data blocks
        for x in range(0, x_width, block_size):
            for y in range(0, y_width, block_size):

                # Lower left coordinate of block (in map units)
                block0_x_min = min([x_min + x, x_max])
                block0_y_min = min([y_min + y, y_max])
                # Upper right coordinate of block (in map units)
                block0_x_max = min([block0_x_min + block_size, x_max])
                block0_y_max = min([block0_y_min + block_size, y_max])

                block_x_max = block0_x_max
                block_x_min = block0_x_min
                block_y_max = block0_y_max
                block_y_min = block0_y_min

                nodes_in_block = []
                for nodeID in nodes:
                    node_x = nodeDict[nodeID]["POINT_X"]
                    node_y = nodeDict[nodeID]["POINT_Y"]
                    if (block0_x_min <= node_x <= block0_x_max and
                            block0_y_min <= node_y <= block0_y_max):
                        nodes_in_block.append([nodeID, node_x, node_y])

                if nodes_in_block:

                    # Minimize the size of the block0 by the true
                    # extent of the nodes in the block
                    node_x_min = min([nodeDict[nodeID]["POINT_X"] for nodeID in nodes_in_block])
                    node_y_min = min([nodeDict[nodeID]["POINT_Y"] for nodeID in nodes_in_block])
                    node_x_max = max([nodeDict[nodeID]["POINT_X"] for nodeID in nodes_in_block])
                    node_y_max = max([nodeDict[nodeID]["POINT_Y"] for nodeID in nodes_in_block])

                    if block0_x_min < node_x_min: block_x_min = node_x_min
                    if block0_x_max > node_x_max: block_x_max = node_x_max
                    if block0_y_min < node_y_min: block_y_min = node_y_min
                    if block0_y_max > node_y_max: block_y_max = node_y_max

                    # Add the block extent for processing and the buffer distance.
                    # Because the node coordinates are used to determine if the node is within the block extent,
                    # a buffer distance is added to ensure each block extent will include all samples around nodes
                    # located at the edge of the block.

                    # order 0 left,      1 bottom,    2 right,     3 top
                    block_extents.append((block_x_min - buffer, block_y_min - buffer,
                                          block_x_max + buffer, block_y_max + buffer))
                    block_nodes.append(nodes_in_block)

        return block_extents, block_nodes
    
        def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
            """converts x/y coordinates to col and row of the array"""
            col_x = int(((easting - block_x_min) - ((easting - block_x_min) % x_cellsize)) / x_cellsize)
            row_y = int(((block_y_max - northing) - ((block_y_max - northing) % y_cellsize)) / y_cellsize)
            return [col_x, row_y]
    
        def sample_raster(block, nodes_in_block, z_raster, cellcoords, con_z_to_m):
    
            if con_z_to_m is not None:
                nodata_to_value = -9999 / con_z_to_m
            else:
                nodata_to_value = -9999
        
            # localize the block extent values
            block_x_min = block[0]
            block_y_min = block[1]
            block_x_max = block[2]
            block_y_max = block[3]
    
            x_cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEX").getOutput(0))
            y_cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEY").getOutput(0))   
    
            # Get the coordinates extent of the input raster
            raster_x_min = float(arcpy.GetRasterProperties_management(z_raster, "LEFT").getOutput(0))
            raster_y_min = float(arcpy.GetRasterProperties_management(z_raster, "BOTTOM").getOutput(0))
            raster_x_max = float(arcpy.GetRasterProperties_management(z_raster, "RIGHT").getOutput(0))
            raster_y_max = float(arcpy.GetRasterProperties_management(z_raster, "TOP").getOutput(0))
    
            # Calculate the block x and y offset from the raster and adjust
            # the block coordinates so they are at the raster cell corners.
            # This is for ESRI's RastertoNumpyArray function which defaults to the adjacent
            # lower left cell
            block_x_min_corner = block_x_min - ((block_x_min - raster_x_min) % x_cellsize)
            block_y_min_corner = block_y_min - ((block_y_min - raster_y_min) % y_cellsize)
            block_x_max_corner = block_x_max + ((raster_x_max - block_x_max) % x_cellsize)
            block_y_max_corner = block_y_max + ((raster_y_max - block_y_max) % y_cellsize)
    
            # calculate the number of cols/rows from the lower left
            ncols = int((block_x_max_corner - block_x_min_corner) / x_cellsize)
            nrows = int((block_y_max_corner - block_y_min_corner) / y_cellsize)
    
            # Construct the array. Note returned array is (row, col) so (y, x)
            try:
                raster_array = arcpy.RasterToNumPyArray(z_raster, arcpy.Point(block_x_min_corner, block_y_min_corner),
                                                        ncols, nrows, nodata_to_value)
            except:
                tbinfo = traceback.format_exc()
                pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
                sys.exit(pymsg)    
    
            # convert array values to meters if needed
            if con_z_to_m is not None:
                raster_array = raster_array * con_z_to_m
    
            z_list = []
            if raster_array.max() > -9999:
                # There is at least one pixel of data
                for node in nodes_in_block:
                    xy = coord_to_array(node[1], node[2], block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize)
                    z_sampleList = []
                    for coord in cellcoords:
                        # Calculate the cell X/Y based on the base coordinate movement
                        cell_x = xy[0] + coord[0]
                        cell_y = xy[1] + coord[1]
                        z_sampleList.append(raster_array[cell_y,cell_x])
                
                    # sample at node:
                    z_node = raster_array[xy[1], xy[0]]
                
                    # Remove no data values (-9999) unless they are all no data
                    if not max(z_sampleList) < -9998:
                        z_sampleList = [z for z in z_sampleList if z > -9999]
                    # Get the lowest elevation        
                    node.append(min(z_sampleList))
                    node.append(z_node)
                    z_list.append(node)
    
            else:
                # No data, add -9999 for elevation and z_node
                for node in nodes_in_block:
                    node.append(-9999)
                    node.append(-9999)
                    z_list.append(node)
        
            return z_list
    
        def from_z_units_to_meters_con(zUnits):
            """Returns the conversion factor to get from the input z
            units to meters"""
        
            try:
                con_z_to_m = float(zUnits)
            except:
                if zUnits == "Meters":
                    con_z_to_m = 1.0 
                elif zUnits == "Feet":
                    con_z_to_m = 0.3048
                else: con_z_to_m = None # The conversion factor will not be used
    
            return con_z_to_m
    
        def from_meters_con(inFeature):
            """Returns the conversion factor to get from
            meters to the spatial units of the input feature class"""
            try:
                con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
            except:
                arcpy.AddError("{0} has a coordinate system that ".format(inFeature)+
                               "is not projected or not recognized. Use a "+
                               "projected coordinate system preferably in "+
                               "linear units of feet or meters.")
                sys.exit("Coordinate system is not projected or not recognized. "+
                         "Use a projected coordinate system, preferably in "+
                         "linear units of feet or meters.")   
            return con_from_m
    
        #enable garbage collection
        gc.enable()

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        searchCells = int(parameters[1].value)
        smooth_flag = bool(parameters[2].value)
        z_raster = parameters[3].valueAsText
        z_units = parameters[4].valueAsText
        block_size = float(parameters[5].value)
        overwrite_data = bool(parameters[6].value)
    
        try:
            arcpy.AddMessage("Step 3: Sample Stream Elevations/Gradient")
            print("Step 3: Sample Stream Elevations/Gradient")
    
            #keeping track of time
            startTime= time.time()
    
            # Check if the node fc exists
            if not arcpy.Exists(nodes_fc):
                arcpy.AddError("\nThis output does not exist: \n" +
                               "{0}\n".format(nodes_fc))
                sys.exit("This output does not exist: \n" +
                         "{0}\n".format(nodes_fc))     
    
            if overwrite_data is True: 
                env.overwriteOutput = True
            else:
                env.overwriteOutput = False
    
            # Determine input point spatial units
            proj = arcpy.Describe(nodes_fc).spatialReference
            proj_ele = arcpy.Describe(z_raster).spatialReference
    
            # Check to make sure the raster and input 
            # points are in the same projection.
            if proj.name != proj_ele.name:
                arcpy.AddError("{0} and {1} do not ".format(nodes_fc,z_raster)+
                               "have the same projection."+
                               "Please reproject your data.")
                sys.exit("Input points and elevation raster do not have the "+
                         "same projection. Please reproject your data.")
    
            if block_size == "#": block_size = 5
    
            # Get the units conversion factor
            con_z_to_m = from_z_units_to_meters_con(z_units)
            con_from_m = from_meters_con(nodes_fc)
    
            # convert block size from km to meters to units of the node fc
            # in the future block size should be estimated based on available memory
            # memorysize = datatypeinbytes*nobands*block_size^2
            # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
            if block_size in ["#", ""]:
                block_size = int(con_from_m * 5000)
            else:
                block_size = int(con_from_m * block_size * 1000)
    
            # Get the elevation raster cell size
            cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEX").getOutput(0))
    
            # calculate the buffer distance (in raster spatial units) to add to 
            # the base bounding box when extracting to an array. The buffer is 
            # equal to the cellsize * the searchcells to make sure the block includes 
            # the surrounding cells at each corner
            buffer = int((searchCells + 1)* cellsize)    
    
            # Make a list of the base x/y coordinate movements
            # from the node origin. These values will be 
            # multiplied by the cell size.
            # searchCells = 0 samples at the node
            # searchCells = 1 cell width around node = 9 cells
            # searchCells = 2 cell widths around node = 25 cells ...
            cell_moves = [i for i in range(searchCells*-1, searchCells+1, 1)]
            cellcoords = list(itertools.product(cell_moves, cell_moves))
    
            # read the data into a nested dictionary
            addFields = ["ELEVATION", "Z_NODE"]
            nodeDict = read_nodes_fc1(nodes_fc, overwrite_data, addFields)
            if len(nodeDict) != 0:
        
                # Get a list of the nodes, sort them
                nodes = list(nodeDict.keys())
                nodes.sort()
                n_nodes = len(nodes)
        
                # Build the block list
                block_extents, block_nodes = create_block_list(nodeDict, nodes, buffer, block_size)
        
                # Iterate through each block, calculate sample coordinates,
                # convert raster to array, sample the raster
                total_samples = 0
        
                for p, block in enumerate(block_extents):
                    nodes_in_block = block_nodes[p]
                    nodes_in_block.sort()
            
                    arcpy.AddMessage("Processing block {0} of {1}".format(p + 1, len(block_extents)))
                    print("Processing block {0} of {1}".format(p + 1, len(block_extents)))
        
                    z_list = sample_raster(block, nodes_in_block, z_raster, cellcoords, con_z_to_m)
            
                    # Update the node fc
                    for row in z_list:
                        nodeDict[row[0]]["ELEVATION"] = row[3]
                        nodeDict[row[0]]["Z_NODE"] = row[4]
                
                    nodes_to_query = [row[0] for row in z_list]
            
                    # Write the elevation data to the TTools point feature class 
                    update_nodes_fc1(nodeDict, nodes_fc, addFields, nodes_to_query)
        
                    total_samples = total_samples + len(z_list)
                    del z_list
                    gc.collect()
        
            else:
                arcpy.AddMessage("The elevation field checked in the input point feature class " +
                      "have existing data. Advancing to gradient processing")
                print("The elevation field checked in the input point feature class " +
                      "have existing data. Advancing to gradient processing")
            del(nodeDict)
    
            # Start on gradients
    
            # read the data into a nested dictionary
            addFields = ["GRADIENT"]
            nodeDict = read_nodes_fc2(nodes_fc, overwrite_data, addFields)    
            nodes_to_query = []
    
            for n, streamID in enumerate(nodeDict):
                arcpy.AddMessage("Calculating gradients stream {0} of {1}".format(n + 1, len(nodeDict)))
                print("Calculating gradients stream {0} of {1}".format(n + 1, len(nodeDict)))
            
                stream_kms = list(nodeDict[streamID].keys())
                stream_kms.sort(reverse=True)
    
                z_list = [nodeDict[streamID][km]["ELEVATION"] for km in stream_kms]
                len_list = [nodeDict[streamID][km]["LENGTH"] for km in stream_kms]
        
                # Calculate Gradient
                gradientList = calculate_gradient(z_list, len_list, smooth_flag)
        
                for i, km in enumerate(stream_kms):
                    nodeDict[streamID][km]["GRADIENT"] = gradientList[i]
                    nodes_to_query.append(nodeDict[streamID][km]["NODE_ID"])
    
            update_nodes_fc2(nodeDict, nodes_fc, addFields, nodes_to_query)
    
            endTime = time.time()
            gc.collect()  
    
            elapsedmin= ceil(((endTime - startTime) / 60)* 10)/10
            mspernode = timedelta(seconds=(endTime - startTime) / n_nodes).microseconds
            arcpy.AddMessage("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))    
            print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))    
    
        # For arctool errors
        except arcpy.ExecuteError:
            msgs = arcpy.GetMessages(2)
            arcpy.AddError(msgs)
            print(msgs)
    
        # For other errors
        except:
            tbinfo = traceback.format_exc()
    
            pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" +str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"
    
            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)
    
            print(pymsg)
            print(msgs)


class MeasureTopographicAngles(object):
    def __init__(self):
        self.label = "Step 4: Measure Topographic Angles"
        self.description = "Calculate the maximum topographic shade angle for each stream node."
        self.canRunInBackground = False

    def getParameterInfo(self):

        param0 = arcpy.Parameter(
            displayName="Stream Nodes Feature Class",
            name="nodes_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Topographic Direction Option",
            name="topo_directions",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param1.filter.type = "ValueList"
        param1.filter.list = [1, 2, 3]
        param1.value = 1

        param2 = arcpy.Parameter(
            displayName="Max Search Distance (km)",
            name="searchDistance_max_km",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        param3 = arcpy.Parameter(
            displayName="Elevation Raster",
            name="z_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        param4 = arcpy.Parameter(
            displayName="Elevation Units",
            name="z_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param4.filter.type = "ValueList"
        param4.filter.list = ["Meters", "Feet"]

        param5 = arcpy.Parameter(
            displayName="Output Topographic Feature Class",
            name="topo_fc",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")
        param5.overwrite = True

        param6 = arcpy.Parameter(
            displayName="Raster Block Size (km)",
            name="block_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        param7 = arcpy.Parameter(
            displayName="Overwrite Existing Data",
            name="overwrite_data",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param7.value = True

        return [param0, param1, param2, param3, param4, param5, param6, param7]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

    # Log input parameters to file
        log_run(self.label, parameters)
        """
        TTools Step 4: Measure Topographic Shade Angles

        This script will take an input point feature (from Step 1) and calculate the maximum topographic 
        shade angle for each stream node in different directions.

        REQUIREMENTS
        TTools steps 1 - 3 must be run before Step 4.
        ESRI ArcGIS Pro w/ Spatial Analyst extension
        Python 3.7+

        INPUT VARIABLES
        0: nodes_fc:
        Path to the TTools point feature class.

        1: topo_directions
        An integer value corresponding to the specific list of azimuth directions to sample (e.g. topo_directions = 1).
        Options listed below.

        1. Heat Source or Washington Department of Ecology's Shade model: [270, 180, 90]
        2. All cardinal and intercardinal directions: [45, 90, 135, 180, 225, 270, 315, 360]
        3. CE-QUAL-W2: [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340]

        2: searchDistance_max_km
        The maximum distance in kilometers to search for the largest topographic shade angle.

        3: z_raster:
        Path and name of the ground elevation raster.

        4: z_units:
        z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is not in feet or meters the elevation
        values must be converted.

        5: topo_fc:
        Path and name of output topographic point feature. This feature identifies the location on the z_raster producing
        the largest topographic shade angle for each node and sample direction.

        6: block_size:
        The x and y size in kilometers for the z_raster blocks pulled into an array. Start with 5 if you aren't sure and reduce
        if there is an error. To increase processing the z_raster is subdivided iteratively into smaller blocks and pulled
        into arrays for processing. Very large block sizes may use a lot of memory and result in an error.

        7: overwrite_data:
        True/False flag if existing data in nodes_fc and topo_fc can be overwritten.

        # OUTPUTS
        0. nodes_fc:
        New fields are added into nodes_fc with the maximum topographic shade angles for each direction at each node.

        1. topo_fc: - point for each x/y location of the maximum elevation.
        New point feature class created with a point at each x/y location that produces the largest topographic shade angle
        for each node and sample direction.

        """
        # Import system modules
        import sys
        import os
        import gc
        import time
        import traceback
        from datetime import timedelta
        import arcpy
        from arcpy import env
        from math import radians, sin, cos, hypot, ceil
        from collections import defaultdict
        import numpy as np

        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------

        # Future Updates
        # update fc after each topo line sample so data isn't lost if the script crashes
        # eliminate arcpy and use gdal for reading/writing feature class data

        # The below are used for debugging. Currently, turned off.
        topo_line_fc = r"C:\workspace\ttools_tests\TTools_py39\jc_test_py39.gdb\topo_line"
        block_fc = r"C:\workspace\ttools_tests\TTools_py39\jc_test_py39.gdb\blocks"
        plot_dir = r"C:\workspace\ttools_tests\TTools_py39\plots"

        def nested_dict():
            """Build a nested dictionary"""
            return defaultdict(nested_dict)


        def read_nodes_fc(nodes_fc, overwrite_data, addFields):
            """Reads the input point feature class and returns the STREAM_ID,
            NODE_ID, and X/Y coordinates as a nested dictionary"""

            arcpy.AddMessage("Reading nodes feature class")
            print("Reading nodes feature class")

            nodeDict = nested_dict()
            incursorFields = ["NODE_ID", "STREAM_ID", "Z_NODE", "SHAPE@X", "SHAPE@Y"]

            # Get a list of existing fields
            existingFields = []
            for f in arcpy.ListFields(nodes_fc):
                existingFields.append(f.name)

                # Check to see if the 1st field exists if yes add it.
            if overwrite_data is False and (addFields[0] in existingFields) is True:
                incursorFields.append(addFields[0])
            else:
                overwrite_data = True

            # Determine input point spatial units
            proj = arcpy.Describe(nodes_fc).spatialReference

            with arcpy.da.SearchCursor(nodes_fc, incursorFields, "", proj) as Inrows:
                if overwrite_data:
                    for row in Inrows:
                        nodeDict[row[0]]["STREAM_ID"] = row[1]
                        nodeDict[row[0]]["Z_NODE"] = row[2]
                        nodeDict[row[0]]["POINT_X"] = row[3]
                        nodeDict[row[0]]["POINT_Y"] = row[4]

                else:
                    for row in Inrows:
                        # if the data is null or zero (0 = default for shapefile),
                        # it is retrieved and will be overwritten.
                        if row[5] is None or row[5] == 0 or row[5] < -9998:
                            nodeDict[row[0]]["STREAM_ID"] = row[1]
                            nodeDict[row[0]]["Z_NODE"] = row[2]
                            nodeDict[row[0]]["POINT_X"] = row[3]
                            nodeDict[row[0]]["POINT_Y"] = row[4]
            if len(nodeDict) == 0:
                sys.exit("The fields checked in the input point feature class " +
                         "have existing data. There is nothing to process. Exiting")

            return (nodeDict)


        def create_topo_line_fc(topo_line, streamID, nodeID, a, topo_line_fc, proj):
            poly_array = arcpy.Array()
            pnt = arcpy.Point()
            cursorfields = ["STREAM_ID", "NODE_ID", "AZIMUTH"]

            # Check to see if the block fc exists, if not create it
            if not arcpy.Exists(topo_line_fc):
                # Create an empty output with the same projection as the input polyline
                arcpy.CreateFeatureclass_management(os.path.dirname(topo_line_fc),
                                                    os.path.basename(topo_line_fc),
                                                    "POLYLINE", "", "DISABLED", "DISABLED",
                                                    proj)

                # Determine Stream ID field properties
                sid_type = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].type
                sid_precision = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].precision
                sid_scale = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].scale
                sid_length = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].length

                # Add attribute fields # TODO add dictionary
                # of field types so they aren't all double
                for f in cursorfields:
                    if f == "STREAM_ID":
                        arcpy.AddField_management(topo_line_fc, f, sid_type,
                                                  sid_precision, sid_scale,
                                                  sid_length, "", "NULLABLE",
                                                  "NON_REQUIRED")
                    else:
                        arcpy.AddField_management(topo_line_fc, f, "DOUBLE", "", "", "",
                                                  "", "NULLABLE", "NON_REQUIRED")

            with arcpy.da.InsertCursor(topo_line_fc, ["SHAPE@"] + cursorfields) as cursor:
                for pnt_x, pnt_y in topo_line:
                    pnt.X = pnt_x
                    pnt.Y = pnt_y
                    poly_array.add(pnt)
                poly = arcpy.Polyline(poly_array)
                cursor.insertRow([poly, streamID, nodeID, a])
                poly_array.removeAll()


        def update_topo_fc(topo_list, topo_fc, nodes_fc, nodes_to_update, overwrite_data, proj):
            """Creates/updates the output topo point feature
            class using the data from the topo list"""

            # Create an empty output with the same projection as the input polyline
            cursorfields = ["POINT_X", "POINT_Y", "STREAM_ID", "NODE_ID",
                            "AZIMUTH", "TOPOANGLE", "TOPO_ELE", "NODE_ELE",
                            "ELE_CHANGE", "TOPODIS", "SEARCHDIS", "NA_SAMPLES"]

            # Check if the output exists and create if not
            if not arcpy.Exists(topo_fc):
                arcpy.CreateFeatureclass_management(os.path.dirname(topo_fc),
                                                    os.path.basename(topo_fc),
                                                    "POINT", "", "DISABLED", "DISABLED", proj)

                # Determine Stream ID field properties
                sid_type = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].type
                sid_precision = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].precision
                sid_scale = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].scale
                sid_length = arcpy.ListFields(nodes_fc, "STREAM_ID")[0].length

                # Add attribute fields # TODO add dictionary of field types
                # so they aren't all double
                for f in cursorfields:
                    if f == "STREAM_ID":
                        arcpy.AddField_management(topo_fc, f, sid_type,
                                                  sid_precision, sid_scale, sid_length,
                                                  "", "NULLABLE", "NON_REQUIRED")

                    else:
                        arcpy.AddField_management(topo_fc, f, "DOUBLE", "", "", "",
                                                  "", "NULLABLE", "NON_REQUIRED")

            if not overwrite_data:
                # Build a query to retrieve existing rows from the nodes
                # that need updating
                whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_update))

                # delete those rows
                with arcpy.da.UpdateCursor(topo_fc, ["NODE_ID"], whereclause) as cursor:
                    for row in cursor:
                        cursor.deleteRow()

            with arcpy.da.InsertCursor(topo_fc, ["SHAPE@X", "SHAPE@Y"] +
                                                cursorfields) as cursor:
                for row in topo_list:
                    cursor.insertRow(row)


        def update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_query):
            """Updates the input point feature class with data from the nodes dictionary"""

            # Build a query to retrieve just the nodes that needs updating
            whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))

            with arcpy.da.UpdateCursor(nodes_fc, ["NODE_ID"] + addFields, whereclause) as cursor:
                for row in cursor:
                    for f, field in enumerate(addFields):
                        nodeID = row[0]
                        row[f + 1] = nodeDict[nodeID][field]
                        cursor.updateRow(row)


        def to_meters_con(inFeature):
            """Returns the conversion factor to get
            from the input spatial units to meters"""
            try:
                con_to_m = arcpy.Describe(inFeature).SpatialReference.metersPerUnit
            except:
                arcpy.AddError("{0} has a coordinate system that is".format(inFeature) +
                               "not projected or not recognized. Use a projected " +
                               "coordinate system preferably in linear units of " +
                               "feet or meters.")
                sys.exit("Coordinate system is not projected or not recognized. " +
                         "Use a projected coordinate system, preferably in" +
                         "linear units of feet or meters.")
            return con_to_m


        def from_meters_con(inFeature):
            """Returns the conversion factor to get from
            meters to the spatial units of the input feature class"""
            try:
                con_from_m = 1 / arcpy.Describe(inFeature).SpatialReference.metersPerUnit
            except:
                arcpy.AddError("{0} has a coordinate system that ".format(inFeature) +
                               "is not projected or not recognized. Use a " +
                               "projected coordinate system preferably in " +
                               "linear units of feet or meters.")
                sys.exit("Coordinate system is not projected or not recognized. " +
                         "Use a projected coordinate system, preferably in " +
                         "linear units of feet or meters.")
            return con_from_m


        def from_z_units_to_meters_con(zUnits):
            """Returns the conversion factor to
            get from the input z units to meters"""

            try:
                con_z_to_m = float(zUnits)
            except:
                if zUnits == "Meters":
                    con_z_to_m = 1.0
                elif zUnits == "Feet":
                    con_z_to_m = 0.3048
                else:
                    con_z_to_m = None  # The conversion factor will not be used

            return con_z_to_m


        def build_search_array(searchDistance_min, searchDistance_max, cellsize, use_skippy):
            """Build a numpy array from the minimum to the max
            search distance by increments of the cellsize."""

            # use next cell over to avoid divide by zero errors
            if searchDistance_min <= 0: searchDistance_min = cellsize

            if use_skippy:
                # This is a modified version of the skippy
                # algorithm from Greg Pelletier. It is not being used but I'm
                # keeping it in here in case someone wants to turn it on.
                # It needs to be fixed so the ncells are adjusted based on distance
                searchDistance = searchDistance_min
                distanceList = [searchDistance_min]

                # ncells = max([int(ceil((block_x_max - block_x_min)/ x_cellsize)), 1])

                ncells = 1
                while not searchDistance > searchDistance_max:
                    if ncells <= 10:
                        searchDistance_m = searchDistance + (cellsize)
                    if 10 < ncells <= 20:
                        searchDistance = searchDistance + (cellsize * 3)
                    if 20 < ncells <= 40:
                        searchDistance = searchDistance + (cellsize * 6)
                    if 40 < ncells <= 50:
                        searchDistance = searchDistance + (cellsize * 12)
                    if 50 < ncells <= 60:
                        searchDistance = searchDistance + (cellsize * 25)
                    if ncells > 60:
                        searchDistance = searchDistance + (cellsize * 50)
                    distanceList.append(searchDistance)
                    ncells = ncells + 1
                distance_array = np.array(distanceList)
            else:
                if (searchDistance_max - searchDistance_min >= cellsize):
                    distance_array = np.arange(searchDistance_min, searchDistance_max, cellsize)
                else:
                    distance_array = np.array([searchDistance_min])
            return distance_array


        def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
            """converts x/y coordinates to col and row of the array"""
            col_x = int(((easting - block_x_min) - ((easting - block_x_min) % x_cellsize)) / x_cellsize)
            row_y = int(((block_y_max - northing) - ((block_y_max - northing) % y_cellsize)) / y_cellsize)

            return [col_x, row_y]

        def plot_it(pts1, pts2, nodeID, a, b, b0, plot_dir):
            """plots the block and topo line"""

            import matplotlib.pyplot as plt

            x1 = [i[0] for i in pts1]
            y1 = [i[1] for i in pts1]

            x2 = [i[0] for i in pts2]
            y2 = [i[1] for i in pts2]

            plt.plot(x1, y1, 'b--', x2, y2, 'r--')
            plt.savefig(r"{0}\node{1}_a{2}_{3}b_{4}b0.png".format(plot_dir, nodeID, a, b, b0))


        def create_block_fc(block_segments, b, block_fc, proj):
            """Creates a poly line feature class of the block segments"""

            poly_array = arcpy.Array()
            pnt = arcpy.Point()
            cursorfields = ["BLOCK", "SEGMENT"]

            # Check to see if the block fc exists, if not create it
            if not arcpy.Exists(block_fc):
                # Create an empty output with the same projection as the input polyline
                arcpy.CreateFeatureclass_management(os.path.dirname(block_fc),
                                                    os.path.basename(block_fc),
                                                    "POLYLINE", "", "DISABLED", "DISABLED",
                                                    proj)

                # Add attribute fields
                for f in cursorfields:
                    arcpy.AddField_management(block_fc, f, "DOUBLE", "", "", "",
                                              "", "NULLABLE", "NON_REQUIRED")

            with arcpy.da.InsertCursor(block_fc, ["SHAPE@"] + cursorfields) as cursor:
                for s, segment in enumerate(block_segments):
                    for pnt_x, pnt_y in segment:
                        pnt.X = pnt_x
                        pnt.Y = pnt_y
                        poly_array.add(pnt)
                    poly = arcpy.Polyline(poly_array)
                    cursor.insertRow([poly, b, s])
                    poly_array.removeAll()


        def create_blocks(NodeDict, block_size, last_azimuth, searchDistance_max):
            """Returns two lists, one containing the coordinate extent
            for each block that will be iteratively extracted to an array
            and the other containing the start and stop distances for each
            topo line that is within the block extent."""

            arcpy.AddMessage("Preparing blocks and topo sampling")
            print("Preparing blocks and topo sampling")

            # Create a dictionary to lookup which nodesIDs can be updated after
            # the block has been sampled.
            blockDict = nested_dict()

            # Get a list of the nodes, sort them
            nodes = list(nodeDict.keys())
            nodes.sort()

            topo_list = []
            x_coord_list = []
            y_coord_list = []

            for nodeID in nodes:
                node_x = nodeDict[nodeID]["POINT_X"]
                node_y = nodeDict[nodeID]["POINT_Y"]
                streamID = nodeDict[nodeID]["STREAM_ID"]
                z_node = nodeDict[nodeID]["Z_NODE"]

                for a in azimuths:
                    # calculate x/y coordinates at max search distance
                    end_x = ((searchDistance_max * sin(radians(a))) + node_x)
                    end_y = ((searchDistance_max * cos(radians(a))) + node_y)

                    x_coord_list.append(end_x)
                    y_coord_list.append(end_y)

                    topo_list.append([nodeID, streamID, a, z_node, node_x, node_y, end_x, end_y])

            # calculate bounding box extent for samples
            x_min = min(x_coord_list)
            x_max = max(x_coord_list)
            y_min = min(y_coord_list)
            y_max = max(y_coord_list)

            x_width = int(x_max - x_min + 1)
            y_width = int(y_max - y_min + 1)

            block_extents = []
            block_samples = []
            b = 0

            # Build blocks
            for x in range(0, x_width, block_size):
                for y in range(0, y_width, block_size):

                    # Lower left coordinate of block (in map units)
                    block_x_min = min([x_min + x, x_max])
                    block_y_min = min([y_min + y, y_max])
                    # Upper right coordinate of block (in map units)
                    block_x_max = min([block_x_min + block_size, x_max])
                    block_y_max = min([block_y_min + block_size, y_max])

                    nodes_to_update = []
                    topo_in_block = []

                    block_segments = (((block_x_min, block_y_max), (block_x_min, block_y_min)),
                                      ((block_x_min, block_y_min), (block_x_max, block_y_min)),
                                      ((block_x_max, block_y_min), (block_x_max, block_y_max)),
                                      ((block_x_max, block_y_max), (block_x_min, block_y_max)))

                    # --------------------------------------------------------
                    # This is an argument for plot_it(). I used it for debugging.
                    # It's a tuple of the x/y coords for each block segment.
                    block_for_plot = ((block_x_min, block_y_max),
                                      (block_x_min, block_y_min),
                                      (block_x_min, block_y_min),
                                      (block_x_max, block_y_min),
                                      (block_x_max, block_y_min),
                                      (block_x_max, block_y_max),
                                      (block_x_max, block_y_max),
                                      (block_x_min, block_y_max))
                    # --------------------------------------------------------

                    # Now start iterating through the topo list to evaluate
                    # if any part of the topo line is in the block extent
                    for nodeID, streamID, a, z_node, node_x, node_y, end_x, end_y in topo_list:

                        # --------------------------------------------------------
                        # This was used for debugging.
                        # It slows the script WAY down.
                        # topo_line = ((node_x, node_y),(end_x, end_y))
                        # plot_it(block_for_plot, topo_line, nodeID, a, b, b0, plot_dir)
                        # create_topo_line_fc(topo_line, streamID, nodeID, a,
                        #                    topo_line_fc, proj)
                        # --------------------------------------------------------

                        contains_node = False
                        contains_end = False
                        last_sample = False

                        # check if the node is inside the block
                        if (block_x_min <= node_x <= block_x_max and
                                block_y_min <= node_y <= block_y_max):
                            contains_node = True

                        # check if the topo line end point is inside the block
                        if (block_x_min <= end_x <= block_x_max and
                                block_y_min <= end_y <= block_y_max):
                            contains_end = True

                        # check if this is the last block to process for this node
                        if a == last_azimuth:
                            if a == 45:
                                # we have to get the coordinate at
                                # the far upper right corner
                                searchDistance_last = (hypot(searchDistance_max, searchDistance_max))
                                last_x = ((searchDistance_last * sin(radians(a))) + node_x)
                                last_y = ((searchDistance_last * cos(radians(a))) + node_y)
                            else:
                                last_x = end_x
                                last_y = end_y

                            if (block_x_min <= last_x <= block_x_max and
                                    block_y_min <= last_y <= block_y_max):
                                last_sample = True

                                # check if the entire topo line segment is inside the block
                        if contains_node and contains_end:
                            block_search_start = 0
                            block_search_end = searchDistance_max

                            topo_in_block.append([nodeID, streamID, a,
                                                  z_node,
                                                  node_x, node_y,
                                                  end_x, end_y,
                                                  block_search_start,
                                                  block_search_end])

                            if last_sample: nodes_to_update.append(nodeID)

                        # check if and where the topo segment cross the block
                        else:
                            distance = []

                            # check each block segment for an intersection
                            for i, block_segment in enumerate(block_segments):

                                intersects, inter1_x, inter1_y, inter2_x, inter2_y = find_intersection(block_segment[0],
                                                                                                       block_segment[1],
                                                                                                       (node_x, node_y),
                                                                                                       (end_x, end_y), True)

                                # if there is an intersection calculate the
                                # distance in fc units from node to
                                # that intersection
                                if intersects:
                                    if a in [0, 180]:
                                        # need to use the y direction for these
                                        distance.append((inter1_y - node_y) / cos(radians(a)))
                                    else:
                                        distance.append((inter1_x - node_x) / sin(radians(a)))

                            if (len(distance) == 1 and
                                    (0 < distance[0] < searchDistance_max) and
                                    (contains_node or contains_end)):
                                # one intersection
                                if contains_node:
                                    # This will be changed from zero to
                                    # the cell size in build_search_array()
                                    # Should probably just change it here
                                    block_search_start = 0
                                    block_search_end = distance[0]
                                elif contains_end:
                                    # end of the topo line
                                    block_search_start = distance[0]
                                    block_search_end = searchDistance_max

                                # part of the topo line is in the block, add it
                                topo_in_block.append([nodeID, streamID, a,
                                                      z_node,
                                                      node_x, node_y,
                                                      end_x, end_y,
                                                      block_search_start,
                                                      block_search_end])

                                if last_sample and not nodeDict[nodeID]["updated"]:
                                    nodes_to_update.append(nodeID)
                                    nodeDict[nodeID]["updated"] = True

                            elif len(distance) > 1:
                                # two intersections, crosses the block
                                # two intersections, end or start on block line
                                # three intersections, collinear over length of block
                                # three intersections, collinear w/ end or start on block line
                                block_search_start = min(i for i in distance if i is not None)
                                block_search_end = max(distance)

                                topo_in_block.append([nodeID, streamID, a,
                                                      z_node,
                                                      node_x, node_y,
                                                      end_x, end_y,
                                                      block_search_start,
                                                      block_search_end])

                                if last_sample and contains_end and not nodeDict[nodeID]["updated"]:
                                    nodes_to_update.append(nodeID)
                                    nodeDict[nodeID]["updated"] = True

                            elif last_sample and not contains_end and not nodeDict[nodeID]["updated"]:
                                nodes_to_update.append(nodeID)
                                nodeDict[nodeID]["updated"] = True

                            del distance[:]

                    if topo_in_block:
                        # order 0 left,      1 bottom,    2 right,     3 top
                        blockDict[b]["extent"] = (block_x_min, block_y_min,
                                                  block_x_max, block_y_max)
                        blockDict[b]["samples"] = topo_in_block
                        blockDict[b]["nodes_to_update"] = nodes_to_update

                        # --------------------------------------------------
                        # Creates a feature class of the blocks.
                        # create_block_fc(block_segments, b, block_fc, proj)
                        # --------------------------------------------------

                    b = b + 1
            return blockDict


        def find_intersection(a, b, c, d, check_collinear=True):
            """Calculates 2D intersection coordinates of segments a-b and c-d.
            a,b,c,d are tuples in the form of (x,y). If checking for collinearity
            then segments that overlap or touch but do not cross will be
            considered an intersection. Returns the min and max points of
            intersection for overlap and the same points for intersections."""

            # a-b = block segment
            # c-d = topo line

            Dx_Cx = d[0] - c[0]
            Ay_Cy = a[1] - c[1]
            Dy_Cy = d[1] - c[1]
            Ax_Cx = a[0] - c[0]
            Bx_Ax = b[0] - a[0]
            By_Ay = b[1] - a[1]

            Cx_Bx = c[0] - b[0]
            Dx_Ax = d[0] - a[0]
            Cy_By = c[1] - b[1]
            Dy_Ay = d[1] - a[1]

            numerator_a = Dx_Cx * Ay_Cy - Dy_Cy * Ax_Cx
            numerator_b = Bx_Ax * Ay_Cy - By_Ay * Ax_Cx
            denominator = Dy_Cy * Bx_Ax - Dx_Cx * By_Ay

            if (check_collinear and numerator_a == 0 and numerator_b == 0 and denominator == 0):
                # Lines are collinear

                # if the signs are different the
                # segments have some overlap
                overlap_x = (Cx_Bx < 0) != (Dx_Ax < 0)
                overlap_y = (Cy_By < 0) != (Dy_Ay < 0)

                point_overlap = (a == d or b == c)

                # if any are True there is an intersection
                if (overlap_x or overlap_y or point_overlap):
                    # There is overlap
                    x = sorted((a[0], b[0], c[0], d[0]))
                    y = sorted((a[1], b[1], c[1], d[1]))

                    # min
                    ixa = x[1]
                    iya = y[1]

                    # max
                    ixb = x[2]
                    iyb = y[2]

                    return True, ixa, iya, ixb, iyb
                return False, None, None, None, None

            if denominator == 0:
                # Lines are parallel or
                # maybe collinear if check_collinear=False

                return False, None, None, None, None

            u_a = numerator_a / denominator
            u_b = numerator_b / denominator

            if (u_a >= 0) and (u_a <= 1) and (u_b >= 0) and (u_b <= 1):
                # segments intersect
                ixa = a[0] + Bx_Ax * u_a
                iya = a[1] + By_Ay * u_a

                ixb = c[0] + Dx_Cx * u_b
                iyb = c[1] + Dy_Cy * u_b

                return True, ixa, iyb, ixa, iyb
            return False, None, None, None, None


        def get_topo_angles(nodeDict, block_extent, block_samples, z_raster, azimuthdisdict, searchDistance_max_m, con_z_to_m):
            """This gets the maximum topographic angle and other information for
            each topo line within the block. The data is saved to the nodeDict
            as a list."""

            nodata_to_value = -9999 / con_z_to_m

            x_cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEX").getOutput(0))
            y_cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEY").getOutput(0))

            # localize the block extent values
            block_x_min = block_extent[0]
            block_y_min = block_extent[1]
            block_x_max = block_extent[2]
            block_y_max = block_extent[3]

            # Get the coordinates extent of the input raster
            # this could be in main so it doesn't have to run each time
            raster_x_min = float(arcpy.GetRasterProperties_management(z_raster, "LEFT").getOutput(0))
            raster_y_min = float(arcpy.GetRasterProperties_management(z_raster, "BOTTOM").getOutput(0))
            raster_x_max = float(arcpy.GetRasterProperties_management(z_raster, "RIGHT").getOutput(0))
            raster_y_max = float(arcpy.GetRasterProperties_management(z_raster, "TOP").getOutput(0))

            # Calculate the block x and y offset from the raster and adjust
            # the block coordinates so they are at the raster cell corners.
            # This is for ESRI's RastertoNumpyArray function which defaults to the adjacent
            # lower left cell
            block_x_min_corner = block_x_min - ((block_x_min - raster_x_min) % x_cellsize)
            block_y_min_corner = block_y_min - ((block_y_min - raster_y_min) % y_cellsize)
            block_x_max_corner = block_x_max + ((raster_x_max - block_x_max) % x_cellsize)
            block_y_max_corner = block_y_max + ((raster_y_max - block_y_max) % y_cellsize)

            # calculate the number of cols/ros from the lower left
            ncols = int((block_x_max_corner - block_x_min_corner) / x_cellsize)
            nrows = int((block_y_max_corner - block_y_min_corner) / y_cellsize)

            # Construct the array. Note returned array is (row, col) so (y, x)
            try:
                z_array = arcpy.RasterToNumPyArray(z_raster, arcpy.Point(block_x_min_corner, block_y_min_corner),
                                                   ncols, nrows, nodata_to_value)
            except:
                tbinfo = traceback.format_exc()
                pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"
                sys.exit(pymsg)

                # convert array values to meters if needed
            z_array = z_array * con_z_to_m

            topo_samples = []
            if z_array.max() > -9999:
                # There is at least one pixel of data
                for (nodeID, streamID, a, z_node,
                     node_x, node_y, end_x, end_y,
                     block_search_start, block_search_end) in block_samples:

                    # create list of distance movements along the topo line
                    # from the block edges
                    # this is in units of the fc
                    cellsize = azimuthdisdict[a]
                    distance_array = build_search_array(block_search_start,
                                                        block_search_end,
                                                        cellsize,
                                                        use_skippy=False)
                    z_topo_list = []
                    for distance in distance_array:
                        pt_x = ((distance * sin(radians(a))) + node_x)
                        pt_y = ((distance * cos(radians(a))) + node_y)

                        xy = coord_to_array(pt_x, pt_y, block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize)
                        z_topo_list.append(z_array[xy[1], xy[0]])

                    # convert distances to meters
                    distance_array_m = distance_array * con_to_m
                    z_topo_array = np.array(z_topo_list)
                    # Calculate the topo angles along the topo line
                    angle_array = np.degrees(np.arctan((z_topo_array - z_node) / distance_array_m))
                    # remove the off raster samples
                    naindex = np.where(z_topo_array < -9998)
                    for x in naindex[0]: angle_array[x] = -9999
                    
                    # Find the max topo angle (in degrees)
                    topoAngle = angle_array.max()
##
                    
                    # array index at the max topo angle
                    arryindex = np.where(angle_array == topoAngle)[0][0]
                    z_topo = z_topo_array[arryindex]
                    # elevation change between topo angle location and node elevation
                    z_change = z_topo - z_node
                    # distance from the node to topo angle location in units of fc
                    topoAngleDistance = distance_array[arryindex]
                    topoAngle_x = (topoAngleDistance * sin(radians(a))) + node_x
                    topoAngle_y = (topoAngleDistance * cos(radians(a))) + node_y
                    topoAngleDistance_m = topoAngleDistance * con_to_m
                    off_rastersamples = (z_topo_array < -9998).sum()

                    topo_samples.append([topoAngle_x, topoAngle_y,
                                         topoAngle_x, topoAngle_y,
                                         streamID, nodeID, a,
                                         topoAngle, z_topo,
                                         z_node, z_change,
                                         topoAngleDistance_m,
                                         searchDistance_max_m,
                                         off_rastersamples])

            return topo_samples


        # enable garbage collection
        gc.enable()

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        topo_directions = int(parameters[1].value)
        searchDistance_max_km = float(parameters[2].value)
        z_raster = parameters[3].valueAsText
        z_units = parameters[4].valueAsText
        topo_fc = parameters[5].valueAsText
        block_size = float(parameters[6].value)
        overwrite_data = bool(parameters[7].value)

        try:
            arcpy.AddMessage("Step 4: Measure Topographic Angles")
            print("Step 4: Measure Topographic Angles")

            # keeping track of time
            startTime = time.time()

            # Check if the node fc exists
            if not arcpy.Exists(nodes_fc):
                arcpy.AddError("\nThis output does not exist: \n" +
                               "{0}\n".format(nodes_fc))
                sys.exit("This output does not exist: \n" +
                         "{0}\n".format(nodes_fc))

                # Check if the topo output exists and delete if needed
            if arcpy.Exists(topo_fc) and overwrite_data:
                arcpy.Delete_management(topo_fc)

            # Check if the block fc exists and delete or throw an error
            if arcpy.Exists(block_fc):
                if overwrite_data:
                    arcpy.Delete_management(block_fc)
                else:
                    arcpy.AddMessage("\nThis output already exists: \n" +
                                     "{0}\n".format(block_fc) +
                                     "overwrite data = False. New data will be " +
                                     "appended to the existing feature class.")
                    print("This output already exists: \n" +
                          "{0}\n".format(block_fc) +
                          "overwrite data = False. New data will be " +
                          "appended to the existing feature class.")

            if overwrite_data:
                env.overwriteOutput = True
            else:
                env.overwriteOutput = False

                # Determine input point spatial units
            proj = arcpy.Describe(nodes_fc).spatialReference
            proj_ele = arcpy.Describe(z_raster).spatialReference

            # Check to make sure the raster and input
            # points are in the same projection.
            if proj.name != proj_ele.name:
                arcpy.AddError("{0} and {1} do not have ".format(nodes_fc, z_raster) +
                               "the same projection. Please reproject your data.")
                sys.exit("Input points and elevation raster do not have the " +
                         "same projection. Please reproject your data.")

            # Get the units conversion factors
            con_z_to_m = from_z_units_to_meters_con(z_units)
            con_to_m = to_meters_con(nodes_fc)
            con_from_m = from_meters_con(nodes_fc)
            searchDistance_max_m = searchDistance_max_km * 1000  # in meters

            # search distance in units of the fc
            searchDistance_max = int(con_from_m * searchDistance_max_m)

            # convert block size from km to meters to units of the node fc
            # in the future block size should be estimated based on available memory
            # memorysize = datatypeinbytes*nobands*block_size^2
            # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
            if block_size in ["#", ""]:
                block_size = int(con_from_m * 5000)
            else:
                block_size = int(con_from_m * block_size * 1000)

                # Get the elevation raster cell size in units of the raster
            x_cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEX").getOutput(0))
            y_cellsize = float(arcpy.GetRasterProperties_management(z_raster, "CELLSIZEY").getOutput(0))

            if topo_directions == 1:
                azimuths = [270, 180, 90]
                last_azimuth = 90

                azimuthdict = {90: "TOPO_E", 180: "TOPO_S", 270: "TOPO_W"}

            elif topo_directions == 2:  # All directions
#                azimuths = [45, 90, 135, 180, 225, 270, 315, 365]   #!MSF- Line commented out. North angle of 365 is incorrect, should be 360 or 0. 360 was erroring out so used 0 instead to point North
                azimuths = [45, 90, 135, 180, 225, 270, 315, 0]
                last_azimuth = 45

                # azimuthdict = {45: "TOPO_NE", 90: "TOPO_E", 135: "TOPO_SE",
                #                180: "TOPO_S", 225: "TOPO_SW", 270: "TOPO_W",
                #                315: "TOPO_NW", 365: "TOPO_N"}                     #!MSF- Line commented out. North angle of 365 is incorrect and was causing model to crash. Updated to 360 or 0. 360 was erroring out so used 0 instead to point North
                azimuthdict = {45: "TOPO_NE", 90: "TOPO_E", 135: "TOPO_SE",
                               180: "TOPO_S", 225: "TOPO_SW", 270: "TOPO_W",
                               315: "TOPO_NW", 0: "TOPO_N"}                     

            elif topo_directions == 3:
                azimuths = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340]
                last_azimuth = 20

                azimuthdict = {0: "TOPO01", 20: "TOPO02", 40: "TOPO03", 60: "TOPO04", 80: "TOPO05", 100: "TOPO06", 120: "TOPO07",
                               140: "TOPO08", 160: "TOPO09", 180: "TOPO10", 200: "TOPO11", 220: "TOPO12", 240: "TOPO13",
                               260: "TOPO14", 280: "TOPO15", 300: "TOPO16", 320: "TOPO17", 340: "TOPO18"}

            azimuthdisdict = dict.fromkeys(azimuths, None)

            for a in azimuthdisdict:
                if a <= 45:
                    azimuthdisdict[a] = abs(y_cellsize / cos(radians(a)))
                elif 45 < a <= 135:
                    azimuthdisdict[a] = abs(x_cellsize / sin(radians(a)))
                elif 135 < a <= 225:
                    azimuthdisdict[a] = abs(y_cellsize / cos(radians(a)))
                elif 225 < a <= 315:
                    azimuthdisdict[a] = abs(x_cellsize / sin(radians(a)))
                elif 315 < a <= 360:
                    azimuthdisdict[a] = abs(y_cellsize / cos(radians(a)))

            # Build the topo field names
            addFields = [azimuthdict[a] for a in azimuths]

            # Get a list of existing fields
            existingFields = []
            for f in arcpy.ListFields(nodes_fc):
                existingFields.append(f.name)

                # Check to see if the field exists and add it if not
            for f in addFields:
                if (f in existingFields) is False:
                    arcpy.AddField_management(nodes_fc, f, "DOUBLE", "", "", "",
                                              "", "NULLABLE", "NON_REQUIRED")

            # Read the feature class data into a nested dictionary
            nodeDict = read_nodes_fc(nodes_fc, overwrite_data, addFields)

            # Build the blockDict
            blockDict = create_blocks(nodeDict, block_size, last_azimuth,
                                      searchDistance_max)

            # Iterate through each block
            total_samples = 0
            blockIDs = list(blockDict.keys())
            blockIDs.sort()

            for p, blockID in enumerate(blockIDs):
                block_extent = blockDict[blockID]["extent"]
                block_samples = blockDict[blockID]["samples"]
                block_samples.sort()

                arcpy.AddMessage("Processing block {0} of {1}".format(p + 1, len(blockIDs)))
                print("Processing block {0} of {1}".format(p + 1, len(blockIDs)))

                # calculate coordinates along the
                # portion of the topo line in the block,
                # convert raster to array, sample the raster
                # calculate the topo angles and other info
                topo_samples = get_topo_angles(nodeDict, block_extent, block_samples,
                                               z_raster, azimuthdisdict,
                                               searchDistance_max, con_z_to_m)
                if topo_samples:
                    # Update the nodeDict
                    for sample in topo_samples:
                        nodeID = sample[5]
                        a = sample[6]
                        topoAngle = sample[7]

                        # Create a key to hold the topo list info for this block
                        topo_key = azimuthdict[a] + "_list"

                        if azimuthdict[a] in nodeDict[nodeID]:
                            if nodeDict[nodeID][azimuthdict[a]] < topoAngle:
                                nodeDict[nodeID][azimuthdict[a]] = topoAngle
                                nodeDict[nodeID][topo_key] = sample

                        else:
                            nodeDict[nodeID][azimuthdict[a]] = topoAngle
                            nodeDict[nodeID][topo_key] = sample

                    del topo_samples

                # Check if any nodes can be updated in the node and topo fc
                if blockDict[blockID]["nodes_to_update"]:
                    nodes_to_update = blockDict[blockID]["nodes_to_update"]

                    # Write the topo data to the TTools point feature class
                    update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_update)

                    # Build/add to the output topo feature class
                    topo_list = []
                    for nodeID in nodes_to_update:
                        for field in addFields:
                            topo_key = field + "_list"
                            topo_list.append(nodeDict[nodeID][topo_key])
                            # delete some data
                            nodeDict[nodeID].pop(field, None)
                            nodeDict[nodeID].pop(topo_key, None)
                            nodeDict[nodeID].pop("updated", None)
                    update_topo_fc(topo_list, topo_fc, nodes_fc,
                                   nodes_to_update, overwrite_data, proj)

                    del topo_list
                    gc.collect()

            endTime = time.time()
            elapsedmin = ceil(((endTime - startTime) / 60) * 10) / 10
            mspernode = timedelta(seconds=(endTime - startTime) / len(nodeDict.keys())).microseconds
            print("Process Complete in {0} minutes. {1} microseconds per node".format(elapsedmin, mspernode))
            arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per node" % (elapsedmin, mspernode))

            # === BEGIN TTOOLS Edit (MSF - Sep 14, 2025): Write radians alongside degrees for CE-QUAL-W2 (option 3) ===
            # This post-processing keeps the original degree fields (e.g., TOPO01) and
            # creates/updates corresponding radian fields (e.g., TOPO01_RAD) when option 3 is selected.
            # No changes to internal calculations; this only writes additional fields to the attribute table.
            if topo_directions == 3:
                try:
                    arcpy.AddMessage("Creating/updating *_RAD fields for CE-QUAL-W2 (radians).")
                    # Cache existing field names once
                    existing_fields = [f.name for f in arcpy.ListFields(nodes_fc)]
                    for field in addFields:
                        rad_field = f"{field}_RAD"
                        if rad_field not in existing_fields:
                            arcpy.AddField_management(nodes_fc, rad_field, "DOUBLE", "", "", "",
                                                      "", "NULLABLE", "NON_REQUIRED")
                            existing_fields.append(rad_field)

                        # Populate radians = degrees * pi / 180
                        with arcpy.da.UpdateCursor(nodes_fc, [field, rad_field]) as cursor:
                            for row in cursor:
                                deg_val = row[0]
                                if deg_val is not None:
                                    row[1] = deg_val * (np.pi / 180.0)
                                else:
                                    # Preserve existing rad value if any, or leave as None
                                    pass
                                cursor.updateRow(row)
                    arcpy.AddMessage("Finished writing *_RAD fields.")
                except Exception as _rad_ex:
                    # Non-fatal: radians fields are extra outputs
                    arcpy.AddWarning(f"Failed to write radians fields: {_rad_ex}")
            # === END TTOOLS EDIT (MSF - Sep 14, 2025) ===


        # For arctool errors
        except arcpy.ExecuteError:
            msgs = arcpy.GetMessages(2)
            arcpy.AddError(msgs)
            print(msgs)

        # For other errors
        except:
            tbinfo = traceback.format_exc()

            pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)

            print(pymsg)
            print(msgs)


class SampleLandcoverStartPattern(object):
    def __init__(self):
        self.label = "Step 5A: Sample Landcover - Star Pattern"
        self.description = "Calculate the maximum topographic shade angle for each stream node."
        self.canRunInBackground = False

    def getParameterInfo(self):
    
        param0 = arcpy.Parameter(
            displayName="Stream Nodes Feature Class",
            name="nodes_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input") # -MSF used DEFeatureClass instead of GPFeatureLayer to explicitly capture path. Made global update to use DE just in case.
    
        param1 = arcpy.Parameter(
            displayName="Number of Transects per Node",
            name="trans_count",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
    
        param2 = arcpy.Parameter(
            displayName="Number of Samples per Transect",
            name="transsample_count",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
    
        param3 = arcpy.Parameter(
            displayName="Distance between samples (m)",
            name="transsample_distance",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
 # Changed caption from "Sample Zone Only" to "Use Zone Sampling Method" - MSF 9/29/2025 
        param4 = arcpy.Parameter(
            displayName="Use Zone Sampling Method",
            name="zone_sample",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param4.value = False
    
        param5 = arcpy.Parameter(
            displayName="HeatSource 8 Output Format",
            name="heatsource8",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param5.value = False
        
 # Made param6 invisible from the user interface by changing the parametertype to "Derived". This was a special custom feature added and is not used   - MSF 9/29/2025    
        param6 = arcpy.Parameter(
            displayName="Add SampleID Field for Code Lookup",
            name="sampleID_for_code",
            datatype="GPBoolean",
            parameterType="Derived",
            direction="Input")
        param6.value = False
    
        param7 = arcpy.Parameter(
            displayName="Landcover Raster",
            name="lc_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
    
        param8 = arcpy.Parameter(
            displayName="Landcover Units",
            name="lc_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param8.filter.type = "ValueList"
        param8.filter.list = ["Meters", "Feet", "None"]
# param9 through param12 are Canopy related optional parameters    
        param9 = arcpy.Parameter(
            displayName="Canopy Cover Type",
            name="canopy_data",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")
        param9.filter.type = "ValueList"
        param9.filter.list = ["CanopyCover", "LAI", "None"]
        param9.value = "None"
  
        param10 = arcpy.Parameter(
            displayName="Canopy Raster",
            name="canopy_raster",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
    
        param11 = arcpy.Parameter(
            displayName="K Raster",
            name="k_raster",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
    
        param12 = arcpy.Parameter(
            displayName="Overhanging Vegetation Raster",
            name="oh_raster",
            datatype="GPRasterLayer",
            parameterType="Optional",
            direction="Input")
    
        param13 = arcpy.Parameter(
            displayName="Elevation Raster",
            name="z_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")
    
        param14 = arcpy.Parameter(
            displayName="Elevation Units",
            name="z_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param14.filter.type = "ValueList"
        param14.filter.list = ["Meters", "Feet"]
    
        param15 = arcpy.Parameter(
            displayName="Output Landcover Point Feature Class",
            name="lc_point_fc",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")
    
        param16 = arcpy.Parameter(
            displayName="Raster Block Size (km)",
            name="block_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
    
        param17 = arcpy.Parameter(
            displayName="Overwrite Existing Data",
            name="overwrite_data",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param17.value = True
            
    
        return [param0, param1, param2, param3, param4, param5, param6, param7,
                param8, param9, param10, param11, param12, param13, param14,
                param15, param16, param17]
            
    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
    # Log input parameters to file
        log_run(self.label, parameters)
        #!/usr/bin/python

        """
        TTools Step 5: Sample Landcover - Star Pattern

        The star pattern sampling method will take an input point feature (from Step 1) and sample input landcover rasters
        along transects in any number of azimuth directions oriented outward from the stream node. The numer of transects,
        sample points along each transect, and the distance between samples is user defined. The measured distance along each
        transect starts at the stream node.

        The star pattern sampling method can be used to develop landcover inputs for Heat Source 7 - 9.

        REQUIREMENTS
        TTools steps 1 - 3 must be run before Step 5.
        ESRI ArcGIS Pro w/ Spatial Analyst extension
        Python 3.7+

        INPUT VARIABLES
        0: nodes_fc:
        Path to the TTools point feature class.

        1: trans_count:
        Number of transects per node. The degrees of separation between each transect is equal to 360 / trans_count. If
        trans_count = 8, the heading of the first transect is 45 degrees (Northeast), the second transect is 90 degrees
        (West), and the eighth is 360 degrees (North). Transects are oriented clockwise from north. trans_count is ignored if
        heatsource8 = True.

        2: transsample_count:
        Number of samples per transect. This number DOES NOT include the sample at the stream node.

        3: transsample_distance:
        The distance between transect samples measured in meters.

        4: zone_sample:
        Boolean (True/False) flag to indicate if the sample should represent a zone and be centered relative to
        the transsample_distance as measured from the stream node. If this is True the distance from the stream
        node to the sample location along each transect is equal to the transsample_distance * (sample number - 0.5).
        This should be True if using heat source 7, heat source 8.0.1 - 8.0.5, or heat source 8.0.7 -  8.0.8 when the model is
        configured to use the zone method. This should be False if using heat source 9 or when using heat source 8.0.7 -  8.0.8
        and the model is configured to use the point method.

        5. heatsource8:
        Boolean (True/False) flag to indicate if the star pattern with 7 transects should be used. Heat source
        version 7 and 8 use 7 transects around each stream node in the following directions: Northeast, East, Southeast, South,
        Southwest, West, Northwest.

        6: sampleID_for_code
        Boolean (True/False) flag to indicate if the sample id should be used as the landcover value. If True the lc_raster
        will be ignored and only the z_raster will be sampled. Generally this should be set to False.

        7: lc_raster:
        Path and name of the land cover code, height, or elevation raster.

        8: lc_units:
        z units of the lc_raster (aka units of height or elevation). Use "Feet", "Meters", or "None" if the lc_raster values are
        codes and do not represent elevation or height units.

        9: canopy_data:
        The input canopy data type being sampled in canopy_raster. Can be one of "CanopyCover", "LAI", or "None" if not sampling
        a canopy raster. If canopy_data = None, the variables canopy_raster, k_raster, and oh_raster are ignored.

        10: canopy_raster
        Path and name of the canopy or LAI raster. Ignored if canopy_data = "None".

        11: k_raster
        Path and name of the k coefficient raster. There must be a k_raster if canopy_data = "LAI". Ignored if
        canopy_data = "None".

        12: oh_raster
        Path and name of the vegetation overhang raster. Ignored if canopy_data = "None".

        13: z_raster:
        Path and name of the ground elevation raster.

        14: z_units:
        z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is not in feet or meters the elevation
        values must be converted.

        15: lc_point_fc:
        Path and name of output sample point feature file.

        16: block_size:
        The x and y size in kilometers for each raster block pulled into an array. Start with 5 if you aren't sure and reduce
        if there is an error. To increase processing speed rasters are subdivided iteratively into smaller blocks and pulled
        into arrays for processing. Very large block sizes may use a lot of memory and result in an error.

        17: overwrite_data:
        True/False flag if existing data in nodes_fc and lc_point_fc can be overwritten.

        OUTPUTS
        0. nodes_fc:
        New fields are added into nodes_fc with the Landcover and elevation values for each transect sample

        1. lc_point_fc:
        New point feature class created with a point at each x/y sample and the sample raster values

        """
        # Import system modules
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

        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------

        # General scripts steps include:
        # 1. open nodes fc. iterate and read info from node fc into a dict

        # 2. create a list with the x/y and related info for each lc sample
        #    calculate the extent bounding box for the entire dataset

        # 3. create a list holding the bounding box coords for each block iteration

        # 4. loop through each block
            #- sample the raster for all samples in the block
            #- update the node dict
            #- update the node fc
            #- update the lc sample fc
            #- continue to the next block

        # Future Updates
        # -Change the node dict so the node is the primary key
        # -Build the block list based on the nodes and then build point list iteratively
        # instead of building them into one huge list. The huge list results
        # in a memory error for large areas
        # -Include stream sample in transect count (True/False)
        # -Eliminate arcpy and use gdal for reading/writing feature class data

        env.overwriteOutput = True

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
            # Grabs last field because often the first field, emergent, is zero
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

        def update_lc_point_fc(lc_point_list, type, lc_point_fc, nodes_fc,
                               nodes_in_block, overwrite_data, proj):
            """Creates/updates the output landcover sample point feature
            class using the data from the landcover point list"""

            cursorfields = ["POINT_X","POINT_Y"] + ["STREAM_ID","NODE_ID", "SAMPLE_ID",
                                                    "TRANS_DIR", "TRANSECT",
                                                    "SAMPLE", "KEY"] + type
    
            # Check if the output exists and create if not
            if not arcpy.Exists(lc_point_fc):
                arcpy.CreateFeatureclass_management(os.path.dirname(lc_point_fc),
                                                    os.path.basename(lc_point_fc),
                                                    "POINT","","DISABLED","DISABLED",proj)
        
                # Determine Stream ID field properties
                sid_type = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].type
                sid_precision = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].precision
                sid_scale = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].scale
                sid_length = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].length
        
                typeDict = {"POINT_X": "DOUBLE", 
                            "POINT_Y": "DOUBLE", 
                            "NODE_ID": "LONG",
                            "SAMPLE_ID": "LONG",
                            "TRANS_DIR": "DOUBLE",
                            "TRANSECT": "SHORT",
                            "SAMPLE": "SHORT",
                            "KEY": "TEXT",}
                
                for t in type:
                    typeDict[t] = "DOUBLE"
    
                # Add attribute fields
                for f in cursorfields:
                    if f == "STREAM_ID":
                        arcpy.AddField_management(lc_point_fc, f, sid_type,
                                                  sid_precision, sid_scale, sid_length,
                                                  "", "NULLABLE", "NON_REQUIRED")
                    else:
                        arcpy.AddField_management(lc_point_fc, f, typeDict[f], "",
                                                  "", "", "", "NULLABLE", "NON_REQUIRED")                 
    
            if not overwrite_data:
                # Build a query to retrieve existing rows from the nodes
                # that need updating
                whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_in_block))        
        
                # delete those rows
                with arcpy.da.UpdateCursor(lc_point_fc,["NODE_ID"], whereclause) as cursor:  
                    for row in cursor:
                        cursor.deleteRow()     

            with arcpy.da.InsertCursor(lc_point_fc, ["SHAPE@X","SHAPE@Y"] +
                                       cursorfields) as cursor:
                for row in lc_point_list:
                    cursor.insertRow(row)

        def setup_lcdata_headers(transsample_count, trans_count, canopy_data):
            """Generates a list of the landcover data file
            column header names and data types"""
    
            type = ["ELE"]
     
            #Use LAI methods   
            if canopy_data == "LAI":
                type = type + ["LAI","k","OH"]
    
            #Use Canopy Cover methods    
            if canopy_data == "CanopyCover":
                type = type + ["CAN","OH"]
    
            lcheaders = []
            otherheaders = []
     
            dirs = ["T{0}".format(x) for x in range(1, trans_count + 1)]

            zones = range(1,int(transsample_count)+1)
    
            # Concatenate the type, dir, and zone and order in the correct way
    
            for d, dir in enumerate(dirs):
                for z, zone in enumerate(zones):
                    if d==0 and z==0:
                        lcheaders.append("LC_T0_S0") # add emergent
                        lcheaders.append("LC_{0}_S{1}".format(dir, zone))
                    else:
                        lcheaders.append("LC_{0}_S{1}".format(dir, zone))
    
            for t in type:
                for d, dir in enumerate(dirs):
                    for z, zone in enumerate(zones):
                        if t !="ELE" and d==0 and z==0:
                            otherheaders.append(t+"_T0_S0") # add emergent
                            otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))
                        else:
                            otherheaders.append("{0}_{1}_S{2}".format(t, dir, zone))    

            return lcheaders, otherheaders

        def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
            """converts x/y coordinates to col and row of the array"""
            col_x = int(((easting - block_x_min) - ((easting - block_x_min) % x_cellsize)) / x_cellsize)
            row_y = int(((block_y_max - northing) - ((block_y_max - northing) % y_cellsize)) / y_cellsize)
    
            return [col_x, row_y]

        def create_lc_point_list(nodeDict, nodes_in_block, dirs, zones, transsample_distance, zone_sample):
            """This builds a unique long form list of information for all the
            landcover samples in the block. This list is used to
            create/update the output feature class."""
    
            lc_point_list = []
            numDirs = len(dirs)
            numZones = len(zones)
            zonesPerNode = (numDirs * numZones) + 1

            if zone_sample:
                adjust = 0.5
            else:
                adjust = 0.0

            for nodeID in nodes_in_block:
                origin_x = nodeDict[nodeID]["POINT_X"]
                origin_y = nodeDict[nodeID]["POINT_Y"]
                streamID = nodeDict[nodeID]["STREAM_ID"]
                sampleID = nodeID * zonesPerNode
        
                # This is the emergent/stream sample
                lc_point_list.append([origin_x, origin_y, origin_x, origin_y,
                                      streamID, nodeID, sampleID,
                                      0, 0, 0, "T0_S0"])
    
                for d, dir in enumerate(dirs):
                    for zone in zones:
                        # Calculate the x and y coordinate of the 
                        # landcover sample location
                        pt_x = ((zone - adjust) * transsample_distance * con_from_m *
                                sin(radians(dir))) + origin_x
                        pt_y = ((zone - adjust) * transsample_distance * con_from_m *
                                cos(radians(dir))) + origin_y
                
                        key = 'T{0}_S{1}'.format(d+1, zone)
                        tran_num = '{:{}{}}'.format(d+1, 0, 3)
                        samp_num = '{:{}{}}'.format(zone, 0, 2)
                        sampleID = (nodeID * zonesPerNode) + (d * numZones) + zone
        
                        # Add to the list          
                        lc_point_list.append([pt_x, pt_y, pt_x, pt_y,
                                              streamID, nodeID, sampleID,
                                              dir, d+1, zone, key])    
     
            return lc_point_list

        def create_block_list(nodes, block_size):
            """Returns two lists, one containing the coordinate extent
            for each block that will be iteratively extracted to an array
            and the other containing node IDs within each block extent."""
    
            arcpy.AddMessage("Calculating block extents")    
            print("Calculating block extents")    
            x_coord_list = [nodeDict[nodeID]["POINT_X"] for nodeID in nodes]
            y_coord_list = [nodeDict[nodeID]["POINT_Y"] for nodeID in nodes]
    
            # calculate the buffer distance (in raster spatial units) to add to 
            # the base bounding box when extracting to an array. The buffer is 
            # equal to the sample distance + 1 to make sure the block includes 
            # all the landcover samples for each node.
            buffer = int((transsample_count + 1) * transsample_distance * con_from_m)
    
            # calculate bounding box extent for samples
            x_min = min(x_coord_list)
            x_max = max(x_coord_list)
            y_min = min(y_coord_list) 
            y_max = max(y_coord_list)
    
            x_width = int(x_max - x_min + 1)
            y_width = int(y_max - y_min + 1)
    
            block_extents = []
            block_nodes = []
      
            # Build data blocks
            for x in range(0, x_width, block_size):
                for y in range(0, y_width, block_size):

                    # Lower left coordinate of block (in map units)
                    block0_x_min = min([x_min + x, x_max])
                    block0_y_min = min([y_min + y, y_max])
                    # Upper right coordinate of block (in map units)
                    block0_x_max = min([block0_x_min + block_size, x_max])
                    block0_y_max = min([block0_y_min + block_size, y_max])
            
                    block_x_min = block0_x_max
                    block_x_max = block0_x_min
                    block_y_min = block0_y_max
                    block_y_max = block0_y_min
            
                    nodes_in_block = []
                    for nodeID in nodes:
                        node_x = nodeDict[nodeID]["POINT_X"]
                        node_y = nodeDict[nodeID]["POINT_Y"]
                        if (block0_x_min <= node_x <= block0_x_max and
                            block0_y_min <= node_y <= block0_y_max):
                    
                            nodes_in_block.append(nodeID)
                    
                            # Minimize the size of the block0 by the true 
                            # extent of the nodes in the block
                            if block_x_min > node_x: block_x_min = node_x
                            if block_x_max < node_x: block_x_max = node_x
                            if block_y_min > node_y: block_y_min = node_y
                            if block_y_max < node_y: block_y_max = node_y
            
                    if nodes_in_block:
                        # add the block extent for processing
                        # order 0 left,      1 bottom,    2 right,     3 top
                        block_extents.append((block_x_min - buffer, block_y_min - buffer,
                                              block_x_max + buffer, block_y_max + buffer))           
                        block_nodes.append(nodes_in_block)
    
            return block_extents, block_nodes
    
        def sample_raster(block, lc_point_list, raster, con):
#-----MSF - 9/15/2025
            # Normalize raster input so dataset path is always used (works for GPRasterLayer or dataset)
            try:
                raster_src = arcpy.Describe(raster).catalogPath
            except:
                raster_src = raster
#-------------------    
            if con is not None:
                nodata_to_value = -9999 / con
            else:
                nodata_to_value = -9999
        
            # localize the block extent values
            block_x_min = block[0]
            block_y_min = block[1]
            block_x_max = block[2]
            block_y_max = block[3]
# updated raster to raster_src - MSF 9/15/2025    
            x_cellsize = float(arcpy.GetRasterProperties_management(raster_src, "CELLSIZEX").getOutput(0))
            y_cellsize = float(arcpy.GetRasterProperties_management(raster_src, "CELLSIZEY").getOutput(0)) 

            # Get the coordinate extent of the input raster
            raster_x_min = float(arcpy.GetRasterProperties_management(raster_src, "LEFT").getOutput(0))
            raster_y_min = float(arcpy.GetRasterProperties_management(raster_src, "BOTTOM").getOutput(0))
            raster_x_max = float(arcpy.GetRasterProperties_management(raster_src, "RIGHT").getOutput(0))
            raster_y_max = float(arcpy.GetRasterProperties_management(raster_src, "TOP").getOutput(0))

            # Calculate the block x and y offset from the raster and adjust
            # the block coordinates so they are at the raster cell corners.
            # This is for ESRI's RastertoNumpyArray function which defaults to the adjacent
            # lower left cell
            block_x_min_corner = block_x_min - ((block_x_min - raster_x_min) % x_cellsize)
            block_y_min_corner = block_y_min - ((block_y_min - raster_y_min) % y_cellsize)
            block_x_max_corner = block_x_max + ((raster_x_max - block_x_max) % x_cellsize)
            block_y_max_corner = block_y_max + ((raster_y_max - block_y_max) % y_cellsize)
    
            # calculate the number of cols/rows from the lower left
            ncols = int((block_x_max_corner - block_x_min_corner) / x_cellsize)
            nrows = int((block_y_max_corner - block_y_min_corner) / y_cellsize)
    
            # Construct the array. Note returned array is (row, col) so (y, x)
            try:
                # updated raster to raster_src - MSF 9/15/2025
                raster_array = arcpy.RasterToNumPyArray(raster_src, arcpy.Point(block_x_min_corner, block_y_min_corner),
                                                        ncols, nrows, nodata_to_value)
            except:
                tbinfo = traceback.format_exc()
                pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
                sys.exit(pymsg)    
    
            # convert array values to meters if needed
            if con is not None:
                raster_array = raster_array * con
    
            lc_point_list_new = []
            if raster_array.max() > -9999:
                # There is at least one pixel of data
                for point in lc_point_list:
                    xy = coord_to_array(point[0], point[1], block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize)
                    point.append(raster_array[xy[1], xy[0]])
                    lc_point_list_new.append(point)
            else:
                # No data, add -9999
                for point in lc_point_list:
                    point.append(-9999)
                    lc_point_list_new.append(point)
            return lc_point_list_new            

        def update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_query):
            """Updates the input point feature class with data from the
            nodes dictionary"""

            # Build a query to retrieve just the nodes that needs updating
            whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))

            with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:  
                for row in cursor:
                    for f, field in enumerate(addFields):
                        nodeID =row[0]
                        row[f+1] = nodeDict[nodeID][field]
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
            """Returns the conversion factor to get from
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

        #enable garbage collection
        gc.enable()

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        trans_count = int(parameters[1].value)
        transsample_count = int(parameters[2].value)
        transsample_distance = float(parameters[3].value)
        zone_sample = bool(parameters[4].value)
        heatsource8 = bool(parameters[5].value)
        sampleID_for_code = bool(parameters[6].value)
        lc_raster = parameters[7].valueAsText
        lc_units = parameters[8].valueAsText
        canopy_data = parameters[9].valueAsText
        canopy_raster = parameters[10].valueAsText if parameters[10].value else None
        k_raster = parameters[11].valueAsText if parameters[11].value else None
        oh_raster = parameters[12].valueAsText if parameters[12].value else None
        z_raster = parameters[13].valueAsText
        z_units = parameters[14].valueAsText
        lc_point_fc = parameters[15].valueAsText
        block_size = float(parameters[16].value)
        overwrite_data = bool(parameters[17].value)

        try:
    
            arcpy.AddMessage("Step 5: Sample Landcover - Star Pattern, Point Method")
            print("Step 5: Sample Landcover - Star Pattern, Point Method")
    
            #keeping track of time
            startTime= time.time()
    
            # Check if the node fc exists
            if not arcpy.Exists(nodes_fc):
                arcpy.AddError("\nThis output does not exist: \n" +
                               "{0}\n".format(nodes_fc))
                sys.exit("This output does not exist: \n" +
                         "{0}\n".format(nodes_fc))
    
            # Check if the lc point fc exists and delete if needed
            if arcpy.Exists(lc_point_fc) and overwrite_data:
                arcpy.Delete_management(lc_point_fc)
    
    
            # Determine input spatial units and set unit conversion factors
            proj = arcpy.Describe(nodes_fc).SpatialReference
            proj_ele = arcpy.Describe(z_raster).spatialReference
            proj_lc = arcpy.Describe(lc_raster).spatialReference
    
            con_from_m = from_meters_con(nodes_fc)
            con_lc_to_m = from_z_units_to_meters_con(lc_units)
            con_z_to_m = from_z_units_to_meters_con(z_units)
    
            # convert block size from km to meters to units of the node fc
            # in the future block size should be estimated based on available memory
            # memorysize = datatypeinbytes*nobands*block_size^2
            # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
            if block_size in ["#", "", None]:
                block_size = int(con_from_m * 5000)
            else:
                block_size = int(con_from_m * block_size * 1000)
    
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
            if sampleID_for_code:
                rasterDict["LC"] = None

            if canopy_data == "LAI":
                # Use LAI methods
                rasterDict["LAI"] = canopy_raster
                rasterDict["k"] = k_raster
                rasterDict["OH"] = oh_raster

            if canopy_data == "CanopyCover":
                # Use Canopy Cover methods
                rasterDict["CAN"] = canopy_raster
                rasterDict["OH"] = oh_raster

            # flag indicating the model should use the heat source 8 methods 
            # (same as 8 directions but no north)
            if heatsource8:
                dirs = [45,90,135,180,225,270,315]
                trans_count = 7
            else:        
                dirs = [x * 360.0 / trans_count for x in range(1,trans_count+ 1)]

            zones = range(1,int(transsample_count+1))

            lcheaders, otherheaders = setup_lcdata_headers(transsample_count, trans_count, canopy_data)
    
            addFields = lcheaders + otherheaders
    
            # Get a list of existing fields
            existingFields = [f.name for f in arcpy.ListFields(nodes_fc)] 
        
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
            nodes = list(nodeDict.keys())
            nodes.sort()
   
            # Build the block list
            block_extents, block_nodes = create_block_list(nodes, block_size)
    
            # Iterate through each block, calculate sample coordinates,
            # convert raster to array, sample the raster
            total_samples = 0
            for p, block in enumerate(block_extents):
                nodes_in_block = block_nodes[p]
                nodes_in_block.sort()
                arcpy.AddMessage("Processing block {0} of {1}".format(p + 1, len(block_extents)))
                print("Processing block {0} of {1}".format(p + 1, len(block_extents)))
        
                # build the landcover sample list
                lc_point_list = create_lc_point_list(nodeDict, nodes_in_block,
                                                    dirs, zones, transsample_distance, zone_sample)
        
                for t, (type, raster) in enumerate(rasterDict.items()):
                    if raster is None:
                        for i in range(0, len(lc_point_list)):
                            lc_point_list[i].append(-9999)
                    else: 
                        if raster == z_raster:
                            con = con_z_to_m
                        elif raster == lc_raster:
                            con = con_lc_to_m
                        else:
                            con = None  
            
                        lc_point_list = sample_raster(block, lc_point_list, raster, con)
        
                    # Update the node dict
                    if (sampleID_for_code and type == "LC"):
                        for row in lc_point_list:
                            key = "{0}_{1}".format(type, row[10])
                            nodeDict[row[5]][key] = row[6]
            
                    else:
                        for row in lc_point_list:
                            key = "{0}_{1}".format(type, row[10])
                            nodeDict[row[5]][key] = row[11 + t]
        
                # Write the landcover data to the TTools point feature class 
                update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_in_block)
        
                # Build the output point feature class using the data         
                update_lc_point_fc(lc_point_list=lc_point_list, type=list(rasterDict.keys()), lc_point_fc=lc_point_fc,
                                   nodes_fc=nodes_fc, nodes_in_block=nodes_in_block, overwrite_data=overwrite_data, proj=proj)
    
                total_samples = total_samples + len(lc_point_list)
                del lc_point_list
                gc.collect()
   
            endTime = time.time()
    
            elapsedmin = ceil(((endTime - startTime) / 60)* 10)/10
            mspersample = timedelta(seconds=(endTime - startTime) /
                                    total_samples).microseconds
            arcpy.AddMessage("Process Complete in {0} minutes. {1} microseconds per sample".format(elapsedmin, mspersample))
            print("Process Complete in {0} minutes. {1} microseconds per sample".format(elapsedmin, mspersample))
            #arcpy.AddMessage("Process Complete in %s minutes. %s microseconds per sample" % (elapsedmin, mspersample))

        # For arctool errors
        except arcpy.ExecuteError:
            msgs = arcpy.GetMessages(2)
            arcpy.AddError(msgs)
            print(msgs)

        # For other errors
        except:
            tbinfo = traceback.format_exc()

            pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" +str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)

            print(pymsg)
            print(msgs)

import arcpy
import os
import sys

class SampleLandcoverOrthogonalMethod(object):
    def __init__(self):
        self.label = "Step 5B: Sample Landcover - Orthogonal Method"
        self.description = "Sample landcover values along orthogonal transects from stream banks or nodes."
        self.canRunInBackground = False

    def getParameterInfo(self):
        param0 = arcpy.Parameter(
            displayName="Stream Nodes Feature Class",
            name = "nodes_fc",
            datatype ="GPFeatureLayer",
            parameterType="Required",
            direction = "Input")
        param1 = arcpy.Parameter(
            displayName="Start at Stream Bank", 
            name ="start_bank", 
            datatype ="GPBoolean", 
            parameterType="Required", 
            direction ="Input")
        param1.value = True    
        param2 = arcpy.Parameter(
            displayName="Samples per Transect", 
            name = "transsample_count", 
            datatype ="GPLong", 
            parameterType="Required", 
            direction = "Input")
        param3 = arcpy.Parameter(
            displayName="Distance Between Samples (m)", 
            name = "transsample_distance", 
            datatype ="GPDouble", 
            parameterType="Required", 
            direction = "Input")
        param4 = arcpy.Parameter(
            displayName="Landcover Raster", 
            name = "lc_raster", 
            datatype ="GPRasterLayer", 
            parameterType="Required", 
            direction = "Input")
        param5 = arcpy.Parameter(
            displayName="Landcover Units", 
            name = "lc_units", 
            datatype ="GPString", 
            parameterType="Required", 
            direction = "Input")
        param5.filter.type = "ValueList"; 
        param5.filter.list = ["Meters", "Feet", "None"]
        param6 = arcpy.Parameter(
            displayName="Elevation Raster", 
            name = "z_raster", 
            datatype ="GPRasterLayer", 
            parameterType="Required", 
            direction = "Input")
        param7 = arcpy.Parameter(
            displayName="Elevation Units", 
            name = "z_units", 
            datatype ="GPString", 
            parameterType="Required", 
            direction = "Input")
        param7.filter.type = "ValueList";
        param7.filter.list = ["Meters", "Feet"]
        param8 = arcpy.Parameter(
            displayName="Output Landcover Points", 
            name = "lc_point_fc", 
            datatype ="DEFeatureClass", 
            parameterType="Required", 
            direction = "Output")
        param9 = arcpy.Parameter(
            displayName="Raster Block Size (km)", 
            name = "block_size", 
            datatype ="GPDouble", 
            parameterType="Required", 
            direction = "Input")
        param10 = arcpy.Parameter(
            displayName="Overwrite Existing Data", 
            name = "overwrite_data", 
            datatype ="GPBoolean", 
            parameterType="Required", 
            direction = "Input")
        param10.value = True    
           

        return [param0, param1, param2, param3, param4, param5, param6,
                param7, param8, param9, param10]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
    # Log input parameters to file
        log_run(self.label, parameters)
        #!/usr/bin/python

        """
        TTools Step 5: Sample Landcover - Orthogonal Method

        The orthogonal sampling method will take an input point
        feature (from Step 1) and sample an input landcover raster along two
        transects that are orthogonal to the stream aspect (left and right looking downstream).
        The numer of sample points along each transect and distance between samples
        is user defined. The transect can start at the stream node or at
        the right and left banks.

        The Orthogonal sampling method can be used to develop inputs for Heat Source 6,
        Washington Department of Ecology's Shade model, and the shade file for CE-QUAL-W2.

        REQUIREMENTS
        TTools steps 1 - 3 must be run before Step 5.
        ESRI ArcGIS Pro w/ Spatial Analyst extension
        Python 3.7+

        INPUT VARIABLES
        0: nodes_fc:
        path to the TTools point feature class.

        1: start_bank:
        boolean (True/False) to indicate if the transect should start at the stream bank. If False the transect will start at
        the stream node. Normally set to True.

        2: transsample_count:
        Number of samples per transect. The number DOES NOT include the sample at the stream node.

        3: transsample_distance:
        The distance between transect samples (meters).

        5: lc_raster:
        Path and name of the land cover code, height, or elevation raster.

        6: lc_units:
        z units of the lc_raster (aka units of height or elevation). Use "Feet", "Meters", or None if the lc_raster values are
        codes and do not represent elevation or height units.

        7: z_raster:
        Path and name of the ground elevation raster.

        8: z_units:
        z_raster ground elevation units. Either "Feet" or "Meters". If the z unit is not in feet or meters the elevation values
        must be converted.

        9: lc_point_fc:
        Path and name of output sample point feature file.

        10: block_size:
        The x and y size in kilometers for each raster block pulled into an array. Start with 5 if you aren't sure and reduce
        if there is an error. To increase processing speed rasters are subdivided iteratively into smaller blocks and pulled
        into arrays for processing. Very large block sizes may use a lot of memory and result in an error.

        11: overwrite_data:
        True/False flag if existing data in nodes_fc and lc_point_fc can be overwritten.

        OUTPUTS
        0. nodes_fc:
        New fields are added into nodes_fc with the Landcover and elevation values for each transect sample

        1. lc_point_fc:
        New point feature class created with a point at each x/y sample and the sample raster values

        """
        # Import system modules
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

        # ----------------------------------------------------------------------
        # ----------------------------------------------------------------------

        # General script steps include:
        # 1. open nodes fc. iterate and read info from node fc into a dict

        # 2. create a list with the x/y and related info for each lc sample
        #    calculate the extent bounding box for the entire dataset

        # 3. create a list holding the bounding box coords for each block iteration

        # 4. loop through each block
            #- sample the raster for all samples in the block
            #- update the node dict
            #- update the node fc
            #- update the lc sample fc
            #- continue to the next block

        # Future Updates
        # -Change the node dict so the node is the primary key
        # -Build the block list based on the nodes and then build point list iteratively instead of building them into
        #   one huge list. The huge list results in a memory error for large areas
        # -Eliminate arcpy and use gdal for reading/writing feature class data

        env.overwriteOutput = True

        def nested_dict():
            """Build a nested dictionary"""
            return defaultdict(nested_dict)

        def read_nodes_fc(nodes_fc, overwrite_data, addFields):
            """Reads the input point feature class and returns the STREAM_ID,
            NODE_ID, and X/Y coordinates as a nested dictionary"""
            nodeDict = nested_dict()
            incursorFields = ["NODE_ID", "STREAM_ID", "STREAM_KM", "ASPECT", "LEFT", "RIGHT", "SHAPE@X", "SHAPE@Y"]

            # Get a list of existing fields
            existingFields = [f.name for f in arcpy.ListFields(nodes_fc)] 

            # Check to see if the last field exists if yes add it. 
            # Grabs last field because often the first field, emergent, is zero
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
                        nodeDict[row[0]]["ASPECT"] = row[3]
                        nodeDict[row[0]]["LEFT"] = row[4]
                        nodeDict[row[0]]["RIGHT"] = row[5]
                        nodeDict[row[0]]["POINT_X"] = row[6]
                        nodeDict[row[0]]["POINT_Y"] = row[7]
                else:
                    for row in Inrows:
                        # Is the data null or zero, if yes grab it.
                        if row[8] is None or row[8] == 0 or row[8] < -9998:
                            nodeDict[row[0]]["STREAM_ID"] = row[1]
                            nodeDict[row[0]]["STREAM_KM"] = row[2]
                            nodeDict[row[0]]["ASPECT"] = row[3]
                            nodeDict[row[0]]["LEFT"] = row[4]
                            nodeDict[row[0]]["RIGHT"] = row[5]
                            nodeDict[row[0]]["POINT_X"] = row[6]
                            nodeDict[row[0]]["POINT_Y"] = row[7]
    
            if len(nodeDict) == 0:
                sys.exit("The fields checked in the input point feature class " +
                         "have existing data. There is nothing to process. Exiting")
              
            return nodeDict

        def update_lc_point_fc(lc_point_list, type, lc_point_fc, nodes_fc,
                               nodes_in_block, overwrite_data, proj):
            """Creates/updates the output landcover sample point feature
            class using the data from the landcover point list"""
            #print("Exporting data to land cover sample feature class")
    
            cursorfields = ["POINT_X","POINT_Y"] +["STREAM_ID","NODE_ID", "SAMPLE_ID", 
                                                   "TRANS_DIR", "TRANSECT",
                                                   "SAMPLE", "KEY"] + type    
    
            # Check if the output exists and create if not
            if not arcpy.Exists(lc_point_fc):
                arcpy.CreateFeatureclass_management(os.path.dirname(lc_point_fc),
                                                    os.path.basename(lc_point_fc),
                                                    "POINT","","DISABLED","DISABLED",proj)
        
                # Determine Stream ID field properties
                sid_type = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].type
                sid_precision = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].precision
                sid_scale = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].scale
                sid_length = arcpy.ListFields(nodes_fc,"STREAM_ID")[0].length
        
                typeDict = {"POINT_X": "DOUBLE", 
                            "POINT_Y": "DOUBLE", 
                            "NODE_ID": "LONG",
                            "SAMPLE_ID": "LONG",
                            "TRANS_DIR": "DOUBLE",
                            "TRANSECT": "SHORT",
                            "SAMPLE": "SHORT",
                            "KEY": "TEXT",}
                
                for t in type:
                    typeDict[t] = "DOUBLE"
    
                # Add attribute fields
                for f in cursorfields:
                    if f == "STREAM_ID":
                        arcpy.AddField_management(lc_point_fc, f, sid_type,
                                                  sid_precision, sid_scale, sid_length,
                                                  "", "NULLABLE", "NON_REQUIRED")
                    else:
                        arcpy.AddField_management(lc_point_fc, f, typeDict[f], "",
                                                  "", "", "", "NULLABLE", "NON_REQUIRED")                 
    
            if not overwrite_data:
                # Build a query to retrieve existing rows from the nodes
                # that need updating
                whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_in_block))        
        
                # delete those rows
                with arcpy.da.UpdateCursor(lc_point_fc,["NODE_ID"], whereclause) as cursor:  
                    for row in cursor:
                        cursor.deleteRow()     

            with arcpy.da.InsertCursor(lc_point_fc, ["SHAPE@X","SHAPE@Y"] +
                                       cursorfields) as cursor:
                for row in lc_point_list:
                    cursor.insertRow(row)

        def setup_lcdata_headers(transsample_count):
            """Generates a list of the landcover data file column header names and data types"""
    
            type = ["ELE"]

            lcheaders = []
            otherheaders = []
     
            dirs = ["L", "R"]

            # Concatenate the type, dir, and sample and order in the correct way
            for d, dir in enumerate(dirs):
                for s, sample in enumerate(range(1, int(transsample_count) + 1)):
                    if d==0 and s==0:
                        lcheaders.append("LC_T0_S0") # add emergent
                        lcheaders.append("LC_{0}_S{1}".format(dir, sample))
                    else:
                        lcheaders.append("LC_{0}_S{1}".format(dir, sample))
    
            for t in type:
                for d, dir in enumerate(dirs):
                    for s, sample in enumerate(range(1, int(transsample_count) + 1)):
                        if t !="ELE" and d==0 and s==0:
                            otherheaders.append(t+"T0_S0") # add emergent
                            otherheaders.append("{0}_{1}_S{2}".format(t, dir, sample))
                        else:
                            otherheaders.append("{0}_{1}_S{2}".format(t, dir, sample))

            return lcheaders, otherheaders

        def coord_to_array(easting, northing, block_x_min, block_y_max, x_cellsize, y_cellsize):
            """converts x/y coordinates to col and row of the array"""
            col_x = int(((easting - block_x_min) - ((easting - block_x_min) % x_cellsize)) / x_cellsize)
            row_y = int(((block_y_max - northing) - ((block_y_max - northing) % y_cellsize)) / y_cellsize)
    
            return [col_x, row_y]

        def create_lc_point_list(nodeDict, nodes_in_block, transsample_count, transsample_distance, start_bank):
            """This builds a unique long form list of information for all the
            landcover samples in the block. This list is used to
            create/update the output feature class."""
    
            lc_point_list = []
            samplesPerNode = (2 * transsample_count) + 1

            if start_bank:
                # need to reduce sample number by 1 so first sample starts at bank
                bx = 1
            else:
                bx = 0

            for nodeID in nodes_in_block:
                origin_x = nodeDict[nodeID]["POINT_X"]
                origin_y = nodeDict[nodeID]["POINT_Y"]
                aspect = nodeDict[nodeID]["ASPECT"]
                streamID = nodeDict[nodeID]["STREAM_ID"]
                sampleID = nodeID * transsample_count

                # Determine the offset from the node to the start of the first transect.
                # Measured in meters here and converted to map units farther down
                if start_bank:
                    offset_l = nodeDict[nodeID]["LEFT"]
                    offset_r = nodeDict[nodeID]["RIGHT"]
                else:
                    offset_l = 0
                    offset_r = 0

                # calculate right and left transect directions in degrees
                dir_left = aspect - 90
                dir_right = aspect + 90

                if dir_left < 0:
                     dir_left = dir_left + 360

                if dir_right > 360:
                    dir_right = dir_right - 360

                dirs = [dir_left, dir_right]

                # This is the emergent/stream sample
                lc_point_list.append([origin_x, origin_y, origin_x, origin_y,
                                      streamID, nodeID, sampleID,
                                      0, 0, 0, "T0_S0"])

                for d, dir in enumerate(dirs):

                    if d == 0:
                        offset = offset_l
                        prefix = "L"
                    if d == 1:
                        offset = offset_r
                        prefix = "R"

                    for sample in range(1, int(transsample_count + 1)):
                        # Calculate the x and y coordinate of the 
                        # landcover sample location
                        node_to_sample_dis = (offset + ((sample - bx) * transsample_distance)) * con_from_m
                        pt_x = (node_to_sample_dis * sin(radians(dir))) + origin_x
                        pt_y = (node_to_sample_dis * cos(radians(dir))) + origin_y
                
                        key = '{0}_S{1}'.format(prefix, sample)
                        sampleID = (nodeID * samplesPerNode) + (d * transsample_count) + sample
        
                        # Add to the list          
                        lc_point_list.append([pt_x, pt_y, pt_x, pt_y,
                                              streamID, nodeID, sampleID,
                                              dir, d+1, sample, key])
     
            return lc_point_list

        def create_block_list(nodes, block_size):
            """Returns two lists, one containing the coordinate extent
            for each block that will be iteratively extracted to an array
            and the other containing node IDs within each block extent."""
    
            arcpy.AddMessage("Calculating block extents")    
            print("Calculating block extents")    
            x_coord_list = [nodeDict[nodeID]["POINT_X"] for nodeID in nodes]
            y_coord_list = [nodeDict[nodeID]["POINT_Y"] for nodeID in nodes]
    
            # calculate the buffer distance (in raster spatial units) to add to 
            # the base bounding box when extracting to an array. The buffer is 
            # equal to the sample distance + 1 to make sure the block includes 
            # all the landcover samples for each node.
            buffer_dis = int((transsample_count + 1) * transsample_distance * con_from_m)
    
            # calculate bounding box extent for samples
            x_min = min(x_coord_list)
            x_max = max(x_coord_list)
            y_min = min(y_coord_list) 
            y_max = max(y_coord_list)
    
            x_width = int(x_max - x_min + 1)
            y_width = int(y_max - y_min + 1)
    
            block_extents = []
            block_nodes = []
      
            # Build data blocks
            for x in range(0, x_width, block_size):
                for y in range(0, y_width, block_size):

                    # Lower left coordinate of block (in map units)
                    block0_x_min = min([x_min + x, x_max])
                    block0_y_min = min([y_min + y, y_max])
                    # Upper right coordinate of block (in map units)
                    block0_x_max = min([block0_x_min + block_size, x_max])
                    block0_y_max = min([block0_y_min + block_size, y_max])
            
                    block_x_min = block0_x_max
                    block_x_max = block0_x_min
                    block_y_min = block0_y_max
                    block_y_max = block0_y_min
            
                    nodes_in_block = []
                    for nodeID in nodes:
                        node_x = nodeDict[nodeID]["POINT_X"]
                        node_y = nodeDict[nodeID]["POINT_Y"]
                        if (block0_x_min <= node_x <= block0_x_max and
                            block0_y_min <= node_y <= block0_y_max):
                    
                            nodes_in_block.append(nodeID)
                    
                            # Minimize the size of the block0 by the true 
                            # extent of the nodes in the block
                            if block_x_min > node_x: block_x_min = node_x
                            if block_x_max < node_x: block_x_max = node_x
                            if block_y_min > node_y: block_y_min = node_y
                            if block_y_max < node_y: block_y_max = node_y
            
                    if nodes_in_block:
                        # add the block extent for processing
                        # order 0 left,      1 bottom,    2 right,     3 top
                        block_extents.append((block_x_min - buffer_dis, block_y_min - buffer_dis,
                                              block_x_max + buffer_dis, block_y_max + buffer_dis))           
                        block_nodes.append(nodes_in_block)
    
            return block_extents, block_nodes
    
        def sample_raster(block, lc_point_list, raster, con):
#-----MSF - 9/15/2025
            # Normalize raster input so dataset path is always used (works for GPRasterLayer or dataset)
            try:
                raster_src = arcpy.Describe(raster).catalogPath
            except:
                raster_src = raster
#-------------------    
    
            if con is not None:
                nodata_to_value = -9999 / con
            else:
                nodata_to_value = -9999
        
            # localize the block extent values
            block_x_min = block[0]
            block_y_min = block[1]
            block_x_max = block[2]
            block_y_max = block[3]
    
            x_cellsize = float(arcpy.GetRasterProperties_management(raster_src, "CELLSIZEX").getOutput(0))
            y_cellsize = float(arcpy.GetRasterProperties_management(raster_src, "CELLSIZEY").getOutput(0)) 

            # Get the coordinate extent of the input raster
            raster_x_min = float(arcpy.GetRasterProperties_management(raster_src, "LEFT").getOutput(0))
            raster_y_min = float(arcpy.GetRasterProperties_management(raster_src, "BOTTOM").getOutput(0))
            raster_x_max = float(arcpy.GetRasterProperties_management(raster_src, "RIGHT").getOutput(0))
            raster_y_max = float(arcpy.GetRasterProperties_management(raster_src, "TOP").getOutput(0))

            # Calculate the block x and y offset from the raster and adjust
            # the block coordinates so they are at the raster cell corners.
            # This is for ESRI's RastertoNumpyArray function which defaults to the adjacent
            # lower left cell
            block_x_min_corner = block_x_min - ((block_x_min - raster_x_min) % x_cellsize)
            block_y_min_corner = block_y_min - ((block_y_min - raster_y_min) % y_cellsize)
            block_x_max_corner = block_x_max + ((raster_x_max - block_x_max) % x_cellsize)
            block_y_max_corner = block_y_max + ((raster_y_max - block_y_max) % y_cellsize)
    
            # calculate the number of cols/rows from the lower left
            ncols = int((block_x_max_corner - block_x_min_corner) / x_cellsize)
            nrows = int((block_y_max_corner - block_y_min_corner) / y_cellsize)
    
            # Construct the array. Note returned array is (row, col) so (y, x)
            try:
# updated raster to raster_src - MSF 9/15/2025    
                raster_array = arcpy.RasterToNumPyArray(raster_src, arcpy.Point(block_x_min_corner, block_y_min_corner),
                                                        ncols, nrows, nodata_to_value)
            except:
                tbinfo = traceback.format_exc()
                pymsg = tbinfo + "\nError Info:\n" + "\nNot enough memory. Reduce the block size"       
                sys.exit(pymsg)    
    
            # convert array values to meters if needed
            if con is not None:
                raster_array = raster_array * con
    
            lc_point_list_new = []
            if raster_array.max() > -9999:
                # There is at least one pixel of data
                for point in lc_point_list:
                    xy = coord_to_array(point[0], point[1], block_x_min_corner, block_y_max_corner, x_cellsize, y_cellsize)
                    point.append(raster_array[xy[1], xy[0]])
                    lc_point_list_new.append(point)
            else:
                # No data, add -9999
                for point in lc_point_list:
                    point.append(-9999)
                    lc_point_list_new.append(point)
            return lc_point_list_new            

        def update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_to_query):
            """Updates the input point feature class with data from the
            nodes dictionary"""

            # Build a query to retrieve just the nodes that needs updating
            whereclause = """{0} IN ({1})""".format("NODE_ID", ','.join(str(i) for i in nodes_to_query))

            with arcpy.da.UpdateCursor(nodes_fc,["NODE_ID"] + addFields, whereclause) as cursor:  
                for row in cursor:
                    for f, field in enumerate(addFields):
                        nodeID =row[0]
                        row[f+1] = nodeDict[nodeID][field]
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
            """Returns the conversion factor to get from
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

        # enable garbage collection
        gc.enable()

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        start_bank = bool(parameters[1].value)
        transsample_count = int(parameters[2].value)
        transsample_distance = float(parameters[3].value)
        lc_raster = parameters[4].valueAsText
        lc_units = parameters[5].valueAsText
        z_raster = parameters[6].valueAsText
        z_units = parameters[7].valueAsText
        lc_point_fc = parameters[8].valueAsText
        block_size = float(parameters[9].value)
        overwrite_data = bool(parameters[10].value)

        try:
            arcpy.AddMessage("Step 5: Sample Landcover - Orthogonal Method")
            print("Step 5: Sample Landcover - Orthogonal Method")
    
            # keeping track of time
            startTime= time.time()
    
            # Check if the node fc exists
            if not arcpy.Exists(nodes_fc):
                arcpy.AddError("\nThis output does not exist: \n" +
                               "{0}\n".format(nodes_fc))
                sys.exit("This output does not exist: \n" +
                         "{0}\n".format(nodes_fc))
    
            # Check if the lc point fc exists and delete if needed
            if arcpy.Exists(lc_point_fc) and overwrite_data:
                arcpy.Delete_management(lc_point_fc)
    
    
            # Determine input spatial units and set unit conversion factors
            proj = arcpy.Describe(nodes_fc).SpatialReference
            proj_ele = arcpy.Describe(z_raster).spatialReference
            proj_lc = arcpy.Describe(lc_raster).spatialReference
    
            con_from_m = from_meters_con(nodes_fc)
            con_lc_to_m = from_z_units_to_meters_con(lc_units)
            con_z_to_m = from_z_units_to_meters_con(z_units)
    
            # convert block size from km to meters to units of the node fc
            # in the future block size should be estimated based on available memory
            # memorysize = datatypeinbytes*nobands*block_size^2
            # block_size = int(sqrt(memorysize/datatypeinbytes*nobands))
            if block_size in ["#", "", None]:
                block_size = int(con_from_m * 5000)
            else:
                block_size = int(con_from_m * block_size * 1000)
    
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

            lcheaders, otherheaders = setup_lcdata_headers(transsample_count)
    
            addFields = lcheaders + otherheaders
    
            # Get a list of existing fields
            existingFields = [f.name for f in arcpy.ListFields(nodes_fc)] 
        
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
            nodes = list(nodeDict.keys())
            nodes.sort()
   
            # Build the block list
            block_extents, block_nodes = create_block_list(nodes, block_size)
    
            # Iterate through each block, calculate sample coordinates,
            # convert raster to array, sample the raster
            total_samples = 0
            for p, block in enumerate(block_extents):
                nodes_in_block = block_nodes[p]
                nodes_in_block.sort()
                arcpy.AddMessage("Processing block {0} of {1}".format(p + 1, len(block_extents)))
                print("Processing block {0} of {1}".format(p + 1, len(block_extents)))
        
                # build the landcover sample list
                lc_point_list = create_lc_point_list(nodeDict, nodes_in_block, transsample_count, transsample_distance, start_bank)
        
                for t, (type, raster) in enumerate(rasterDict.items()):
                    if raster is None:
                        for i in range(0, len(lc_point_list)):
                            lc_point_list[i].append(-9999)
                    else: 
                        if raster == z_raster:
                            con = con_z_to_m
                        elif raster == lc_raster:
                            con = con_lc_to_m
                        else:
                            con = None  
            
                        lc_point_list = sample_raster(block, lc_point_list, raster, con)
        
                    # Update the node dict
                    for row in lc_point_list:
                        key = "{0}_{1}".format(type, row[10])
                        nodeDict[row[5]][key] = row[11 + t]
        
                # Write the landcover data to the TTools point feature class 
                update_nodes_fc(nodeDict, nodes_fc, addFields, nodes_in_block)
        
                # Build the output point feature class using the data         
                update_lc_point_fc(lc_point_list, list(rasterDict.keys()), lc_point_fc,
                                   nodes_fc, nodes_in_block, overwrite_data, proj)
    
                total_samples = total_samples + len(lc_point_list)
                del lc_point_list
                gc.collect()
    
            endTime = time.time()
    
            elapsedmin = ceil(((endTime - startTime) / 60)* 10)/10
            mspersample = timedelta(seconds=(endTime - startTime) /
                                    total_samples).microseconds
            arcpy.AddMessage("Process Complete in {0} minutes. {1} microseconds per sample".format(elapsedmin, mspersample))
            print("Process Complete in {0} minutes. {1} microseconds per sample".format(elapsedmin, mspersample))

        # For arctool errors
        except arcpy.ExecuteError:
            msgs = arcpy.GetMessages(2)
            arcpy.AddError(msgs)
            print(msgs)

        # For other errors
        except:
            tbinfo = traceback.format_exc()

            pymsg = "PYTHON ERRORS:\n" + tbinfo + "\nError Info:\n" + str(sys.exc_info()[1])
            msgs = "ArcPy ERRORS:\n" + arcpy.GetMessages(2) + "\n"

            arcpy.AddError(pymsg)
            arcpy.AddError(msgs)

            print(pymsg)
            print(msgs)

            
