#!/usr/bin/env python3
"""
TTools

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
import traceback

class Toolbox:
    def __init__(self):
        self.label = "TTools"
        self.alias = "TTools"
        self.tools = [CreateStreamNodes,
                      MeasureChannelWidth,
                      SampleElevationGradient,
                      MeasureTopographicAngles,
                      SampleLandcoverStarPattern,
                      SampleLandcoverOrthogonalMethod]


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

# Variable assignment based on user input - MSF
        streamline_fc = parameters[0].valueAsText
        sid_field = parameters[1].valueAsText
        node_dx = float(parameters[2].value)
        cont_stream_km = bool(parameters[3].value)
        checkDirection = bool(parameters[4].value)
        z_raster = parameters[5].valueAsText if parameters[5].value else None
        nodes_fc = parameters[6].valueAsText

        try:
            from ttools import step1

            step1(streamline_fc, sid_field, node_dx, cont_stream_km,
                  nodes_fc, checkDirection=checkDirection, z_raster=z_raster)

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

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        rb_fc = parameters[1].valueAsText
        lb_fc = parameters[2].valueAsText
        overwrite_data = bool(parameters[3].value)

        try:
            from ttools import step2

            step2(nodes_fc, rb_fc, lb_fc, overwrite_data=overwrite_data)

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
        param5.value = 10

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

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        searchCells = int(parameters[1].value)
        smooth_flag = bool(parameters[2].value)
        z_raster = parameters[3].valueAsText
        z_units = parameters[4].valueAsText
        block_size = float(parameters[5].value)
        overwrite_data = bool(parameters[6].value)

        try:
            from ttools import step3

            step3(nodes_fc, searchCells, smooth_flag, z_raster, z_units,
                  block_size=block_size, overwrite_data=overwrite_data)

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
            displayName="Output Topographic Point Feature Class",
            name="topo_fc",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        param6 = arcpy.Parameter(
            displayName="Raster Block Size (km)",
            name="block_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param6.value = 10

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
            from ttools import step4

            step4(nodes_fc, topo_directions, searchDistance_max_km,
                  z_raster, z_units, topo_fc,
                  block_size=block_size, overwrite_data=overwrite_data)

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


class SampleLandcoverStarPattern(object):
    def __init__(self):
        self.label = "Step 5A: Sample Landcover - Star Pattern"
        self.description = "Sample landcover raster values along transects oriented outward from the stream node in a star pattern."
        self.canRunInBackground = False

    def getParameterInfo(self):

        param0 = arcpy.Parameter(
            displayName="Stream Nodes Feature Class",
            name="nodes_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")

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

        param6 = arcpy.Parameter(
            displayName="Landcover Raster",
            name="lc_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        param7 = arcpy.Parameter(
            displayName="Landcover Units",
            name="lc_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param7.filter.type = "ValueList"
        param7.filter.list = ["Meters", "Feet", "None"]

        param8 = arcpy.Parameter(
            displayName="Elevation Raster",
            name="z_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input")

        param9 = arcpy.Parameter(
            displayName="Elevation Units",
            name="z_units",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param9.filter.type = "ValueList"
        param9.filter.list = ["Meters", "Feet"]

        param10 = arcpy.Parameter(
            displayName="Output Landcover Point Feature Class",
            name="lc_point_fc",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        param11 = arcpy.Parameter(
            displayName="Raster Block Size (km)",
            name="block_size",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param11.value = 10

        param12 = arcpy.Parameter(
            displayName="Overwrite Existing Data",
            name="overwrite_data",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param12.value = True


        return [param0, param1, param2, param3, param4, param5, param6,
                param7, param8, param9, param10, param11, param12]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):
    # Log input parameters to file
        log_run(self.label, parameters)

# Variable assignment based on user input - MSF
        nodes_fc = parameters[0].valueAsText
        trans_count = int(parameters[1].value)
        transsample_count = int(parameters[2].value)
        transsample_distance = float(parameters[3].value)
        zone_sample = bool(parameters[4].value)
        heatsource8 = bool(parameters[5].value)
        lc_raster = parameters[6].valueAsText
        lc_units = parameters[7].valueAsText
        z_raster = parameters[8].valueAsText
        z_units = parameters[9].valueAsText
        lc_point_fc = parameters[10].valueAsText
        block_size = float(parameters[11].value)
        overwrite_data = bool(parameters[12].value)

        try:
            from ttools import step5_star

            step5_star(nodes_fc, trans_count, transsample_count,
                       transsample_distance, zone_sample, heatsource8,
                       lc_raster, lc_units, z_raster, z_units,
                       lc_point_fc, block_size, overwrite_data)

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

class SampleLandcoverOrthogonalMethod(object):
    def __init__(self):
        self.label = "Step 5B: Sample Landcover - Orthogonal Method"
        self.description = "Sample landcover raster values along orthogonal transects from stream banks or nodes."
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
        param9.value = 10
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
            from ttools import step5_ortho

            step5_ortho(nodes_fc, start_bank, transsample_count,
                        transsample_distance,
                        lc_raster, lc_units, z_raster, z_units,
                        lc_point_fc, block_size, overwrite_data)

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
