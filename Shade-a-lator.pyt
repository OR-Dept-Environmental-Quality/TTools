

import os
import datetime
from pathlib import Path

def log_run(tool_name, params, log_dir=None):
    """
    Every time shade-a-lator is executed log_run will:
    Append a line to a .log file.
    The log filename includes the current date.
    The log line includes the timestamp and parameter values.
    """

    # Note, any None, empty string, or "#" is turned into "<null>". Blank or no entries cause a crash

    # Directory to store logs
    if log_dir is None:
        # Same folder as this script
        log_dir = os.path.join(os.path.dirname(__file__), "logs")

    os.makedirs(log_dir, exist_ok=True)

    # Log filename with date
    today_str = datetime.datetime.now().strftime("%Y-%m-%d")
    log_file = os.path.join(log_dir, f"shadealator_log_{today_str}.log")

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

            param_values.append(f"{p.name} = {val}")
        except Exception as e:
            param_values.append(f"<error: {e}>")

    log_entry = f"{timestamp} | Shade-a-lator: {tool_name} | " + " | ".join(param_values) + "\n"

    # Append to log file
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(log_entry)

# ------------------------------------------------------------

import sys
import traceback

ttools_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(ttools_dir))
hs_exe = Path(ttools_dir, "shadealator", "bin", "hs.exe")

from ttools.geo_package import set_geo_package
from ttools.utils import set_message_method, set_warning_method

set_geo_package(use_arcpy=True)
set_message_method(arcpy.AddMessage)
set_warning_method(arcpy.AddWarning)

class Toolbox:
    def __init__(self):
        self.label = "Shade-a-lator"
        self.alias = "shadealator"
        self.tools = [SetupSolarModel, RunSolarModel]

class SetupSolarModel:

    def __init__(self):
        self.label = "1. Setup Heat Source 9 Solar Model"
        self.description = "Create Heat Source version 9 solar model control files and model inputs from TTools nodes. Requires TTools steps 1 - 5 be run first."

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            displayName="Stream Nodes Feature Class",
            name="nodes_fc",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input")
        param0.filter.list = ["Point"]

        param1 = arcpy.Parameter(
            displayName="Model Grouping Field",
            name="model_grouping_field",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        param1.parameterDependencies = [param0.name]

        param2 = arcpy.Parameter(
            displayName="Spacing Between Nodes (meters)",
            name="longsample",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        param3 = arcpy.Parameter(
            displayName="Number of Transects per Node",
            name="trans_count",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        param4 = arcpy.Parameter(
            displayName="Number of Samples per Transect",
            name="transsample_count",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        param5 = arcpy.Parameter(
            displayName="Distance Between Transect Samples (meters)",
            name="transsample_distance",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")

        param6 = arcpy.Parameter(
            displayName="Land Cover Codes Input File",
            name="lccodes",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        param6.filter.list = ["csv"]

        param7 = arcpy.Parameter(
            displayName="Canopy Data Type",
            name="canopy_data",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param7.filter.type = "ValueList"
        param7.filter.list = ["LAI", "CanopyCover"]

        param8 = arcpy.Parameter(
            displayName="Use Heat Source 8 Land Cover Methods",
            name="heatsource8",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param8.value = False

        param9 = arcpy.Parameter(
            displayName="Land Cover Sample Method",
            name="lcsampmethod",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param9.filter.type = "ValueList"
        param9.filter.list = ["point", "zone"]
        param9.value = "point"

        param10 = arcpy.Parameter(
            displayName="Account For Emergent Vegetation Shading",
            name="emergent",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param10.value = False

        param11 = arcpy.Parameter(
            displayName="Modeling Start Date (yyyy-mm-dd)",
            name="modelstart",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        param12 = arcpy.Parameter(
            displayName="Modeling End Date (yyyy-mm-dd)",
            name="modelend",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        param13 = arcpy.Parameter(
            displayName="Time Offset From UTC (hours)",
            name="offset",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")

        param14 = arcpy.Parameter(
            displayName="Model Time Step (minutes)",
            name="dt",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param14.value = 1

        param15 = arcpy.Parameter(
            displayName="Output Time Step (minutes)",
            name="outputdt",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param15.value = 60

        param16 = arcpy.Parameter(
            displayName="Model Directory",
            name="model_dir",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        param17 = arcpy.Parameter(
            displayName="Input Directory",
            name="inputdir",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param17.value = 'inputs'

        param18 = arcpy.Parameter(
            displayName="Output Directory",
            name="outputdir",
            datatype="GPString",
            parameterType="Required",
            direction="Input")
        param18.value = 'outputs'

        param19 = arcpy.Parameter(
            displayName="CSV Mode",
            name="csv_mode",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param19.value = False

        param20 = arcpy.Parameter(
            displayName="Overwrite Existing Setup Files",
            name="overwrite",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param20.value = True

        param21 = arcpy.Parameter(
            displayName="Model Description/User Notes",
            name="usertxt",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        param22 = arcpy.Parameter(
            displayName="Simulation Name",
            name="name",
            datatype="GPString",
            parameterType="Optional",
            direction="Input")

        return [param0, param1, param2, param3, param4, param5, param6,
                param7, param8, param9, param10, param11, param12, param13,
                param14, param15, param16, param17, param18, param19, param20, param21, param22]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        nodes_fc = parameters[0].valueAsText
        model_grouping_field = parameters[1].valueAsText
        longsample = int(parameters[2].value)
        trans_count = int(parameters[3].value)
        transsample_count = int(parameters[4].value)
        transsample_distance = float(parameters[5].value)
        lccodes = parameters[6].valueAsText
        canopy_data = parameters[7].valueAsText
        heatsource8 = bool(parameters[8].value)
        lcsampmethod = parameters[9].valueAsText
        emergent = bool(parameters[10].value)
        modelstart = parameters[11].valueAsText
        modelend = parameters[12].valueAsText
        offset = int(parameters[13].value)
        dt = float(parameters[14].value)
        outputdt = float(parameters[15].value)
        model_dir = parameters[16].valueAsText
        inputdir = parameters[17].valueAsText
        outputdir = parameters[18].valueAsText
        csv_mode = bool(parameters[19].value)
        overwrite = bool(parameters[20].value)
        usertxt = parameters[21].valueAsText
        name = parameters[22].valueAsText

        log_dir = Path(model_dir).resolve()
        log_run(self.label, parameters, log_dir)

        try:

            if not hs_exe.is_file():
                msg = "Heat Source was not found at {0}. Include hs.exe in the shadealator/bin folder before running this tool.".format(hs_exe)
                raise FileNotFoundError(msg)

            from shadealator import hs
            hs.setup.solar(nodes_fc, model_grouping_field, longsample, trans_count, transsample_count, transsample_distance,
                           lccodes, canopy_data, heatsource8, lcsampmethod, emergent, modelstart, modelend,
                           offset, dt, outputdt, model_dir, inputdir, outputdir, csv_mode,
                           overwrite, usertxt, name)

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

class RunSolarModel:

    def __init__(self):
        self.label = "2. Run Heat Source 9 Solar Model"
        self.description = "Run Heat Source solar models from a model list CSV or Excel file and update their statuses."

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            displayName="Model List",
            name="models_file",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")
        param0.filter.list = ["csv", "xlsx"]

        param1 = arcpy.Parameter(
            displayName="Overwrite Completed Models",
            name="overwrite",
            datatype="GPBoolean",
            parameterType="Required",
            direction="Input")
        param1.value = False

        param2 = arcpy.Parameter(
            displayName="Concurrent Models To Run",
            name="concurrent_models",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param2.value = 1

        return [param0, param1, param2]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def execute(self, parameters, messages):

        # Variable assignment based on user input
        models_file = parameters[0].valueAsText
        overwrite = bool(parameters[1].value)
        concurrent_models = int(parameters[2].value)

        log_dir = Path(models_file).resolve().parent
        log_run(self.label, parameters, log_dir)

        try:
            if not hs_exe.is_file():
                msg = "Heat Source was not found at {0}. Include hs.exe in the shadealator/bin folder before running this tool.".format(hs_exe)
                raise FileNotFoundError(msg)

            from shadealator import hs
            hs.run_solar(models_file, overwrite, concurrent_models)

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