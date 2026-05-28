TTools 
======
Current Version: TTools 9.2

## About

TTools is a collection of python scripts used in conjunction with ESRI ArcGIS Pro 3.x to sample and assemble 
stream channel and land cover geospatial data for input into the [Heat Source][1] model, Washington Department 
of Ecology's [Shade][2] model, or the shade input file for [CE-QUAL-W2][3].

Requires
ESRI ArcGIS Pro w/ Spatial Analyst extension
Python 3.10+
Numpy 1.26+

Previous TTools extensions for ESRI ArcView 3, ArcGIS 9 - 10.8 can be found here: http://www.oregon.gov/deq/wq/tmdls/Pages/TMDLs-Tools.aspx

[1]: https://github.com/DEQrmichie/heatsource-9
[2]: https://ecology.wa.gov/Research-Data/Data-resources/Models-spreadsheets/Modeling-the-environment/Models-tools-for-TMDLs
[3]: http://www.ce.pdx.edu/w2/

## File Structure

The contents of `TTools/` can be saved anywhere on your computer but the file structure and file names must be maintained, especially if used as a toolbox directly within ArcPro. 

```text
TTools/
├── TTools.pyt
├── TTools.pyt.xml
├── TTools.CreateStreamNodes.pyt.xml
├── TTools.MeasureChannelWidth.pyt.xml
├── TTools.MeasureTopographicAngles.pyt.xml
├── TTools.SampleElevationGradient.pyt.xml
├── TTools.SampleLandcoverOrthogonalMethod.pyt.xml
├── TTools.SampleLandcoverStarPattern.pyt.xml
└── ttools/
    ├── __init__.py
    ├── __version__.py
    ├── step1.py
    ├── step2.py
    ├── step3.py
    ├── step4.py
    ├── step5.py
    └── utils.py
```

## Add TTools to ArcPro

Before starting the tool, the Arc GIS Pro TTools Toolbox requires the following: 
1. Toolbox Files. Ensure the following files are present in the same directory:
   * `TTools.pyt` (the main toolbox file)
   * `TTools.pyt.xml` (the associated XML metadata file)
2. Metadata Help Files. Each tool within the toolbox requires its own metadata XML help file. These should also reside in the same directory as the `TTools.pyt` file. The required files are:
   * `TTools.CreateStreamNodes.pyt.xml`
   * `TTools.MeasureChannelWidth.pyt.xml`
   * `TTools.SampleElevationGradient.pyt.xml`
   * `TTools.MeasureTopographicAngles.pyt.xml`
   * `TTools.SampleLandcoverStarPattern.pyt.xml`
   * `TTools.SampleLandcoverOrthogonalMethod.pyt.xml`
   
   These metadata files provide help documentation accessible within ArcGIS Pro. To view the help, hover the cursor next to the start of each field name in a tool and click the blue icon that appears.

3.	Required Extension. The ArcGIS Pro Spatial Analyst extension must be enabled to use the toolbox.
4. Load the ArcToolbox in ArcGIS Pro. Open or create a new ArcGIS Pro project from the Catalog pane. If it is not visible, go to: View → Catalog Pane. This will show the ArcGIS Pro Catalog Pane and show Toolboxes
5. Add the toolbox to ArcGIS Pro. From the main menu ribbon in ArcGIS Pro, right-click Toolbox → Add Toolbox → Browse to the toolbox file (`TTools.pyt`)  → Select the Toolbox file `TTools.pyt` → Click OK 
6. From the Catalog pane, click on the arrow next to Toolboxes to show all available Toolboxes. 
7. Click on the arrow next to the `TTools.pyt` toolbox in the Catalog to display the tools available in the toolbox.

See the TTools user guide or [ESRI documentation][4] for additional instructions.

[4]:(https://doc.esri.com/en/arcgis-pro/latest/help/projects/connect-to-a-toolbox.html#A9E)

## Use TTools from a python script

TTools can be used from a python script. There is no python package install required as long as 
1. The lower case `ttools` folder shown in the file structure diagram is saved somewhere on your computer. 
2. The script points to the `ttools` folder as shown in the example below
3. There is a valid ArcPro license with spatial analyst extension.
4. The python environment being used can import arcpy.

Documentation for each TTools step can be found in the user guide, in each python file e.g. `step1.py`, or using (`help(step1)`).

```python
# This example is also saved to `TTools_Example_Script.py` located in the repository root.

import os
import sys
import arcpy

# ttools_package_folder points to the ttools package folder (lowercase).
ttools_package_folder = r"C:\Workspace\TTools\ttools"

# Working geodatabase. This is where the input and output features are saved
run_gdb = r"C:\Workspace\TTools_JC_test_features_gdb\JohnsonCreek.gdb"

# Input features
streamline_fc = os.path.join(run_gdb, "jc_thalweg")
rb_fc = os.path.join(run_gdb, "jc_rightbank")
lb_fc = os.path.join(run_gdb, "jc_leftbank")

# Output features
nodes_fc = os.path.join(run_gdb, "jc_nodes")
topo_samples_fc = os.path.join(run_gdb, "jc_topo_samples")
lc_samples_fc = os.path.join(run_gdb, "jc_lc_samples")

# Input rasters
raster_dir = r"C:\Workspace\TTools_JC_test_rasters"
z_raster = os.path.join(raster_dir, "jcw_be_m_mosaic.tif")
lc_raster = os.path.join(raster_dir, "jcw_vght_m_mosaic.tif")

# ------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(ttools_package_folder))

from ttools import step1, step2, step3, step4, step5_ortho, step5_star

# ------------------------------------------------------------------------
# TTools Step 1: Create Stream Nodes
step1(streamline_fc=streamline_fc,
      sid_field="NAME",
      node_dx=50,
      cont_stream_km=False,
      nodes_fc=nodes_fc,
      checkDirection=True,
      z_raster=z_raster)

# ------------------------------------------------------------------------
# TTools Step 2: Measure Channel Widths
step2(nodes_fc=nodes_fc,
      rb_fc=rb_fc,
      lb_fc=lb_fc,
      overwrite_data=True)

# ------------------------------------------------------------------------
# TTools Step 3: Sample Elevation and Gradient
step3(nodes_fc=nodes_fc,
      searchCells=2,
      smooth_flag=True,
      z_raster=z_raster,
      z_units="Meters",
      block_size=10,
      overwrite_data=True)
# ------------------------------------------------------------------------
# TTools Step 4: Measure Topographic Shade Angles
step4(nodes_fc=nodes_fc,
      topo_directions=2,
      searchDistance_max_km=10,
      z_raster=z_raster,
      z_units="Meters",
      topo_fc=topo_samples_fc,
      block_size=10,
      overwrite_data=True)

# ------------------------------------------------------------------------
# TTools Step 5: Sample Landcover - Star Pattern
step5_star(nodes_fc=nodes_fc,
           trans_count=8,
           transsample_count=5,
           transsample_distance=8,
           zone_sample=False,
           heatsource8=False,
           lc_raster=lc_raster,
           lc_units="Meters",
           z_raster=z_raster,
           z_units="Meters",
           lc_point_fc=lc_samples_fc,
           block_size=10,
           overwrite_data=True)
```