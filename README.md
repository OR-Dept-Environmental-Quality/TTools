TTools 
======
Current Version: TTools 9.3

## About

TTools is a collection of python scripts used to sample and assemble stream channel
and land cover geospatial data for input into the [Heat Source][1] model, Washington Department
of Ecology's [Shade][2] model, or the shade input file for [CE-QUAL-W2][3].

Requires
Python 3.11+

TTools can optionally be used in conjunction with ESRI ArcGIS Pro 3.3+, although it is not required.

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
    ├── arcpy_methods.py
    ├── geo_package.py
    ├── open_methods.py
    ├── step1.py
    ├── step2.py
    ├── step3.py
    ├── step4.py
    ├── step5.py
    └── utils.py
```

## Download TTools

1. Download the [TTools package][4].
2. Unzip the file and save the unpacked files somewhere on your computer. Keep all the files in the uppercase `TTools` directory together and organized the same as shown in the file structure diagram.

[4]: https://github.com/OR-Dept-Environmental-Quality/TTools/archive/refs/heads/master.zip

## Add TTools to ArcPro

1. Enable the ArcGIS Pro Spatial Analyst extension.
2. Load the Toolbox in ArcGIS Pro. Open or create a new ArcGIS Pro project from the Catalog pane. If it is not visible, go to: View → Catalog Pane. This will show the ArcGIS Pro Catalog Pane and show Toolboxes
3. Add the toolbox to ArcGIS Pro. From the main menu ribbon in ArcGIS Pro, right-click Toolbox → Add Toolbox → Browse to the toolbox file (`TTools.pyt`)  → Select the Toolbox file `TTools.pyt` → Click OK
4. From the Catalog pane, click on the arrow next to Toolboxes to show all available Toolboxes.
5. Click on the arrow next to the `TTools.pyt` toolbox in the Catalog to display the tools available in the toolbox.

The XML files shown in the file structure diagram contain the help documentation for each tool. The help information is accessible within ArcGIS Pro. To view the help, hover the cursor next to the start of each field name in a tool and click the blue icon that appears.

See the TTools user guide or [ESRI documentation][5] for additional instructions.

[5]: https://doc.esri.com/en/arcgis-pro/latest/help/projects/connect-to-a-toolbox.html#A9E

## Use TTools from a Python script

TTools can be run from a Python script instead of using the ArcPro toolbox. There are two ways to run scripts:

1. Open source methods, the default. Open source python packages are used for TTools geoprocessing instead of ESRI's arcpy. The relevant open source geospatial packages must be installed including `geopandas`, `rasterio`, `fiona`, `shapely`, `pyproj`, `pandas`, and `numpy`.

2. ArcPy methods. This uses ESRI's arcpy package for geospatial processing. There must be a valid ArcPro license with spatial analyst extension. The Python environment running the script must also be able to import arcpy. ESRI recommends running stand-along scripts that need arcpy by calling [propy.bat or proenv.bat][6]. Point the script to the `ttools` folder as shown in the example below, import ttools, and set the geopackage to `use_arcpy=True`.

Documentation for each TTools step can be found in the user guide, in each ttools Python file e.g. `step1.py`, or using (`help(step1)`).

[6]: https://pro.arcgis.com/en/pro-app/3.6/arcpy/get-started/using-conda-with-arcgis-pro.htm

```python
# This example is also saved to `TTools_Example_Script.py` located in the repository root.

import os
import sys

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
# Add TTools to python path
sys.path.insert(0, os.path.dirname(ttools_package_folder))

# Import ttools
from ttools.geo_package import set_geo_package
from ttools import step1, step2, step3, step4, step5_ortho, step5_star

# Set arcpy as the geoprocessing package, 
# use_arcpy=False if using open source.
set_geo_package(use_arcpy=True)

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
# Use topo_directions=1 for Heat Source 6-8 and WA Ecology shade model
# Use topo_directions=2 for Heat Source 9
# Use topo_directions=3 for CE-QUAL-W2
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
# README --
# Normally step 5 is run once using either the star pattern or orthogonal method.
# Both methods are shown in this script as an example.
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

# ------------------------------------------------------------------------
# TTools Step 5: Sample Landcover - Orthogonal Method
# README --
# Normally step 5 is run once using either the star pattern or orthogonal method.
# Both methods are shown here as an example.
step5_ortho(nodes_fc=nodes_fc,
            start_bank=True,
            transsample_count=9,
            transsample_distance=8,
            lc_raster=lc_raster,
            lc_units="Meters",
            z_raster=z_raster,
            z_units="Meters",
            lc_point_fc=lc_samples_fc,
            block_size=10,
            overwrite_data=True)
```