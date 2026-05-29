"""
If you are using TTools from within the ESRI ArcPro toolbox, follow the TTools ArcPro instructions.
You don't need to use this script.

If you want to run ttools form a python script, this script shows an example of how to do that.

There is no install required as long as the ttools folder is saved somewhere on your computer, and you are using
the python environment with arcpy installed.
"""

import os
import sys
import arcpy

# ttools_package_folder points to the ttools package folder (lower case).
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
# This allows python to import ttools.
# Don't edit these lines.
sys.path.insert(0, os.path.dirname(ttools_package_folder))

from ttools import step1, step2, step3, step4, step5_ortho, step5_star

# ------------------------------------------------------------------------
# TTools Step 1: Create Stream Nodes
step1(
    streamline_fc=streamline_fc,
    sid_field="NAME",
    node_dx=50,
    cont_stream_km=False,
    nodes_fc=nodes_fc,
    checkDirection=True,
    z_raster=z_raster,
)

# ------------------------------------------------------------------------
# TTools Step 2: Measure Channel Widths
step2(
    nodes_fc=nodes_fc,
    rb_fc=rb_fc,
    lb_fc=lb_fc,
    overwrite_data=True,
)

# ------------------------------------------------------------------------
# TTools Step 3: Sample Elevation and Gradient
step3(
    nodes_fc=nodes_fc,
    searchCells=2,
    smooth_flag=True,
    z_raster=z_raster,
    z_units="Meters",
    block_size=10,
    overwrite_data=True,
)
# ------------------------------------------------------------------------
# TTools Step 4: Measure Topographic Shade Angles
step4(
    nodes_fc=nodes_fc,
    topo_directions=2,
    searchDistance_max_km=10,
    z_raster=z_raster,
    z_units="Meters",
    topo_fc=topo_samples_fc,
    block_size=10,
    overwrite_data=True,
)

# ------------------------------------------------------------------------
# TTools Step 5: Sample Landcover - Star Pattern
# Normally step 5 is run once using either the star pattern or orthogonal method.
# Both are shown in this script as an example.
step5_star(
    nodes_fc=nodes_fc,
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
    overwrite_data=True,
)

# ------------------------------------------------------------------------
# TTools Step 5: Sample Landcover - Orthogonal Method
# Normally step 5 is run once using either the star pattern or orthogonal method.
# Both are shown in this script as an example.
step5_ortho(
    nodes_fc=nodes_fc,
    start_bank=True,
    transsample_count=9,
    transsample_distance=8,
    lc_raster=lc_raster,
    lc_units="Meters",
    z_raster=z_raster,
    z_units="Meters",
    lc_point_fc=lc_samples_fc,
    block_size=10,
    overwrite_data=True,
)
