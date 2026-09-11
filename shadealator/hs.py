"""Setup and run Heat Source solar models using the TTools nodes feature class"""

import csv
import subprocess
import sys
from pathlib import Path
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait

from shadealator.utils import read_to_list, write_input

from ttools.geo_package import read_fc
from ttools.utils import message

# Name of Heat Source excel sheets
# copied from heat source v 9.1.0 at commit 763f4e0
sheetnames = {
    "controlfile": "Control Settings",
    "bcfile": "Boundary Conditions",
    "accretionfile": "Accretion Flow",
    "morphfile": "Morphology Data",
    "lcdatafile": "Land Cover Data",
    "lccodefile": "Land Cover Codes",
    "metsitesfile": "Meteorological Sites",
    "tribsitesfile": "Tributary Sites",
    "metfiles": "Meteorological Data",
    "tribfiles": "Tributary Data",
}

def setup_solar(nodes_fc, model_grouping_field, longsample, trans_count, transsample_count, transsample_distance,
                lccodes, canopy_data, heatsource8, lcsampmethod, emergent, modelstart, modelend,
                offset, dt, outputdt, model_dir, inputdir="inputs", outputdir="outputs", csv_mode=False,
                overwrite=False, usertxt=None, name=None):

    """Create Heat Source solar model from TTools nodes"""

    hs_exe = Path(Path(__file__).resolve().parent, "bin", "hs.exe")
    model_dir = Path(model_dir).expanduser().resolve()
    ext = ".csv" if csv_mode else ".xlsx"
    models_file = Path(model_dir, "model_list" + ext)

    nodeDict = read_fc(nodes_fc)
    groups = {}

    for nodeID, data in nodeDict.items():
        streamID = data[model_grouping_field]
        groups[streamID] = groups.get(streamID, 0) + 1

    folders = set()

    for streamID, node_count in groups.items():
        folder = str(streamID)

        # TODO Some error checking for names that can't be used as folder names

        folders.add(folder.casefold())

    control_file = "HeatSourceControl" + ext
    lccodes_data = read_to_list(Path(lccodes).expanduser().resolve())
    control_files = []
    model_rows = []

    for streamID, node_count in groups.items():
        msg = "Setting up solar model for {0}".format(streamID)
        message(msg)

        group_nodes = {nodeID: data for nodeID, data in nodeDict.items()
                       if data[model_grouping_field] == streamID}
        node_ids = sorted(group_nodes, reverse=True)

        run_dir = Path(model_dir, str(streamID))
        run_inputdir = Path(run_dir, inputdir).resolve()
        run_outputdir = Path(run_dir, outputdir).resolve()
        run_dir.mkdir(parents=True, exist_ok=True)
        run_inputdir.mkdir(parents=True, exist_ok=True)
        run_outputdir.mkdir(parents=True, exist_ok=True)

        # This makes a unique list of every file in the run and input dir so an existing file won't
        # get overwritten if overwrite = False
        existing_files = set(run_dir.iterdir()) | set(run_inputdir.iterdir())
        length = (node_count - 1) * longsample / 1000

        # setup the control file
        control_path = Path(run_dir, control_file)
        control_cmd = [
            str(hs_exe), "setup", "-cf", "-md", str(run_dir),
            "-set", "inputdir={0}".format(run_inputdir),
            "-set", "outputdir={0}".format(run_outputdir),
            "-set", "length={0}".format(length),
            "-set", "longsample={0}".format(longsample),
            "-set", "dx={0}".format(longsample),
            "-set", "outputkm=all",
            "-set", "datastart={0}".format(modelstart),
            "-set", "modelstart={0}".format(modelstart),
            "-set", "modelend={0}".format(modelend),
            "-set", "dataend={0}".format(modelend),
            "-set", "flushdays=0",
            "-set", "offset={0}".format(offset),
            "-set", "dt={0}".format(dt),
            "-set", "outputdt={0}".format(outputdt),
            "-set", "metsites=1",
            "-set", "tribsites=0",
            "-set", "morphfile=morphology{0}".format(ext),
            "-set", "lcdatafile=landcover{0}".format(ext),
            "-set", "lccodefile=lccodes{0}".format(ext),
            "-set", "trans_count={0}".format(trans_count),
            "-set", "transsample_count={0}".format(transsample_count),
            "-set", "transsample_distance={0}".format(transsample_distance),
            "-set", "emergent={0}".format(emergent),
            "-set", "canopy_data={0}".format(canopy_data),
            "-set", "lcsampmethod={0}".format(lcsampmethod),
            "-set", "heatsource8={0}".format(heatsource8),
            "-set", "usertxt={0}".format("" if usertxt is None else usertxt),
            "-set", "name={0}".format("" if name is None else name),
        ]

        if csv_mode:
            control_cmd.append("-csv")

        # Setup the model input files
        input_cmd = [str(hs_exe), "setup", "-mi", "-md", str(run_dir)]

        for cmd in (control_cmd, input_cmd):
            if overwrite:
                cmd.append("-o")
            try:
                result = subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
            except subprocess.CalledProcessError as exc:
                msg = "Heat Source setup failed for model {0!r}.\nModel directory: {1}\n{2}".format(
                    streamID, run_dir, exc.stderr)
                raise RuntimeError(msg) from exc
            sys.stderr.write(result.stderr)

        lookup_rows = []

        input_files = (("morphfile", "morphology" + ext), ("lcdatafile", "landcover" + ext))

        for file_key, filename in input_files:
            path = Path(run_inputdir, filename)
            data = read_to_list(path=path, sheetname=sheetnames[file_key])
            headers = data[0]
            rows = data[1:]

            if path in existing_files and not overwrite:
                if file_key == "morphfile":
                    lookup_rows = [[row[0], row[1], row[2], str(run_outputdir)] for row in rows]
                continue


            if file_key == "lcdatafile" and not heatsource8:
                topo_fields = ("TOPO_NE", "TOPO_E", "TOPO_SE", "TOPO_S", "TOPO_SW", "TOPO_W", "TOPO_NW", "TOPO_N")
                node_fields = group_nodes[node_ids[0]]

                # This looks for either topo3 or topo8 fields
                if not all(field in node_fields for field in topo_fields):
                    headers = headers[:5] + ["TOPO_W", "TOPO_S", "TOPO_E"] + headers[13:]
                    rows = [row[:5] + [None, None, None] + row[13:] for row in rows]

            fields = headers[3:]
            if file_key == "lcdatafile" and heatsource8:
                for i, field in enumerate(fields):
                    if field.startswith(("LC_", "ELE_")) and field != "LC_T0_S0":
                        prefix, tran, sample = field.split("_")
                        tran = ("NE", "E", "SE", "S", "SW", "W", "NW").index(tran) + 1
                        fields[i] = "{0}_T{1}_{2}".format(prefix, tran, sample)

            for row, nodeID in zip(rows, node_ids):
                row[0], row[1] = streamID, nodeID
                for i, field in enumerate(fields, start=3):
                    if file_key == "lcdatafile" or field == "ELEVATION":
                        row[i] = group_nodes[nodeID][field]
                if file_key == "morphfile":
                    lookup_rows.append([streamID, nodeID, row[2], str(run_outputdir)])

            write_input(path=path, headers=headers, rows=rows, sheetname=sheetnames[file_key], csv_mode=csv_mode)

        path = Path(run_dir, "HeatSource_Met_Sites" + ext)
        if overwrite or path not in existing_files:
            data = read_to_list(path, sheetname=sheetnames["metsitesfile"])
            data[1][2] = lookup_rows[0][2]

            write_input(path=path, headers=data[0], rows=data[1:],
                        sheetname=sheetnames["metsitesfile"], csv_mode=csv_mode)

        path = Path(run_inputdir, "metsite1" + ext)
        if overwrite or path not in existing_files:
            data = read_to_list(path=path, sheetname=sheetnames["metfiles"])

            for row in data[1:]:
                row[1] = 0

            write_input(path=path, headers=data[0], rows=data[1:],
                        sheetname=sheetnames["metfiles"], csv_mode=csv_mode)

        path = Path(run_inputdir, "lccodes" + ext)
        if overwrite or path not in existing_files:
            write_input(path=path, headers=lccodes_data[0], rows=lccodes_data[1:],
                        sheetname=sheetnames["lccodefile"], csv_mode=csv_mode)

        path = Path(run_outputdir, "node_lookup.csv")
        if overwrite or not path.exists():
            with path.open("w", newline="", encoding="utf-8") as output:
                writer = csv.writer(output)
                writer.writerow(["STREAM_ID", "NODE_ID", "STREAM_KM", "OUTPUTDIR"])
                writer.writerows(lookup_rows)

        control_files.append(str(control_path))
        model_rows.append([streamID, str(run_dir), ""])

        del group_nodes

    if overwrite or not models_file.exists():
        write_input(path=models_file, headers=["stream_id", "model_dir", "status"], rows=model_rows,
                    sheetname="Model List", csv_mode=csv_mode)

    return control_files


def run_solar(models_file, overwrite=False, concurrent_models=1):
    """Run Heat Source solar models"""

    hs_exe = Path(Path(__file__).resolve().parent, "bin", "hs.exe")
    models_file = Path(models_file).expanduser().resolve()
    csv_mode = models_file.suffix.lower() == ".csv"
    model_list = read_to_list(models_file, sheetname="Model List")
    model_list_headers = model_list[0]
    model_list_rows = model_list[1:]

    run_models = iter(model_row for model_row in model_list_rows if overwrite or model_row[2] != "complete")
    error = None
    with ThreadPoolExecutor(max_workers=concurrent_models) as executor:
        pending = {}
        while True:
            if error is None:
                for _ in range(concurrent_models - len(pending)):
                    model_row = next(run_models, None)
                    if model_row is None:
                        break
                    streamID, model_dir, status = model_row
                    msg = "Running solar model for {0}".format(streamID)
                    message(msg)
                    model_row[2] = "running"
                    write_input(path=models_file, headers=model_list_headers,
                                rows=model_list_rows, sheetname="Model List", csv_mode=csv_mode)
                    cmd = [str(hs_exe), "run", "-s", "-md", model_dir]
                    future = executor.submit(subprocess.run, cmd, check=True,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE, text=True)
                    pending[future] = model_row
            if not pending:
                break

            complete, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in complete:
                model_row = pending.pop(future)
                try:
                    result = future.result()
                except Exception as exc:
                    model_row[2] = "failed"
                    if error is None:
                        error = exc
                        error_model = model_row
                else:
                    model_row[2] = "complete"
                    msg = result.stdout
                    message(msg)
                    sys.stderr.write(result.stderr)
            write_input(path=models_file, headers=model_list_headers,
                        rows=model_list_rows, sheetname="Model List", csv_mode=csv_mode)

    if error is not None:
        if isinstance(error, subprocess.CalledProcessError):
            msg = "Heat Source {0!r} failed.\nModel directory: {1}\n{2}".format(error_model[0],
                                                                                error_model[1],
                                                                                error.stderr)
            raise RuntimeError(msg) from error
        raise error
