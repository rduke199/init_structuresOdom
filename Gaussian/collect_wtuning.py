import os
import re
import json


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def runtime_from_log(logfile):
    """
    Collects runtime in core hours from a logfile
    """
    time_patt = re.compile(r"\d+\.d+|\d+")
    time_data = []
    with open(logfile, "r") as f:
        line = f.readline()
        while line != "":
            if re.match("Job cpu time", line.strip()):
                time_data.extend(time_patt.findall(line))
            line = f.readline()
    if time_data:
        time_data = [float(time) for time in time_data]
        runtime = (time_data[0] * 86400 + time_data[1] * 3600 + time_data[2] * 60 + time_data[3]) / 3600
    else:
        runtime = 0
    return round(runtime, 3)


def make_json(mol_dir, mol_name, json_file):
    """
    For a molecule directory, mol_dir, this finds the molecule name, omegas
    value, and runtiime and places those into  json_file.
    """
    try:
        # Gather info from log file in mol_dir
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('output.log')][0]
        log_path = os.path.join(mol_dir, log_fn)
        runtime = runtime_from_log(log_path)
        with open(log_path, 'r') as fn:
            w_data = fn.readlines()[-2].split()[1]
            w = "0{}".format(w_data.split('.')[1])
        # Write json
        json_data = {"molecule_name": mol_name, "omega": w, "runtime": runtime}
        write_json(json_data, json_file)
        print("Data collected for {}".format(mol_name))

    except:
        print("Error. Omega value NOT collected for {}. wtuning may not have finished!".format(mol_name))


def write_master_json(json_dir, master_json_file):
    """
    Write a master data file from all json files in json_dir
    """
    data = {}
    for f in os.listdir(json_dir):
        fpath = os.path.join(json_dir, f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
            data[mol_info["molecule_name"]] = mol_info["omega"]
    write_json(data, master_json_file)


def main():
    # Establish paths that will be used
    home = os.getcwd()
    opt_runs_path = os.path.join(home, 'wtuning/')
    out_home = os.path.join(home, 'omegas')
    if not os.path.isdir(out_home): os.mkdir(out_home)
    master_json_file = os.path.join(out_home, "master_omegas.json")

    # Iterate through mol directories and gather information
    mols = [m for m in os.listdir(opt_runs_path) if m.startswith('mol')]
    for mol in mols:
        molpath = os.path.join(opt_runs_path, mol)
        json_file = "{}/omega_{}.json".format(out_home, mol)
        make_json(molpath, mol, json_file)

    write_master_json(out_home, master_json_file)


if __name__ == "__main__":
    main()
