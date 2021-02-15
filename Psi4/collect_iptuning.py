import os
import re
import json


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def gather_ipfitting_omega(out_file_path):
    molecule_name = out_file_path.split('/')[-1]
    with open(out_file_path) as fn:
        out_data = fn.readlines()
        normal_line = out_data[-5]
        if re.search("IP Fitting Converged", normal_line):
            omega_line = out_data[-4]
            omega = omega_line.split(" ")[-1]
            return omega
        else:
            raise UserWarning("Termination error in ip fitting for {}".format(molecule_name))


def make_json_file(mol_dir, mol_name, json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    try:
        out_fn = [x for x in os.listdir(mol_dir) if x.endswith('ip_fitting.dat')][0]
        out_path = os.path.join(mol_dir, out_fn)
        omega = gather_ipfitting_omega(out_path)
        json_data = {
            "molecule_name": mol_name,
            "omega": omega,
            }
        write_json(json_data, json_file)
        print("Data collected for {}".format(mol_name))
    except UserWarning:
        print("Termination error in ip fitting for {}.".format(mol_name))
    except IndexError:
        print("Omega line in out file not found. IP fitting for {} may not have begun.".format(mol_name))


def write_master_json(json_dir, master_json_file):
    """
    Write a master data file from all json files in json_dir
    """
    data = {}
    mols = [m for m in os.listdir(json_dir) if m.startswith("omega")]
    for f in mols:
        fpath = os.path.join(json_dir, f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
            data[mol_info["molecule_name"]] = mol_info["omega"]
    write_json(data, master_json_file)


def main():
    # Establish paths that will be used
    home = os.getcwd()
    ipf_runs_path = os.path.join(home, 'ipfitting/')
    out_home = os.path.join(home, 'omegas')
    if not os.path.isdir(out_home): os.mkdir(out_home)
    master_json_file = os.path.join(out_home, "master_omegas.json")

    # Iterate through mol directories and gather information
    mols = [m for m in os.listdir(ipf_runs_path) if m.startswith('mol')]
    for mol in mols:
        molpath = os.path.join(ipf_runs_path, mol)
        json_file = "{}/omega_{}.json".format(out_home, mol)
        make_json_file(molpath, mol, json_file)

    write_master_json(out_home, master_json_file)


if __name__ == "__main__":
    main()
