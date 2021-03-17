import os
import re
import json


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def check_normal_optimization(out_file_path):
    with open(out_file_path) as fn:
        data = fn.readlines()
    for line in data:
        if re.search('Optimization is complete', line):
            return True
    print("Optimization did not finish for {}".format(out_file_path))
    return False


def gather_energy(out_file_path):
    with open(out_file_path) as fn:
        data = fn.readlines()
    for line in data:
        if re.search('Final energy is', line):
            energy = line.split(" ")[-1]
            return energy


def make_json_file(mol_dir, mol_name, json_file_path):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    try:
        out_fn = [x for x in os.listdir(mol_dir) if x.endswith('geometry_optimization.dat')][0]
        out_path = os.path.join(mol_dir, out_fn)
        if check_normal_optimization(out_path):
            energy = gather_energy(out_path)
            json_data = {
                "molecule_name": mol_name,
                "energy": energy,
            }
            write_json(json_data, json_file_path)
            print("Data collected for {}".format(mol_name))
    except IndexError:
        print("Error. No geometry_optimization.dat file found for {}".format(mol_name))


def write_master_json(json_dir, master_json_file_path):
    data = {'molecules': []}
    for f in os.listdir(json_dir):
        fpath = os.path.join(json_dir, f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
        data['molecules'].append(mol_info)
    write_json(data, master_json_file_path)


home = os.getcwd()
opt_runs_path = os.path.join(home, 'optimization/')
out_home = os.path.join(home, 'jsons_opt')
if not os.path.isdir(out_home): os.mkdir(out_home)
master_json_file = os.path.join(out_home, "master.json")

mols = [m for m in os.listdir(opt_runs_path) if m.startswith('mol')]
for mol in mols:
    molpath = os.path.join(opt_runs_path, mol)
    json_file = "{}/{}.json".format(out_home, mol)
    make_json_file(molpath, mol, json_file)

    cat1path = os.path.join(molpath, 'cation1/')
    json1_file = "{}/{}cat1.json".format(out_home, mol)
    make_json_file(cat1path, mol + 'cat1', json1_file)

    cat2path = os.path.join(molpath, 'cation2/')
    json2_file = "{}/{}cat2.json".format(out_home, mol)
    make_json_file(cat2path, mol + 'cat2', json2_file)

write_master_json(out_home, master_json_file)
