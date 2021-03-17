import os
import re
import json
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianOutput

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

def write_json(data, filename):
    with open(filename,'w') as f:
        json.dump(data, f, indent=2)

def MakeJSON(mol_dir,json_file,master_json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    mol_name = str(mol_dir.split('/')[-1])
    try:
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('.log')][0]
        log_path = os.path.join(mol_dir,log_fn)
        mol = GaussianOutput(log_path)
        if mol.properly_terminated:
            num_electrons = mol.electrons[0]
            eigens = list(mol.eigenvalues.values())[0]
            homo = eigens[num_electrons - 1] * 27.2114
            lumo = eigens[num_electrons] * 27.2114
            homo_lumo_gap = lumo - homo
            energy = mol.final_energy * 27.2114
            structure = mol.final_structure
            mol_xyz = structure.to(fmt="xyz")
            runtime = runtime_from_log(log_path)
            functional = mol.functional
            basis_set = mol.basis_set
            normal = 1
        else:
            normal = 0
            print("Error. Calculation did not properly terminate for {}".format(mol_name))
            print(mol_dir)
    except:
        normal = 0
        print("Error. Data NOT collected for {}. Energy calculations may not have finished!".format(mol_name))
        print(mol_dir)

    try:
        catpath = os.path.join(mol_dir, 'cation/')
        log_fn = [x for x in os.listdir(catpath) if x.endswith('.log')][0]
        log_path = os.path.join(catpath,log_fn)
        cat = GaussianOutput(log_path)
        if cat.properly_terminated:
            cat_energy = cat.final_energy * 27.2114
            cnormal = 1
        else:
            cnormal = 0
            print("Error. Calculation did not properly terminate for {} cation with {} and {}".format(mol_name,functional,basis_set))
            print(mol_dir)
    except:
        cnormal = 0
        print("Error. Data NOT collected for {} cation. Energy calculations for cation may not have finished!".format(mol_name))
        print(mol_dir)
    if normal == 1 and cnormal == 1:
        json_data = {
        "molecule_name" : mol_name,
        "energy" : energy,
        "cation_energy" : cat_energy,
        "runtime" : runtime,
        "functional" : functional,
        "basis_set" : basis_set,
        "homo" : homo,
        "lumo" : lumo,
        "structure" : mol_xyz,
        "homo_lumo_gap" : homo_lumo_gap
        }

        write_json(json_data, json_file)
        print("Data collected for {} with {} and {}".format(mol_name,functional,basis_set))

def WriteMaster(json_dir,master_json_file):
    data = {}
    data['molecules'] = []
    for f in os.listdir(json_dir):
        fpath = os.path.join(json_dir,f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
        data['molecules'].append(mol_info)
    write_json(data, master_json_file)

home = os.getcwd()
opt_runs_path = os.path.join(home,'opt_runs/')
out_home = os.path.join(home,'jsons/')
funcitonals = ['LC-wHPBE','wB97XD']
basis_sets = ['cc-pVDZ', '6-31+G**', 'TZVP', 'aug-cc-pVTZ']
master_json_file =os.path.join(out_home,"master.json")

for f in funcitonals:
    fpath = os.path.join(opt_runs_path,f)
    for b in basis_sets:
        bpath = os.path.join(fpath,b)
        for mol in os.listdir(bpath):
            molpath = os.path.join(bpath, mol)
            json_file = "{}/{}_{}_{}.json".format(out_home,mol,f,b)
            MakeJSON(molpath,json_file,master_json_file)
WriteMaster(out_home,master_json_file)
