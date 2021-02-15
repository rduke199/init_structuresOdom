import os
import re
import json


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def gather_energy(out_file_path):
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
        energy = gather_energy(out_path)
        json_data = {
            "molecule_name": mol_name,
            "energy": energy,
            }
        write_json(json_data, json_file)
        print("Data collected for {}".format(mol_name))
    except UserWarning:
        print("Error. Data NOT collected for {}. Energy calculations may not have finished!".format(mol_name))


def MakeJSON(mol_dir,mol_name,json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    try:
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('.log')][0]
        log_path = os.path.join(mol_dir,log_fn)
        mol = GaussianOutput(log_path)
        if mol.properly_terminated:
            structure = mol.final_structure
            xyz = structure.to(fmt="xyz")
            charge = mol.charge
            energy = mol.final_energy * 27.2114
            num_electrons = mol.electrons[0]
            eigens = list(mol.eigenvalues.values())[0]
            homo = eigens[num_electrons - 1] * 27.2114
            lumo = eigens[num_electrons] * 27.2114
            homo_lumo_gap = lumo - homo
            runtime = runtime_from_log(log_path)
            normal = 1
        else:
            normal = 0
            print("Error. Calculation did not properly terminate for {}".format(mol_name))
    except:
        normal = 0
        print("Error. Data NOT collected for {}. Energy calculations may not have finished!".format(mol_name))
    if normal == 1:
        json_data = {
        "molecule_name" : mol_name,
        "charge" : charge,
        "structure" : xyz,
        "energy" : energy,
        "homo" : homo,
        "lumo" : lumo,
        "homo_lumo_gap" : homo_lumo_gap,
        "runtime" : runtime
        }

        write_json(json_data, json_file)
        print("Data collected for {}".format(mol_name))


def write_master_json(json_dir, master_json_file):
    data = {'molecules': []}
    for f in os.listdir(json_dir):
        fpath = os.path.join(json_dir, f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
        data['molecules'].append(mol_info)
    write_json(data, master_json_file)

home = os.getcwd()
dft_runs_path = os.path.join(home,'dft_Nomega/')
out_home = os.path.join(home,'jsons')
if not os.path.isdir(out_home): os.mkdir(out_home)
master_json_file =os.path.join(out_home,"master.json")

mols = [m for m in os.listdir(dft_runs_path) if m.startswith('mol')]
for mol in mols:
    molpath = os.path.join(dft_runs_path, mol)
    json_file = "{}/{}.json".format(out_home,mol)
    MakeJSON(molpath,mol,json_file)

    cat1path = os.path.join(molpath,'cation1/')
    json1_file = "{}/{}cat1.json".format(out_home,mol)
    MakeJSON(cat1path,mol+'cat1',json1_file)

    cat2path = os.path.join(molpath,'cation2/')
    json2_file = "{}/{}cat2.json".format(out_home,mol)
    MakeJSON(cat2path,mol+'cat2',json2_file)

WriteMaster(out_home,master_json_file)
