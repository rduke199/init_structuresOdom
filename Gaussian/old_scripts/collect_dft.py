import os
from pymatgen.io.gaussian import GaussianOutput
from functions import write_json, runtime_from_log, write_master_json, find_ground_charge


def make_json(mol_dir, mol_name, json_file):
    """
    For a molecule directory, mol_dir, this finds the energy, xyz, charge, homo,
    lumo, homo_lumo_gap, and runtime and places those into  json_file.
    """
    # Gather info from log file in mol_dir
    try:
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('.log')][0]
        log_path = os.path.join(mol_dir, log_fn)
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
    # write json file if data was correctly gathered
    if normal == 1:
        json_data = {
            "molecule_name": mol_name,
            "charge": charge,
            "structure": xyz,
            "energy": energy,
            "homo": homo,
            "lumo": lumo,
            "homo_lumo_gap": homo_lumo_gap,
            "runtime": runtime
        }

        write_json(json_data, json_file)
        print("Data collected for {}".format(mol_name))


def main():
    # Establish paths that will be used
    home = os.getcwd()
    opt_runs_path = os.path.join(home, 'opt_Nomega/')
    out_home = os.path.join(home, 'jsons')
    if not os.path.isdir(out_home): os.mkdir(out_home)
    master_json_file = os.path.join(out_home, "master.json")

    # Iterate through mol directories and gather information
    mols = [m for m in os.listdir(opt_runs_path) if m.startswith('mol')]
    for mol in mols:
        # molpath = os.path.join(opt_runs_path, mol)
        # json_file = "{}/{}.json".format(out_home, mol)
        # make_json(molpath, mol, json_file)
        #
        # cat1path = os.path.join(molpath, 'cation1/')
        # json1_file = "{}/{}cat1.json".format(out_home, mol)
        # make_json(cat1path, mol + 'cat1', json1_file)
        #
        # cat2path = os.path.join(molpath, 'cation2/')
        # json2_file = "{}/{}cat2.json".format(out_home, mol)
        # make_json(cat2path, mol + 'cat2', json2_file)
        #
        mol_name = mol.split('.')[0]
        ground_charge = find_ground_charge(mol_name)
        charges = {'neutral': ground_charge, 'cation1': ground_charge + 1, 'cation2': ground_charge + 2}
        for name in charges.keys():
            try:
                mol_dir_path = os.path.join(opt_runs_path, mol, name+'/')
                json2_file = "{}/{}_{}.json".format(out_home, mol, name)
                make_json(mol_dir_path, mol_name + '_' + name, json2_file)
            except UserWarning:
                pass

    write_master_json(out_home, master_json_file)


if __name__ == "__main__":
    main()
