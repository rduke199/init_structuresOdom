import os
from pymatgen.io.gaussian import GaussianOutput
from functions import write_json, runtime_from_log, write_master_json, find_ground_charge


def make_json(mol_dir, mol_name, json_file):
    """
    For a molecule directory, mol_dir, this finds the molecule name, frequencies
    value, and runtiime and places those into  json_file.
    """
    try:
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('.log')][0]
        log_path = os.path.join(mol_dir, log_fn)
        mol = GaussianOutput(log_path)
        runtime = runtime_from_log(log_path)
        if mol.properly_terminated:
            frequencies = [f["frequency"] for f in mol.frequencies[0]]
            charge = mol.charge
        else:
            print("Error. Calculation did not properly terminate for {}".format(mol_name))

        # Write json
        json_data = {"molecule_name": mol_name, "frequencies": frequencies, "charge": charge, "runtime": runtime}
        write_json(json_data, json_file)
        print("Data collected for {}".format(mol_name))

    except:
        print("Error. Frequencies NOT collected for {}. Frequency calculations may not have finished!".format(mol_name))


def main():
    # Establish paths that will be used
    home = os.getcwd()
    opt_runs_path = os.path.join(home, 'opt_omega/')
    out_home = os.path.join(home, 'jsons_freq')
    if not os.path.isdir(out_home): os.mkdir(out_home)
    master_json_file = os.path.join(out_home, "master_omegas.json")

    # Iterate through mol directories and gather information
    mols = [m for m in os.listdir(opt_runs_path) if m.startswith('mol')]
    for mol in mols:
        # molpath = os.path.join(opt_runs_path, mol)
        # json_file = "{}/omega_{}.json".format(out_home, mol)
        # make_json(molpath, mol, json_file)
        mol_name = mol.split('.')[0]
        ground_charge = find_ground_charge(mol_name)
        charges = {'neutral': ground_charge, 'cation1': ground_charge + 1, 'cation2': ground_charge + 2}
        for name in charges.keys():
            try:
                mol_dir_path = os.path.join(opt_runs_path, mol, name + '/')
                json2_file = "{}/{}_{}.json".format(out_home, mol, name)
                make_json(mol_dir_path, mol_name + '_' + name, json2_file)
            except UserWarning:
                pass

    write_master_json(out_home, master_json_file, prop='frequencies')


if __name__ == "__main__":
    main()
