import os
from pymatgen.io.gaussian import GaussianOutput
from functions import write_json, runtime_from_log, write_master_json, find_ground_charge


def make_json(mol_dir, mol_name, json_file):
    """
    For a molecule directory, mol_dir, this finds the molecule name, frequencies
    value, and runtiime and places those into  json_file.
    """
    log_files = [x for x in os.listdir(mol_dir) if x.endswith('.log')]
    if len(log_files) > 0:
        log_fn = log_files[0]
        log_path = os.path.join(mol_dir, log_fn)
        mol = GaussianOutput(log_path)
        runtime = runtime_from_log(log_path)
        if mol.properly_terminated:
            charge = mol.charge
            try:
                excitations = mol.read_excitation_energies()
            except IndexError:
                excitations = None


            # Write json
            json_data = {"molecule_name": mol_name, "excitations": excitations, "charge": charge, "runtime": runtime}
            write_json(json_data, json_file)
            print("Data collected for {}".format(mol_name))
        else:
            print("Error. Calculation did not properly terminate for {}".format(mol_name))
    else:
        print("No log file for ", mol_name)


def main():
    # Establish paths that will be used
    home = os.getcwd()
    opt_runs_path = os.path.join(home, 'opt_omega/')
    out_home = os.path.join(home, 'jsons_tddft')
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
                mol_dir_path = os.path.join(opt_runs_path, mol, name, 'tddft/')
                json2_file = "{}/{}_{}.json".format(out_home, mol, name)
                make_json(mol_dir_path, mol_name + '_' + name, json2_file)
            except UserWarning:
                pass

    write_master_json(out_home, master_json_file, prop='excitations')


if __name__ == "__main__":
    main()
