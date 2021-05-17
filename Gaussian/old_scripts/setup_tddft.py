from Gaussian.old_scripts.functions import *


def main():
    home = os.getcwd()

    calc_dir = os.path.join(home, 'opt_omega')
    omega_file = os.path.join(home, "omegas/master_omegas_gaus.json")
    with open(omega_file, 'r') as fn:
        omega_dict = json.load(fn)

    mols = [m for m in os.listdir(calc_dir) if m.startswith('mol')]
    for mol in mols:
        mol_name = mol.split('.')[0]
        ground_charge = find_ground_charge(mol_name)
        omega = omega_dict[mol]

        charges = {'neutral': ground_charge, 'cation1': ground_charge+1, 'cation2': ground_charge+2}
        for ch_name in charges.keys():
            try:
                in_file = os.path.join(calc_dir, mol, ch_name, 'opt.log')
                out_path = get_out_paths(mol, calc_dir, cation=ch_name, subtype='tddft')
                setup_a_folder(mol_name, in_file, out_path, calc_dir, omega=omega, charge=charges[ch_name],
                               calculation='tddft')
                get_run_folders(out_path, calc_dir)
            except:
                pass


if __name__ == "__main__":
    main()
