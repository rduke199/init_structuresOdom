from functions import *


def main():
    home = os.getcwd()
    calc_dir = os.path.join(home, 'opt_omega')
    omega_file = os.path.join(home, 'omegas/master_omegas_gaus.json')
    with open(omega_file, 'r') as fn:
        omega_dict = json.load(fn)

    mols = [m for m in os.listdir(calc_dir) if m.startswith('mol')]
    for mol in mols:
        mol_name = mol.split('.')[0]
        ground_charge = find_ground_charge(mol_name)
        omega = omega_dict[mol]

        charges = {'neutral': ground_charge, 'cation1': ground_charge+1, 'cation2': ground_charge+2}
        for charge_type in charges.keys():
            try:
                mol_path = os.path.join(calc_dir, mol_name, charge_type)
                opt_log_file = os.path.join(mol_path, mol_name + '.log')
                setup_a_folder(mol_name, opt_log_file, mol_path, calc_dir, omega=omega, charge=charges[charge_type])
                # get_run_folders(mol_path, calc_dir)
            except UserWarning:
                pass


if __name__ == "__main__":
    main()
