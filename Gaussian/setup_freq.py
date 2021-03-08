import os
import json
from functions import setup_a_folder, find_ground_charge, get_in_out_paths


def main():
    home = os.getcwd()
    in_dir = os.path.join(home, 'dft_omega')
    freq_dir = os.path.join(home, 'freq/')
    omega_file = os.path.join(home, 'omegas/master_omegas_gaus.json')
    with open(omega_file, 'r') as fn:
        omega_dict = json.load(fn)

    mols = [m for m in os.listdir(in_dir) if m.startswith('mol')]
    for mol in mols:
        mol_name = mol.split('.')[0]
        ground_charge = find_ground_charge(mol_name)
        omega = omega_dict[mol]

        charges = {'neutral': ground_charge, 'cation1': ground_charge+1, 'cation2': ground_charge+2}
        for name in charges.keys():
            try:
                in_file, out_path = get_in_out_paths(mol, in_dir, freq_dir, cation=name)
                setup_a_folder(mol_name, in_file, out_path, freq_dir, omega=omega, charge=charges[name], calculation='freq')
            except UserWarning:
                pass


if __name__ == "__main__":
    main()
