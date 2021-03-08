import os
import json
from functions import setup_a_folder, find_ground_charge, get_in_out_paths


def main():
    home = os.getcwd()

    # # NO OMEGA
    # in_file_path = os.path.join(home, 'xyz/')
    # dft_dir = os.path.join(home, 'dft_Nomega/')

    # OMEGA
    in_dir = os.path.join(home, 'wtuning/')
    dft_dir = os.path.join(home, 'dft_omega/')
    omega_file = os.path.join(home, 'omegas/master_omegas_gaus.json')
    with open(omega_file, 'r') as fn:
        omega_dict = json.load(fn)

    mols = [m for m in os.listdir(in_dir) if m.startswith('mol')]
    for mol in mols:
        mol_name = mol.split('.')[0]
        ground_charge = find_ground_charge(mol_name)
        if omega_dict is not None:
            file_name = [f for f in os.listdir(os.path.join(in_dir, mol)) if f.endswith('opt_0.log')][0]
            in_file = os.path.join(in_dir, mol, file_name)
            omega = omega_dict[mol]
        else:
            omega = None
            in_file = os.path.join(in_dir, mol)

        charges = {'neutral': ground_charge, 'cation1': ground_charge+1, 'cation2': ground_charge+2}
        for name in charges.keys():
            try:
                _, out_path = get_in_out_paths(mol, in_dir, dft_dir, cation=name)
                setup_a_folder(mol_name, in_file, out_path, dft_dir, omega=omega, charge=charges[name])
            except UserWarning:
                pass


if __name__ == "__main__":
    main()
