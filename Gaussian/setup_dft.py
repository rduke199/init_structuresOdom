import os
import re
import json
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput


def xyz_to_gjf(xyz_fn, out_dir, functional, basis_set, omega=None, charge=0):
    """
    Convert an individually inputted xyz file to a gjf Gaussian input file
    """
    mol = Molecule.from_file(xyz_fn)
    mol_name = xyz_fn.split('/')[-1].split('.')[0]
    if omega is None:
        gau = GaussianInput(mol=mol, charge=charge, functional=functional, basis_set=basis_set,
                            route_parameters={'opt': ''},
                            link0_parameters={'%mem': '5GB', '%chk': '{}.chk'.format(mol_name)})
    else:
        gau = GaussianInput(mol=mol, charge=charge, functional=functional, basis_set=basis_set,
                            route_parameters={'iop(3/107={}, 3/108={})'.format(omega, omega): '', 'opt': ''},
                            link0_parameters={'%mem': '5GB', '%chk': '{}.chk'.format(mol_name)})
    gjf_file = gau.write_file('{}/{}.gjf'.format(out_dir, mol_name))
    return gjf_file


def get_run_folders(molecule_dir, out_dir, nflag=''):
    mol_name = molecule_dir.split('/')[-1].split('.')[0]
    log_files = [x for x in os.listdir(molecule_dir) if x.endswith('.log') and not x.endswith('deg.log')]
    gjf_files = [x for x in os.listdir(molecule_dir) if x.endswith('.gjf')]
    normal = None
    if len(log_files) == 1:
        with open(os.path.join(molecule_dir, log_files[0])) as fn:
            log_fn = fn.readlines()
        for line in log_fn:
            if re.search(' Optimized Parameters', line):
                last_line = log_fn[-1]
                if re.search('Normal termination', last_line):
                    normal = 1
                    break
            if re.search('Error termination', line):
                normal = 0
                print("Error in termination for {}".format(mol_name))
    else:
        normal = 0
        print("{} has {} log file(s).".format(mol_name, len(log_files)))
    if normal is None:
        if len(gjf_files) == 0:
            print("Error. {} does not contain a gjf file.".format(mol_name))
        else:
            print("{} may still be running".format(mol_name))
    if normal == 0:
        txt_file = os.path.join(out_dir, 'folders{}_to_run.txt'.format(nflag))
        with open(txt_file, 'a+') as fn:
            fn.write(molecule_dir + '\n')


def setup_a_folder(mol_name, molpath, xyz_file, dft_runs_path, omega=None, functional='LC-wHPBE', basis_set='TZVP',
                   charge=0):
    if not os.path.isdir(molpath): os.mkdir(molpath)
    try:
        # xyz_to_gjf(xyz_file, molpath, functional, basis_set, omega=omega, charge=charge)
        get_run_folders(molpath, dft_runs_path)  # , nflag=str(charge))
        print("Done setting up {} charge{} with {}/{}!".format(mol_name, charge, functional, basis_set))
    except:
        print("Error. Calculation for {} charge{} with {}/{} was not set up.".format(mol_name, charge, functional,
                                                                                     basis_set))


def main():
    home = os.getcwd()
    xyz_path = os.path.join(home, 'xyz/')
    dft_runs_path = os.path.join(home, 'dft_omega/')
    omega_file = os.path.join(home, 'omegas/master_omegas.json')
    with open(omega_file, 'r') as fn:
        omega_dict = json.load(fn)

    for mol in os.listdir(xyz_path):
        mol_name = mol.split('.')[0]
        xyz_file = os.path.join(xyz_path, mol)
        omega = omega_dict[mol[:-4]]
        # omega = None

        if 'cp' in mol_name:
            if 'whol' in mol_name:
                molpath = os.path.join(dft_runs_path, mol_name)
                setup_a_folder(mol_name, molpath, xyz_file, dft_runs_path, omega=omega, charge=2)
                cat1path = os.path.join(molpath, 'cation1/')
                setup_a_folder(mol_name, cat1path, xyz_file, dft_runs_path, omega=omega, charge=3)
                cat2path = os.path.join(molpath, 'cation2/')
                setup_a_folder(mol_name, cat2path, xyz_file, dft_runs_path, omega=omega, charge=4)
            else:
                molpath = os.path.join(dft_runs_path, mol_name)
                setup_a_folder(mol_name, molpath, xyz_file, dft_runs_path, omega=omega, charge=1)
                cat1path = os.path.join(molpath, 'cation1/')
                setup_a_folder(mol_name, cat1path, xyz_file, dft_runs_path, omega=omega, charge=2)
                cat2path = os.path.join(molpath, 'cation2/')
                setup_a_folder(mol_name, cat2path, xyz_file, dft_runs_path, omega=omega, charge=3)
        else:
            molpath = os.path.join(dft_runs_path, mol_name)
            setup_a_folder(mol_name, molpath, xyz_file, dft_runs_path, omega=omega, charge=0)
            cat1path = os.path.join(molpath, 'cation1/')
            setup_a_folder(mol_name, cat1path, xyz_file, dft_runs_path, omega=omega, charge=1)
            cat2path = os.path.join(molpath, 'cation2/')
            setup_a_folder(mol_name, cat2path, xyz_file, dft_runs_path, omega=omega, charge=2)


if __name__ == "__main__":
    main()
