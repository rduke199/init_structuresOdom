import os
import re
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput


def generate_gjf(in_fn, out_dir, functional, basis_set, calculation, omega=None, charge=0):
    """
    Convert an individually inputted xyz file to a gjf Gaussian input file
    """
    mol = Molecule.from_file(in_fn)
    mol_name = in_fn.split('/')[-1][:14]
    if omega is None:
        gau = GaussianInput(mol=mol, charge=charge, functional=functional, basis_set=basis_set,
                            route_parameters={calculation: '', 'SCF': '(MaxCycle=250)'},
                            link0_parameters={'%mem': '5GB', '%chk': '{}.chk'.format(mol_name)})
    else:
        gau = GaussianInput(mol=mol, charge=charge, functional=functional, basis_set=basis_set,
                            route_parameters={'iop(3/107={}, 3/108={})'.format(omega, omega): '', calculation: '',
                                              'SCF': '(MaxCycle=250)'},
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


def find_ground_charge(mol_name):
    if 'cp' in mol_name:
        if 'whol' in mol_name:
            ground_charge = 2
        else:
            ground_charge = 1
    else:
        ground_charge = 0
    return ground_charge


def get_in_out_paths(mol, in_dir, out_dir, cation='', log_end='.log'):
    mol_name = mol.split('.')[0]
    try:
        file_name = [f for f in os.listdir(os.path.join(in_dir, mol, cation)) if f.endswith(log_end)][0]
        in_file = os.path.join(in_dir, mol_name, cation, file_name)

        mol_out_path = os.path.join(out_dir, mol_name)
        if not os.path.isdir(mol_out_path): os.mkdir(mol_out_path)
        out_path = os.path.join(mol_out_path, cation)
        if not os.path.isdir(out_path): os.mkdir(out_path)

        return in_file, out_path
    except IndexError:
        raise UserWarning('No opt_0.log file found for {}'.format(mol_name))


def get_out_paths(mol, out_dir, cation=''):
    mol_name = mol.split('.')[0]
    mol_out_path = os.path.join(out_dir, mol_name)
    if not os.path.isdir(mol_out_path): os.mkdir(mol_out_path)
    out_path = os.path.join(mol_out_path, cation)
    if not os.path.isdir(out_path): os.mkdir(out_path)
    return out_path


def setup_a_folder(mol_name, in_file, out_path, runs_folder, omega=None, functional='LC-wHPBE', basis_set='TZVP',
                   charge=0, calculation='opt'):
    if not os.path.isdir(out_path): os.mkdir(out_path)
    try:
        # generate_gjf(in_file, out_path, functional, basis_set, calculation, omega=omega, charge=charge)
        get_run_folders(out_path, runs_folder)  # , nflag=str(charge))
        print("Done setting up {} charge{} with {}/{}!".format(mol_name, charge, functional, basis_set))
    except:
        print("Error. Calculation for {} charge{} with {}/{} was not set up.".format(mol_name, charge, functional,
                                                                                     basis_set))
