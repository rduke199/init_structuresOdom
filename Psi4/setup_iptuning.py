import os
import re
from pymatgen.core import Molecule


def write_xyz(in_file, out_dir):
    """
    Write xyz file.
    """
    mol = Molecule.from_file(in_file)
    fout = os.path.join(out_dir, in_file.split('/')[-1])
    mol.to(filename=fout)


def write_ip_fitting(out_dir, charge=0, multiplicity=1):
    """
    Write python setup file for ipfitting job
    """
    fout = os.path.join(out_dir, 'ipfitting.py')
    with open(fout, 'w+') as wt:
        wt.write("""import json
from sys import argv
import psi4
from psi4.driver.frac import ip_fitting

mol_file = argv[1]
molecule_name = (mol_file.split('/')[-1]).split('.')[0]
with open(mol_file,'r') as f:
    mol = psi4.core.Molecule.from_string(f.read(), dtype='xyz')
mol.set_molecular_charge(""" + str(charge) + """)
mol.set_multiplicity(""" + str(multiplicity) + """)

psi4.set_memory('2 GB')
psi4.set_num_threads(2)
psi4.set_output_file(molecule_name + '_ip_fitting.dat', False)
psi4.set_options({'basis': 'def2-TZVP','scf_type': 'df'})
omega = ip_fitting('LRC-wPBEH', 0.1, 2.0,molecule=mol)
json_data = {"molecule_name" : molecule_name, "omega" : omega}
json_file = ("omega_{}.txt".format(molecule_name))
with open(json_file,'w') as f:
    json.dump(json_data, f, indent=2)""")


def gather_ipfitting_omega(out_file_path):
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


def get_run_folders(molecule_dir, out_dir, nflag=''):
    mol_name = molecule_dir.split('/')[-1].split('.')[0]
    out_files = [x for x in os.listdir(molecule_dir) if x.endswith('_ip_fitting.dat')]
    xyz_files = [x for x in os.listdir(molecule_dir) if x.endswith('.xyz')]
    if len(out_files) == 1:
        out_path = os.path.join(molecule_dir, out_files[0])
        with open(out_path) as fn:
            out_data = fn.readlines()
            normal_line = out_data[-5]
        if re.search("IP Fitting Converged", normal_line):
            normal = True
        else:
            normal = False
            print("Error. {} job did not have a normal termination.".format(mol_name))
    else:
        normal = False
        print("{} has {} ip_fitting.dat file(s).".format(mol_name, len(out_files)))
    if len(xyz_files) == 0:
        normal = False
        print("Error. {} does not contain a xyz file.".format(mol_name))
    else:
        pass
    if not normal:
        txt_file = os.path.join(out_dir, 'folders{}_to_run.txt'.format(nflag))
        with open(txt_file, "a+") as fn:
            fn.write(molecule_dir + "\n")


def implement_setup(molpath, xyz_file, ipfitting_path, charge=0, multiplicity=1):
    mol_name = (molpath.split('/')[-1]).split('.')[0]
    if not os.path.isdir(molpath): os.mkdir(molpath)
    try:
        # write_xyz(xyz_file, molpath)
        write_ip_fitting(molpath, charge, multiplicity)
        get_run_folders(molpath, ipfitting_path)
        print("Done setting up ip fitting for {}.".format(mol_name))
    except:
        print("Error setting up ip fitting for {}!".format(mol_name))


def main():
    home = os.getcwd()
    xyz_path = os.path.join(home, 'xyz/')
    ipfitting_path = os.path.join(home, 'ipfitting/')

    for mol in os.listdir(xyz_path):
        mol_name = mol.split('.')[0]
        xyz_file = os.path.join(xyz_path, mol)
        molpath = os.path.join(ipfitting_path, mol_name)

        if 'cp' in mol_name:
            if 'whol' in mol_name:
                implement_setup(molpath, xyz_file, ipfitting_path, charge=2, multiplicity=1)
            else:
                implement_setup(molpath, xyz_file, ipfitting_path, charge=1, multiplicity=2)
        else:
            implement_setup(molpath, xyz_file, ipfitting_path, charge=0, multiplicity=1)


if __name__ == "__main__":
    main()
