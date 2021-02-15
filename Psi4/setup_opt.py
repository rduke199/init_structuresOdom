import os
import re
import json
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

def WriteXYZ(in_file, out_dir):
    """
    Write xyz file.
    """
    mol = Molecule.from_file(in_file)
    fout = os.path.join(out_dir, in_file.split('/')[-1])
    mol.to(filename=fout)

def WriteOpt(out_dir, omega, charge=0, multiplicity=1):
    """
    Write python setup file for optimization job
    """
    fout = os.path.join(out_dir, 'optimize.py')
    with open(fout,'w+') as wt:
        wt.write("""import json
import psi4
from sys import argv

mol_file = argv[1]
molecule_name = (mol_file.split('/')[-1]).split('.')[0]
molecule_dir = '/'.join(mol_file.split('/')[:-1])
with open(mol_file,'r') as mol:
    mol = psi4.core.Molecule.from_string(mol.read(), dtype='xyz')
mol.set_molecular_charge("""+str(charge)+""")
mol.set_multiplicity("""+str(multiplicity)+""")

psi4.set_memory('2 GB')
psi4.set_num_threads(2)
psi4.set_module_options('alpha',{'DFT_OMEGA':"""+str(omega)+"""})
psi4.set_output_file(molecule_name + '_geometry_optimization.dat', False)
psi4.set_options({'basis': 'def2-TZVP'})
final_energy = psi4.optimize('LRC-wPBEH', molecule=mol)
mol.save_xyz_file(molecule_name + '_geometry_final.xyz',False)

json_data = {molecule_name: final_energy}
json_file = ("{}/energy_{}.txt".format(molecule_dir,molecule_name))
with open(json_file,'w') as f:
    json.dump(json_data, f, indent=2)""")

def GetRunFolders(molecule_dir, out_dir, nflag=''):
    mol_name = molecule_dir.split('/')[-1].split('.')[0]
    out_files = [x for x in os.listdir(molecule_dir) if x.endswith('_ip_fitting.dat')]
    xyz_files = [x for x in os.listdir(molecule_dir) if x.endswith('.xyz')]
    if len(out_files) == 1:
        log_path = os.path.join(molecule_dir,out_files[0])
        # if mol.properly_terminated:
        #     normal = True
        # else:
        #     normal = False
        #     print("Error. {} job did not have a normal termination.".format(mol_name))
    else:
        normal = False
        print("{} has {} log file(s).".format(mol_name, len(out_files)))
    if len(xyz_files) == 0:
        normal = False
        print("Error. {} does not contain a xyz file.".format(mol_name))
    else:
        pass
    if normal == False:
        txt_file = os.path.join(out_dir,'folders{}_to_run.txt'.format(nflag))
        with open(txt_file, "a+") as fn:
            fn.write(molecule_dir+"\n")

def ImplementSetup(molpath, xyz_file, opt_runs_path, omega, multiplicity=1, charge=0):
    mol_name = (molpath.split('/')[-1]).split('.')[0]
    if not os.path.isdir(molpath): os.mkdir(molpath)
    WriteXYZ(xyz_file, molpath)
    WriteOpt(molpath, omega, charge, multiplicity)
    GetRunFolders(molpath, opt_runs_path)#, nflag=str(charge))
    print("Done setting up {} charge{}!".format(mol_name, charge))


home = os.getcwd()
xyz_path = os.path.join(home, 'xyz/')
opt_runs_path = os.path.join(home,'optimization/')
omega_file = os.path.join(home, 'omegas/master_omegas.json')
with open(omega_file, 'r') as fn:
    omega_dict = json.load(fn)

for mol in os.listdir(xyz_path):
    mol_name = mol.split('.')[0]
    xyz_file = os.path.join(xyz_path, mol)
    try:
        omega = omega_dict[mol[:-4]]
    except(KeyError):
        print("No omega value for {}. Optimization files not set up.".format(mol[:-4]))
        continue


    if 'cp' in mol_name:
        if 'whol' in mol_name:
            molpath = os.path.join(opt_runs_path, mol_name)
            ImplementSetup(molpath, xyz_file, opt_runs_path, omega, multiplicity=1, charge=2)
            cat1path = os.path.join(molpath, 'cation1/')
            ImplementSetup(molpath,cat1path, xyz_file, opt_runs_path, omega, multiplicity=2, charge=3)
            cat2path = os.path.join(molpath, 'cation2/')
            ImplementSetup(cat2path, xyz_file, opt_runs_path, omega, multiplicity=1, charge=4)
        else:
            molpath = os.path.join(opt_runs_path, mol_name)
            ImplementSetup(molpath, xyz_file, opt_runs_path, omega, multiplicity=2, charge=1)
            cat1path = os.path.join(molpath, 'cation1/')
            ImplementSetup(cat1path, xyz_file, opt_runs_path, omega, multiplicity=1, charge=2)
            cat2path = os.path.join(molpath, 'cation2/')
            ImplementSetup(cat2path, xyz_file, opt_runs_path, omega, multiplicity=2, charge=3)
    else:
        molpath = os.path.join(opt_runs_path, mol_name)
        ImplementSetup(molpath, xyz_file, opt_runs_path, omega, multiplicity=1, charge=0)
        cat1path = os.path.join(molpath, 'cation1/')
        ImplementSetup(cat1path, xyz_file, opt_runs_path, omega, multiplicity=2, charge=1)
        cat2path = os.path.join(molpath, 'cation2/')
        ImplementSetup(cat2path, xyz_file, opt_runs_path, omega, multiplicity=1, charge=2)
