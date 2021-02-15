import os
import re
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput

def XYZtoGJF(xyz_fn, out_dir, functional, basis_set, charge=0):
    """
    Convert an individually imputted xyz file to a gjf Gaussian input file
    """
    mol = Molecule.from_file(xyz_fn)
    mol_name = xyz_fn.split('/')[-1].split('.')[0]
    gau = GaussianInput(mol=mol,charge=charge,functional=functional,basis_set=basis_set,  route_parameters={"opt":""}, link0_parameters={'%mem':'5GB','%chk':'{}.chk'.format(mol_name)})
    gjf_file = gau.write_file('{}/{}.gjf'.format(out_dir,mol_name))
    return gjf_file

def GetRunFolders(molecule_dir, out_dir, nflag='_'):
    mol_name = molecule_dir.split('/')[-1].split('.')[0]
    log_files = [x for x in os.listdir(molecule_dir) if x.endswith('.log') and not x.endswith('deg.log')]
    gjf_files = [x for x in os.listdir(molecule_dir) if x.endswith('.gjf')]
    normal = 0
    if len(log_files) == 1:
        with open(os.path.join(molecule_dir,log_files[0])) as fn:
            log_fn = fn.readlines()
            last_line = log_fn[-1]
            if re.search('Normal termination', last_line):
                for line in log_fn:
                    if re.search(' Optimized Parameters',line):
                        normal = 1
                        break
                    else:
                        continue
                if normal == 0:
                    print("Error. {} job had a normal termination but did not have optimized parameters.".format(mol_name))
            else:
                print("Error. {} job did not have a normal termination.".format(mol_name))
    if len(gjf_files) == 0:
        normal = 2
        print("Error. {} does not contain a gjf file.".format(mol_name))
    if normal == 0:
        txt_file = os.path.join(out_dir,'folders{}_to_run.txt'.format(nflag))
        with open(txt_file, "a+") as fn:
            fn.write(molecule_dir+"\n")

home = os.getcwd()
xyz_path = os.path.join(home, 'xyz/')
dft_runs_path = os.path.join(home,'dft_runs/')

funcitonals = ['LC-wHPBE','wB97XD']
basis_sets = ['cc-pVDZ', '6-31+G**', 'TZVP', 'aug-cc-pVTZ']

for f in funcitonals:
    fpath = os.path.join(dft_runs_path,f)
    if not os.path.isdir(fpath): os.mkdir(fpath)
    for b in basis_sets:
        bpath = os.path.join(fpath,b)
        if not os.path.isdir(bpath): os.mkdir(bpath)
        # cprop = [m for m in os.listdir(xyz_path) if m.startswith('cprop')] #
        # for mol in cprop: #
        for mol in os.listdir(xyz_path):
            mol_name = mol.split('.')[0]
            molpath = os.path.join(bpath, mol_name)
            if not os.path.isdir(molpath): os.mkdir(molpath)
            xyz_file = os.path.join(xyz_path, mol)
            try:
                # XYZtoGJF(xyz_file, molpath, f, b, charge=0)
                print("Done setting up {} with {}/{}!".format(mol_name, f, b))
            except:
                print("Error. Calculation for {} with {}/{} was not set up.".format(mol_name, f, b))
            GetRunFolders(molpath, dft_runs_path,nflag=f[0:2])

            catpath = os.path.join(molpath, 'cation/')
            if not os.path.isdir(catpath): os.mkdir(catpath)
            try:
                # XYZtoGJF(xyz_file, catpath, f, b, charge=1)
                print("Done setting up {} cation with {}/{}!".format(mol_name, f, b))
            except:
                print("Error. Calculation for {} cation with {}/{} was not set up.".format(mol_name, f, b))
            GetRunFolders(catpath, dft_runs_path,nflag=str(f[0:2]+'cat'))
