import os
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianOutput


def write_xyz(in_file, out_dir):
    """
    Write xyz file.
    """
    mol = Molecule.from_file(in_file)
    fout = os.path.join(out_dir, in_file.split('/')[-1])
    mol.to(filename=fout)


def write_wtuning(out_dir, charge=0):
    """
    Write python setup file for wtuning job
    """
    fout = os.path.join(out_dir, 'wtuning.py')
    with open(fout, 'w+') as wt:
        wt.write("""from sys import argv
from ocelot.task.wtuning import WtuningJob
from pymatgen.core.structure import Molecule

fin = argv[1]
name = (fin.split('/')[-1]).split('.')[0]
pymol = Molecule.from_file(fin)
pymol.perturb(0.05)
job = WtuningJob(func='LC-wHPBE', basis='TZVP',name=name,nproc=8,mem=30,n_charge="""+str(charge)+""",n_spin=1,scheme='Jh')
job.mol = pymol
job.mol = job.geo_opt()
job.wtuning_cycle(max_cycles=0)""")


def get_run_folders(molecule_dir, out_dir, nflag=''):
    mol_name = molecule_dir.split('/')[-1].split('.')[0]
    log_files = [x for x in os.listdir(molecule_dir) if x.endswith('opt_0.log')]
    xyz_files = [x for x in os.listdir(molecule_dir) if x.endswith('.xyz')]
    if len(log_files) == 1:
        log_path = os.path.join(molecule_dir, log_files[0])
        mol = GaussianOutput(log_path)
        if mol.properly_terminated:
            normal = 1
        else:
            normal = 0
            print(log_path)
            print("Error. {} job did not have a normal termination.".format(mol_name))
    else:
        normal = 0
        print("{} has {} log file(s).".format(mol_name, len(log_files)))
    if len(xyz_files) == 0:
        normal = 2
        print("Error. {} does not contain a xyz file.".format(mol_name))
    if normal == 0:
        txt_file = os.path.join(out_dir,'folders{}_to_run.txt'.format(nflag))
        with open(txt_file, "a+") as fn:
            fn.write(molecule_dir+"\n")


def setup_a_folder(molpath, xyz_file, wtuning_path, charge=0):
    mol_name = (molpath.split('/')[-1]).split('.')[0]
    if not os.path.isdir(molpath): os.mkdir(molpath)
    try:
        write_xyz(xyz_file, molpath)
        write_wtuning(molpath, charge)
        get_run_folders(molpath, wtuning_path)
        print("Done setting up wtuning for {}.".format(mol_name))
    except:
        print("Error setting up wtuning for {}!".format(mol_name))


def implement_setup(home):
    xyz_path = os.path.join(home, 'xyz/')
    wtuning_path = os.path.join(home, 'wtuning/')
    if not os.path.isdir(wtuning_path): os.mkdir(wtuning_path)

    for mol in os.listdir(xyz_path):
        mol_name = mol.split('.')[0]
        xyz_file = os.path.join(xyz_path, mol)
        molpath = os.path.join(wtuning_path, mol_name)

        if 'cp' in mol_name:
            if 'whol' in mol_name:
                setup_a_folder(molpath, xyz_file, wtuning_path, charge=2)
            else:
                setup_a_folder(molpath, xyz_file, wtuning_path, charge=1)
        else:
            setup_a_folder(molpath, xyz_file, wtuning_path, charge=0)


def main():
    home = os.getcwd()
    implement_setup(home)


if __name__ == "__main__":
    main()
