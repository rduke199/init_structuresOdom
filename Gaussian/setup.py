import os
import re
from sys import argv
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianOutput, GaussianInput


class SetUpGaus:
    """
    For a molecule this sets up calculations
    """

    def __init__(self, structure_fn, calc_dir, ground_charge=0):
        self.structure_fn = structure_fn
        self.ground_charge = ground_charge
        self.charges_dict = {'neutral': self.ground_charge, 'cation1': self.ground_charge + 1,
                             'cation2': self.ground_charge + 2}

        self.calc_dir = calc_dir
        self.mol_dir = os.path.join(self.calc_dir, self.mol_name)
        self.wtuning_dir = os.path.join(self.mol_dir, 'wtuning/')
        self.wtuning_log = os.path.join(self.wtuning_dir, self.mol_name+'_opt_0.log')

        if not os.path.isdir(self.mol_dir):
            os.mkdir(self.mol_dir)
        if not os.path.isdir(self.wtuning_dir):
            os.mkdir(self.wtuning_dir)

    def get_mol_calc_path(self, cation_name):
        return os.path.join(self.mol_dir, cation_name)

    def get_log(self, cation_name, calc_type):
        return os.path.join(self.mol_dir, cation_name, calc_type+'.log')

    @property
    def mol_name(self):
        """
        Get molecule name from file name
        """
        return self.structure_fn.split('/')[-1].split('.')[0]

    @property
    def omega(self):
        output_path = os.path.join(self.wtuning_dir, 'output.log')
        try:
            with open(output_path, 'r') as fn:
                w_data = fn.readlines()[-2].split()[1]
            return "0{}".format(w_data.split('.')[1])
        except FileNotFoundError:
            return None

    @staticmethod
    def mol_object(mol_fn):
        """
        Get pymatgen object in preparation for further processing
        Returns:
            A pymatgen object
        """
        try:
            return Molecule.from_file(mol_fn)
        except FileNotFoundError:
            print("File does not exist: {}".format(mol_fn))

    @staticmethod
    def write_wtuning(out_dir, charge=0):
        """
        Write python setup file for wtuning job
        """
        fout = os.path.join(out_dir, 'old_scripts/wtuning.py')
        with open(fout, 'w+') as wt:
            wt.write("""from sys import argv
from ocelot.task.wtuning import WtuningJob
from pymatgen.core.structure import Molecule

fin = argv[1]
name = (fin.split('/')[-1]).split('.')[0]
pymol = Molecule.from_file(fin)
pymol.perturb(0.05)
job = WtuningJob(func='LC-wHPBE', basis='TZVP',name=name,nproc=8,mem=30,n_charge=""" + str(charge) + """,n_spin=1,scheme='Jh')
job.mol = pymol
job.mol = job.geo_opt()
job.wtuning_cycle(max_cycles=0)""")

    @staticmethod
    def generate_gjf(in_fn, out_dir, functional='LC-wHPBE', basis_set='TZVP', charge=0, calculation='opt', omega=None, oldchk=None):
        # Convert an individually inputted xyz file to a gjf Gaussian input file
        mol = Molecule.from_file(in_fn)
        mol.perturb(0.1)
        mol_name = in_fn.split('/')[-1][:14]
        link0_parameters = {'%mem': '5GB', '%chk': '{}.chk'.format(calculation)}
        route_parameters = {calculation: '', 'SCF': '(MaxCycle=250)', 'Int': '(Grid=SuperFine)'}
        if calculation == 'tddft':
            route_parameters = {'TD(NStates=5, 50-50)': ''}
        if omega is not None:
            route_parameters["iop(3/107={}, 3/108={})".format(omega, omega)] = ''
        if oldchk:
            link0_parameters['%oldchk'] = '{}.chk'.format(oldchk)
            route_parameters['Geom'] = 'AllCheck'
        gau = GaussianInput(mol=mol, charge=charge, functional=functional, basis_set=basis_set,
                            route_parameters=route_parameters,
                            link0_parameters=link0_parameters)
        gjf_file = gau.write_file('{}/{}.gjf'.format(out_dir, calculation))
        return gjf_file

    def check_if_run_finished(self, calc_type, cation_name=''):
        txt_file = os.path.join(self.calc_dir, 'folders_to_run_{}.txt'.format(calc_type))
        if calc_type == 'wtuning':
            log_path, dir_path = self.wtuning_log, self.wtuning_dir
        else:
            log_path, dir_path = self.get_log(cation_name, calc_type), self.get_mol_calc_path(cation_name)

        if not os.path.isfile(log_path):
            print("No log file for {} calculation for {} {} ".format(calc_type, self.mol_name, cation_name))
        elif calc_type == 'opt':
            with open(log_path) as fn:
                log_fn = fn.readlines()
            for line in log_fn:
                if re.search(' Optimized Parameters', line):
                    last_line = log_fn[-1]
                    if re.search('Normal termination', last_line):
                        return True
            print("Error. {} {} optimization did NOT terminated normally".format(self.mol_name, cation_name))
        else:
            if GaussianOutput(log_path).properly_terminated:
                print("{} {} {} terminated normally".format(calc_type, self.mol_name, cation_name))
                return True
            print("Error. {} {} {} did NOT terminated normally".format(calc_type, self.mol_name, cation_name))

        with open(txt_file, "a+") as fn:
            fn.write(dir_path + "\n")
        return False

    def setup_wtuning(self):
        xyz_fn = os.path.join(self.wtuning_dir, self.mol_name + '.xyz')
        self.mol_object(self.structure_fn).to(filename=xyz_fn)
        self.write_wtuning(self.wtuning_dir, self.ground_charge)
        self.check_if_run_finished('wtuning')

    def setup_general(self, cation_name, calc_type, **kwargs):
        out_path = self.get_mol_calc_path(cation_name)
        structure = self.structure_fn if calc_type == 'opt' else self.get_log(cation_name, 'opt')
        self.generate_gjf(structure, out_path, charge=self.charges_dict[cation_name], omega=self.omega,
                          calculation=calc_type, **kwargs)
        self.check_if_run_finished(calc_type, cation_name=cation_name)
        print("{} Setup Completed for {} {}".format(calc_type, self.mol_name, cation_name))


def setup_for_mol(xyz_file, calc_type, calc_dir, ground_charge=0):
    molecule = SetUpGaus(xyz_file, calc_dir)
    if calc_type == "wtuning":
        molecule.setup_wtuning()
    else:
        charges = {'neutral': ground_charge, 'cation1': ground_charge + 1, 'cation2': ground_charge + 2}
        for charge_name in charges.keys():
            molecule.setup_general(charge_name, calc_type)


def get_all_runfiles(calc_type, ground_charge=0):
    # Establish paths that will be used
    home = os.getcwd()
    xyz_dir = os.path.join(home, 'xyz/')
    calc_dir = os.path.join(home, 'opt/')

    # Iterate through mol directories and gather information
    mols = [m for m in os.listdir(xyz_dir) if m.startswith('mol')]
    for mol in mols:
        xyz_file = os.path.join(xyz_dir, mol)
        molecule = SetUpGaus(xyz_file, calc_dir)
        if calc_type == "wtuning":
            molecule.check_if_run_finished(calc_type)
        else:
            charges = {'neutral': ground_charge, 'cation1': ground_charge + 1, 'cation2': ground_charge + 2}
            for charge_name in charges.keys():
                molecule.check_if_run_finished(calc_type, cation_name=charge_name)

def main():
    # Establish paths that will be used
    home = os.getcwd()
    xyz_dir = os.path.join(home, 'xyz/')
    calc_dir = os.path.join(home, 'opt/')
    calc_type = argv[1]

    # Iterate through mol directories and gather information
    mols = [m for m in os.listdir(xyz_dir) if m.startswith('mol')]
    for mol in mols:
        mol_structure_path = os.path.join(xyz_dir, mol)
        setup_for_mol(mol_structure_path, calc_type, calc_dir)


if __name__ == "__main__":
    main()
