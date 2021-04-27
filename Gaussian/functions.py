import os
import re
import json
from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput


def generate_gjf(in_fn, out_dir, functional='LC-wHPBE', basis_set='TZVP', charge=0, calculation='opt', omega=None,
                 oldchk=None):
    """
    Convert an individually inputted xyz file to a gjf Gaussian input file
    """
    mol = Molecule.from_file(in_fn)
    mol.perturb(0.1)
    mol_name = in_fn.split('/')[-1][:14]
    link0_parameters = {'%mem': '5GB', '%chk': '{}.chk'.format(calculation)}
    if calculation == 'tddft':
        route_parameters = {'TD(NStates=5, 50-50)': ''}
    else:
        route_parameters = {calculation: '', 'SCF': '(MaxCycle=250)', 'Int': '(Grid=Ultrafine)'}
    if omega is not None:
        route_parameters["iop(3/107={}, 3/108={})".format(omega, omega)] = ''
    else:
        pass
    if oldchk:
        link0_parameters['%oldchk'] = '{}.chk'.format(oldchk)
        route_parameters['Geom'] = 'AllCheck'
    else:
        pass
    gau = GaussianInput(mol=mol, charge=charge, functional=functional, basis_set=basis_set,
                        route_parameters=route_parameters,
                        link0_parameters=link0_parameters)
    gjf_file = gau.write_file('{}/{}.gjf'.format(out_dir, calculation))
    return gjf_file


def get_run_folders(molecule_dir, out_dir, nflag=''):
    mol_name = molecule_dir.split('/')[-1].split('.')[0]
    opt_log = os.path.join(molecule_dir, 'opt.log')
    opt_gjf = os.path.join(molecule_dir, 'opt.gjf')
    freq_gjf = os.path.join(molecule_dir, 'freq.gjf')
    tddft_gjf = os.path.join(molecule_dir, 'tddft.gjf')
    if os.path.isfile(tddft_gjf):
        txt_file = os.path.join(out_dir, 'folders_to_run_tddft.txt')
        with open(txt_file, 'a+') as fn:
            fn.write(molecule_dir + '\n')
        return None
    normal = None
    if os.path.isfile(opt_log):
        with open(opt_log) as fn:
            log_fn = fn.readlines()
        for line in log_fn:
            if re.search(' Optimized Parameters', line):
                last_line = log_fn[-1]
                if re.search('Normal termination', last_line):
                    if os.path.isfile(freq_gjf):
                        txt_file = os.path.join(out_dir, 'folders_to_run_freq.txt')
                        with open(txt_file, 'a+') as fn:
                            fn.write(molecule_dir + '\n')
                    else:
                        print("{} is finished optimizing, but there is not freq.gjf file.".format(mol_name))
                    return None
            if re.search('Error termination', line):
                normal = 0
                txt_file = os.path.join(out_dir, 'folders_to_run{}.txt'.format(nflag))
                with open(txt_file, 'a+') as fn:
                    fn.write(molecule_dir + '\n')
                print("Error in termination for {}".format(mol_name))
                return None

    else:
        print("{} does not have an opt.log file.")
    if os.path.isfile(opt_gjf):
        print("{} may still be running".format(mol_name))
    else:
        print("Error. {} does not contain a gjf file.".format(mol_name))


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


def get_out_paths(mol, out_dir, cation='', subtype=''):
    mol_name = mol.split('.')[0]
    mol_out_path = os.path.join(out_dir, mol_name)
    if not os.path.isdir(mol_out_path): os.mkdir(mol_out_path)
    out_path = os.path.join(mol_out_path, cation, subtype)
    if not os.path.isdir(out_path): os.mkdir(out_path)
    return out_path


def setup_a_folder(mol_name, in_file, out_path, runs_folder, charge=0, calculation='opt', **kwargs):
    if not os.path.isdir(out_path): os.mkdir(out_path)
    # try:
    generate_gjf(in_file, out_path, charge=charge, calculation=calculation, **kwargs)
    # txt_file = os.path.join(runs_folder, 'folders_to_run.txt')
    # with open(txt_file, 'a+') as fn:
    #     fn.write(out_path + '\n')
    print("Done setting up {} charge{}!".format(mol_name, charge))
    # except:
    #     print("Error. Calculation for {} charge{} was not set up.".format(mol_name, charge))


def write_json(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)


def runtime_from_log(logfile):
    """
    Collects runtime in core hours from a logfile
    """
    time_patt = re.compile(r"\d+\.d+|\d+")
    time_data = []
    with open(logfile, "r") as f:
        line = f.readline()
        while line != "":
            if re.match("Job cpu time", line.strip()):
                time_data.extend(time_patt.findall(line))
            line = f.readline()
    if time_data:
        time_data = [float(time) for time in time_data]
        runtime = (time_data[0] * 86400 + time_data[1] * 3600 + time_data[2] * 60 + time_data[3]) / 3600
    else:
        runtime = 0
    return round(runtime, 3)


def write_master_json(json_dir, master_json_file, prop=None):
    """
    Write a master data file from all json files in json_dir
    """
    data = {}
    json_files = [f for f in os.listdir(json_dir) if f.startswith("mol")]
    for f in json_files:
        fpath = os.path.join(json_dir, f)
        with open(fpath) as fn:
            mol_info = json.load(fn)
        mol_name = mol_info["molecule_name"]
        if prop is None:
            data[mol_name] = mol_info
        else:
            prop_info = mol_info[prop]
            data[mol_name] = prop_info
    write_json(data, master_json_file)
