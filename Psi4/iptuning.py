from sys import argv
import json
import psi4
from psi4.driver.frac import ip_fitting
from psi4.core.Molecule import from_string

mol_file = argv[1]
molecule_name = (mol_file.split('/')[-1]).split('.')[0]
with open(mol_file,'r') as f:
    mol = from_string(f.read(), dtype='xyz')
mol.set_molecular_charge(0)
mol.set_multiplicity(1)

psi4.set_memory('2 GB')
psi4.set_num_threads(2)
psi4.set_output_file(molecule_name + '_ip_fitting.dat', False)
psi4.set_options({'basis': 'def2-TZVP'})
omega = ip_fitting('LRC-wPBEH', 0.1, 2.0,molecule=mol)
json_data = {"molecule_name" : molecule_name, "omega" : omega}
json_file = ("omega_{}.txt".format(molecule_name))
with open(json_file,'w') as f:
    json.dump(json_data, f, indent=2)
