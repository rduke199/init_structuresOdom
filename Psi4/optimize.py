import psi4
from sys import argv

mol_file = argv[1]
molecule_name = (mol_file.split('/')[-1]).split('.')[0]
molecule_dir = '/'.join(mol_file.split('/')[:-1])
with open(mol_file,'r') as mol:
    mol = psi4.core.Molecule.from_string(mol.read(), dtype='xyz')
mol.set_molecular_charge(charge) ##input
mol.set_multiplicity(multiplicity) ##input

psi4.set_memory('2 GB')
psi4.set_num_threads(2)
psi4.set_module_options('alpha',{'DFT_OMEGA':omega}) ##input
psi4.set_output_file(molecule_name + '_geometry_optimization.dat', False)
psi4.set_options({'basis': 'def2-TZVP'})
final_energy = psi4.optimize('LRC-wPBEH', molecule=mol)
mol.save_xyz_file(molecule_name + '_geometry_final.xyz',False)

json_data = {"molecule_name" : molecule_name, "final_energy" : final_energy}
json_file = ("{}/energy_{}.txt".format(molecule_dir,molecule_name))
with open(json_file,'w') as f:
    json.dump(json_data, f, indent=2)
