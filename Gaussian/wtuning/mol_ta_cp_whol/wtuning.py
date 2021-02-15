from sys import argv
from ocelot.task.wtuning import WtuningJob
from pymatgen.core.structure import Molecule

fin = argv[1]
name = (fin.split('/')[-1]).split('.')[0]
pymol = Molecule.from_file(fin)
pymol.perturb(0.05)
job = WtuningJob(func='LC-wHPBE', basis='TZVP',name=name,nproc=8,mem=30,n_charge=2,n_spin=1,scheme='Jh')
job.mol = pymol
job.mol = job.geo_opt()
job.wtuning_cycle(max_cycles=0)