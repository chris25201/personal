from pylab import *
from ase.io import read,write
from thermo.gpumd.preproc import add_basis, repeat
from thermo.gpumd.io import create_basis, create_kpoints, ase_atoms_to_gpumd

B=read('POSCAR')
B.wrap()
B1 = B.repeat([70,16,8])
#num_atoms = len(B)
#B.unitcell = list(range(num_atoms))
#B.basis = list(range(num_atoms))
write("model.xyz", B1)

