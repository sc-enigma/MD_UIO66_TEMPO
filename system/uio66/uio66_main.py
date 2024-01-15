import sys
import numpy as np
import pickle

sys.path.append('../../../MD_ZIF7_TEMPO/system/components')
from atom import Atom, mol2_to_atoms, count_atoms, remove_non_bonded_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, compose_itp_files
from pbc import apply_pbc, calculate_lattice_vectors

from uio66_utils import define_uio66_atom_types, define_uio66_atom_names, remove_Zr_Zr_bond
from uio66_params import get_uio66_params

# STEP 1. Read data
# Set lower bound in Mercury calculate packing = 0.0
a     = 20.7465
b     = 20.7465
c     = 20.7465
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi

bounds_a, bounds_b, bounds_c = [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]

vec_a, vec_b, vec_c = calculate_lattice_vectors(a, b, c, alpha, beta, gamma)
'''atoms = mol2_to_atoms(read_mol2_file('__uio66_source.mol2'))

# STEP 2. Apply periodic boundary conditions
atoms = apply_pbc(atoms, a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
atoms = remove_Zr_Zr_bond(atoms)
with open('__tmp/atoms_uio66.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)'''

with open('__tmp/atoms_uio66.pickle', 'rb') as handle:
    atoms = pickle.load(handle)

# STEP 3. Define atom types and names
atoms = remove_non_bonded_atoms(atoms)
atoms = define_uio66_atom_types(atoms)
atoms = define_uio66_atom_names(atoms)

# STEP 4. Write .gro file for large pore structure
with open('__tmp/atoms_uio66.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)
write_gro_file(atoms, 'uio66.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)

# STEP 5. Write .itp files
mass, charge, bond_params, angle_params, dihedral_params = get_uio66_params()
write_atoms(atoms, charge, mass, 'UIO', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds.itp')
write_angles(atoms, angle_params, 'angles.itp')
write_dihedrals(atoms, dihedral_params, 'dihedrals.itp')
compose_itp_files(['atomtypes.itp',
                   'moleculetype.itp', 'atoms.itp', 'bonds.itp', 'angles.itp', 'dihedrals.itp'], 'uio66.itp')








