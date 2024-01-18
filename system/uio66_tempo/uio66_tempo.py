import sys
import numpy as np
import pickle

sys.path.append('../../../MD_ZIF7_TEMPO/system/components')
sys.path.append('../uio66')
from atom import Atom, mol2_to_atoms, count_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from pbc import calculate_lattice_vectors

atoms_tempo = mol2_to_atoms(read_mol2_file('../tempo/tempo.mol2'))
for atom_idx in range(len(atoms_tempo)):
    atoms_tempo[atom_idx].resid_idx = 1
    atoms_tempo[atom_idx].resid = 'TMP'
    atoms_tempo[atom_idx].atom_idx = atom_idx
 
def put_tempo_in_lattice(atoms_uio66, atoms_tempo):
    # Calculate cell center
    cell_center = np.zeros(3)
    for atom_idx in [1704, 2611]:
        cell_center += atoms_uio66[atom_idx].r / 2.0
    cell_center -= np.array([2.0, 1.0, 2.0])

    # Calculate TEMPO center
    tempo_center = np.zeros(3)
    for atom in atoms_tempo:
        tempo_center += atom.r / len(atoms_tempo)
        
    # Translate TEMPO atoms inside cell
    for atom_idx in range(len(atoms_tempo)):
        atoms_tempo[atom_idx].r += cell_center - tempo_center

    # Calculate minimal distance from TEMPO to cell
    min_dist = 1e10
    for atom_tempo in atoms_tempo:
        for atom_uio66 in atoms_uio66:
            dist = pow(sum(np.power(atom_tempo.r - atom_uio66.r, 2)), 0.5)
            min_dist = min(min_dist, dist)
    print(min_dist)
    return atoms_uio66.copy(), atoms_tempo.copy()

a     = 20.7465
b     = 20.7465
c     = 20.7465
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi

bounds_a, bounds_b, bounds_c = [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]

# Load data
with open('../uio66/__tmp/atoms_uio66.pickle', 'rb') as handle:
    atoms_uio66 = pickle.load(handle)

# Save .gro file
atoms_uio66, atoms_tempo = put_tempo_in_lattice(atoms_uio66, atoms_tempo)
with open('__tmp/atoms_uio66.pickle', 'wb') as handle:
    pickle.dump(atoms_uio66, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('__tmp/atoms_tempo.pickle', 'wb') as handle:
    pickle.dump(atoms_tempo, handle, protocol=pickle.HIGHEST_PROTOCOL)
atoms_uio66_tempo = atoms_uio66.copy()
atoms_uio66_tempo.extend(atoms_tempo.copy())
write_gro_file(atoms_uio66_tempo, 'uio66_tempo.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)

