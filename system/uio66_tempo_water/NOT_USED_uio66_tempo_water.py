import sys
import numpy as np
import pickle
from random import random

sys.path.append('../../../MD_ZIF7_TEMPO/system/components')
sys.path.append('../uio66')
from atom import Atom, mol2_to_atoms, count_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from pbc import calculate_lattice_vectors

a     = 20.7465
b     = 20.7465
c     = 20.7465
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi

bounds_a, bounds_b, bounds_c = [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]

def select_loc(atoms, mol, num_iter=1000):
    r_selected_best = np.zeros(3)
    dist_best = 0.0
    for iter in range(num_iter):
        r_selected = np.array([random() * 41.493, random() * 41.493, random() * 41.493])
        for i in range(len(atoms)):
            dist = 1e10
            for j in range(len(mol)):
                delta = np.abs(atoms[i].r - (mol[j].r + r_selected))
                delta[0] = min(delta[0], 41.493 - delta[0])
                delta[1] = min(delta[1], 41.493 - delta[1])
                delta[2] = min(delta[2], 41.493 - delta[2])
                dist = min(dist, pow(sum(np.power(delta, 2)), 0.5))
        if dist > dist_best:
            dist_best = dist
            r_selected_best = r_selected
        if dist_best > 20:
            break
    return r_selected

# Load data
for water_count in [100, 500, 1000]:
    with open('../uio66_tempo/__tmp/atoms_uio66.pickle', 'rb') as handle:
        atoms_uio66 = pickle.load(handle)
    with open('../uio66_tempo/__tmp/atoms_tempo.pickle', 'rb') as handle:
        atoms_tempo = pickle.load(handle)

    atoms = atoms_uio66.copy()
    atoms.extend(atoms_tempo.copy())

    for w_i in range(water_count):
        print(f'atom {w_i} / {water_count}')
        water = []
        water.append(Atom('OW' , [0.0,    0.0,  0.0]))
        water.append(Atom('HW1', [0.41,  -0.58, 0.65]))
        water.append(Atom('HW2', [-0.93, -0.08, 0.17]))
        water.append(Atom('MW' , [-0.06, -0.08, 0.1]))
        
        if w_i == 0:
            shift = 0.5 * (atoms[379].r + atoms[3514].r)
        if w_i == 1:
            shift = 0.5 * (atoms[434].r + atoms[3518].r)
        if w_i == 2:
            shift = 0.5 * (atoms[1704].r + atoms[3518].r)
        if w_i > 3:
            shift = select_loc(atoms, water)
            
        for atom_idx in range(4):
            water[atom_idx].r += shift
            water[atom_idx].resid_idx = w_i
            water[atom_idx].resid = 'SOL'
            water[atom_idx].atom_idx = len(atoms) + atom_idx
        
        atoms.extend(water.copy())
    
    write_gro_file(atoms, f'uio66_tempo_{water_count}water.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)

