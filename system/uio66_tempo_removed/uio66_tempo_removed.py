import sys
import numpy as np
import pickle

sys.path.append('../../../MD_ZIF7_TEMPO/system/components')
sys.path.append('../uio66')

from atom import Atom, mol2_to_atoms, count_atoms, remove_atoms, remove_non_bonded_atoms, copy_atom
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, compose_itp_files
from pbc import apply_pbc, calculate_lattice_vectors

from uio66_utils import define_uio66_atom_types, define_uio66_atom_names, remove_Zr_Zr_bond, remove_periodic_bonds
from uio66_params import get_uio66_params

def find_loc(atoms, idx):
    r_best = np.zeros(3)
    dist_best = 0.0
    for theta in np.linspace(0, np.pi, 10):
        for phi in np.linspace(0, 2.0 * np.pi, 10):
            r = np.copy(atoms[idx].r)
            r[0] += 1.25 * np.sin(theta) * np.cos(phi)
            r[1] += 1.25 * np.sin(theta) * np.sin(phi)
            r[2] += 1.25 * np.cos(theta)
            
            dist = 1e10
            for atom_idx in range(len(atoms)):
                if atom_idx == idx:
                    continue
                delta = pow(sum(np.power(np.copy(atoms[atom_idx].r) - np.copy(r), 2)), 0.5)
                dist = min(dist, delta)
            if dist > dist_best:
                dist_best = dist
                r_best = np.copy(r)
    return r_best

def add_H(atoms, idx, resid=''):
    l = len(atoms)
    atoms.append(Atom('HR', np.copy(find_loc(atoms, idx))))
    atoms[idx].adjacency.append(l)
    atoms[l].adjacency.append(idx)
    
    atoms[l].resid_idx = 1
    atoms[l].resid = resid
    atoms[l].atom_idx = l + 1
    atoms[l].atom_type = 'repl_H'
    
    return atoms

def add_OH(atoms, idx, resid=''):
    l = len(atoms)
    atoms.append(Atom('OR', np.copy(find_loc(atoms, idx))))
    atoms[idx].adjacency.append(l)
    atoms[l].adjacency.append(idx)
    atoms.append(Atom('HR', np.copy(find_loc(atoms, l))))
    atoms[l].adjacency.append(l + 1)
    atoms[l + 1].adjacency.append(l)
    
    atoms[l].resid_idx = 1
    atoms[l + 1].resid_idx = 1
    atoms[l].resid = resid
    atoms[l + 1].resid = resid
    atoms[l].atom_idx = l + 1
    atoms[l + 1].atom_idx = l + 2
    atoms[l].atom_type = 'repl_O'
    atoms[l + 1].atom_type = 'repl_H'
    
    return atoms
    
def add_OH2(atoms, idx, resid = ''):
    l = len(atoms)
    atoms.append(Atom('OR', np.copy(find_loc(atoms, idx))))
    atoms[idx].adjacency.append(l)
    atoms[l].adjacency.append(idx)
    atoms.append(Atom('HR', np.copy(find_loc(atoms, l))))
    atoms[l].adjacency.append(l + 1)
    atoms[l + 1].adjacency.append(l)
    atoms.append(Atom('HR', np.copy(find_loc(atoms, l))))
    atoms[l].adjacency.append(l + 2)
    atoms[l + 2].adjacency.append(l)
    
    atoms[l].resid_idx = 1
    atoms[l + 1].resid_idx = 1
    atoms[l + 2].resid_idx = 1
    atoms[l].resid = resid
    atoms[l + 1].resid = resid
    atoms[l + 2].resid = resid
    atoms[l].atom_idx = l + 1
    atoms[l + 1].atom_idx = l + 2
    atoms[l + 2].atom_idx = l + 3
    atoms[l].atom_type = 'repl_O'
    atoms[l + 1].atom_type = 'repl_H'
    atoms[l + 2].atom_type = 'repl_H'
    
    return atoms

def create_linker(atoms, linker_ids, linker_cnt):
    if linker_cnt == 1:
        resid = 'TRA'
    if linker_cnt == 2:
        resid = 'TRB'
    if linker_cnt == 3:
        resid = 'TRC'
    shift = {}
    atoms_linker = [copy_atom(atoms[linker_id]) for linker_id in linker_ids]
    for atom_idx in range(len(linker_ids)):
        shift[linker_ids[atom_idx]] = atom_idx
    for atom_idx in range(len(atoms_linker)):
        for adj_idx in range(len(atoms_linker[atom_idx].adjacency)):
            adj = atoms_linker[atom_idx].adjacency[adj_idx]
            if adj in shift.keys():
                atoms_linker[atom_idx].adjacency[adj_idx] = shift[adj]
            else:
                atoms_linker[atom_idx].adjacency[adj_idx] = -1
        while -1 in atoms_linker[atom_idx].adjacency:
            atoms_linker[atom_idx].adjacency.remove(-1)
    atoms_linker = add_H(atoms_linker, 0)
    atoms_linker = add_H(atoms_linker, 15)
    for atom_idx in range(len(atoms_linker)):
        atoms_linker[atom_idx].atom_idx = atom_idx + 1
        atoms_linker[atom_idx].resid = resid
        atoms_linker[atom_idx].resid_idx = linker_cnt
    with open(f'__tmp/atoms_linker{linker_cnt}.pickle', 'wb') as handle:
        pickle.dump(atoms_linker, handle, protocol=pickle.HIGHEST_PROTOCOL)
    write_gro_file(atoms_linker, f'__tmp/atoms_linker{linker_cnt}.gro', 41.493, 41.493, 41.493)
    
    mass, charge, bond_params, angle_params, dihedral_params = get_uio66_params()
    write_atoms(atoms_linker, charge, mass, resid, f'__tmp/linker_atoms{linker_cnt}.itp')
    write_bonds(atoms_linker, bond_params,         f'__tmp/linker_bonds{linker_cnt}.itp')
    write_angles(atoms_linker, angle_params,       f'__tmp/linker_angles{linker_cnt}.itp')
    write_dihedrals(atoms_linker, dihedral_params, f'__tmp/linker_dihedrals{linker_cnt}.itp')
    # '../uio66/atomtypes.itp'
    compose_itp_files([f'linker_moleculetype{linker_cnt}.itp',
                       f'__tmp/linker_atoms{linker_cnt}.itp', f'__tmp/linker_bonds{linker_cnt}.itp', f'__tmp/linker_angles{linker_cnt}.itp', f'__tmp/linker_dihedrals{linker_cnt}.itp'],
                       f'linker{linker_cnt}.itp')
    return atoms_linker
    
                   
# Load data
a     = 20.7465
b     = 20.7465
c     = 20.7465
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 90.00 / 180 * np.pi

bounds_a, bounds_b, bounds_c = [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]

with open('../uio66_tempo/__tmp/atoms_uio66.pickle', 'rb') as handle:
    atoms_uio66 = pickle.load(handle)
with open('../uio66_tempo/__tmp/atoms_tempo.pickle', 'rb') as handle:
    atoms_tempo = pickle.load(handle)

# Create linkers
linker1_ids = [399, 307, 312, 308, 309, 400, 401, 310, 426, 2681, 2677, 431, 2678, 2679, 2680, 2676]
linker2_ids = [344, 348, 378, 379, 345, 346, 380, 347, 418, 1645, 421, 1656, 1644, 1646, 1655, 1643]
linker3_ids = [382, 320, 321, 325, 322, 383, 323, 384, 434, 438, 3434, 3391, 3390, 3393, 3389, 3433]
zr_cluster_ids = [317, 373, 324, 381, 336, 377, 311, 369, 319, 365, 343, 306, 331, 354] # 6 Zr, 8 O3

atoms_linker = []
linker_cnt = 0
for linker_ids in [linker1_ids, linker2_ids, linker3_ids]:
    linker_cnt += 1
    atoms_linker.append(create_linker(atoms_uio66, linker_ids, linker_cnt))

# main cycle
for removed_count in [1, 2, 3]:
    for remove_zr_cluster in [False, True]:
        print(f'removed_count = {removed_count} remove_zr_cluster = {remove_zr_cluster}')
        
        # Load data
        with open('../uio66_tempo/__tmp/atoms_uio66.pickle', 'rb') as handle:
            atoms_uio66 = pickle.load(handle)
        with open('../uio66_tempo/__tmp/atoms_tempo.pickle', 'rb') as handle:
            atoms_tempo = pickle.load(handle)
        
        # Add H / OH / OH2
        ids_to_remove = []
        if remove_zr_cluster:
            atoms_uio66 = add_H(atoms_uio66, 307, 'UIO')
            atoms_uio66 = add_H(atoms_uio66, 344, 'UIO')
            atoms_uio66 = add_H(atoms_uio66, 320, 'UIO')
            atoms_uio66 = add_H(atoms_uio66, 385, 'UIO')
            atoms_uio66 = add_H(atoms_uio66, 394, 'UIO')
            atoms_uio66 = add_H(atoms_uio66, 374, 'UIO')
            ids_to_remove.extend(zr_cluster_ids)
        if removed_count > 0:
           atoms_uio66 = add_OH(atoms_uio66, 306,   'UIO')
           atoms_uio66 = add_OH2(atoms_uio66, 319,  'UIO')
           atoms_uio66 = add_OH(atoms_uio66, 2643,  'UIO')
           atoms_uio66 = add_OH2(atoms_uio66, 2611, 'UIO')
           ids_to_remove.extend(linker1_ids)
        if removed_count > 1:
            atoms_uio66 = add_OH(atoms_uio66, 343,  'UIO')
            atoms_uio66 = add_OH2(atoms_uio66, 306,  'UIO')
            atoms_uio66 = add_OH(atoms_uio66, 3388,  'UIO')
            atoms_uio66 = add_OH2(atoms_uio66, 3407, 'UIO')
            ids_to_remove.extend(linker2_ids)
        if removed_count > 2:
            atoms_uio66 = add_OH(atoms_uio66, 319,  'UIO')
            atoms_uio66 = add_OH2(atoms_uio66, 343, 'UIO')
            atoms_uio66 = add_OH(atoms_uio66, 1651,  'UIO')
            atoms_uio66 = add_OH2(atoms_uio66, 1642, 'UIO')
            ids_to_remove.extend(linker3_ids)
        
        # Remove linkers / Zr cluster
        atoms_uio66 = remove_atoms(atoms_uio66, ids_to_remove)
        atoms_uio66 = remove_non_bonded_atoms(atoms_uio66)
        
        # Write .itp files
        mass, charge, bond_params, angle_params, dihedral_params = get_uio66_params()
        if remove_zr_cluster:
            cluster_fl = 'cl'
        else:
            cluster_fl = ''
        write_atoms(atoms_uio66, charge, mass, 'UIO', f'__tmp/uio66_atoms{removed_count}{cluster_fl}.itp')
        write_bonds(atoms_uio66, bond_params,         f'__tmp/uio66_bonds{removed_count}{cluster_fl}.itp')
        write_angles(atoms_uio66, angle_params,       f'__tmp/uio66_angles{removed_count}{cluster_fl}.itp')
        write_dihedrals(atoms_uio66, dihedral_params, f'__tmp/uio66_dihedrals{removed_count}{cluster_fl}.itp') 
        compose_itp_files(['../uio66/atomtypes.itp', '../uio66/moleculetype.itp',
                           f'__tmp/uio66_atoms{removed_count}{cluster_fl}.itp', f'__tmp/uio66_bonds{removed_count}{cluster_fl}.itp', f'__tmp/uio66_angles{removed_count}{cluster_fl}.itp', f'__tmp/uio66_dihedrals{removed_count}{cluster_fl}.itp'],
                           f'uio66_{removed_count}{cluster_fl}.itp')
        
        # Write .gro file
        atoms_uio66_tempo_removed = atoms_uio66.copy()
        atoms_uio66_tempo_removed.extend(atoms_tempo.copy())
        if removed_count > 0:
            atoms_uio66_tempo_removed.extend(atoms_linker[0].copy())
        if removed_count > 1:
            atoms_uio66_tempo_removed.extend(atoms_linker[1].copy())
        if removed_count > 2:
            atoms_uio66_tempo_removed.extend(atoms_linker[2].copy())
        write_gro_file(atoms_uio66_tempo_removed, f'uio66_{removed_count}{cluster_fl}.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
        

