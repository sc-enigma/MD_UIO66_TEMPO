import sys
import numpy as np
import pickle

sys.path.append('../../../MD_ZIF7_TEMPO/system/components')

from atom import Atom
from pbc import calculate_lattice_vectors

def define_uio66_atom_types(atoms):
    def getNeighbours(atom_idx):
        return [atoms[adj_idx] for adj_idx in atoms[atom_idx].adjacency]
    
    def getNeighbourElems(atom_idx):
        return np.sort([atom.name[0] for atom in getNeighbours(atom_idx)])
    
    def getElem(atom_idx):
        return atoms[atom_idx].name[0]
        
    # define Zr: flexFF_Zr
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'Z':
            atoms[atom_idx].atom_type = 'flexFF_Zr'
            
    # define H: flexFF_H1
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'H':
            atoms[atom_idx].atom_type = 'flexFF_H1'
  
    # define O1: flexFF_O1
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'O' and atoms[atom_idx].atom_type == 'unk_type':
            if (np.array_equal(getNeighbourElems(atom_idx), ['C', 'Z'])):
                atoms[atom_idx].atom_type = 'flexFF_O1'
    
    # define O3: flexFF_O3
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'O' and atoms[atom_idx].atom_type == 'unk_type':
            if (np.array_equal(getNeighbourElems(atom_idx), ['Z', 'Z', 'Z'])):
                atoms[atom_idx].atom_type = 'flexFF_O3'
    
    # define C3: flexFF_C3
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'C' and atoms[atom_idx].atom_type == 'unk_type':
            if (np.array_equal(getNeighbourElems(atom_idx), ['C', 'C', 'H'])):
                atoms[atom_idx].atom_type = 'flexFF_C3'
    
    # define C2: flexFF_C2
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'C' and atoms[atom_idx].atom_type == 'unk_type':
            if (np.array_equal(getNeighbourElems(atom_idx), ['C', 'C', 'C'])):
                atoms[atom_idx].atom_type = 'flexFF_C2'
            
    # define C1: flexFF_C1
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'C' and atoms[atom_idx].atom_type == 'unk_type':
            if (np.array_equal(getNeighbourElems(atom_idx), ['C', 'O', 'O'])):
                atoms[atom_idx].atom_type = 'flexFF_C1'
            
    # check that all atom types defined
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].atom_type == 'unk_type':
            # print('ERROR')            
            # print(getElem(atom_idx), getNeighbourElems(atom_idx))
            atoms[atom_idx].atom_type = 'flexFF_C1'
            
    return atoms

def define_uio66_atom_names(atoms):  
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].name = atoms[atom_idx].atom_type.replace('flexFF_', '')
        atoms[atom_idx].resid_idx = 1
        atoms[atom_idx].resid = 'UIO'
        atoms[atom_idx].atom_idx = atom_idx
        
    return atoms

def remove_Zr_Zr_bond(atoms):
    for atom_idx in range(len(atoms)):
        if not 'Zr' in atoms[atom_idx].name:
            continue
        adj_cnt = 0
        while adj_cnt < len(atoms[atom_idx].adjacency):
            adj_idx = atoms[atom_idx].adjacency[adj_cnt]
            if 'Zr' in atoms[adj_idx].name:
                del(atoms[atom_idx].adjacency[adj_cnt])
            else:
                adj_cnt += 1

    return atoms

def remove_periodic_bonds(atoms):
    for atom_idx in range(len(atoms)):
            adj_cnt = 0
            while adj_cnt < len(atoms[atom_idx].adjacency):
                adj_idx = atoms[atom_idx].adjacency[adj_cnt]
                dist = pow(sum(np.power(atoms[atom_idx].r - atoms[adj_idx].r, 2)), 0.5)
                if dist > 10.0:
                    del(atoms[atom_idx].adjacency[adj_cnt])
                else:
                    adj_cnt += 1
    return atoms
