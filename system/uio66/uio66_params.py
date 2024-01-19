def get_uio66_params():
    # [atom_type] = m
    mass = {}
    mass['flexFF_Zr'] = 91.224
    mass['flexFF_C1'] = 12.011
    mass['flexFF_C2'] = 12.011
    mass['flexFF_C3'] = 12.011
    mass['flexFF_O1'] = 16.000
    mass['flexFF_O3'] = 16.000
    mass['flexFF_H1'] =  1.008
    mass['repl_H']    =  1.008
    mass['repl_O']    = 16.000
 
    # [atom_type] = q
    charge = {}
    charge['flexFF_Zr'] = 1.968
    charge['flexFF_C1'] = 0.630
    charge['flexFF_C2'] = -0.082
    charge['flexFF_C3'] = -0.065
    charge['flexFF_O1'] = -0.533
    charge['flexFF_O3'] = -0.902
    charge['flexFF_H1'] =  0.133
    charge['repl_H']    =  0.4
    charge['repl_O']    = -0.8

    # [atom_type-atom_type] = [funct, r0, k]
    bond_params = {}
    bond_params['flexFF_Zr-flexFF_O3'] = [1, 0.2098, 107733.800]
    bond_params['flexFF_Zr-flexFF_O1'] = [1, 0.2232, 287290.200]
    bond_params['flexFF_C1-flexFF_O1'] = [1, 0.1273, 451872.000]
    bond_params['flexFF_C1-flexFF_C2'] = [1, 0.1487, 293928.300]
    bond_params['flexFF_C2-flexFF_C3'] = [1, 0.1393, 401664.000]
    bond_params['flexFF_C3-flexFF_C3'] = [1, 0.1393, 401664.000]
    bond_params['flexFF_C3-flexFF_H1'] = [1, 0.1080, 304106.800]
    
    bond_params['flexFF_Zr-repl_O'] = [1, 0.2098, 107733.800]
    bond_params['flexFF_O1-repl_H'] = [1, 0.0960, 462750.000]
    bond_params['repl_O-repl_H']    = [1, 0.0960, 462750.000]
    
    # [atom_type-atom_type-atom_type] = [funct, angle, k]
    angle_params = {}
    angle_params['flexFF_O3-flexFF_Zr-flexFF_O3']  = [1,  71.100, 115.776]
    angle_params['flexFF_O3-flexFF_Zr-flexFF_O3']  = [1, 100.400, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O1']  = [1,  76.000, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O1']  = [1,  81.900, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O1']  = [1,  69.800, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O1']  = [1,  74.400, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O3']  = [1,  85.700, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O3']  = [1, 124.300, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O3']  = [1, 127.600, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O3']  = [1, 143.500, 115.776]
    angle_params['flexFF_O1-flexFF_Zr-flexFF_O3']  = [1, 150.500, 115.776]
    angle_params['flexFF_Zr-flexFF_O1-flexFF_C1']  = [1, 135.800, 231.637]
    angle_params['flexFF_O1-flexFF_C1-flexFF_O1']  = [1, 125.000,1213.360]
    angle_params['flexFF_O1-flexFF_C1-flexFF_C2']  = [1, 117.300, 456.013]
    angle_params['flexFF_C1-flexFF_C2-flexFF_C3']  = [1, 120.000, 290.201]
    angle_params['flexFF_C3-flexFF_C2-flexFF_C3']  = [1, 120.000, 753.120]
    angle_params['flexFF_C2-flexFF_C3-flexFF_C3']  = [1, 120.000, 753.120]
    angle_params['flexFF_C2-flexFF_C3-flexFF_H1']  = [1, 120.000, 309.616]
    angle_params['flexFF_C3-flexFF_C3-flexFF_H1']  = [1, 120.000, 309.616]
    angle_params['flexFF_Zr-repl_O-repl_H']        = [1, 108.500, 372.560]

    # [atom_type-atom_type-atom_type-atom_type] = [funct, angle, k, n]        - periodic
    # [atom_type-atom_type-atom_type-atom_type] = [funct, c1, c2, c3, c4, c5] - fourier
    dihedral_params = {}
    dihedral_params['flexFF_Zr-flexFF_O1-flexFF_C1-flexFF_C2'] = [9, 180.000, 86.837, 2]
    dihedral_params['flexFF_O1-flexFF_C1-flexFF_C2-flexFF_C3'] = [9, 180.000, 10.460, 2]
    dihedral_params['flexFF_C1-flexFF_C2-flexFF_C3-flexFF_C3'] = [9, 180.000, 12.552, 2]
    dihedral_params['flexFF_C1-flexFF_C2-flexFF_C3-flexFF_H1'] = [9, 180.000, 12.552, 2]
    dihedral_params['flexFF_C2-flexFF_C3-flexFF_C3-flexFF_C2'] = [9, 180.000, 12.552, 2]
    dihedral_params['flexFF_C2-flexFF_C3-flexFF_C3-flexFF_H1'] = [9, 180.000, 12.552, 2]
    dihedral_params['flexFF_C3-flexFF_C3-flexFF_C2-flexFF_C3'] = [9, 180.000, 12.552, 2]
    dihedral_params['flexFF_H1-flexFF_C3-flexFF_C2-flexFF_C3'] = [9, 180.000, 12.552, 2]
    dihedral_params['flexFF_H1-flexFF_C3-flexFF_C3-flexFF_H1'] = [9, 180.000, 12.552, 2]
    dihedral_params['flexFF_C2-flexFF_C1-flexFF_O1-flexFF_O1'] = [9, 180.000, 41.840, 2]
    dihedral_params['flexFF_C1-flexFF_C2-flexFF_C3-flexFF_C3'] = [9, 180.000, 41.840, 2]
    dihedral_params['flexFF_C2-flexFF_C3-flexFF_C3-flexFF_H1'] = [9, 180.000,  1.548, 2]

    return mass, charge, bond_params, angle_params, dihedral_params