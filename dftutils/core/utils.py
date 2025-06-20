import sys
import os
import numpy as np

from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from scipy.linalg import polar

import matplotlib.pyplot as plt

def use_matplotlib_style():
    plt.style.use(os.path.join(os.path.dirname(__file__), 'dftutils.mplstyle'))

def format_numeric_folder(root, i):
    return os.path.join(root, f"{i:02d}" if i <= 9 else f"{i}")

def folders_from_path(root):
    min_f, max_f = sys.maxsize, 0
    for _, dirs, _ in os.walk(root):
        numeric_dirs = [int(d) for d in dirs if d.isnumeric()]
        if numeric_dirs:
            min_f, max_f = min(numeric_dirs), max(numeric_dirs)
        break

    return [format_numeric_folder(root, i) for i in range(min_f, max_f + 1)]

def outcars_from_folders(folders):
    return [Outcar(os.path.join(folder, 'OUTCAR')) for folder in folders]

def structures_from_folders(folders):
    return [Structure.from_file(os.path.join(folder, 'CONTCAR') if os.path.exists(os.path.join(folder, 'CONTCAR')) else os.path.join(folder, 'POSCAR')) for folder in folders]

def outcars_and_structures_from_path(root):
    folders = folders_from_path(root)
    return outcars_from_folders(folders), structures_from_folders(folders)

def distance_between_structures(a, b):
    return sum(si.distance(sj) for si, sj in zip(a.sites, b.sites))

def match_indices_from_structs(sta, stb):
    pairings = []

    for si, sa in enumerate(sta.sites):
        clj = -1
        mind = 10000000
        for sj, sb in enumerate(stb.sites):
            d = sa.distance(sb)
            if d <= mind:
                mind = d
                clj = sj
        pairings.append((si, clj))

    stc = stb.copy()

    for si, sb in enumerate(stb.sites):
        stb.sites[pairings[si][0]] = stc.sites[pairings[si][1]]

def match_indices_from_paths(path_a, path_b):
    sta = Structure.from_file(path_a)
    stb = Structure.from_file(path_b)

    match_indices_from_structs(sta, stb)

    sta.to(path_a + "_sorted", fmt='poscar')
    stb.to(path_b + "_sorted", fmt='poscar')

def interp_from_structures(structures: list[Structure], 
                           n: int):
    interp_structures = []
    ts = np.linspace(0, 1, n)*(len(structures)-1)
    for t in ts:
        tmin = int(np.floor(t))
        tmax = int(np.ceil(t))
        ti = 0.0
        if tmin != tmax:
            ti = (t-float(tmin))/float(tmax-tmin)
        elif tmin == tmax and tmin == len(structures) - 1 and tmax == len(structures) - 1:
            ti = 1.0
        else:
            ti = 0.0
        
        s0 = structures[tmin]
        s1 = structures[tmax]
        
        start_coords = np.array(s0.frac_coords)
        end_coords = np.array(s1.frac_coords)

        if (end_coords - start_coords >= 0.8).any():
            end_coords[(end_coords-start_coords >= 0.8)] = - (1.0 - end_coords[(end_coords-start_coords >= 0.8)])

        if (end_coords - start_coords <= -0.8).any():
            end_coords[(end_coords-start_coords <= -0.8)] = (1.0 - end_coords[(end_coords-start_coords <= -0.8)])

        vec = ti * (end_coords - start_coords)
        
        _u, p = polar(np.dot(s1.lattice.matrix.T, np.linalg.inv(s0.lattice.matrix.T)))
        lvec = ti * (p - np.identity(3))
        lstart = s0.lattice.matrix.T
        
        l_a = np.dot(np.identity(3) + lvec, lstart).T
        lattice = Lattice(l_a)
        frac_coords = start_coords + vec

        new_s = s0.copy()
        new_s.lattice = lattice
        for s, fc in zip(new_s.sites, frac_coords):
            s.frac_coords = fc
        interp_structures.append(new_s)


    return interp_structures
