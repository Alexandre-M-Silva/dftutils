import sys
import os
import numpy as np

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

def format_numeric_folder(root, i):
    if i <= 9:
        return os.path.join(root, "{:02d}".format(i))
    else:
        return os.path.join(root, "{:d}".format(i))

def folders_from_path(root):
    # TODO: This assumes paths are numbered in a numeric manner, make this smarter.
    min_f = sys.maxsize
    max_f = 0
    for subdir, dirs, files in os.walk(root):
        for subdir in dirs:
            if subdir.isnumeric():
                min_f = min(min_f, int(subdir))
                max_f = max(max_f, int(subdir))
        break

    folders = [format_numeric_folder(root, i) for i in range(min_f, max_f+1)]
    
    return folders

def outcars_from_folders(folders):
    return [Outcar(os.path.join(folder, 'OUTCAR')) for folder in folders]
    
def structures_from_folders(folders):
    folders_with_car = [os.path.join(folder, 'CONTCAR') if os.path.exists(os.path.join(folder, 'CONTCAR')) else os.path.join(folder, 'POSCAR') for folder in folders]    
    return [Structure.from_file(folder) for folder in folders_with_car]

def outcars_and_structures_from_path(root):
    folders = folders_from_path(root)
    outcars = outcars_from_folders(folders)
    structures = structures_from_folders(folders)
    return outcars, structures

def distance_between_structures(a, b):
    dist_cum = 0
    for si, sj in zip(a.sites, b.sites):
        dist_cum += si.distance(sj)
    return dist_cum

def match_structure_indices(path_a, path_b):
    sta = Structure.from_file(path_a)
    stb = Structure.from_file(path_b)

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

    sta.to(path_a + "_sorted", fmt='poscar')
    stb.to(path_b + "_sorted", fmt='poscar')

    
def interp_from_structures(structures, n):
    interp_structures = []
    ts = np.linspace(0, len(structures)-1, n)
    for t in ts:
        tmin = int(np.floor(t))
        tmax = int(np.ceil(t))
        ti = t - float(tmin)

        s = structures[tmin].copy()
        s.lattice = Lattice.from_parameters(structures[tmin].lattice.parameters + ti*(structures[tmax].lattice.parameters - structures[tmin].lattice.parameters))
        s.coords = structures[tmin].coords + ti*(structures[tmax].coords - structures[tmin].coords)
        interp_structures.append(s)

    return interp_structures
