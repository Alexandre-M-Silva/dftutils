import sys
import os
import numpy as np

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.core.structure import Structure

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

    folders = []
    for i in range(min_f, max_f+1):
        if i <= 9:
            folders.append(os.path.join(root, "{:02d}".format(i)))
        else:
            folders.append(os.path.join(root, "{:d}".format(i)))
    
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

def match_structure_indices(path_a, path_b, output_path):
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

    stb.to(output_path, fmt='poscar')