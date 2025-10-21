import sys
import os
import math
import re
import warnings
from typing import List, Dict, IO

import numpy as np
import importlib.resources
import shutil
from pathlib import Path
import pandas as pd

import matplotlib.pyplot as plt

from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from scipy.ndimage import gaussian_filter1d
from scipy.linalg import polar

def use_matplotlib_style():
    with importlib.resources.path('dftutils.core', 'dftutils.mplstyle') as style_path:
        plt.style.use(str(style_path))

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

def get_gaussian_smeared(energies: np.ndarray, 
                         dos: dict[str, np.ndarray],
                         sigma: float):
     
    diff = [energies[idx + 1] - energies[idx] for idx in range(len(energies) - 1)]
    avg_diff = sum(diff) / len(diff)
    return {spin: gaussian_filter1d(dens, sigma / avg_diff) for spin, dens in dos.items()}


def potcar_zvals(potcar: str | IO) -> Dict[str, int]:
    """
    Parse a POTCAR file and return {element: ZVAL}.
    """
    zvals: Dict[str, int] = {}

    if isinstance(potcar, str):
        potcar = open(potcar, "r", errors="ignore")

    lines = potcar.readlines()

    for i, line in enumerate(lines):
        if line.strip().startswith("PAW_"):
            m = re.search(r"PAW_[A-Z]+(?:\s+)([A-Za-z]{1,2})", line)
            if not m:
                continue
            element = m.group(1).capitalize()

            j = i + 1
            while j < len(lines) and lines[j].strip() == "":
                j += 1
            if j >= len(lines):
                continue

            token = lines[j].strip().split()[0]
            try:
                zval = float(token)
            except ValueError:
                m = re.search(r"([0-9]+(?:\.[0-9]+)?)", lines[j])
                if not m:
                    raise ValueError(f"Could not find numeric ZVAL after line {j+1}")
                zval = float(m.group(1))

            zvals[element] = zval

    if not zvals:
        raise ValueError("No ZVALs found in POTCAR.")
    return zvals

def structure_symbols_and_counts(structure: str | Structure) -> Dict[str, int]:
    structure = structure if isinstance(structure, Structure) else Structure.from_file(structure)
    symbols, counts = np.unique([s.species.elements[0].symbol for s in structure.sites], return_counts=True)
    return dict(zip(symbols, counts))

def nelect_from_structure_potcar(structure: str | Structure, potcar: str | IO) -> int:
    symbols_counts = structure_symbols_and_counts(structure)
    zvals = potcar_zvals(potcar)

    missing_in_potcar = [el for el in symbols_counts if el not in zvals]
    extra_in_potcar = [el for el in zvals if el not in symbols_counts]
    if missing_in_potcar:
        raise ValueError(f"Elements {missing_in_potcar} in structure not found in POTCAR!")
    if extra_in_potcar:
        print(f"Warning: Extra elements in POTCAR not in structure: {extra_in_potcar}")

    return sum(zvals[el] * symbols_counts[el] for el in symbols_counts)

def kpoints_count(kpoints: str | Kpoints) -> int:
    kp = kpoints if isinstance(kpoints, Kpoints) else Kpoints.from_file(kpoints)
    style = kp.style.name.lower()
    if style in ("monkhorst", "gamma", "automatic"):
        kx, ky, kz = map(int, kp.kpts[0])
        return kx * ky * kz
    return len(kp.kpts)

def divisors(n: int) -> List[int]:
    a = np.arange(1, int(np.sqrt(n)) + 1)
    divs = a[n % a == 0]
    return np.unique(np.concatenate((divs, n // divs))).tolist()

def estimate_nbands_upper_end(
    nelect: float,
    nodes: int,
    cores_per_node: int,
    upper_end_mult: float = 2.5,
) -> int:
    """
    Estimate NBANDS from NELECT with a large margin
    """
    total_ranks = nodes * cores_per_node
    nbands = int(math.ceil((nelect / 2) * upper_end_mult))
    nbands = int(math.ceil(nbands / total_ranks) * total_ranks)
    return nbands

def estimate_nbands(
    nelect: float,
    nodes: int,
    cores_per_node: int,
    ncore: int,
    margin_frac: float = 0.2,
) -> tuple[int, int, int, float]:
    """
    Estimate NBANDS from NELECT and ensure it's a multiple of total ranks.
    Returns (nbands_adj, total_ranks, ranks_per_node, bands_per_rank)
    """
    if ncore <= 0:
        raise ValueError("ncore must be > 0")

    if cores_per_node % ncore != 0:
        ranks_per_node = cores_per_node // ncore
    else:
        ranks_per_node = cores_per_node // ncore

    total_ranks = nodes * ranks_per_node
    nbands_est = int(math.ceil((nelect / 2) * (1.0 + margin_frac)))
    nbands_adj = int(math.ceil(nbands_est / total_ranks) * total_ranks)
    bands_per_rank = nbands_adj / total_ranks
    return nbands_adj, total_ranks, ranks_per_node, bands_per_rank

