
import sys
import os
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core import Lattice

def format_strain_directory(strain):
    dir = f"{strain[0]:.6E},{strain[1]:.6E},{strain[2]:.6E}"
    dir.replace("+", "p")
    dir.replace("-", "m")
    dir.replace(".", "_")
    return dir

def apply_strain(structure, strain):
    strain_matrix = np.matrix([[strain[0], 0, 0], [0, strain[1], 0], [0, 0, strain[2]]])
    new_lattice = Lattice(np.dot(structure.lattice.matrix.T, strain_matrix).T)
    new_structure = structure.copy()
    new_structure.lattice = new_lattice
    return new_structure