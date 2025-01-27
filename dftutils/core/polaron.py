import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pymatgen.core import Structure
from pymatgen.core.composition import Element
from pymatgen.analysis.local_env import get_neighbors_of_site_with_index

from monty.serialization import loadfn, dumpfn

from dftutils.core.utils import *

class Polaron:
    def __init__(self, reference, structure, sites):
        self.reference = reference
        self.structure = structure
        self.sites = sites

    def from_reference_and_site(reference, site):
        return Polaron(reference, None, [site])
    
    def from_reference_path_and_site(reference_path, site):
        return Polaron(Structure.from_file(reference_path), None, [site])

    def from_dict(d):
        return Polaron(Structure.from_dict(d["reference"]), Structure.from_dict(d["structure"]), d["sites"])
    
    def to_dict(self):
        return {"reference": self.reference.as_dict() if self.reference is not None else None, 
                "structure": self.structure.as_dict() if self.structure is not None else None, 
                "sites": self.sites}

    def from_yaml(path):
        d = loadfn(path)
        return Polaron.from_dict(d)
    
    def from_yaml_and_structure(path, structure_path):
        d = loadfn(path)
        p = Polaron.from_dict(d)
        p.structure = Structure.from_file(structure_path)
        return p
    
    def to_yaml(self, path):
        dumpfn(self.to_dict(), path)

    def _apply_bdm_to_site_index(self, site_index, amount):
        if self.structure is None:
            self.structure = self.reference

        structure = self.reference
        site = structure.sites[site_index]
        cs = get_neighbors_of_site_with_index(structure, site_index, "voronoi", 0.1, 2.0)
        for c in cs:
            if c.species.elements[0] == Element("O"):
                v = site.frac_coords + (c.frac_coords - site.frac_coords) * (1.0 + amount)
                structure.replace(c.index, "O", coords = v)

    def apply_bdm(self, amount):
        for site in self.sites:
            self._apply_bdm_to_site_index(site, amount)

    def write_structure(self, path):
        self.structure.to(path, fmt="poscar")

    def write_reference(self, path):
        self.reference.to(path, fmt="poscar")

    def print_bond_lengths_of_structure_site(self, structure, site_index):
        site = structure.sites[site_index]
        cs = get_neighbors_of_site_with_index(structure, site_index, "voronoi", 0.1, 2.0)
        for c in cs:
            if c.species.elements[0] == Element("O"):
                bond_length = site.distance(structure.sites[c.index])
                print(f"{c.index:3d} - {bond_length:10.2f} Angstrom")

    def print_structure_bond_lengths(self):
        for site_index in self.sites:
            self.print_bond_lengths_of_structure_site(self.structure, site_index)
    
    def print_reference_bond_lengths(self):
        for site_index in self.sites:
            self.print_bond_lengths_of_structure_site(self.reference, site_index)
    
    def print_bond_lengths_comparison(self):
        if self.structure is None:
            self.structure = self.reference

        for site_index in self.sites:
            site1 = self.structure.sites[site_index]
            site2 = self.reference.sites[site_index]
            cs = get_neighbors_of_site_with_index(self.structure, site_index, "voronoi", 0.1, 2.0)
            for c in cs:
                if c.species.elements[0] == Element("O"):
                    bond_length_1 = site1.distance(self.structure.sites[c.index])
                    bond_length_2 = site2.distance(self.reference.sites[c.index])
                    print(f"{c.index:3d}: {bond_length_1:10.4f} - {bond_length_1:10.4f} = {bond_length_1 - bond_length_2:10.4f} Angstrom")
