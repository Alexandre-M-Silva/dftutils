import sys
import os
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

from dftutils.core.matching import StructureMatcher

from dftutils.core.utils import *

def _format_folder_name(i):
    return f"{i:02d}" if i <= 9 else f"{i}"

class Neb:
    def __init__(self, structures):
        self.structures = structures

    def to_path(self,
                path: str):
        folders = [_format_folder_name(i) for i in range(0, len(self.structures)+1)]
        if path is not None:
            folders = [os.path.join(path, f) for f in folders]

        for s, f in zip(self.structures, folders):
            if not os.path.exists(f):
                os.makedirs(f)
            s.to(os.path.join(f, "POSCAR"), fmt="poscar")

    def to_movie(self,
                 path: str):
        final_str = ""
        for i, s in enumerate(self.structures):
            final_str += s.to(fmt="poscar")

        final_path = os.path.join(path, "movie.vasp") if path is not None else "movie.vasp"
        with open(final_path, "w") as f:
            f.write(final_str)

    def from_initial_and_final(initial: Structure | str, 
                               final: Structure | str, 
                               n: int,
                               match: bool = False):
        if isinstance(initial, str):
            initial = Structure.from_file(initial)
        if isinstance(final, str):
            final = Structure.from_file(final)

        if match:
            match_config = {
                "type": "all",
                "selection": None,
                "sort_first": False,
            }
            matcher = StructureMatcher(initial, final, match_config)
            initial, final = matcher.match()

        structures = [initial, final]
        interp_structures = interp_from_structures(structures, n+2)
        return Neb(interp_structures)

    def from_path(path: str):
        folders = folders_from_path(path)
        if len(folders) <= 2:
            raise ValueError("Not enough numerically named folders on directory to extract a pathway.")
        
        structures = []
        for folder in folders:
            if os.path.exists(os.path.join(folder, "CONTCAR")):
                structures.append(Structure.from_file(os.path.join(folder, "CONTCAR")))
            elif os.path.exists(os.path.join(folder, "POSCAR")):
                structures.append(Structure.from_file(os.path.join(folder, "POSCAR")))
            else:
                raise Exception(f"Could not find CONTCAR or POSCAR files in {folder}.")
        return Neb(structures)   

    def interp(self, n: int, interp_type: str = 'linear'):
        self.structures = interpolate_structures(self.structures, n, interpolate_lattices=True, method=interp_type)   
    
class NebData:
    def __init__(self, name: str = None, path: str = None):
        if path == None:
            return
        
        outcars, structures = outcars_and_structures_from_path(path)

        ref_structure = structures[0]
        ref_energy = outcars[0].final_energy

        distances = [0]
        energies = [0]
        for i, v in enumerate(zip(outcars[1:], structures[1:])):
            distances.append(distance_between_structures(v[1], ref_structure))
            energies.append(v[0].final_energy-ref_energy)

        self.data = pd.DataFrame({"Image": range(len(structures)),
                                  "Distances": distances,
                                  "Distances01": distances/distances[-1],
                                  "Energies": energies})

        self.name = name
        if self.name == None:
            self.name = "{:x}".format(pd.util.hash_pandas_object(self.data))
        
class NebPlotter:
    def __init__(self):
        self.nebs = []

    def add_neb(self, neb_data: NebData):
        self.nebs.append(neb_data)

    def plot(self, name='plot', format='png', save=False):
        fig, ax = plt.subplots()
        ax.set_xlabel("Transition coordinate (a.u.)")
        ax.set_ylabel(r"$\Delta E$ (eV)")
        for neb in self.nebs:
            rx = neb.data["Distances01"]
            ry = neb.data["Energies"]

            spline = make_interp_spline(rx, ry)
            x = np.linspace(rx.min(), rx.max(), 100)
            y = spline(x)

            ax.plot(x, y, '-', alpha=1/3)
            ax.plot(neb.data["Distances01"], neb.data["Energies"], 'o', color=plt.gca().lines[-1].get_color(), label=neb.name)
        ax.legend(frameon=False, loc='upper right', bbox_to_anchor=(1.02, 1.0))
        
        if save:
            assert(format == 'png' or format == 'svg')
            fig.savefig(name + '.' + format, format=format, dpi=300, bbox_inches='tight')

