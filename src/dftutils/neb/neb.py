import sys
import os
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

import dftutils.core.utils as utils

class Neb:
    def __init__(self, name=None, path=None):
        if path == None:
            return
        
        outcars, structures = utils.outcars_and_structures_from_path(path)

        ref_structure = structures[0]
        ref_energy = outcars[0].final_energy

        distances = [0]
        energies = [0]
        for i, v in enumerate(zip(outcars[1:], structures[1:])):
            distances.append(utils.distance_between_structures(v[1], ref_structure))
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

    def add_neb(self, neb: Neb):
        self.nebs.append(neb)

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
        ax.legend(frameon=False, loc='upper right', bbox_anchor=(1.02, 1.0))
        
        if save:
            assert(format == 'png' or format == 'svg')
            fig.savefig(name + '.' + format, format=format, dpi=300, bbox_inches='tight')

