import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess as sp
from io import StringIO
from pymatgen.core import Structure
import re

from dftutils.utils import *

eV = 1.60218E-19
Angstrom = 1e-10

def parse_polarization_from_outcar(path):
    pion = None
    pelc = None

    with open(path, 'r') as f:
        for line in f:
            if 'Total electronic dipole moment' in line:
                pelc =  np.array(list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line))))
            elif 'Ionic dipole moment' in line:
                pion = np.array(list(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", line))))\
    
    if isinstance(pion, np.ndarray) and isinstance(pelc, np.ndarray):
        return pelc, pion
    else:
        raise ValueError("Could not find polarization in OUTCAR.")

def polarization_from_outcar_structure(outcar_path, structure_path):
    """
    Obtain polarization magnitude and quanta in uC cm^-2.
    """
    pelc, pion = parse_polarization_from_outcar(outcar_path)

    if isinstance(structure_path, list):
        structure_path = next((sp for sp in structure_path if os.path.exists(sp)), structure_path[0])

    structure = Structure.from_file(structure_path)
    vol = structure.volume
    abc = np.diag(structure.lattice.matrix)

    ptot = pelc + pion
    ptotv = ptot / vol
    quanta = abc / vol
    ratio = ptotv / quanta
    n = np.abs(np.round(ratio))
    dp = ptotv - n * np.sign(ratio) * quanta
    dp *= 100 * eV / (Angstrom ** 2)  # uC/cm^2
    quanta *= 100 * eV / (Angstrom ** 2)

    return dp, quanta
    
def polarization_scatter_from_path(path, bmin=-5, bmax=5):
    """
    Obtain polarization scatter (multiple branches) in uC cm^-2, 
    for a given interval of branches.
    """
    folders = folders_from_path(path)

    ps, qs = zip(*[
        polarization_from_outcar_structure(
            os.path.join(folder, "OUTCAR"),
            [os.path.join(folder, "CONTCAR"), os.path.join(folder, "POSCAR")]
        )
        for folder in folders
    ])

    data = [
        [i, p[0] + b*q[0], p[1] + b*q[1], p[2] + b*q[2]]
        for i, (p, q) in enumerate(zip(ps, qs))
        for b in range(bmax, bmin-1, -1)
    ]

    df = pd.DataFrame(data, columns=['Image', 'Px', 'Py', 'Pz'])

    return df

def polarization_branch(data, axis, start, bias=0):
    n = data.iloc[:, 0].nunique()
    P = np.zeros(n)
    P[0] = start
    for i in range(0, n-1):
        values = data[data["Image"]==i+1].iloc[:, axis+1].values
        dists = np.abs((values - bias) - P[i]) 
        closest_at = np.argmin(dists)
        P[i+1] = values[closest_at]
    
    return P

def polarization_branch_derivative_bias(branch):
    grad = np.gradient(branch)
    return np.sum(grad) / len(grad)

class Polarization:
    def __init__(self, path=None):
        self.path = path
        self.data = polarization_scatter_from_path(path)
        self.branches = None
        self.switch = None

    def _calc_branches_and_switch(self, axis=2):
        if self.branches is None:
            starts = self.data[self.data["Image"] == 0].iloc[:, axis+1]
            self.branches = [polarization_branch(self.data, axis, start) for start in starts]
            branches_bias = [polarization_branch_derivative_bias(branch) for branch in self.branches]
            self.branches = [polarization_branch(self.data, axis, start, bias) for start, bias in zip(starts, branches_bias)]

        if self.switch is None:
            abs_means = [np.abs(np.mean(branch)) for branch in self.branches]
            i = np.argmin(abs_means)
            self.switch = self.branches[i]  

    def get_spontaneous(self, axis=2):
        self._calc_branches_and_switch(axis)
        return 0.5*(self.switch[-1] - self.switch[0])
    
class PolarizationPlotter:
    def __init__(self, pol):
        self.pol = pol

    def from_path(path):
        pol = Polarization(path)
        return PolarizationPlotter(pol)
    
    def from_pol(pol):
        return PolarizationPlotter(pol)
    
    def plot(self, axis=2, save=True):
        use_matplotlib_style()

        fig, ax = plt.subplots()
        ax.set_xlabel('Image')
        ax.set_ylabel(r'$P_s\ (\mu C/cm^2)$')
        ax.scatter(self.pol.data['Image'], self.pol.data.iloc[:, axis+1], s=0.5, color='black')
    
        self.pol._calc_branches_and_switch(axis)
        for b in self.pol.branches:
            ax.plot(b, "o-")

        if save:
            fig.savefig(os.path.join(self.pol.path, 'plot.png'), bbox_inches='tight', pad_inches=0.05)

        return fig, ax


