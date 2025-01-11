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

    datas = []
    for i in range(0, 3):
        data = []
        for p, q in zip(ps, qs):
            branch_data = []
            for n in range(bmax, bmin-1, -1):
                branch_data.append(p[i]+n*q[i])
            data.append(branch_data)
        datas.append(data)

    return [pd.DataFrame(data).transpose() for data in datas]

def polarization_branch(data, start, bias=0):
    n = len(data.columns)
    P = np.zeros(n)
    P[0] = start
    for i in range(0, n-1):
        values = data.iloc[:, i+1].values
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
        
    def get_branches_and_switch(self, axis=2):
        starts = self.data[axis].iloc[:, 0].values
        branches = [polarization_branch(self.data[axis], start) for start in starts]
        branches_bias = [polarization_branch_derivative_bias(branch) for branch in branches]
        branches = [polarization_branch(self.data[axis], start, bias) for start, bias in zip(starts, branches_bias)]

        abs_means = [np.abs(np.mean(branch)) for branch in branches]
        i = np.argmin(abs_means)
        switch = branches[i]

        return branches, switch

    def get_spontaneous(self, switch=None, axis=2):
        if switch is None:
            _, switch = self.get_branches_and_switch(axis)
        
        return 0.5*(switch[-1] - switch[0])

    def plot(self, path, axis=2, raw=False):
        fig, ax = plt.subplots()
        
        df = self.data[axis]
        for series_name, series in df.items():
            y = series
            x = np.full(shape=y.shape, fill_value=float(series_name))
            ax.scatter(x, y, color='black', s=2)

        if not raw:
            _, switch = self.get_branches_and_switch(axis)
            ax.plot(np.arange(0, len(df.columns)), switch, 'o-', color='red')

            Ps = self.get_spontaneous(switch)
            ax.text(0.01, 0.95, rf'Ps = {Ps:6.2f} $\mu C/cm^2$', fontsize=12, color='red',
                    transform = ax.transAxes)

        ax.set_xlabel("Image")
        ax.set_ylabel(r"$P\ (\mu C/cm^2)$")
        fig.tight_layout()
        fig.savefig(path,
                    bbox_inches='tight', 
                    pad_inches=0.05)
    