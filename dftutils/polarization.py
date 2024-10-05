import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import subprocess as sp
from io import StringIO
from pymatgen.core import Structure

from dftutils.utils import *

def pion_vec(file):
    return pd.read_csv(StringIO(sp.run(["grep 'Ionic dipole moment' " + file + " | awk '{print $5,$6,$7}'"], shell=True, capture_output=True).stdout.decode('utf-8').replace('*', '')), header=None, sep=" ", on_bad_lines='skip').values[0]
      
def pelc_vec(file):
    return pd.read_csv(StringIO(sp.run(["grep 'Total electronic dipole moment' " + file + " | awk '{print $6,$7,$8}'"], shell=True, capture_output=True).stdout.decode('utf-8').replace('*', '')), header=None, sep=" ", on_bad_lines='skip').values[0]
    
def polarization_from_path(path, bmin=-5, bmax=5):
    """
    Obtain polarization points (multiple branches) in uC cm^-2, 
    for a given interval of branches.
    """

    eV = 1.60218E-19
    Angstrom = 1e-10

    folders = folders_from_path(path)
    structures = structures_from_folders(folders)

    pdata_e = []
    pdata_i = []
    for folder in folders:
        pdata_e.append(pelc_vec(os.path.join(folder, "OUTCAR")))
        pdata_i.append(pion_vec(os.path.join(folder, "OUTCAR")))
            
    branch_count = bmax - bmin + 1

    df = pd.DataFrame(columns=['Image', 'Px', 'Py', 'Pz'])
    for branch_n in range(bmin, bmax+1):
        for i in range(0, len(structures)):
            structure = structures[i]

            if structure != None:
                vol = structure.volume
                abc = np.array([structure.lattice.matrix[0][0],
                               structure.lattice.matrix[1][1],
                               structure.lattice.matrix[2][2]])
                pelc = np.array(pdata_e[i])
                pion = np.array(pdata_i[i])
                ptot = pelc + pion
                ptotv = ptot/vol
                quanta = abc/vol
                ratio = ptotv / quanta
                n = np.abs(np.round(ratio)) + branch_n
                dp = ptotv-n*np.sign(ratio)*quanta
                dp *= 100*eV/(Angstrom**2) # uC/cm^2
                df.loc[len(df.index)] = [i, dp[0], dp[1], dp[2]]

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
        self.data = polarization_from_path(path)
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


