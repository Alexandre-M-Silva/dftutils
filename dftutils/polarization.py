import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pymatgen.core as mg
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

class PolarizationPlotter:
    def __init__(self, path=None):
        self.path = path
        self.data = polarization_from_path(path)
        self.branches = None
        self.switch = None
    
    def get_spontaneous(self):
        return 0.5*(self.switch[-1] - self.switch[0])
    
    def plot(self, ylim=None, save=True):
        plt.style.use(os.path.join(os.path.dirname(__file__), 'dftutils.mplstyle'))

        fig, ax = plt.subplots()
        ax.set_xlabel('Image')
        ax.set_ylabel(r'$P_s\ (\mu C/cm^2)$')
        ax.scatter(self.data['Image'], self.data.iloc[:, 3], s=0.5, color='black')
    
        if not ylim is None:
            #self.branches = branches_from_polarization(self.data)
            #self.switch = midpoint_branch_from_branches(self.branches)
            ax.set_ylim(ylim)

        if save:
            fig.savefig(os.path.join(self.path, 'plot.png'), bbox_inches='tight', pad_inches=0.05)

        return fig, ax


