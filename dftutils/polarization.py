import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pymatgen.core as mg
import pandas as pd
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.analysis.ferroelectricity.polarization import Polarization
from pymatgen.core import Structure

from dftutils.utils import *

def polarization_from_path(path, bmin=-5, bmax=5):
    """
    Obtain polarization points (multiple branches) in uC cm^-2, 
    for a given interval of branches.
    """

    e = -1.60218E-19
    Angstrom = 1e-10
    eV = 1.60218E-19

    outcars, structures = outcars_and_structures_from_path(path)

    pdata_e = []
    pdata_i = []
    for outcar, structure in zip(outcars, structures):
        if outcar != None and structure != None:
            pdata = Polarization.from_outcars_and_structures([outcar], [structure]).get_pelecs_and_pions()
            pdata_e.append(pdata[0][0])     
            pdata_i.append(pdata[1][0])           
        else:
            pdata_e.append(None)
            pdata_i.append(None)
            
    branch_count = bmax - bmin + 1

    df = pd.DataFrame(columns=['Image', 'Px', 'Py', 'Pz'])
    for branch_n in range(bmin, bmax+1):
        for i in range(0, len(structures)):
            structure = structures[i]

            if structure != None:
                vol = structure.volume
                abc = np.array(structure.lattice.abc)
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

def branch_from_polarization(pol, axis=2, start=0):
    """
    Extracts the polarization switching branch from the polarization data and a starting polarization,
    assuming that its monotonically increasing.
    """
    branch = []
    nimages = pol['Image'].nunique()
    image = 0
    P = start if start != 0 else np.min(pol[pol["Image"] == 0].iloc[:, axis+1].values)
    for i in range(0, nimages):
        ps = np.sort(pol[pol['Image']==i].iloc[:, axis+1].values)
        if (P < ps).any():
            for p in ps:
                if p > P:
                    P = p
                    break
        else:
            break
            
        branch.append(P)
        
    if len(branch) == nimages:
        return np.array(branch)
    else:
        return None

def branches_from_polarization(pol, axis=2):
    """
    Extracts a list of polarization switching branches from the scatter data.
    """
    branches = []

    starts = np.sort(pol[pol["Image"] == 0].iloc[:, axis+1].values)
    start = starts[0]
    while (start <= starts).any():
        branch = branch_from_polarization(pol, axis=axis, start=start)
        if branch is None:
            break
                
        branches.append(branch)
            
        for s in starts:
            if s > start:
                start = s
                break
    
    return branches

def midpoint_branch_from_branches(branches):
    """
    Picks the branch that centers around zero the most.
    """
    midpoints = [0.5*(branch[len(branch)-1] + branch[0]) for branch in branches]
    ith = np.argmin(np.abs(midpoints))
    return branches[ith]

def midpoint_branch_from_polarization(pol, axis=2):
    """
    Picks the branch that centers around zero the most from polarization.
    """
    branches = branches_from_polarization(pol, axis)
    return midpoint_branch_from_branches(branches)

class PolarizationPlotter:
    def __init__(self, path=None):
        self.path = path
        self.data = polarization_from_path(path)
        self.branches = branches_from_polarization(self.data)
        self.switch = midpoint_branch_from_branches(self.branches)
    
    def get_spontaneous(self):
        return 0.5*(self.switch[-1] - self.switch[0])
    
    def plot(self, show_switch=True, save=True):
        plt.style.use(os.path.join(os.path.dirname(__file__), 'dftutils.mplstyle'))

        fig, ax = plt.subplots()
        ax.set_xlabel('Image')
        ax.set_ylabel(r'$P_s\ (\mu C/cm^2)$')
        ax.scatter(self.data['Image'], self.data.iloc[:, 3], s=0.5, color='black')

        if show_switch:
            ax.plot(self.switch, 'o-')

        if save:
            fig.savefig(os.path.join(self.path, 'plot.png'), bbox_inches='tight', pad_inches=0.05)

        return fig, ax


