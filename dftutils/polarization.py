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

def polarization_scatter_from_path(path, bmin=-5, bmax=5):
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

def branch_from_polarization_scatter(pol, axis=2, start=0):
    """
    Extracts the polarization swithcing branch from the scatter data,
    assuming that its monotonically increasing.
    """
    branch = []
    nimages = pol['Image'].nunique()
    image = 0
    P = start if start != 0 else np.min(pol[pol['Image'] == 0])
    for i in range(0, nimages):
        df = pol[pol['Image']==i]
        df.iloc[:, axis+1] = df.iloc[:, axis+1] - P
        df = df[df.iloc[:, axis+1] >= 0]
        if df.iloc[:, axis+1].shape[0] > 1:
            df = df.iloc[df.iloc[:, axis+1].argsort()[:-1]]
        
        P = pol.loc[df.iloc[:, axis+1].index[0]].values[axis+1]
        branch.append(P)    
    return np.array(branch)

def branches_from_polarization_scatter(pol, axis=2):
    """
    Extracts a list of polarization switching branches from the scatter data.
    """
    branches = []

    starts = np.sort(pol[pol['Image'] == 0].values)
    print(starts)
    start = starts[0]
    while (start <= starts).any():
        print(start)
        branches.append(branch_from_polarization_scatter(pol, axis=axis, start=start))
        for s in starts:
            if s > start:
                start = s
    
    return branches

def find_branch_around_zero(branches):
    """
    Picks the branch that centers around zero the most.
    """
    midpoints = [0.5*(branch[len(branch)-1] + branch[0]) for branch in branches]
    ith = np.argmin(np.abs(midpoints))
    return branches[ith]
