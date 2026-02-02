#!/bin/bash
#SBATCH --partition=queue24
#SBATCH --nodes=16.0
#SBATCH --ntasks-per-node=24
#SBATCH --exclusive
#SBATCH --mem=0

module load vasp
mpirun --mca pml ucx /opt/ohpc/pub/apps/vasp/6.5/vasp_std
