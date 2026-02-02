import os
import math
import re
import numpy as np
from typing import List, Dict, IO
import shutil
from pathlib import Path
import pandas as pd

import warnings

from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.inputs import Potcar

from dftutils.core.utils import *

VASP_SETUP_CMD = "module load vasp"
VASP_RUN_CMD = "mpirun --mca pml ucx /opt/ohpc/pub/apps/vasp/6.5/vasp_std"

def recommend_parallelization(
    poscar: str | Structure,
    potcar: str | Potcar,
    kpoints: str | Kpoints,
    max_nodes: int,
    cores_per_node: int,
    allowed_ncore_values: int | List[int] = None,
    min_bands_per_rank: int = 2,
    max_bands_per_rank: int = 32,
    target_bands_per_rank: int = 8,
    margin_frac: float = 0.2,
    top_n: int = 10,
    force_kpar: int | None = None,
    force_nbands: int = None,
    unbalanced_kpar_punishment: bool | str = "auto",
) -> List[Dict]:
    """
    Search for best (nodes, NCORE, KPAR) combos and score them.
    """
    poscar = Structure.from_file(poscar) if isinstance(poscar, str) else poscar
    potcar = Potcar.from_file(potcar) if isinstance(potcar, str) else potcar
    kpoints = Kpoints.from_file(kpoints) if isinstance(kpoints, str) else kpoints

    if isinstance(unbalanced_kpar_punishment, str):
        unbalanced_kpar_punishment = False if len(poscar) < 48 else True

    nelect = nelect_from_structure_potcar(poscar, potcar)
    nkpts = kpoints_count(kpoints)

    if allowed_ncore_values is None:
        allowed_ncore_values = divisors(cores_per_node)
    allowed_ncore_values = [v for v in allowed_ncore_values if v <= cores_per_node] if isinstance(allowed_ncore_values, list) else [allowed_ncore_values]

    candidates = []
    for ncore in allowed_ncore_values:
        for nodes in range(1, max_nodes + 1):
            if force_nbands is not None:
                if cores_per_node % ncore != 0:
                    ranks_per_node = cores_per_node // ncore
                else:
                    ranks_per_node = cores_per_node // ncore

                total_ranks = nodes * ranks_per_node
                nbands = force_nbands
                bands_per_rank = nbands / total_ranks
            else:
                nbands, total_ranks, ranks_per_node, bands_per_rank = estimate_nbands(
                    nelect=nelect,
                    nodes=nodes,
                    cores_per_node=cores_per_node,
                    ncore=ncore,
                    margin_frac=margin_frac,
                )
            nbands_upper_end = estimate_nbands_upper_end(nelect, max_nodes, cores_per_node)
            
            if total_ranks < 1:
                continue

            if force_kpar is not None:
                kp_list = [force_kpar]
            else:
                kp_candidates = {d for d in divisors(total_ranks) if d <= nkpts}
                kp_candidates.update({1, min(nkpts, total_ranks)})
                kp_list = sorted(kp_candidates)

            for kpar in kp_list:
                if total_ranks % kpar != 0:
                    continue
                
                nodes_coef = 1.0 if len(poscar) < 48 else 0.1
                nbands_coef = 0.1 if len(poscar) < 48 else 1.0
                kpar_coef = 1.0 if len(poscar) < 48 else 0.1
                ncore_coef = 0.1 if len(poscar) < 48 else 0.05

                score = {
                    "Score_BandsPerRankTarget": np.abs((bands_per_rank-target_bands_per_rank)/target_bands_per_rank),
                    "Score_BandsPerRankLimits": (1.0 if bands_per_rank < min_bands_per_rank else 0.0) + (1.0 if bands_per_rank < max_bands_per_rank else 0.0),
                    "Score_Nodes": nodes_coef*(float(nodes)/float(max_nodes)),
                    "Score_NumBands": nbands_coef*((float(nbands) / float(nbands_upper_end))),
                    "Score_Kpar":  (1.0/kpar_coef) * (1.0 - float(kpar)/float(nkpts) + (1.0 if (nkpts % kpar != 0 and unbalanced_kpar_punishment) else 0.0)),
                    "Score_Ncore": (1.0/ncore_coef) * (1.0 - float(ncore)/float(cores_per_node)),
                }
                score["Score_Total"] = np.sum(np.array(list(score.values())))
 
                candidate = {
                    "nodes": nodes,
                    "ncore": ncore,
                    "kpar": kpar,
                    "nbands": nbands,
                    "nelect": nelect,
                    "nkpts": nkpts,
                    "ranks_per_node": ranks_per_node,
                    "total_ranks": total_ranks,
                    "bands_per_rank": bands_per_rank,
                    **score
                }
                candidates.append(candidate)                

    candidates = pd.DataFrame(candidates).sort_values(by="Score_Total").reset_index(drop=True)
    candidates["nodes"] = candidates["nodes"].astype('int')
    candidates["ncore"] = candidates["ncore"].astype('int')
    candidates["kpar"] = candidates["kpar"].astype('int')
    candidates["nbands"] = candidates["nbands"].astype('int')
    candidates["nelect"] = candidates["nelect"].astype('int')
    candidates["nkpts"] = candidates["nkpts"].astype('int')
    candidates["ranks_per_node"] = candidates["ranks_per_node"].astype('int')
    candidates["total_ranks"] = candidates["total_ranks"].astype('int')
    if top_n > 0:
        return candidates.head(top_n)

    return candidates

def launch_script_write(path: str, 
                        partition: str, 
                        nodes: int, 
                        cores_per_node: int,
                        vasp_setup_cmd: str = VASP_SETUP_CMD,    
                        vasp_run_cmd: str = VASP_RUN_CMD):
    
    launch_script = f"""#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --nodes={int(nodes)}
#SBATCH --ntasks-per-node={int(cores_per_node)}
#SBATCH --exclusive
#SBATCH --mem=0

{vasp_setup_cmd}
{vasp_run_cmd}
"""

    with open(path, "w") as f:
        f.write(launch_script)

def setup_vasp_run(
    incar: str | Incar,
    kpoints: str | Kpoints,
    poscar: str | Structure,
    potcar: str | Potcar,
    output_dir: str | Path,
    nbands: int | None = None,
    ncore: int | None = None,
    kpar: int | None = None,
    launch_script: str | None = None,
    partition: str | None = None,
    nodes: int | None = None,
    cores_per_node: int | None = None,
    vasp_setup_cmd: str = VASP_SETUP_CMD,
    vasp_run_cmd: str = VASP_RUN_CMD,
    launch_script_name: str = "launch.sh",
):
    incar = incar if isinstance(incar, Incar) else Incar.from_file(incar)
    kpoints = incar if isinstance(kpoints, Kpoints) else Kpoints.from_file(kpoints)
    poscar = poscar if isinstance(poscar, Structure) else Structure.from_file(poscar)
    potcar = potcar if isinstance(potcar, Potcar) else Potcar.from_file(potcar)

    output_dir = Path(output_dir) if isinstance(output_dir, str) else output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if ncore is not None:
        incar["NCORE"] = int(ncore)
    if nbands is not None:
        incar["NBANDS"] = int(nbands)
    if kpar is not None:
        incar["KPAR"] = int(kpar)

    incar_out = output_dir / "INCAR"
    incar.write_file(incar_out)

    kpoints_out = output_dir / "KPOINTS"
    kpoints.write_file(kpoints_out)

    poscar_out = output_dir / "POSCAR"
    poscar.to(poscar_out, fmt="poscar")

    potcar_out = output_dir / "POTCAR"
    potcar.write_file(potcar_out)

    if launch_script is None:
        if partition is None or nodes is None or cores_per_node is None:
            raise Exception("Partition, nodes and cores_per_node must be set if launch_script is None!")
        launch_out = output_dir / launch_script_name
        launch_script_write(launch_out, partition, nodes, cores_per_node, vasp_setup_cmd, vasp_run_cmd)
    else:
        launch_out = output_dir / launch_script_name
        shutil.copy(launch_script, launch_out)

    return (incar, kpoints, poscar, potcar)

def setup_vasp_run_auto(
    incar: str | Incar,
    kpoints: str | Kpoints,
    poscar: str | Structure,
    potcar: str | Potcar,
    output_dir: str | Path,
    partition: str,
    max_nodes: int,
    cores_per_node: int,
    force_kpar: int | None = None,
    force_ncore: int | None = None,
    vasp_setup_cmd: str = VASP_SETUP_CMD,
    vasp_run_cmd: str = VASP_RUN_CMD,
    margin_frac: float = 0.2,
    launch_script_name: str = "launch.sh",
    verbose: bool = False,
):
    recs = recommend_parallelization(
        poscar=poscar,
        kpoints=kpoints,
        potcar=potcar,
        max_nodes=max_nodes,
        cores_per_node=cores_per_node,
        top_n=1,
        margin_frac=margin_frac,
        force_kpar=force_kpar,
        allowed_ncore_values=force_ncore,
    )
    best = recs.loc[0]

    if verbose:
        print(recs.head(0))

    return setup_vasp_run(
        incar=incar,
        poscar=poscar,
        kpoints=kpoints,
        potcar=potcar,
        partition=partition,
        cores_per_node=cores_per_node,
        nodes=int(best["nodes"]),
        nbands=int(best["nbands"]),
        ncore=int(best["ncore"]),
        kpar=int(best["kpar"]),
        vasp_setup_cmd=vasp_setup_cmd,
        vasp_run_cmd=vasp_run_cmd,
        launch_script_name=launch_script_name,
        output_dir=output_dir,
    )