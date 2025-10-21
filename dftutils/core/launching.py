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

from dftutils.core.utils import *

HIGGS_VASP_SETUP_CMD = "module load vasp"
HIGGS_VASP_RUN_CMD = "mpirun --mca pml ucx /opt/ohpc/pub/apps/vasp/6.5/vasp_std"

def recommend_parallelization(
    poscar: str | Structure,
    potcar: str,
    kpoints: str | Kpoints,
    partition_info: Dict[str, int],
    allowed_ncore_values: List[int] = None,
    min_bands_per_rank: int = 4,
    max_bands_per_rank: int = 64,
    target_bands_per_rank: int = 8,
    margin_frac: float = 0.2,
    top_n: int = 8,
    force_kpar: int | None = None,
    force_nbands: int = None,
    unbalanced_kpar_punishment: bool | str = "auto",
) -> List[Dict]:
    """
    Search for best (nodes, NCORE, KPAR) combos and score them.
    """
    poscar = Structure.from_file(poscar) if isinstance(poscar, str) else poscar
    potcar = open(potcar, "r")
    kpoints = Kpoints.from_file(kpoints) if isinstance(kpoints, str) else kpoints

    max_nodes = partition_info["Nodes"]
    cores_per_node = partition_info["CoresPerNode"]

    if isinstance(unbalanced_kpar_punishment, str):
        unbalanced_kpar_punishment = False if len(poscar) < 48 else True

    nelect = nelect_from_structure_potcar(poscar, potcar)
    nkpts = kpoints_count(kpoints)

    if allowed_ncore_values is None:
        allowed_ncore_values = divisors(cores_per_node)
    allowed_ncore_values = [v for v in allowed_ncore_values if v <= cores_per_node]

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

def setup_vasp_run(
    incar_path: str,
    poscar_path: str,
    potcar_path: str,
    kpoints_path: str,
    partition: str,
    nodes: int,
    cores_per_node: int,
    nbands: int,
    ncore: int,
    kpar: int,
    vasp_setup_cmd: str = HIGGS_VASP_SETUP_CMD,
    vasp_run_cmd: str = HIGGS_VASP_RUN_CMD,
    launch_script_name: str = "launch.sh",
    output_dir: str | Path | None = None,
):
    output_dir = Path(output_dir) if output_dir else Path(poscar_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    incar = Incar.from_file(incar_path)
    poscar = Structure.from_file(poscar_path)
    kpoints = Kpoints.from_file(kpoints_path)

    incar["NBANDS"] = nbands
    incar["NCORE"] = ncore
    incar["KPAR"] = kpar
    incar_out = output_dir / "INCAR"
    incar.write_file(incar_out)

    kpoints_out = output_dir / "KPOINTS"
    kpoints.write_file(kpoints_out)

    poscar_out = output_dir / "POSCAR"
    poscar.to(poscar_out, fmt="poscar")

    potcar_out = output_dir / "POTCAR"
    shutil.copy(potcar_path, potcar_out)

    launch_script = f"""#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={cores_per_node}
#SBATCH --exclusive
#SBATCH --mem=0

{vasp_setup_cmd}
{vasp_run_cmd}
"""
    launch_out = output_dir / launch_script_name
    with open(launch_out, "w") as f:
        f.write(launch_script)

def setup_vasp_run_auto(
    incar_path: str,
    poscar_path: str,
    potcar_path: str,
    kpoints_path: str,
    partition_info: Dict[str, int], 
    force_kpar: int | None = None,
    vasp_setup_cmd: str = HIGGS_VASP_SETUP_CMD,
    vasp_run_cmd: str = HIGGS_VASP_RUN_CMD,
    margin_frac: float = 0.2,
    launch_script_name: str = "launch.sh",
    output_dir: str | Path | None = None,
    verbose: bool | None = None,
):
    recs = recommend_parallelization(
        poscar=poscar_path,
        potcar=potcar_path,
        kpoints=kpoints_path,
        max_nodes=partition_info["Nodes"],
        cores_per_node=partition_info["CoresPerNode"],
        top_n=1,
        margin_frac=margin_frac,
        force_kpar=force_kpar,
    )
    best = recs.loc[0]

    if verbose:
        print(recs.head(0))

    setup_vasp_run(
        incar_path=incar_path,
        poscar_path=poscar_path,
        potcar_path=potcar_path,
        kpoints_path=kpoints_path,
        partition=partition_info["Name"],
        cores_per_node=partition_info["CoresPerNode"],
        nodes=best["nodes"],
        nbands=best["nbands"],
        ncore=best["ncore"],
        kpar=best["kpar"],
        vasp_setup_cmd=vasp_setup_cmd,
        vasp_run_cmd=vasp_run_cmd,
        launch_script_name=launch_script_name,
        output_dir=output_dir,
    )

