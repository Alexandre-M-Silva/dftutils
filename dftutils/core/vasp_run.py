import os
import time
import shutil
import subprocess
import re
import numpy as np
import pickle as pkl
import pandas as pd
from enum import Enum
import datetime
import sys
from typing import List, Dict, IO, Any

from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure

from dftutils.core.utils import *
from dftutils.core.launching import *

class RunStatus(Enum):
    NOT_STARTED = 0
    RUNNING = 1
    UNCONVERGED = 2
    CONVERGED = 3
    UNKNOWN = 4
    
class Sentinel:
    def __init__(self):
        self.run_id = 0

    def increment_run_id(self):
        self.run_id += 1

    def get_run_id(self) -> int:
        return self.run_id
    
    def is_running(self) -> bool:
        return True
    
    def begin(self):
        pass

    def end():
        pass

    def write_data(self):
        pass

    def run(self):
        pass

class VASPRun:
    def __init__(self,
                 dir: str | Path, 
                 incar: str | Incar, 
                 kpoints: str | Kpoints, 
                 poscar: str | Structure,
                 potcar: str | Potcar,
                 can_restart: bool,
                 sentinel: Sentinel | None = None,
                 launch_script: str | None = None,
                 partition: str | None = None,
                 max_nodes: int | None = None,
                 cores_per_node: int | None = None,
                 force_kpar: int | None = None,
                 force_ncore: int | None = None,
                 vasp_setup_cmd: str = VASP_SETUP_CMD,
                 vasp_run_cmd: str = VASP_RUN_CMD):
        self.sentinel = sentinel
        self.id = self.sentinel.get_run_id() if sentinel is not None else 0
        self.dir = dir if isinstance(dir, Path) else Path(dir)
        self.incar = incar if isinstance(incar, Incar) else Incar.from_file(incar)
        self.kpoints = incar if isinstance(kpoints, Kpoints) else Kpoints.from_file(kpoints)
        self.poscar = poscar if isinstance(poscar, Structure) else Structure.from_file(poscar)
        self.potcar = potcar if isinstance(potcar, Potcar) else Potcar.from_file(potcar)
        self.can_restart = can_restart
        self.launch_script = launch_script
        self.partition = partition
        self.max_nodes = max_nodes
        self.cores_per_node = cores_per_node
        self.force_kpar = force_kpar
        self.force_ncore = force_ncore
        self.vasp_setup_cmd = vasp_setup_cmd
        self.vasp_run_cmd = vasp_run_cmd

        self.poll_interval = 30
        self.vr = None
        self.status = RunStatus.NOT_STARTED

        self.incars = []
        self.poscars = []

    def start_print(self):
        pass

    def wait_print(self):
        pass

    def end_print(self):
        pass

    def relaunch_print(self):
        pass

    def prepare(self, extra_incar_settings: Dict[str, Any] | None = None) -> bool:
        if self.get_status() == RunStatus.CONVERGED:
            return False
        
        if extra_incar_settings is not None:
            self.incar.update(extra_incar_settings)

        os.makedirs(self.dir, exist_ok=True)
        return True
    
    def write_input_files(self):
        if self.launch_script is not None:
            self.incar, self.kpoints, self.poscar, self.potcar = setup_vasp_run(incar=self.incar,
                                                                                kpoints=self.kpoints,
                                                                                poscar=self.poscar,
                                                                                potcar=self.potcar,
                                                                                output_dir=self.dir,
                                                                                launch_script=self.launch_script)
        else:
            self.incar, self.kpoints, self.poscar, self.potcar = setup_vasp_run_auto(incar=self.incar,
                                                                                     kpoints=self.kpoints,
                                                                                     poscar=self.poscar,
                                                                                     potcar=self.potcar,
                                                                                     output_dir=self.dir,
                                                                                     partition=self.partition,
                                                                                     max_nodes=self.max_nodes,
                                                                                     cores_per_node=self.cores_per_node,
                                                                                     force_kpar=self.force_kpar,
                                                                                     force_ncore=self.force_ncore,
                                                                                     vasp_setup_cmd=self.vasp_setup_cmd,
                                                                                     vasp_run_cmd=self.vasp_run_cmd,
                                                                                     verbose=True)
            
        self.incars.append(self.incar)
        self.poscars.append(self.poscar)

    def launch(self, print: bool = True):
        if print:
            self.start_print()

        self.write_input_files()
        subprocess.run(f"pushd {self.dir} >/dev/null; sbatch launch.sh; popd >/dev/null",
                       shell=True, capture_output=True, text=True)

    def get_status(self) -> RunStatus:
        path = self.dir / "vasprun.xml"
        if not os.path.exists(path):
            return RunStatus.NOT_STARTED
        try:
            self.vr = Vasprun(path, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
            return RunStatus.UNKNOWN
        except Exception:
            return RunStatus.UNKNOWN

    def move_to_backups(self, files: str | list):
        files = files if isinstance(files, list) else [files]
        for f in files:
            path = self.dir / f
            if os.path.exists(path):
                shutil.move(path, f"{path}.old{self.id}")

    def copy_to_backups(self, files: str | list):
        files = files if isinstance(files, list) else [files]
        for f in files:
            path = self.dir / f
            if os.path.exists(path):
                shutil.copy(path, f"{path}.old{self.id}")

    def remove_file(self, files: str | list):
        files = files if isinstance(files, list) else [files]
        for f in files:
            path = self.dir / f
            if os.path.exists(path):
                os.remove(path)

    def relaunch(self):
        contcar_path = os.path.join(self.dir, "CONTCAR")
        if os.path.exists(contcar_path):
            self.poscar = Structure.from_file(contcar_path)
            self.move_to_backups(["vasprun.xml"])
            self.relaunch_print()
            self.launch()
        else:
            log(f"[ERROR] Cannot relaunch {self.dir}, missing CONTCAR.", flush=True)
            exit(-1)

    def converge(self):
        while True:
            self.status = self.get_status()
            if self.status == RunStatus.CONVERGED:
                break
            elif self.status == RunStatus.UNCONVERGED:
                if self.can_restart:
                    self.relaunch()
                    time.sleep(self.poll_interval)
                else:
                    log(f"[ERROR] Cannot relaunch {self.dir}, run is not restartable.", flush=True)
                    exit(-1)
                    break
            else:
                self.wait_print()
                time.sleep(self.poll_interval)

        self.end_print()

    @classmethod
    def from_vasprun_xml(path: str | Path, sentinel: Sentinel | None = None, structure_to_extract: int = -1) -> "VASPRun":
        vr = Vasprun(path, parse_dos=False, parse_eigen=False, parse_projected_eigen=False)
        result = VASPRun(os.path.dirname(path),
                         vr.incar,
                         vr.kpoints,
                         vr.structures[structure_to_extract],
                         vr.get_potcars(True),
                         sentinel)
        return result

class VASPSpeRun(VASPRun):
    def __init__(self,
                 prev_run: str | VASPRun = None,
                 dir: str | Path = None, 
                 incar: str | Incar = None, 
                 kpoints: str | Kpoints = None, 
                 poscar: str | Structure = None,
                 potcar: str | Potcar = None,
                 sentinel: Sentinel = None,
                 launch_script: str | None = None,
                 partition: str | None = None,
                 max_nodes: int | None = None,
                 cores_per_node: int | None = None,
                 force_kpar: int | None = None,
                 force_ncore: int | None = None,
                 vasp_setup_cmd: str = VASP_SETUP_CMD,
                 vasp_run_cmd: str = VASP_RUN_CMD):
        
        self.prev_run = None
        if prev_run is not None:
            self.prev_run = prev_run if isinstance(prev_run, VASPRun) else VASPRun.from_vasprun_xml(prev_run)

        if prev_run is None:
            super().__init__(dir, incar, kpoints, poscar, potcar, sentinel, True, launch_script, partition,
                         max_nodes, cores_per_node, force_kpar, force_ncore, vasp_setup_cmd, vasp_run_cmd)
        else:
            super().__init__(dir, prev_run.incar, prev_run.kpoints, prev_run.poscar, prev_run.potcar, sentinel, True, launch_script, partition,
                         max_nodes, cores_per_node, force_kpar, force_ncore, vasp_setup_cmd, vasp_run_cmd)

    def start_print(self):
        log(f"{"[SPE]":>10s} Starting SPE in {self.dir}", color="cyan")

    def relaunch_print(self):
        log(f"{"[SPE]":>10s} Did not converge in time, relaunching...", color="red")

    def wait_print(self):
        log(f"{"[SPE]":>10s} Waiting...", color="yellow")

    def end_print(self):
        log(f"{"[SPE]":>10s} Ended SPE in {self.dir}", color="green")
    
    def get_status(self):
        if super().get_status() == RunStatus.NOT_STARTED:
            return RunStatus.NOT_STARTED
        
        if self.vr.converged_electronic and \
           self.incar.get("NSW", 0) == 0 and \
           self.incar.get("LEPSILON", 0) == 0 and \
           self.incar.get("LCALCEPS", 0) == 0 and \
           self.incar.get("LCALCPOL", 0) == 0:
            return RunStatus.CONVERGED
            
        return RunStatus.UNKNOWN

    def prepare(self, prev_run: str | VASPRun = None, extra_incar_settings: Dict[str, Any] | None = None):
        if not super().prepare(extra_incar_settings):
            return False

        if self.prev_run is not None:
            self.incar = prev_run.incar
            self.kpoints = prev_run.kpoints
            self.poscar = prev_run.poscar
            self.potcar = prev_run.potcar

        self.incar["ISTART"] = 1
        self.incar["NSW"] = 0
        self.incar["IBRION"] = -1
        self.incar["NELM"] = 200

        return True

    def converge(self):
        super().converge()

class VASPRelaxRun(VASPRun):
    def __init__(self,
                 sentinel: Sentinel,
                 dir: str | Path, 
                 incar: str | Incar, 
                 kpoints: str | Kpoints, 
                 poscar: str | Structure,
                 potcar: str | Potcar,
                 relax_kind: str = "full",
                 launch_script: str | None = None,
                 partition: str | None = None,
                 max_nodes: int | None = None,
                 cores_per_node: int | None = None,
                 force_kpar: int | None = None,
                 force_ncore: int | None = None,
                 vasp_setup_cmd: str = VASP_SETUP_CMD,
                 vasp_run_cmd: str = VASP_RUN_CMD):
        super().__init__(sentinel, dir, incar, kpoints, poscar, potcar, True, launch_script, partition,
                         max_nodes, cores_per_node, force_kpar, force_ncore, vasp_setup_cmd, vasp_run_cmd)
        self.relax_kind = relax_kind
        
    def start_print(self):
        log(f"{"[RELAX]":>10s} Starting relaxation in {self.dir}", color="cyan")

    def relaunch_print(self):
        log(f"{"[RELAX]":>10s} Did not converge in time, relaunching...", color="red")

    def wait_print(self):
        log(f"{"[RELAX]":>10s} Waiting...", color="yellow")

    def end_print(self):
        log(f"{"[RELAX]":>10s} Ended relaxation in {self.dir}", color="green")
    
    def get_status(self):
        if super().get_status() == RunStatus.NOT_STARTED:
            return RunStatus.NOT_STARTED
        
        if self.vr.converged and self.vr.incar.get("NSW", 0) > 0:
            return RunStatus.CONVERGED
            
        nsw = self.incar.get("NSW", 0)
        if nsw > 0: 
            nionic_steps = self.vr.nionic_steps
            if nionic_steps < nsw:
                return RunStatus.RUNNING
            elif nionic_steps >= nsw:
                return RunStatus.UNCONVERGED
            
        return RunStatus.UNKNOWN

    def prepare(self, prev_run):
        if not super().prepare():
            return False

        self.incar['ISIF'] = 2
        self.incar["LWAVE"] = True
        self.incar["ISTART"] = 1
        self.incar["EFOR"] = build_efor_str(self.bec, self.efield, None, self.efield_dir)
        self.write_incar()

        return True

    def converge(self):
        super().converge()
