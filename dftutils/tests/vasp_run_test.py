import unittest

from click.testing import CliRunner

from dftutils.core.vasp_run import VASPRun
from dftutils.core.vasp_run import Sentinel

class RelaxSentinel(Sentinel):
    def __init__(self):
        super().__init__()

class VASPRunTest(unittest.TestCase):
    def test_vasprun(self):
        sentinel = RelaxSentinel()

        run = (
            sentinel=sentinel,
            dir="dftutils/tests/data/vasp_run/test1",
            incar="dftutils/tests/data/vasp_run/INCAR", 
            kpoints="dftutils/tests/data/vasp_run/KPOINTS", 
            poscar="dftutils/tests/data/vasp_run/POSCAR", 
            potcar="dftutils/tests/data/vasp_run/POTCAR", 
            can_restart = True,
            partition="queue24",
            max_nodes=18,
            cores_per_node=24,
        )

        if run.prepare():
            run.launch()
            run.converge()

if __name__ == "__main__":
    unittest.main()