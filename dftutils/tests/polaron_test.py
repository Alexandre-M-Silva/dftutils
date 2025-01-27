import unittest

from click.testing import CliRunner

from dftutils.core.polaron import Polaron

class PolaronTest(unittest.TestCase):
    def test_polaron(self):
        runner = CliRunner()
        
        p = Polaron.from_structure_path_and_site("dftutils/tests/data/polaron/222.vasp", 28)
        p.to_yaml("dftutils/tests/data/polaron/222.yaml")

        p = Polaron.from_yaml("dftutils/tests/data/polaron/222.yaml")
        p.apply_bdm(0.04)
        p.write_structure("dftutils/tests/data/polaron/222_bdm.vasp")
        p.to_yaml("dftutils/tests/data/polaron/222_bdm.yaml")
        
if __name__ == "__main__":
    unittest.main()