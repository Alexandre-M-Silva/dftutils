import unittest

from click.testing import CliRunner

from dftutils.cli import polarization
from dftutils.cli import polarization_scatter

class PolarizationTest(unittest.TestCase):
    def test_polarization(self):
        runner = CliRunner()
        result = runner.invoke(polarization, ["-o", "dftutils/tests/data/polarization-scatter/pol_nc/00/OUTCAR", "-p", "dftutils/tests/data/polarization-scatter/pol_nc/00/CONTCAR"])
        print(result.output)
        self.assertEqual(result.exit_code, 0)

    def test_polarization_scatter(self):
        runner = CliRunner()
        result = runner.invoke(polarization_scatter, ["-p", "dftutils/tests/data/polarization-scatter/pol_nc", "-r"])
        self.assertEqual(result.exit_code, 0)
        
        result = runner.invoke(polarization_scatter, ["-p", "dftutils/tests/data/polarization-scatter/pol_c", "-r"])
        self.assertEqual(result.exit_code, 0)
        
if __name__ == "__main__":
    unittest.main()