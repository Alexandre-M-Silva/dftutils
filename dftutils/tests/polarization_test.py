import unittest

from click.testing import CliRunner

from dftutils.cli.polarization import polarization

class PolarizationTest(unittest.TestCase):
    def test_polarization(self):
        runner = CliRunner()
        result = runner.invoke(polarization, ["-o", "dftutils/tests/data/polarization-scatter/pol_nc/00/OUTCAR", "-c", "dftutils/tests/data/polarization-scatter/pol_nc/00/CONTCAR"])
        print(result.output)
        print(result.exception)
        self.assertEqual(result.exit_code, 0)

        result = runner.invoke(polarization, ["-s", "-p", "dftutils/tests/data/polarization-scatter/pol_nc", "-ed", "-ef"])
        print(result.output)
        print(result.exception)
        self.assertEqual(result.exit_code, 0)
        
        result = runner.invoke(polarization, ["-s", "-p", "dftutils/tests/data/polarization-scatter/pol_c", "-ed", "-ef"])
        print(result.output)
        print(result.exception)
        self.assertEqual(result.exit_code, 0)
        
if __name__ == "__main__":
    unittest.main()