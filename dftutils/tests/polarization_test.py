from click.testing import CliRunner
from dftutils.cli import polarization
import unittest

class PolarizationTest(unittest.TestCase):
    def test_polarization(self):
        runner = CliRunner()
        result = runner.invoke(polarization, ["-o", "dftutils/tests/data/polarization-scatter/00/OUTCAR", "-p", "dftutils/tests/data/polarization-scatter/00/CONTCAR"])
        print(result.output)
        self.assertEqual(result.exit_code, 0)

if __name__ == "__main__":
    unittest.main()