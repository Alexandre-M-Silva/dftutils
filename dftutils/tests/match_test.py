import unittest

from click.testing import CliRunner

from dftutils.cli.match import match

class MatchTest(unittest.TestCase):
    def test_match(self):
        runner = CliRunner()
        result = runner.invoke(match, ["-a", "dftutils/tests/data/match/POSCAR_I", "-b", "dftutils/tests/data/match/POSCAR_F", "-t", "all", "-s"])
        print(result.output)
        print(result.exception)
        self.assertEqual(result.exit_code, 0)

if __name__ == "__main__":
    unittest.main()