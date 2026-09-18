"""Numerical and physical checks independent of the supplied coursework outputs."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
from analysis import cross_section, lifetime, load_data, fit, parameters, profile_grid


class AnalysisChecks(unittest.TestCase):
    def test_pole_normalization(self):
        self.assertAlmostEqual(cross_section(91.2, 91.2, 2.5, 2.0), 2.0)

    def test_lifetime_units(self):
        self.assertAlmostEqual(lifetime(2.5) / 2.6328478278e-25, 1.0, places=8)
        for width in [0, -1, float("nan")]:
            with self.assertRaises(ValueError):
                lifetime(width)

    def test_noiseless_recovery_and_profile(self):
        energy = np.linspace(85.2, 94.8, 120)
        truth = np.array([91.1, 2.4, 1.9])
        data = np.column_stack([energy, cross_section(energy, *truth), np.full(120, .08)])
        result = fit(data)
        np.testing.assert_allclose(parameters(result), truth, rtol=1e-7)
        _, _, grid = profile_grid(data, [truth[0]], [truth[1]])
        self.assertLess(grid[0, 0], 1e-20)

    def test_invalid_rows_audited(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "data.csv"
            path.write_text("E,sigma,error\n91,2,.1\nfail,2,.1\n92,1,0\n93,nan,.1\n85,1,.1\n91,200,.1\n", encoding="utf-8")
            data, audit = load_data([path])
            self.assertEqual(data.shape, (1, 3))
            self.assertEqual(len(audit), 6)
            self.assertEqual(sum(r["status"] == "candidate" for r in audit), 1)

    def test_insufficient_data(self):
        with self.assertRaises(ValueError):
            fit(np.ones((3, 3)))


if __name__ == "__main__":
    unittest.main()
