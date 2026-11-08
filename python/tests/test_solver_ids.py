import unittest

import homlib


class SolverIdsTestCase(unittest.TestCase):
    """homlib.solver_ids: the solvers, all of them or those of a feature type."""

    FEATURE_TYPES = (homlib.FeatureType.Points, homlib.FeatureType.Affine, homlib.FeatureType.Orientation)

    def test_all_solvers(self):
        self.assertEqual(homlib.solver_ids(), list(homlib.SolverId.__members__.values()))

    def test_feature_types(self):
        for feature_type in self.FEATURE_TYPES:
            with self.subTest(feature_type=feature_type.name):
                ids = homlib.solver_ids(feature_type)
                self.assertTrue(ids)
                for solver_id in ids:
                    self.assertEqual(homlib.Solver(solver_id).feature_type, feature_type)
        self.assertIn(homlib.SolverId.DLT, homlib.solver_ids(homlib.FeatureType.Points))
        self.assertIn(homlib.SolverId.BarathVISAPP2016Affine, homlib.solver_ids(homlib.FeatureType.Affine))
        self.assertIn(homlib.SolverId.ValtonenOrnhagICPR2026OrientationOneSided,
                      homlib.solver_ids(homlib.FeatureType.Orientation))

    def test_partition(self):
        """Each solver has exactly one of the feature types, and the order is the enum's."""
        by_type = [s for t in self.FEATURE_TYPES for s in homlib.solver_ids(t)]
        self.assertEqual(sorted(by_type, key=int), sorted(homlib.solver_ids(), key=int))
        for feature_type in self.FEATURE_TYPES:
            ids = homlib.solver_ids(feature_type)
            self.assertEqual(ids, [s for s in homlib.solver_ids() if s in ids])


if __name__ == "__main__":
    unittest.main()
