import unittest
import numpy as np
from .vix_gaussian_stability import balanced_pilot_split,fit_variants
from esn_eq.gaussian_controls import fit_replication_controls


class GaussianPilotStabilityTests(unittest.TestCase):
    def test_balanced_batches_are_disjoint_and_exhaust_training(self):
        a,b=balanced_pilot_split([12,16])
        self.assertEqual(len(a),14);self.assertEqual(len(b),14)
        self.assertEqual(len(set(a)&set(b)),0)
        np.testing.assert_array_equal(np.sort(np.r_[a,b]),np.arange(28))
        self.assertEqual(int((a<12).sum()),6);self.assertEqual(int((b<12).sum()),6)
        with self.assertRaises(ValueError):balanced_pilot_split([12,15])

    def test_legacy_fit_frozen_and_balanced_fits_use_only_their_batch(self):
        rng=np.random.default_rng(29);C=rng.normal(size=(28,9));Y=rng.normal(size=(28,3));legacy=rng.normal(size=(9,3))
        fits=fit_variants(Y,C,legacy,1.);a,b=balanced_pilot_split([12,16])
        np.testing.assert_array_equal(fits['legacy12'],legacy)
        np.testing.assert_allclose(fits['gaussian'],fit_replication_controls(Y,C))
        altered=Y.copy();altered[b]+=100*rng.normal(size=Y[b].shape)
        changed=fit_variants(altered,C,legacy,1.)
        np.testing.assert_array_equal(changed['gaussian_batch_a'],fits['gaussian_batch_a'])
        self.assertGreater(np.max(abs(changed['gaussian_batch_b']-fits['gaussian_batch_b'])),1.)
        with self.assertRaises(ValueError):fit_variants(Y[:27],C[:27],legacy,1.)
