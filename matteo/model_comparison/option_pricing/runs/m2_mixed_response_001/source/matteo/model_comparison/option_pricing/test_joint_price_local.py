import json
from pathlib import Path
import time
import unittest
import numpy as np
import pandas as pd
from .calibration_policy import LOWER, UPPER, MINIMUM_VARIANCE_FLOOR, validate
from .full_readout_fit import SCALE
from .joint_price_local import JointPriceLocal


class JointLocalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = Path(__file__).parent / 'runs/m2_parity_parameter_fit_003'
        s = json.loads((cls.folder/'score_selected.json').read_text())
        q = pd.read_csv(cls.folder/'prices_canonical_selected.csv')
        o = pd.read_csv(cls.folder/'prices_otm_selected.csv')
        f = pd.read_csv(cls.folder/'forwards_selected.csv')
        with np.load(cls.folder/'jacobian_selected.npz') as j:
            cls.anchor = dict(prices=q.model_price_U.to_numpy(), forwards=f.model_forward.to_numpy(),
                              price_jacobian=j['price'], forward_jacobian=j['forward'])
        cls.p = np.array(s['parameters'])
        cls.score = s
        cls.model = JointPriceLocal(cls.p, cls.anchor, q, o, f)
        cls.native_quotes = o

    def test_five_percent_policy_and_all_eleven_free(self):
        p = self.p.copy();p[10] = MINIMUM_VARIANCE_FLOOR
        validate(p)
        self.assertTrue(np.all(UPPER > LOWER))
        self.assertAlmostEqual(np.sqrt(LOWER[10]), .05)
        self.assertAlmostEqual(np.sqrt(UPPER[10]), .06)
        p[10] = .049 ** 2
        with self.assertRaises(ValueError):validate(p)
        p[10] = UPPER[10];validate(p)
        ratio = np.exp(np.linspace(-50, 30, 100))
        for xi in p[7:10]:
            self.assertTrue(np.all(p[10] + (xi-p[10])*ratio >= p[10]))

    def test_saved_full_anchor_replay_both_objectives(self):
        r = self.model.evaluate(self.p)
        np.testing.assert_allclose(r['criteria'][:2], [self.score['canonical_loss'], self.score['original_otm_loss']], rtol=1e-12)
        np.testing.assert_allclose(r['native_prices'], self.native_quotes.model_price_U, rtol=1e-12, atol=1e-12)
        self.assertEqual(r['simulator_calls'], 0)

    def test_all_eleven_gradient_columns_both_conventions(self):
        p = self.p + SCALE*np.array([.002, -.001, .001, .002, -.001, .001, .001, -.001, .001, .002, .001])
        r = self.model.evaluate(p)
        for j in range(11):
            h = 2e-7*SCALE[j];a=p.copy();b=p.copy();a[j]+=h;b[j]-=h
            fd=(self.model.evaluate(a)['criteria']-self.model.evaluate(b)['criteria'])/(2*h)
            np.testing.assert_allclose(fd, r['criteria_gradients'][:, j], rtol=3e-5, atol=2e-4)

    def test_prediction_is_affine_and_native_parity_uses_model_forward(self):
        d = SCALE*np.linspace(-.002,.002,11)
        a=self.model.evaluate(self.p+d);b=self.model.evaluate(self.p-d)
        np.testing.assert_allclose((a['prices']+b['prices'])/2, self.anchor['prices'], atol=1e-12)
        for day,g in self.native_quotes[self.native_quotes.instrument.eq('VIX')].groupby('expiry_days'):
            j=list(self.model.loss.forwards.expiry_days).index(day)
            for i in g[g.option_type.eq('put')].index:
                self.assertAlmostEqual(a['prices'][i]-a['native_prices'][i], a['forwards'][j]-g.loc[i,'strike'])

    def test_local_search_never_calls_simulator_and_obeys_new_bounds(self):
        from .joint_local_search import propose
        r=propose(self.model,radius=.02)
        self.assertEqual(r['full_pricer_calls'],0)
        self.assertIsNotNone(r['parameters'])
        validate(r['parameters'])
        self.assertLess(r['prediction']['criteria'][0], self.score['canonical_loss'])
        self.assertLess(r['prediction']['criteria'][3], self.score['original_otm_components']['VIX']+self.score['original_otm_components']['forwards'])


if __name__ == '__main__':unittest.main()
