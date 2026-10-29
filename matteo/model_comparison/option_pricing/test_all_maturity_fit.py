"""Coverage, split promotion and objective correspondence for feasibility fits."""
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

from .all_maturity_core import populations,nearest_same_side,InSampleLoss,SPX_DAYS,EXCLUDED_DAYS
from esn_eq.calibration_loss import JointPriceObjective


class AllMaturityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=pd.read_csv(Path(__file__).resolve().parents[3]/'data/processed/options/2026-02-19/quotes_calibration.csv')
        cls.f=pd.read_csv(Path(__file__).resolve().parents[3]/'data/processed/options/2026-02-19/vix_forwards.csv')

    def test_all_seven_dates_and_all_original_splits_are_calibration_inputs(self):
        original=self.raw.copy(deep=True);spx,vix,m,a=populations(self.raw)
        pd.testing.assert_frame_equal(self.raw,original)
        self.assertEqual(set(spx.expiry_days),set(SPX_DAYS));self.assertFalse(set(m.expiry_days)&set(EXCLUDED_DAYS))
        self.assertEqual(len(spx),1294);self.assertEqual(len(vix),217)
        self.assertEqual(set(spx.original_split),{'train','validation','final_test','diagnostic','final_maturity_test'})
        self.assertTrue(spx.split.eq('in_sample_fit').all());self.assertTrue(vix.split.eq('in_sample_fit').all())
        np.testing.assert_allclose(m.scaled_moneyness,abs(np.log(m.K_over_F))/(m.atm_proxy*np.sqrt(m.T_years)),rtol=1e-14)
        np.testing.assert_array_equal(m.priority_selected,m.scaled_moneyness<=3)
        self.assertEqual(list(a.selected),[277,298,228,141,167,112,71])

    def test_frozen_control_donors_reproduce_exact_existing_keys_and_never_switch_side(self):
        spx,_,_,_=populations(self.raw);donors=spx[spx.expiry_days==29].reset_index(drop=True)
        np.testing.assert_array_equal(nearest_same_side(donors,donors),np.arange(len(donors)))
        target=donors.iloc[[4,9,-1]].copy();target['K_over_F']*=1.0003
        ix=nearest_same_side(target,donors)
        np.testing.assert_array_equal(target.option_type,donors.iloc[ix].option_type)
        with self.assertRaises(ValueError):nearest_same_side(target,donors[donors.option_type=='put'].reset_index(drop=True))

    def test_explicit_in_sample_loss_agrees_with_production_for_nonnegative_prices(self):
        spx,vix,_,_=populations(self.raw);q=pd.concat([spx,vix],ignore_index=True)
        loss=InSampleLoss(q,self.f);production=JointPriceObjective(q,self.f,split='in_sample_fit')
        price=q.price_mid_U.to_numpy()*1.2;forward=self.f.forward_vix_points.to_numpy()+.7;xi=np.array([.05,.06,.04])
        actual=loss.evaluate(price,forward,xi)
        expected=production.evaluate(price,forward,curve_penalty=.1*np.sum(np.diff(xi)**2))
        self.assertAlmostEqual(actual['total'],expected['total'],places=9)
        for instrument in ['SPX','VIX']:
            self.assertAlmostEqual(loss.weights[q.instrument==instrument].sum(),1.,places=14)
        with self.assertRaises(ValueError):JointPriceObjective(q,self.f)


if __name__=='__main__':unittest.main()
