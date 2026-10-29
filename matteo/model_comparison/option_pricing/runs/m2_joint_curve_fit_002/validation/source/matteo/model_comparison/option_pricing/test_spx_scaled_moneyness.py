"""Maturity scaling and held-out isolation for central-band coverage."""
import unittest
import numpy as np
import pandas as pd
from matteo.preprocessing.options import black76_forward
from .spx_scaled_moneyness import scaled_training_quotes


class ScaledMoneynessTests(unittest.TestCase):
    def quotes(self):
        rows=[]
        for expiry,T in [(1,.04),(2,.16)]:
            for K in [100.1,105.]:
                price=black76_forward(100.,K,T,.2,'call')
                rows.append(dict(instrument='SPX',split='train',expiry=expiry,expiry_days=int(T*365),
                    T_years=T,forward=100.,strike=K,option_type='call',source_row_ids=f'{expiry}:{K}',
                    price_bid_U=price,price_mid_U=price,price_ask_U=price))
        return pd.DataFrame(rows)

    def test_maturity_scaling_and_price_unit_invariance(self):
        quotes=self.quotes();scaled,anchors=scaled_training_quotes(quotes)
        np.testing.assert_allclose(anchors.iv_mid,.2,atol=1e-11)
        z=scaled[scaled.strike==105.].sort_values('T_years').z_mid.to_numpy()
        self.assertAlmostEqual(z[0]/z[1],2.)
        self.assertGreater(z[0],1.);self.assertLess(z[1],1.)
        enlarged=quotes.copy()
        enlarged[['forward','strike','price_bid_U','price_mid_U','price_ask_U']]*=10
        other,_=scaled_training_quotes(enlarged)
        np.testing.assert_allclose(scaled.z_mid,other.z_mid,atol=1e-11)

    def test_held_out_quotes_cannot_change_anchor_or_membership(self):
        quotes=self.quotes();expected,anchors=scaled_training_quotes(quotes)
        held=quotes.iloc[[0]].copy();held['split']='final_test';held['strike']=100.
        held[['price_bid_U','price_mid_U','price_ask_U']]=np.nan
        actual,actual_anchors=scaled_training_quotes(pd.concat([quotes,held],ignore_index=True))
        pd.testing.assert_frame_equal(expected.reset_index(drop=True),actual.reset_index(drop=True))
        pd.testing.assert_frame_equal(anchors,actual_anchors)
