"""Market-only side/band panel selection and incomplete-grid rejection."""
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from .spx_central_panel import choose_panel,MEMBERSHIP,record_trial
from .spx_innovation_study import train


class CentralPanelTests(unittest.TestCase):
    def test_market_panel_coverage_and_selection_ignores_prices(self):
        quotes=pd.read_csv(MEMBERSHIP);panel=choose_panel(quotes)
        self.assertEqual(len(panel),24)
        self.assertEqual(set(panel.expiry_days),{29,57,211,302})
        self.assertEqual(panel.groupby(['expiry_days','option_type','band']).size().tolist(),[1]*24)
        changed=quotes.copy();changed[['price_mid_U','price_bid_U','price_ask_U']]=np.nan
        other=choose_panel(changed)
        pd.testing.assert_frame_equal(panel[['expiry','strike','option_type','band']],other[['expiry','strike','option_type','band']])
        self.assertTrue((panel.z_mid<=3).all())
        with self.assertRaises(ValueError):choose_panel(quotes[~((quotes.expiry_days==29)&(quotes.option_type=='put')&(quotes.z_mid>=2)&(quotes.z_mid<=3))])

    def test_price_bookkeeping_and_negative_estimates_preserved(self):
        quotes=choose_panel(pd.read_csv(MEMBERSHIP)).iloc[:6].copy()
        rng=np.random.default_rng(123);y=np.abs(rng.normal(size=(32,6)));b=rng.normal(size=(32,6,5));v=rng.normal(size=(32,6,3))
        fits=train(y,b,v);fits['augmented']=np.full((6,8),100.)
        data=dict(g3_y=y,g3_base=b,g3_vol=v)
        x,rows,paired=record_trial(data,fits,quotes,0)
        self.assertEqual(x.shape,(32,6,8));self.assertEqual(len(rows),48);self.assertEqual(len(paired),6)
        self.assertTrue(any(row['estimate']<0 for row in rows))
        np.testing.assert_allclose(x[:,:,0],y*quotes.forward.to_numpy()[None,:])
