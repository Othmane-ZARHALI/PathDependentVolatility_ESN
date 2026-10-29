"""Restore a saved full true-price/Jacobian anchor without resimulating it."""
import json
import numpy as np
import pandas as pd


def load_result(folder,name):
    s=json.loads((folder/f'score_{name}.json').read_text());q=pd.read_csv(folder/f'prices_canonical_{name}.csv');o=pd.read_csv(folder/f'prices_otm_{name}.csv');f=pd.read_csv(folder/f'forwards_{name}.csv')
    j=np.load(folder/f'jacobian_{name}.npz');a=np.load(folder/f'paths_{name}.npz')
    return dict(total=s['canonical_loss'],components=s['canonical_components'],original_otm_total=s['original_otm_loss'],original_otm_components=s['original_otm_components'],
                prices=q.model_price_U.to_numpy(),errors=q.mc_se_U.to_numpy(),forwards=f.model_forward.to_numpy(),forward_errors=f.mc_se.to_numpy(),
                native_prices=o.model_price_U.to_numpy(),native_errors=o.mc_se_U.to_numpy(),actual_paths=q.actual_paths.to_numpy(),effective_paths=q.effective_paths.to_numpy(),
                price_jacobian=j['price'],forward_jacobian=j['forward'],jacobian=j['residual'],residuals=j['residuals'],raw={k:a[k] for k in a.files},
                normalizer=s['normalizer'],allocations=s['allocations'],moments=s['moments'],derivatives_computed=True)
