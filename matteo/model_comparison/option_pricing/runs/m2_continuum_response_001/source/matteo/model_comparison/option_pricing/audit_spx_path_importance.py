"""Descriptive SPX importance-pilot summary; no convergence certification."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def audit(output):
    e = pd.read_csv(output/'estimates.csv')
    c = pd.read_csv(output/'paired_changes.csv')
    d = pd.read_csv(output/'diagnostics.csv')
    fine = e[(e.steps == 256) & (e.normalizer_order == 48)]
    result = dict(description='Dependent, small-sample descriptive screens; zero-hit SE passes excluded', methods={})
    for method, g in fine.groupby('method'):
        h = g[g.estimator == 'hedging']
        raw = g[g.estimator == 'raw']
        result['methods'][method] = dict(contract_trials=len(h),
            positive_hit_screens=int(((h.positive_payoffs > 0) & (h.se <= h.scale/8)).sum()),
            zero_hit_contract_trials=int((h.positive_payoffs == 0).sum()), negative_adjusted_prices=int((h.price < 0).sum()),
            median_hedging_se_over_scale=float((h.se/h.scale).median()),
            median_raw_se_over_scale=float((raw.se/raw.scale).median()))
    h = fine[fine.estimator == 'hedging'].pivot(index=['trial','contract'], columns='method', values='se')
    valid = (h > 0).all(axis=1)
    gain = (h.loc[valid,'target']/h.loc[valid,'importance'])**2
    result['variance_ratio_target_over_importance'] = dict(valid_contract_trials=int(valid.sum()),
        excluded_zero_se=int((~valid).sum()), median=float(gain.median()), min=float(gain.min()), max=float(gain.max()),
        improved=int((gain > 1).sum()))
    result['refinements'] = []
    for (method, kind), g in c[c.estimator == 'hedging'].groupby(['method','change_type']):
        ratio = g.change.abs()/g.scale
        worst = g.iloc[int(np.argmax(ratio.to_numpy()))]
        result['refinements'].append(dict(method=method, change_type=kind,
            median_abs_change_over_scale=float(ratio.median()), max_abs_change_over_scale=float(ratio.max()),
            observed_changes_above_scale_eighth=int((ratio > .125).sum()),
            worst=dict(trial=int(worst.trial), contract=int(worst.contract), change=float(worst.change), se=float(worst.se), scale=float(worst.scale))))
    result['monitors'] = d.to_dict(orient='records')
    (output/'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'monitors'}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    audit(p.parse_args().output)
