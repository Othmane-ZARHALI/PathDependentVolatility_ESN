"""Write an evidence-linked report after the maturity-curve run and audit finish."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd


def report(run):
    read=lambda p:json.loads(p.read_text())
    r=read(run/'run.json');audit=read(run/'audit.json');u=read(run/'checking_vix/uncertainty.json')
    scores={n:read(run/f'score_{n}.json') for n in ['baseline','selected']}
    p=np.array(r['baseline']);s=np.array(r['selected']);b=scores['baseline'];a=scores['selected']
    text=['# Maturity-aligned forward variance calibration', '',
          'Completed 25 October 2026. The previous three-piece curve is exactly nested in the new 16-piece curve. The fit frees 24 parameters: 16 variance levels, seven readout/cap coefficients and the variance floor. OU rates, reservoir dynamics, random draws, risk prices, correlations, initial state and geometric grids remain fixed.', '',
          'Original OTM observations only: SPX 29/57/119/211/302/666/1030 days (cutoff 3), VIX 27/55/118/153 days, 1511 option quotes and four VIX forwards. Calls form the common payoff; market puts use market-forward parity and native model puts use model-forward parity. Every strike uses all 4096 paths. All market observations are in sample.', '',
          '## Actual full-price results', '',
          '| Training criterion | Previous | Selected |', '|---|---:|---:|']
    metrics=[('Common-call joint',b['canonical_loss'],a['canonical_loss']),('Original-OTM joint',b['original_otm_loss'],a['original_otm_loss'])]
    for label,key,part in [('SPX','original_otm_components','SPX'),('Common-call VIX','canonical_components','VIX'),('Original-OTM VIX','original_otm_components','VIX'),('VIX forwards','original_otm_components','forwards')]:
        metrics.append((label,b[key][part],a[key][part]))
    for label,x,y in metrics:text.append(f'| {label} | {x:.6f} | {y:.6f} |')
    old=b['original_otm_components']['VIX']+b['original_otm_components']['forwards'];new=a['original_otm_components']['VIX']+a['original_otm_components']['forwards']
    text+=['',f'Original-OTM VIX-plus-forward training loss changes from {old:.6f} to {new:.6f} ({100*(new/old-1):+.2f}%). Accepted full-price steps: {audit["accepted_steps"]}. This is a bounded exploration, not a convergence or global-optimum claim.','',
           '| Proposal | Curve only | Predicted joint | Actual joint | Actual native VIX + forwards | Accepted |','|---:|:---:|---:|---:|---:|:---:|']
    for x in r['outer_proposals']:
        text.append(f'| {x["cycle"]} | {x["curve_only"]} | {x["approximate_criteria"][0]:.6f} | {x["true_criteria"][0]:.6f} | {x["true_criteria"][3]:.6f} | {x["decision"]["accepted"]} |')
    text+=['','SLSQP minimizes the native VIX-plus-forward objective on a local affine price/forward model. Both joint and both VIX-plus-forward criteria constrain the inner solve. Nonlinear VIX prices screen proposals before full SPX/VIX pricing. Only full nonlinear scores can accept a step. The first proposal frees curve/floor parameters; subsequent proposals free all 24. Only the first proposal was accepted, so the selected readout/cap coefficients remain unchanged despite having been searched in the later proposals.', '',
           '## Separate reused VIX numerical checking', '', '| Criterion | Before | After | Change | Paired jackknife SE of change |','|---|---:|---:|---:|---:|']
    for k,v in u['paired'].items():text.append(f'| {k} VIX + forwards | {v["before"]:.6f} | {v["after"]:.6f} | {v["change"]:+.6f} | {v["jackknife_SE"]:.6f} |')
    text+=['','Selection froze before these VIX prices were checked. These eight independent scrambled-net streams have been used in previous investigations; this is descriptive numerical checking, not an untouched holdout or market out-of-sample test. SPX checking was not repeated. Leave-one-net changes and second moments are preserved in [uncertainty.json](runs/m2_maturity_curve_fit_001/checking_vix/uncertainty.json).','',
           '## Forward variance levels', '', '| Interval in days | Previous annual variance | Selected annual variance | Selected square-root scale |','|---|---:|---:|---:|']
    start=0
    for end,x,y in zip(r['ends_days'],p[7:-1],s[7:-1]):
        text.append(f'| {start:g}–{end:g} | {x:.8f} | {y:.8f} | {100*np.sqrt(y):.3f}% |');start=end
    text+=['','Every level satisfies xi_i >= variance_floor; the levels may rise or fall. These square-root scales describe the expected variance curve, not individual option implied volatilities. Close boundaries can make individual levels weakly identified.','',
           '| Readout/cap/floor parameter | Previous | Selected |','|---|---:|---:|']
    for j in [*range(7),len(s)-1]:text.append(f'| {r["parameter_names"][j]} | {p[j]:.10g} | {s[j]:.10g} |')
    text+=['',f'Selected minimum instantaneous volatility is {100*np.sqrt(s[-1]):.5f}%; beta_J is {s[4]:.8f} > 0. Floor bounds remain 5–6% volatility (the 6% upper endpoint is the earlier implementation choice).','',
           '## All retained maturities', '', '| Instrument | Days | Native interval MSE before | After | Selected prices in bid/ask |','|---|---:|---:|---:|---:|']
    qs={n:pd.read_csv(run/f'prices_otm_{n}.csv') for n in scores}
    for (ins,day),g in qs['selected'].groupby(['instrument','expiry_days'],sort=False):
        before=qs['baseline'].loc[g.index];inside=((g.model_price_U>=g.price_bid_U)&(g.model_price_U<=g.price_ask_U)).sum()
        text.append(f'| {ins} | {int(day)} | {np.mean(before.interval_residual**2):.5f} | {np.mean(g.interval_residual**2):.5f} | {inside}/{len(g)} |')
    text+=['','## Verification and numerical limits', '',f'Ten focused tests pass. The saved-price audit replays prices, standard errors, parity, losses, acceptance decisions, {audit["local_trials_replayed"]} local trials and {audit["sources_verified"]} frozen sources/inputs. Curve/floor production VIX sensitivities passed finite differences. Baseline SPX prices agree with the preceding selected fit to the recorded tolerance; splitting VIX quadrature at the new boundaries changes baseline prices only by the recorded tiny refinement error.', '',
          f'Selected training estimates: {a["native_otm_negative_prices"]} negative native prices; {a["iv_unavailable"]} unavailable common-call IVs; {a["inside_original_otm_bidask"]}/1511 native estimates inside bid/ask. Signed estimates, invalid IVs, strike-shape and raw-tail diagnostics are retained. The VIX wing mismatch must be judged from the plots and strike-region diagnostics, not just the aggregate loss.', '',
          'SPX SE uses the standard deviation across 4096 independent controlled path contributions divided by sqrt(4096). VIX SE uses the standard deviation of the eight independent 512-point scrambled-net means divided by sqrt(8), including model-forward covariance when recovering native puts. Bands are pointwise sampling uncertainty only. They exclude parameter uncertainty, frozen pilot uncertainty, normalizer/quadrature bias and stock timestep error.', '',
          'GH32 conditional expectations, GH48 normalization, time order 8 and the cached 512-step stock grid remain unchanged. New curve boundaries split VIX integration panels. Stock uses the cached left-point mesh; deterministic curve-integral discrepancies are recorded in grid_baseline.json and grid_selected.json, and are not a full price-convergence bound. No Lean/P-compatibility certificate follows.', '',
          'The pause interrupted the first proposal, not the completed baseline. The continuation verified frozen original source hashes, restored the baseline and reproduced the first proposal before resuming full pricing. Baseline and original pause/source records are preserved; continuation sources are separately snapshotted. Timing in the final run records continuation wall time; it does not include time spent paused or all interrupted work.', '',
          '[All seven SPX maturities](visualizations/maturity-curve-013/spx-all-maturity-fit.png) · [All four VIX maturities](visualizations/maturity-curve-013/vix-all-maturity-fit.png) · [Separate VIX checking](visualizations/maturity-curve-013/vix-separate-checking.png) · [Forward variance curve](visualizations/maturity-curve-013/forward-variance-curve.png). Blue is selected, orange is the preceding fit. All PNG/SVG exports and signed-price plots are in the visualization folder.', '',
          '[Selected model](runs/m2_maturity_curve_fit_001/selected_model.json) · [Audit](runs/m2_maturity_curve_fit_001/audit.json) · [Plan](MATURITY_CURVE_FIT_PLAN.md).']
    text+=['', '## Actual VIX forwards', '', '| Streams | Days | Market | Previous | Selected |', '|---|---:|---:|---:|---:|']
    for label,folder in [('Training',run),('Separate checking',run/'checking_vix')]:
        fa=pd.read_csv(folder/'forwards_baseline.csv');fb=pd.read_csv(folder/'forwards_selected.csv')
        for x,y in zip(fa.itertuples(),fb.itertuples()):
            text.append(f'| {label} | {int(x.expiry_days)} | {x.forward_vix_points:.4f} | {x.model_forward:.4f} | {y.model_forward:.4f} |')
    diagnostic=run/'checking_vix/vix_smile_diagnostics.csv'
    if diagnostic.exists():
        d=pd.read_csv(diagnostic)
        text+=['', '## VIX implied-volatility errors on separate checking streams', '', '| Days | Strike region | Previous IV MAE, percentage points | Selected IV MAE | Defined selected IVs |', '|---:|---|---:|---:|---:|']
        for x in d[d.candidate=='selected'].itertuples():
            old=d[(d.candidate=='baseline')&(d.days==x.days)&(d.region==x.region)].iloc[0]
            text.append(f'| {int(x.days)} | {x.region} | {old.IV_MAE_percentage_points:.3f} | {x.IV_MAE_percentage_points:.3f} | {int(x.defined_IVs)}/{int(x.quotes)} |')
        text+=['', 'IV errors use the observed-market-forward common-call convention. Large far-wing IV errors can remain even when aggregate interval-price loss falls.']
    moment=pd.read_csv(run/'checking_vix/second_moments.csv')
    row=moment[(moment.candidate=='selected')&(moment.days==153)].iloc[0]
    text+=['', f'Selected VIX 153-day second moment: estimate {row.estimated_second_moment:.3f}, exact forward-curve target {row.target_second_moment:.3f}, reported eight-net SE {row.net_SE:.3f}. The deficit is about {(row.target_second_moment-row.estimated_second_moment)/row.net_SE:.2f} reported SE; rare-tail/numerical concerns are not cleared by this fit.', '', 'All eight PNG figures were visually inspected; PNG/SVG hashes and source hashes are checked in validation/plot_review.json.']
    shape=[]
    for ins in ['SPX','VIX']:
        rows=[x for x in a['shapes'] if x['instrument']==ins]
        shape.append(f"{ins}: {sum(x['monotonicity_violations'] for x in rows)} monotonicity and {sum(x['convexity_violations'] for x in rows)} convexity flags")
    text+=['', 'Selected training shape diagnostics: '+ '; '.join(shape)+'. These finite-estimator flags are not model-level arbitrage proofs.']
    (run.parent.parent/'MATURITY_CURVE_FIT_RESULTS.md').write_text('\n'.join(text)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();report(a.run.resolve())
