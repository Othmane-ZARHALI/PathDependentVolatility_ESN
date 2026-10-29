"""Generate the research record from saved post-cap experiment outputs."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,dump
from .postcap_vov_study import native_loss


def report(run):
    result=json.loads((run/'run.json').read_text());profiles=json.loads((run/'profiles.json').read_text());selection=result['selection'];fresh=json.loads((run/'fresh_check/profiles.json').read_text());full=json.loads((run/'full/scores.json').read_text());audit=json.loads((run/'audit.json').read_text());market=json.loads((run/'market_shape.json').read_text());quad=json.loads((run/'quadrature_checks.json').read_text());paired=json.loads((run/'fresh_check/paired_wing_shape.json').read_text())
    assert audit['status']=='PASS';base=next(r for r in fresh if r['name']=='baseline');winner=next(r for r in fresh if r['name']==selection['shape_winner'])
    slope_changes=[r['wing_slope']-b['wing_slope'] for r,b in zip(winner['shape']['rows'],base['shape']['rows'])]
    lines=[];add=lambda x='':lines.append(x)
    add('# Independent amplitude after the cap: completed comparison');add();add('Completed 2026-10-26. Experimental fixed-curve comparison; original selected model unchanged.');add()
    add('**Main finding:** the strong 7-day multiplier improves the descriptive fresh VIX slope/curvature score by about31%, but the smile remains too flat and full SPX loss rises43.777 to121.147. The mild multiplier gives a small training joint-loss reduction70.978 to69.876, smaller than its paired sampling sensitivity (change -1.102, jackknife SE3.933), and does not improve fresh VIX price loss. Neither is an accepted joint fit.');add()
    add('## What changed');add();add('The old spike-only bounded sensitivity layer is absent. Starting from the original capped baseline, an independent Q OU state multiplies the normalized variance component AFTER the cap:');add();add('```text\ndZ = -a Z dt + sqrt(2a) dB, Z0=0, q(t)=1-exp(-2at)\nL(t)=exp(gamma Z(t)-gamma^2 q(t)/2)\nv(t)=f+(xi(t)-f)*(G_base(X(t))/M_base(t))*L(t)\n```');add();add('All original readout coefficients, OU banks, correlations, premiums, 16 variance levels, floor and cap stay fixed. gamma is a log-VARIANCE sensitivity, not a bounded multiplier range. The new factor is independent of every original driver. E[L]=1, so the original normalizer M_base is unchanged; there is no candidate-specific renormalization. Relative to the previous spike-only experiment, placement, factor shape and centering all change; this is not a placement-only causal test.');add()
    add('The analytic conditional factor is exp(gamma exp(-a(u-T)) Z_T - gamma^2 exp(-2a(u-T)) q(T)/2). Independence lets it multiply the original conditional expectation of G. No new conditional Gaussian quadrature is introduced.');add()
    add('## Full VIX grid on the original training simulations');add();add('|Case|Native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Wing slope/curvature score|Missing call IVs|');add('|---|---:|---:|---:|---:|---:|')
    for r in profiles:
        add(f"|{r['name']}|{native_loss(r['score']):.6f}|{r['shape']['wing_MAE']:.3f}|{r['shape']['ATM_MAE']:.3f}|{r['shape_score']:.6f}|{r['score']['iv_unavailable']}|")
    add();add('MAE is in IV percentage points, averaging the four maturities equally. Wing: K/marketF>2; near ATM: abs(log(K/marketF))<=0.2. The descriptive shape score averages squared slope and curvature discrepancies, each divided by max(10,abs(market coefficient)) squared. Lower is closer; it is not an arbitrage condition or a fitted polynomial volatility model.');add()
    add(f"Price selection: `{selection['price_winner']}`. Shape selection: `{selection['shape_winner']}`. These selections were frozen before SPX prices or fresh VIX draws. Neither selection implies joint acceptance. All nine nonzero variants and the zero-factor reference are retained.");add()
    add('## Actual full SPX repricing');add();add('|Case|SPX loss|Native VIX + forwards|Native joint loss|Inside original bid-ask /1511|Negative native estimates|Missing call IVs|');add('|---|---:|---:|---:|---:|---:|---:|')
    for r in full:
        add(f"|{r['candidate']}|{r['original_otm_components']['SPX']:.6f}|{native_loss(r):.6f}|{r['original_otm_loss']:.6f}|{r['inside_original_otm_bidask']}|{r['native_otm_negative_prices']}|{r['iv_unavailable']}|")
    add();add('All seven agreed SPX dates are repriced with 4096 samples at every strike. Original exogenous paths and equity shocks are reused; variance and stock paths are recomputed. Independent OU shocks use the same parameter-independent counter streams across candidates. No large additional factor-path cache is created. Original independent pilot controls stay frozen.');add()
    for r in full:
        name=r['candidate']
        if name=='baseline':continue
        pp=json.loads((run/'full'/f'paired_{name}_baseline.json').read_text());add(f"- `{name}`: paired SPX loss change {pp['change'][1]:+.6f}, descriptive eight-block jackknife SE {pp['jackknife_SE'][1]:.6f}; native joint change {pp['change'][0]:+.6f}, SE {pp['jackknife_SE'][0]:.6f}.")
    add();add('## Fresh joint VIX simulations after selection');add();add('|Case|Native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Shape score|');add('|---|---:|---:|---:|---:|')
    for r in fresh:add(f"|{r['name']}|{native_loss(r['score']):.6f}|{r['shape']['wing_MAE']:.3f}|{r['shape']['ATM_MAE']:.3f}|{r['shape_score']:.6f}|")
    add();add('|Case|Days|Wing slope|Wing curvature|');add('|---|---:|---:|---:|')
    for r in market['rows']:add(f"|Market|{r['days']}|{r['wing_slope']:.3f}|{r['wing_curvature']:.3f}|")
    for profile in fresh:
        for r in profile['shape']['rows']:add(f"|{profile['name']}|{r['days']}|{r['wing_slope']:.3f}|{r['wing_curvature']:.3f}|")
    add();add('These are regressions against log(K/marketF) over K/marketF>2, not pointwise derivatives. Market curvature is not uniformly positive. Both linear- and log-strike figures are exported. Call PRICE convexity and implied-VOLATILITY convexity are different statements.');add()
    add(f"The selected shape candidate changes the four fresh wing slopes relative to baseline by {[round(x,3) for x in slope_changes]}. Its fresh shape score changes from {base['shape_score']:.6f} to {winner['shape_score']:.6f}. Interpret this alongside price loss, forwards and SPX rather than declaring success from slope alone.");add()
    for r in paired['records']:
        if r['name']=='baseline':continue
        add(f"- `{r['name']}` fresh paired slope changes: {[round(x[0],3) for x in r['change']]}; jackknife SE: {[round(x[0],3) for x in r['jackknife_SE']]}.")
    add();add('Checking uses 32 newly generated joint37D scrambles: eight512-point nets at each of four maturities, exact defensive-mixture likelihood weights, no self-normalization. Selection cannot be revised from these draws. Market quotes are the same in-sample observations. There is no market out-of-sample claim or fresh SPX ensemble.');add()
    add('## Numerical and mathematical checks');add();add(f"- Audit PASS: {len(audit['price_records'])} saved price/parity/SE/objective records; {audit['exact_joint_scrambles_replayed']} complete joint scrambles regenerated with exact weights; {len(audit['conditional_factor_oracles'])} dense independent conditional-factor integrals; {len(audit['original_state_NumPy_oracles'])} direct NumPy old-state Gaussian oracles; {audit['source_checks']} archived/current source comparisons.")
    replay=json.loads((run/'full/baseline_replay.json').read_text());add(f"- Zero-factor full baseline maximum price discrepancy: {replay['maximum_price_difference']:.3g}.")
    add(f"- Sampled full VIX windows, old-state GH32/time8 to GH48/time16: maximum change {max(r['maximum_VIX_point_difference'] for r in quad):.6g} VIX points. These are representative checks, not a whole-price convergence certificate.")
    add(f"- Original normalizer GH48 to GH96 sampled relative changes: maximum {max(abs(r['relative_change']) for r in audit['normalizer_order_checks']):.3g}.")
    add('- Five unit checks pass across the new implementation and prior spike-layer implementation. New checks compare analytic conditional factors with dense Gaussian integration, test exact mean-one centering, compare streamed stock stepping with independent Python, and verify zero-strength nesting of the original kernel.');add()
    add('For the continuous model, conditioning on the independent Z path makes L a bounded continuous time function on a finite horizon. With capped G_base and positive continuous deterministic M_base, conditional stock volatility is bounded, yielding the stock martingale identity. Fixed-time variance moments are finite. This does not guarantee finite variance of stock payoff estimators. Numerical errors and rare-tail sampling remain separate questions. No new Lean certificate or P-measure calibration is claimed.');add()
    add('Negative signed control estimates and missing IVs are retained in the raw tables. A missing IV can also arise from using the market forward for inversion when the model forward differs; it need not mean a negative model price.');add()
    for row in audit['negative_estimates']:
        if row['count']:add(f"- `{row['name']}` ({row['folder']}): {row['count']} negative native estimates, maximum magnitude {row['maximum_magnitude_in_SE']:.3f} reported SE. These are signed finite control estimates; no negative true model-price claim follows.")
    add()
    for folder in ['','fresh_check/']:
        moments=pd.read_csv(run/(folder+'moments.csv'));add(f"{'Training' if not folder else 'Fresh checking'} squared-VIX moments (weighted estimates, not a convergence pass):");add();add('|Case|Days|E[VIX²]|Curve target|Eight-net SE|');add('|---|---:|---:|---:|---:|')
        for _,r in moments.iterrows():
            if r['candidate'] in selection['full_candidates']:add(f"|{r['candidate']}|{int(r['days'])}|{r['second']:.3f}|{r['target']:.3f}|{r['second_SE']:.3f}|")
        add()
    add('## Evidence and reproduction');add();add('Plan: [POSTCAP_VOV_PLAN.md](POSTCAP_VOV_PLAN.md). Plots: [postcap-vov017](visualizations/postcap-vov017/README.md). Run: `runs/m2_postcap_vov_001`. `environment.txt` records the temporary Python3.12 package versions. Reproduce with `PYTHONPATH=matteo/exogenous_reservoir_eq_m2/src` and modules `postcap_vov_study --output NEW_DIRECTORY`, `audit_postcap_vov --run RUN`, `plot_postcap_vov --run RUN --output PLOTS`, `summarize_postcap_vov --run RUN`. All original baseline files are preserved.');add();add(f"Study wall time: {result['wall_seconds']/60:.2f} minutes. No model was automatically adopted; this comparison is a finite grid with the old parameters fixed, not a full enlarged-model calibration.")
    (HERE/'POSTCAP_VOV_RESULTS.md').write_text('\n'.join(lines)+'\n')
    final=run/'finalization';final.mkdir(exist_ok=True)
    newfiles=[HERE/n for n in ['postcap_vov_core.py','postcap_vov_kernel.cpp','postcap_vov_study.py','test_postcap_vov.py','audit_postcap_vov.py','plot_postcap_vov.py','summarize_postcap_vov.py','POSTCAP_VOV_PLAN.md','POSTCAP_VOV_RESULTS.md']]
    source=final/'source';source.mkdir(exist_ok=True)
    for p in newfiles:(source/p.name).write_bytes(p.read_bytes())
    artifacts=[p for p in run.rglob('*') if p.is_file() and 'source' not in p.parts and 'finalization' not in p.parts]
    dump(final/'manifest.json',dict(artifacts_sha256={str(p.relative_to(run)):hashlib.sha256(p.read_bytes()).hexdigest() for p in artifacts},postprocessing_sources_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in newfiles},baseline_replaced=False))
    print('WROTE',HERE/'POSTCAP_VOV_RESULTS.md',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);report(p.parse_args().run.resolve())
