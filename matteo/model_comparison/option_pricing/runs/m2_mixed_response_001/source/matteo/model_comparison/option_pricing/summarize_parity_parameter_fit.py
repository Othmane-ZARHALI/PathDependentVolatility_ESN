"""Readable before/after joint calibration and independent numerical check."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import NAMES


def summarize(run,out):
    r=json.loads((run/'run.json').read_text());audit=json.loads((run/'audit.json').read_text());uncertainty=json.loads((run/'checking/validation/block_uncertainty.json').read_text())
    initial=np.array(r['baseline']);selected=np.array(r['selected']);lo=np.array(r['bounds_lower']);hi=np.array(r['bounds_upper'])
    accepted=sum(x['decision']['accepted'] for x in r['outer_proposals'])
    lines=['# Eleven-parameter joint recalibration with 4096 paths','',
           'Completed bounded calibration attempt. All eleven agreed parameters were free; OU reservoir, Q prices, initial state and geometric grids remain fixed.',
           '',f'Run: `{run.name}`. Start at the previous fitted positive-shift vector, not the historical zero-shift reference. '
           'Keep all 1294 SPX cutoff 3 OTM observations at seven maturities, 217 VIX OTM observations at four maturities and four VIX forwards. '
           'Every strike uses 4096 paths; VIX uses eight independent 512-point nets. GH32/GH48/time8/stock512 remain unchanged.',
           '', '## Results on identical inputs','',
           'Loss measures squared price error outside bid/ask, divided by the existing quote scale, then averaged equally over maturities within SPX and VIX. A price inside bid/ask contributes zero. Forward error is scaled by 0.5 VIX points. These are price scores, not percentages or IV errors.', '',
           '| Ensemble / convention | Before | After |','|---|---:|---:|']
    for role,key in [('Training','training_scores'),('Independent checking','checking_scores')]:
        a,b=r[key]
        for label,field in [('Common-call loss','canonical_loss'),('Original-OTM loss','original_otm_loss')]:
            lines.append(f'| {role}: {label} | {a[field]:.6f} | {b[field]:.6f} |')
        for label,field in [('SPX component','SPX'),('VIX-option component','VIX'),('Forward component','forwards')]:
            lines.append(f"| {role}: {label} | {a['canonical_components'][field]:.6f} | {b['canonical_components'][field]:.6f} |")
        lines.append(f"| {role}: original-OTM VIX-option component | {a['original_otm_components']['VIX']:.6f} | {b['original_otm_components']['VIX']:.6f} |")
        lines.append(f"| {role}: original-OTM VIX plus forwards | {a['original_otm_components']['VIX']+a['original_otm_components']['forwards']:.6f} | {b['original_otm_components']['VIX']+b['original_otm_components']['forwards']:.6f} |")
    lines+=['','The checking streams are NEW and independent of training, with 4096 paths per SPX expiry and 8x512 VIX points. '
            'Parameters were frozen before checking prices were read. Both candidates use identical checking draws and fixed controls. '
            'All market observations still enter calibration: this is numerical checking, not market out-of-sample performance.',
            '', 'The first prototype (run001) was stopped after its first true proposal: common-call joint loss fell from98.119 to79.216, but native OTM VIX-plus-forward loss worsened by9.060. The guard rejected it. Run003 aligns its inner constraints with that safeguard; the rejected vector is not an adopted calibration.',
            '', '## Parameters','', '| Parameter | Before | After | Search bounds | Near bound |','|---|---:|---:|---|---|']
    for n,a,b,l,h in zip(NAMES,initial,selected,lo,hi):
        tag='lower' if (b-l)/(h-l)<.01 else ('upper' if (h-b)/(h-l)<.01 else '')
        lines.append(f'| {n} | {a:.9g} | {b:.9g} | [{l:.6g}, {h:.6g}] | {tag} |')
    lines+=['','## Independent VIX forwards','', '| Days | Market | Before | After | After SE | Error after |','|---:|---:|---:|---:|---:|---:|']
    before=pd.read_csv(run/'checking/forwards_baseline.csv');after=pd.read_csv(run/'checking/forwards_selected.csv')
    for a,b in zip(before.itertuples(),after.itertuples()):lines.append(f'| {b.expiry_days} | {b.forward_vix_points:.6f} | {a.model_forward:.6f} | {b.model_forward:.6f} | {b.mc_se:.6f} | {b.error_points:+.6f} |')
    lines+=['','## Every retained maturity','', '| Instrument | Days | Checking common-call MSE before | After | Native OTM MSE before | After | Native quotes in spread after |','|---|---:|---:|---:|---:|---:|---:|']
    q=[pd.read_csv(run/'checking'/f'prices_canonical_{n}.csv') for n in ['baseline','selected']]
    o=[pd.read_csv(run/'checking'/f'prices_otm_{n}.csv') for n in ['baseline','selected']]
    rows=[]
    for (ins,day),g in q[0].groupby(['instrument','expiry_days']):
        t=q[1][(q[1].instrument==ins)&(q[1].expiry_days==day)];native=[x[(x.instrument==ins)&(x.expiry_days==day)] for x in o]
        vals=[float(np.mean(x.interval_residual**2)) for x in [g,t,*native]];inside=int(native[1].interval_residual.eq(0).sum())
        lines.append(f'| {ins} | {day} | {vals[0]:.4f} | {vals[1]:.4f} | {vals[2]:.4f} | {vals[3]:.4f} | {inside}/{len(t)} |')
        rows.append(dict(instrument=ins,days=int(day),canonical_before=vals[0],canonical_after=vals[1],native_before=vals[2],native_after=vals[3],inside_after=inside))
    pd.DataFrame(rows).to_csv(run/'checking/expiry_summary.csv',index=False)
    diagnostic=pd.read_csv(run/'checking/vix_smile_diagnostics.csv')
    lines+=['','## High-strike VIX smile still misses','',
            '| Days | IV MAE before, percentage points | After |', '|---:|---:|---:|']
    for day in sorted(diagnostic.days.unique()):
        t=diagnostic[(diagnostic.days==day)&diagnostic.region.eq('K/F > 2')].set_index('candidate')
        lines.append(f"| {day} | {t.loc['baseline','IV_MAE_percentage_points']:.4f} | {t.loc['selected','IV_MAE_percentage_points']:.4f} |")
    lines+=['','The price objective is not an IV-error objective. Large IV discrepancies can correspond to comparatively small prices in the wings. '
            'The lower aggregate price loss does not clear the visible smile misfit, and the high-strike errors above must not be hidden.',
            '', '## Independent paired VIX loss sensitivity','',
            '| Criterion | Before | After | Change | Eight-block jackknife SE of change |',
            '|---|---:|---:|---:|---:|']
    for metric in ['canonical_VIX_plus_forwards','native_VIX_plus_forwards']:
        u=uncertainty['paired_components'][metric]
        lines.append(f"| {metric} | {u['before']:.6f} | {u['after']:.6f} | {u['change']:+.6f} | {u['jackknife_SE']:.6f} |")
    lines+=['','These paired changes use the same checking paths for both candidates and frozen controls. The eight-block diagnostic is descriptive; it excludes parameter and deterministic quadrature uncertainty.',
            '', '## What the optimizer did','',
            'The objective keeps the previous family weights, quote scales, interval loss and penalties. '
            'It prices one call payoff at all strikes, converting only the existing OTM market inputs. '
            'Native model puts use the MODEL forward; no observed-forward substitution is used to hide a VIX mismatch.',
            '', 'At each outer anchor, calculate full joint 4096-path prices and analytic derivatives. '
            'Inside a bounded parameter box, use first-order SPX prices and actual nonlinear 4096-point VIX prices/forwards. '
            'A VIX-focused SLSQP solve proposes the next vector, prioritizing original-OTM VIX-option loss plus forwards. The warm proposal comes from the preceding average-VIX local solve. Both local VIX and joint criteria are constrained to improve, with margins for SPX approximation error. A normalized 1e-4 tolerance applies only to these extra reserve margins; the actual-price improvement tests remain strict. Reprice the proposal with the FULL SPX/VIX simulator before accepting it. '
            'An acceptance requires lower true joint and VIX-plus-forward losses in BOTH common-call and original-OTM conventions. '
            'This checked local approximation reduces expensive SPX evaluations; no surrogate score is represented as an actual fit.',
            '', f'{len(r["outer_proposals"])} full outer proposal(s), {accepted} accepted. Two subsequent local cycles reached their 20-evaluation limits without a reserve-feasible proposal. Maximum 3 proposals and 75 optimization minutes were prespecified. '
            'All inner/outer evaluations and rejection reasons are retained. Termination at the declared budget is not convergence, global optimality or parameter identification.',
            '', '## Numerical checking and limitations','',
            f"Saved prices/SEs/parity/both objectives and every acceptance decision replay; {audit['sources_verified']} source/input files verify. "
            'Production baseline replay and derivative correspondence screens pass. '
            'Five new tests check cached VIX equivalence, local objective/native criteria derivatives and true-price acceptance rules. The unchanged true pricer passed production VIX/full-price directional checks in retained run001; run002 adds a production native/canonical-criteria derivative screen. Run003 restores the verified full baseline Jacobian checkpoint and reuses those screens.',
            '', 'SPX SE is sample SD/sqrt(4096). VIX SE is SD of eight complete independent scrambled-net prices/sqrt(8). '
            'They are pointwise sampling estimates, conditional on frozen pilot controls, excluding pilot/parameter uncertainty and deterministic bias. '
            'No clipping or post-price/IV smoothing is applied. Invalid IVs and signed prices are retained in the diagnostic tables.',
            '',f"Independent checking paired loss change: {uncertainty['paired_change']['full']:+.6f}; descriptive eight-block jackknife SE: {uncertainty['paired_change']['jackknife_SE']:.6f}. "
            'Inspect leave-one-block changes and largest stock contributions; this is not a price-convergence certificate.',
            '', '| Checking candidate | Inside native OTM bid/ask | Negative native prices | Unavailable common-call IVs | SPX order failures | SPX convexity failures | VIX order / convexity |',
            '|---|---:|---:|---:|---:|---:|---|']
    for s in r['checking_scores']:
        sp=[x for x in s['shapes'] if x['instrument']=='SPX'];vx=[x for x in s['shapes'] if x['instrument']=='VIX']
        lines.append(f"| {s['candidate']} | {s['inside_original_otm_bidask']}/1511 | {s['native_otm_negative_prices']} | {s['iv_unavailable']} | {sum(x['monotonicity_violations'] for x in sp)} | {sum(x['convexity_violations'] for x in sp)} | {sum(x['monotonicity_violations'] for x in vx)} / {sum(x['convexity_violations'] for x in vx)} |")
    selected_q=q[1];gaps=selected_q[selected_q.iv_unavailable].groupby(['instrument','expiry_days']).size().to_dict()
    lines+=['',f'Independent selected IV gaps by instrument/expiry: {gaps}. Eight selected SPX signed estimates are negative at 57/302 days; '
            'two VIX55 positive calls remain below observed-forward intrinsic. These fresh-stream failures show that common-call parity and4096 paths do not establish estimator convergence. '
            'The previous fixed-parameter ensemble had no negative selected native prices; that fact did not carry over to independent draws. No clipping is applied.']
    zeros={n:sum(x['zero_hit_contracts'] for x in audit['checking'][n]['tails']) for n in ['baseline','selected']}
    lines+=['',f'Checking zero-raw-payoff-hit contract counts: {zeros}. Zero empirical SE/hits cannot clear unseen tail contributions. '
            'The model fit can remain outside bid/ask even when sampling error decreases. '
            'No new P compatibility assessment or Lean certificate is asserted. The selected model remains exploratory until its actual residuals and numerical limitations are accepted.',
            '', 'Plots: [all seven SPX maturities](visualizations/parity-parameter-fit-011/spx-all-maturity-fit.png), '
            '[all four VIX maturities](visualizations/parity-parameter-fit-011/vix-all-maturity-fit.png), '
            '[VIX forwards](visualizations/parity-parameter-fit-011/forwards-and-path-budgets.png). '
            'All plots compare old/new parameters on identical NEW checking draws. Signed errors and SVG exports are included.',
            '', f'Selected configuration: [selected_model.json](runs/{run.name}/selected_model.json); full records and checking data under '
            f'`runs/{run.name}`. New checking cache/control directories are separate from training; previous runs and plots remain untouched.', '']
    out.write_text('\n'.join(lines));print('REPORT',out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();summarize(a.run.resolve(),a.output.resolve())
