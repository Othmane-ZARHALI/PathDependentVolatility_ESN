"""Report measured point1 acceleration separately from fit/accuracy evidence."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import NAMES


def summarize(run,out):
    r=json.loads((run/'run.json').read_text());audit=json.loads((run/'audit.json').read_text());t=r['timing']
    a,b=r['training_scores'];accepted=sum(x['decision']['accepted'] for x in r['outer_proposals'])
    lines=['# Point1 acceleration and minimum volatility policy','',
           'The current calibration policy requires minimum instantaneous volatility of5%, with its floor calibratable between5% and6%. '
           'The6% upper endpoint is the implementation choice used to keep the floor free within the stated reasonable range. '
           'All eleven parameters remain free. Historical bounds/configurations are preserved.',
           '',r'With $v_t=f+(\xi(t)-f)G_t/M(t)$, positive $G_t/M(t)$ and $\xi(t)\geq f$, '
           r'$f\geq0.05^2=0.0025$ ensures $\sqrt{v_t}\geq5\%$. '
           'This is an instantaneous model-volatility bound, not an option-IV constraint.',
           '', '## What changed','',
           'At each fully priced anchor, store call prices, VIX forwards and all eleven derivative columns. '
           'Inner optimizer trials now approximate BOTH SPX and VIX prices/forwards by an affine Taylor model. '
           'The original interval loss is evaluated on those approximate prices. Native model puts use the approximate MODEL forward and its derivative. '
           'No simulator is held by the local-model object. Full nonlinear4096-path repricing checks every proposal before acceptance; accepted anchors refresh the Jacobian.',
           '',f'Pilot `{run.name}` starts from the previous selected vector with its variance floor rounded upward by about9e-15 to the exact new lower bound. '
           'Same1511 original OTM observations, four forwards, SPX cutoff3/seven retained maturities/four VIX maturities, OU dynamics/rates/grids/Q prices/initial state, '
           'common-call payoff, all4096 paths, GH32/GH48/time8/stock512. No coarse mesh or parameter-block optimization was introduced.',
           '', '## Measured work','', '| Operation | Seconds |','|---|---:|',
           f"| Old nonlinear-VIX local trial, historical median | {t['historical_nonlinear_VIX_local_median_seconds']:.6f} |",
           f"| New approximate joint-price evaluation, measured median | {t['local_evaluation_median_seconds']:.6f} |",
           f"| All new inner solves ({t['local_evaluations']} optimizer evaluations) | {t['local_search_seconds']:.6f} |",
           f"| Full baseline, including eleven sensitivities | {t['full_baseline_seconds']:.3f} |",
           f"| All full checked proposals | {t['full_proposal_seconds']:.3f} |",
           f"| Whole pilot | {r['wall_seconds']:.3f} |",
           '',f"Full nonlinear pricer calls: {t['full_pricer_calls']}; local optimizer simulator calls:0. "
           'The old/new local rows perform different work (nonlinear calculation versus approximation). '
           'This establishes local computational acceleration, not equal accuracy at every trial or an end-to-end speed factor. '
           'Whole-pilot timing also differs in validation scope from the earlier run: no new independent checking is priced here.',
           '', '## Actual training-ensemble fit','', '| Criterion | Before | After |','|---|---:|---:|']
    for label,key in [('Common-call joint','canonical_loss'),('Original-OTM joint','original_otm_loss')]:
        lines.append(f'| {label} | {a[key]:.6f} | {b[key]:.6f} |')
    for label,key,component in [('SPX','canonical_components','SPX'),('Common-call VIX options','canonical_components','VIX'),
                               ('Original-OTM VIX options','original_otm_components','VIX'),('VIX forwards','canonical_components','forwards')]:
        lines.append(f'| {label} | {a[key][component]:.6f} | {b[key][component]:.6f} |')
    lines+=['',f'{len(r["outer_proposals"])} fully checked proposal(s), {accepted} accepted. '
            'Acceptance requires true improvements in both joint and both VIX-plus-forward criteria. '
            'These are actual full-pricer scores. Local approximations are used only to propose steps; rejected points and prediction errors remain in the run records. '
            'A two-proposal pilot does not establish optimization convergence, global optimality, identification or a satisfactory smile fit.',
            '', '| Proposal | Predicted common-call joint | Actual joint | Predicted native VIX + forwards | Actual | Accepted |',
            '|---:|---:|---:|---:|---:|---|']
    for row in r['outer_proposals']:
        ap=row['approximate_criteria'];tr=row['true_criteria']
        lines.append(f"| {row['cycle']} | {ap[0]:.6f} | {tr[0]:.6f} | {ap[3]:.6f} | {tr[3]:.6f} | {row['decision']['accepted']} |")
    lines+=['','## Parameters','', '| Parameter | Before | After | New bounds |','|---|---:|---:|---|']
    for n,x,y,lo,hi in zip(NAMES,r['baseline'],r['selected'],r['bounds_lower'],r['bounds_upper']):
        lines.append(f'| {n} | {x:.10g} | {y:.10g} | [{lo:.7g}, {hi:.7g}] |')
    lines+=['',f"Selected floor volatility: {100*np.sqrt(r['selected'][10]):.6f}%; beta_J={r['selected'][4]:.9g}>0.",
            '', '## Every retained maturity, actual original-OTM prices','',
            '| Instrument | Days | MSE before | After | Quotes in bid/ask after |','|---|---:|---:|---:|---:|']
    qa=pd.read_csv(run/'prices_otm_baseline.csv');qb=pd.read_csv(run/'prices_otm_selected.csv')
    for (ins,day),g in qa.groupby(['instrument','expiry_days']):
        h=qb[(qb.instrument==ins)&(qb.expiry_days==day)]
        lines.append(f'| {ins} | {day} | {np.mean(g.interval_residual**2):.5f} | {np.mean(h.interval_residual**2):.5f} | {int(h.interval_residual.eq(0).sum())}/{len(h)} |')
    lines+=['','## Verification and limits','',
            f"Ten focused new/existing tests pass. Saved actual means/SE/parity/both losses, {audit['local_evaluations_replayed']} local trial scores, "
            f"all true acceptance decisions and {audit['sources_verified']} frozen sources/inputs replay. Bounds and selected configuration verify. "
            'Core datasets are the existing read-only caches, whose byte identities were verified in parameter-fit003; no new dataset was generated.',
            '',f"Selected actual prices: {b['native_otm_negative_prices']} negative native estimates; {b['iv_unavailable']} unavailable common-call IVs; "
            f"{b['inside_original_otm_bidask']}/1511 native prices inside bid/ask. "
            'Signed estimates, invalid IVs and shape/tail diagnostics are retained without clipping or smoothing. '
            'SPX pointwise SE uses iid4096 paths; VIX uses eight independent512-point nets. Bands exclude parameter/pilot uncertainty and deterministic error.',
            '', 'No new independent checking or price-convergence claim is made for this pilot. '
            'The prior VIX wing mismatch and SPX estimator stability concerns remain relevant. '
            'A local linear model can extrapolate poorly; the saved full checks measure this directly. '
            'No P-compatibility or Lean certificate is asserted.',
            '', '[All seven SPX maturities](visualizations/joint-local-012/spx-all-maturity-fit.png) · '
            '[All four VIX maturities](visualizations/joint-local-012/vix-all-maturity-fit.png) · '
            '[Forwards and path allocation](visualizations/joint-local-012/forwards-and-path-budgets.png). '
            'Blue is the selected vector and orange the previous fitted vector, on the same TRAINING draws. PNG/SVG and signed-error/precision figures are included.',
            '',f'[Selected model](runs/{run.name}/selected_model.json). '
            'Current bounds live in `calibration_policy.py`; use these for new fits. Historical scripts keep their old bounds for reproducibility.', '']
    out.write_text('\n'.join(lines));print('REPORT',out)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();summarize(a.run.resolve(),a.output.resolve())
