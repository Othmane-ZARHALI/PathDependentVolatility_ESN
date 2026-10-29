"""Save the protected-fit empirical conclusion and auditable restart state."""
import hashlib
import json
import platform
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from .full_readout_fit import HERE,dump
from .two_speed_precision import FREE


def spaced_prose(text):
    """Separate prose numbers while preserving math, identifiers and links."""
    protected=r'(`[^`]*`|\[[^\]]*\]\([^)]+\)|\$\$[\s\S]*?\$\$|\$[^\n$]*\$|\bM2\b|\bGH48\b)'
    chunks=re.split(protected,text)
    for j in range(0,len(chunks),2):
        chunks[j]=re.sub(r'(?<=[A-Za-z])(?=\d)|(?<=\d)(?=[A-Za-z])',' ',chunks[j])
    return ''.join(chunks)


def finalize(run):
    plots=HERE/'visualizations/protected-curve028';fresh=run/'fresh_check';audit=json.loads((run/'audit.json').read_text());assert audit['status'].startswith('implementation_checks_PASS');scores=json.loads((fresh/'scores.json').read_text());metrics=pd.read_csv(plots/'SPX_metrics.csv');curve=pd.read_csv(plots/'curve_levels.csv');fitted=json.loads((run/'run.json').read_text());base=next(r for r in scores['summaries'] if r['model']=='baseline' and r['steps']==128);selected=next(r for r in scores['summaries'] if r['model']=='selected' and r['steps']==128);biv=float(metrics[metrics.model=='baseline'].IV_MAE_pp.mean());siv=float(metrics[metrics.model=='selected'].IV_MAE_pp.mean());gain=100*(1-selected['loss']/base['loss']);pair=next(r for r in scores['paired'] if r['steps']==128);interval=pair['paired_batch_bootstrap_loss_delta_interval'];grid=json.loads((plots/'grid_summary.json').read_text());inputs={}
    # Audit the actual training primitives and binary controls against immutable
    # manifests that predate this experiment; do not silently refresh a hash.
    for root,filter_fn in [(HERE/'runs/m2_two_speed_joint_cache_001',lambda n:n.startswith('search/')),(HERE/'runs/m2_two_speed_precision_conditional_002',lambda n:n.startswith('search/') and n.endswith('_reference.npz') or n=='conditional.dylib'),(HERE/'runs/m2_two_speed_precision_002',lambda n:n=='precision.dylib')]:
        manifest=json.loads((root/'artifact_manifest.json').read_text());manifest=manifest.get('files',manifest)
        for name,expected in manifest.items():
            if not filter_fn(name):continue
            f=root/name;actual=hashlib.sha256(f.read_bytes()).hexdigest();assert actual==expected,str(f);inputs[str(f)]=actual
    dump(run/'input_integrity_audit.json',dict(status='PASS',verified_inputs=len(inputs),input_hashes=inputs,meaning='Training inputs and reused control binaries checked against pre-existing immutable manifests. These are integrity checks, not independent price estimates.'))
    final=run/'finalization';final.mkdir(exist_ok=True);(final/'source').mkdir(exist_ok=True);(final/'logs').mkdir(exist_ok=True)
    table='| Days | Baseline price loss | Protected price loss | Baseline IV MAE (pp) | Protected IV MAE (pp) |\n|---:|---:|---:|---:|---:|\n'
    for day in sorted(metrics.day.unique()):
        a=metrics[(metrics.day==day)&(metrics.model=='baseline')].iloc[0];b=metrics[(metrics.day==day)&(metrics.model=='selected')].iloc[0];table+=f'| {day} | {a.loss:.2f} | {b.loss:.2f} | {a.IV_MAE_pp:.3f} | {b.IV_MAE_pp:.3f} |\n'
    levels='| Index | Days | Baseline xi | Protected xi |\n|---:|---:|---:|---:|\n'
    for r in curve[curve.free].itertuples():levels+=f'| {r.index} | {r.start_days:g}–{r.end_days:g} | {r.baseline_xi:.8f} | {r.candidate_xi:.8f} |\n'
    report=rf'''# Protected nine-level SPX fit: results and durable conclusion

28 October 2026. The user authorized this experiment and explicitly asked to save its conclusion so it is not rediscovered in a later conversation. The complete run is `runs/m2_protected_curve_001`. See [the experiment plan](PROTECTED_CURVE_PLAN.md).

**The nine-level fit improves the specified SPX price objective without sacrificing VIX, but does not solve SPX smile shape. It is not a satisfactory joint calibration.** On fresh128-step checking, price loss falls by{gain:.1f}%, while equal-maturity average native implied-volatility error increases from{biv:.3f} to{siv:.3f}percentage points. Both facts must be retained when describing the result.

## Exact experimental scope

The existing two-speed variance model remains

$$v_t=f+(\xi_t-f)[(1-w)R_t+w\{{(1-r)A_t^s+rA_t^f\}}].$$

The original M2 dynamics, readout, normalizer, positive spike shift, risk premia, full initial state, both independent stress factors, stress weights and sensitivities, half-lives and floor are fixed at the successful VIX-only values. In particular, w=0.922838 and the volatility floor remains12%. Only xi indices0,5,6,9,11,12,13,14,15 are free, with the original floor-to0.2variance bounds. All seven VIX-sensitive xi levels stay fixed. No smoothing, new factor, correlation, jump or maturity-specific rule was added.

All1294primary SPX quotes at29,57,119,211,302,666,1030days enter the unchanged native bid–ask-scaled interval objective, with equal maturity weights. All1751eligible quotes are checked and retained in CSVs. All217VIX native quotes and four forwards remain protected. Signed estimates and undefined IVs are retained; there is no clipping or smoothing.

## Why VIX stays unchanged

VIX at expiry T depends on conditional variance over[T,T+30days]. The free curve intervals lie outside all four VIX windows. They do not enter the exact VIX parameter map, and none of the exogenous variance-factor dynamics depends on the stock or on xi. Therefore these changes leave the VIX distribution itself unchanged. This is stronger than merely preserving four forward means.

The audit checks exact array equality of baseline/candidate VIX inputs and every other fixed coordinate. `VIX_protection/invariance.json` records this and hashes inherited original price/forward artifacts. The original checking VIX option loss remains0.050081, with153/217quotes inside bid–ask. The previously measured checking forward errors remain as previously measured: no new population-forward equality certificate is claimed, and inherited VIX CSVs are not presented as new Monte Carlo estimates.

## Fitting and genuinely fresh checking

Fitting uses512existing full-M2 search paths and4096cheap reference paths, with the terminal-Black-delta auxiliary estimator from the precision study and analytic derivatives. This is an unbiased estimator of the existing discrete model, with sampling and numerical errors assessed separately. It is cheaper than the one-step conditional-value estimator; they must not be conflated. All nine native-price derivative directions pass finite differences. Two bounded starts were tested, each with at most40least-squares evaluations. Both terminated by ftol. Training loss changes from{fitted['baseline_loss']:.3f} to{fitted['selected_loss']:.3f}. This is a bounded local search, not a global fit or uniqueness result.

Selection was frozen before checking. Fresh checking generates512new target-Q paths (256antithetic pairs in eight batches), including new M2, equity and both stress innovations. Each path carries both candidate models and all seven maturities. The stock64/128grids are coupled through the same exact fine OU/stress/equity paths. All curve knots and maturity events are inserted exactly. Independent reference checking uses eight new scrambled OU-bridge nets of4096points per expiry, shared between candidates. The old expensive checking cache is not reused. This is independent simulation checking at the same market snapshot, not market out-of-sample testing.

The12orthogonal OU modes have zero C/Y/F readout projection and no feedback into retained dynamics, so they are analytically marginalized. The exact retained24-mode joint OU/equity Gaussian law is preserved and audited, including its Q drift and covariance. This is not a reduced three-projection Markov approximation. State paths are streamed in batches, not stored as an expanded full fine-state cache.

## Price loss versus implied-volatility error

At64steps/tradingday, fresh price loss changes5964.131→3895.393. At128steps it changes{base['loss']:.3f}→{selected['loss']:.3f}. The paired128-step loss difference is{pair['loss_delta_selected_minus_baseline']:.3f}; its approximate paired batch-bootstrap95%interval is[{interval[0]:.3f},{interval[1]:.3f}]. It excludes zero, supporting improvement of this price metric beyond the observed sampling noise. It excludes fitting variability, normalizer error and stock-grid bias; it is not a general price-error bound.

{table}

The price objective and IV MAE measure different things. An option's price sensitivity to implied volatility varies with strike and maturity. Reducing the chosen bid–ask-scaled price loss can therefore worsen an average IV error. The fit mostly lowers the SPX curves: it helps the price objective, but the left wing becomes too low and the right wing remains too high relative to the market. The curves retain a much more U-shaped profile than the observed downward skew.

At128steps only{selected['inside']}/1294primary quotes are inside bid–ask, versus{base['inside']}for the baseline. All fresh checking native prices are nonnegative in this sample; signed estimates were retained. Median candidate sampling SE is{selected['median_SE_scales']:.3f}price-scale units. Thus this is clearly not a final quote-level calibration despite the statistically supported loss improvement.

[All seven SPX smiles](visualizations/protected-curve028/SPX-protected-smiles.png), [the four inherited VIX smiles](visualizations/protected-curve028/VIX-protected-smiles.png), and [the variance-curve changes](visualizations/protected-curve028/protected-variance-curve.png). Model IVs use the existing percentage-valued inversion routine, and invalid values are flagged in the all-quote CSVs.

## What changed in the curve

{levels}

Several short or intermediate off-window levels approach the frozen floor. Those boundary/near-boundary values are numerical outcomes of this restricted fit, not uniquely identified economic estimates. The original model and this exploratory candidate remain separate.

## Numerical checks and limits

The implementation audit passes: all9price derivative directions; independent scalar/NumPy reconstruction of the fresh stock kernel (maximum error{audit['independent_fresh_kernel_numpy_replay_error']:.3g}); all8new batch-seed replays; the complete retained Gaussian covariance recurrence over91040fine steps (maximum absolute covariance error{audit['retained_OU_covariance_semigroup_max_error']:.3g}); Q projected means; source fingerprints; input-manifest checks; exact protected-parameter invariance; and positive spike shift. Normalization still uses the frozen GH48/interpolated normalizer, with its earlier sampled interpolation audits. No new global GH-convergence certificate is supplied.

The broad improvement survives the64→128grid change. Nevertheless the candidate's median absolute paired price change is{grid[1]['median_absolute_grid_change_scales']:.3f}price-scale units, maximum{grid[1]['max_absolute_grid_change_scales']:.3f}; finer-grid convergence is not established. Sampling SEs use antithetic pairs and independent reference-net means. Loss intervals resample the eight expensive and eight cheap batch means separately, with common resamples across candidates. Eight batches give an approximate uncertainty assessment, not a guaranteed tail bound.

## Durable conclusion: do not repeat this experiment as if it were untested

1. Freeing off-window forward-variance levels can improve the original SPX price objective while preserving the successful VIX distribution.
2. In these two bounded starts it does not produce the required SPX downward skew. Price loss improves while average native-IV error worsens. Do not summarize the experiment only by the lower price loss.
3. Further curve-only iterations should have a new concrete hypothesis or new evidence. Repeating the same unrestricted or nine-level optimization is not the recommended next step.
4. A plausible explanation is that the dominant stress variance is independent of the stock shock. It supplies volatility variation without directly linking stress increases to stock falls. This is a mechanism hypothesis, not an established theorem or proof that the entire fixed-factor family is impossible to calibrate.
5. The next proposed structural test is one global stock/stress leverage parameter, with no added stress speed. A constructive driver rotation can preserve the joint variance-factor law and hence VIX; see [STOCK_DRIVER_ROTATION_PROPOSAL.md](STOCK_DRIVER_ROTATION_PROPOSAL.md). That proposal is not implemented here, changes the old stock/M2 correlations by construction, and requires new joint-transition and pricing-control checks.

The original VIX-only artifact and production baseline are unchanged. The candidate is `experimental_model_protected.json`, explicitly not adopted. No new Lean theorem or global impossibility claim is made. Current continuation notes link directly to this result so a future conversation can resume without repeating it.
'''
    (HERE/'PROTECTED_CURVE_RESULTS.md').write_text(spaced_prose(report))
    conclusion=f'''# Research decision: variance levels versus SPX skew

Saved28October2026 at the user's explicit request. Evidence: [full result](PROTECTED_CURVE_RESULTS.md), [plan](PROTECTED_CURVE_PLAN.md), run `runs/m2_protected_curve_001`, [SPX plot](visualizations/protected-curve028/SPX-protected-smiles.png).

**Remember both findings:** the protected nine-level fit lowers fresh SPX native price loss by{gain:.1f}% while leaving the VIX distribution unchanged, but equal-maturity IV MAE worsens{biv:.3f}→{siv:.3f}percentage points and the SPX downward skew remains missing. It is not a satisfactory joint calibration.

Only xi indices0,5,6,9,11,12,13,14,15 moved. Both stress factors, their weights/speeds/sensitivities, the floor, old M2 and seven VIX-sensitive levels stayed fixed. All seven SPX maturities/1294primary quotes were fitted; all1751eligible quotes checked. All217VIX quotes/four forwards remain protected. Two bounded starts; training selection before genuinely fresh512-path paired Q validation. The128-step paired loss-improvement interval excludes zero, but final sampling/grid/normalizer accuracy remains unresolved.

Do not spend another iteration rediscovering that lowering the variance curve can improve the price objective while worsening left-wing IVs. The next proposed mechanism is stock/stress leverage, not more stress speeds or another blind curve fit. [A driver-rotation construction](STOCK_DRIVER_ROTATION_PROPOSAL.md) can preserve the variance-factor joint law and VIX. It remains a proposal: no correlation parameter has been implemented or fitted.

This is a bounded empirical finding at one market snapshot, not proof that every possible curve, continuous stress measure or model in the family is incapable of joint calibration. New data or a clearly different hypothesis can justify revisiting it. Preserve price and IV metrics, numerical uncertainty and the successful VIX baseline when continuing.
'''
    (HERE/'PROTECTED_CURVE_CONCLUSION.md').write_text(spaced_prose(conclusion))
    checkpoint=f'''## 28 October 2026 — protected nine-level fit completed; durable conclusion saved

User explicitly asked to save this conclusion. Read [PROTECTED_CURVE_CONCLUSION.md](PROTECTED_CURVE_CONCLUSION.md) and [PROTECTED_CURVE_RESULTS.md](PROTECTED_CURVE_RESULTS.md) before repeating a curve-only experiment. Run `runs/m2_protected_curve_001`, figures `visualizations/protected-curve028`. Only nine off-window xi levels changed; full VIX law and original successful VIX artifact preserved. All1294primary/seven SPX maturities fitted, all1751eligible checked; all217VIX quotes/four forwards retained. Revised numerical policy explicitly authorized by user: all-quote precision screen is final acceptance, not a prerequisite for a bounded exploratory fit.

Training price loss6406.258→3665.656. Genuinely fresh512-path Q check,256antithetic pairs/eightbatches, coupled64/128stock grids and new reference scrambles:128-step price loss5979.323→3895.071({gain:.1f}% decrease), paired approximate95%loss-difference interval[{interval[0]:.1f},{interval[1]:.1f}]. BUT equal-maturity native-IV MAE{biv:.3f}→{siv:.3f}pp, inside bid–ask11→8/1294; downward skew still missing. This is not a satisfactory joint calibration. Full-M2/equity/stress checking draws are new, not a reused old-state cache; same market snapshot, not market OOS. Exact24-mode retained Gaussian marginal used;12irrelevant orthogonal modes marginalized, not a three-state Markov approximation. Independent fresh-kernel, Gaussian-law,9gradient,8seed, input/source and VIX-invariance checks PASS. Final sampling/grid/normalizer gates unresolved. Candidate saved separately, not adopted; no Lean/infeasibility theorem.

Next proposed structural experiment: one global stock/stress leverage parameter, no new stress speed. [STOCK_DRIVER_ROTATION_PROPOSAL.md](STOCK_DRIVER_ROTATION_PROPOSAL.md) derives a valid rotation preserving the joint variance-factor law and hence VIX; it also weakens old stock/M2 correlations and invalidates unchanged reuse of the current independent-stress pricing controls. Proposal only, not implemented or calibrated. Do not silently treat a lower price loss as improved IV shape or repeat the same curve-only search without a new hypothesis. Earlier precision note's optimizer prohibition was superseded by the later explicit user authorization recorded in PROTECTED_CURVE_PLAN.md.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        f=HERE/name;old=f.read_text()
        if old.startswith('## 28 October 2026 — protected nine-level fit completed; durable conclusion saved'):
            next_header=old.find('\n## ',1);assert next_header>=0;old=old[next_header+1:]
        f.write_text(spaced_prose(checkpoint)+old)
    proposal=HERE/'STOCK_DRIVER_ROTATION_PROPOSAL.md';proposal.write_text(spaced_prose(proposal.read_text()))
    fitted.update(status='protected exploratory fit and fresh paired validation completed; not satisfactory joint calibration',fresh_check_summary=scores,saved_conclusion='PROTECTED_CURVE_CONCLUSION.md',VIX_parameters_changed=False,production_model_changed=False,adopted=False);dump(run/'run.json',fitted)
    dump(final/'environment.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()))
    for name in ['protected_curve_fit.py','protected_curve_check.py','protected_curve_fresh.cpp','audit_protected_curve.py','plot_protected_curve.py','finalize_protected_curve.py','PROTECTED_CURVE_PLAN.md','PROTECTED_CURVE_RESULTS.md','PROTECTED_CURVE_CONCLUSION.md','STOCK_DRIVER_ROTATION_PROPOSAL.md','NEXT_STEPS.md','CONTINUE_FITTING_PROMPT.md','RESEARCH_LEDGER.md']:shutil.copy2(HERE/name,final/'source'/name)
    for f in Path('/private/tmp').glob('m2-protected-curve-*.log'):shutil.copy2(f,final/'logs'/f.name)
    manifest={}
    for folder in [run,plots]:
        for f in sorted(folder.rglob('*')):
            if f.is_file() and f!=run/'artifact_manifest.json':manifest[str(f.relative_to(HERE))]=hashlib.sha256(f.read_bytes()).hexdigest()
    dump(run/'artifact_manifest.json',dict(files=manifest,root=str(HERE),count=len(manifest)))
    print(json.dumps(dict(archived_artifacts=len(manifest),verified_prior_inputs=len(inputs),fresh_price_loss_reduction_percent=gain,IV_MAE_before=biv,IV_MAE_after=siv,conclusion=str(HERE/'PROTECTED_CURVE_CONCLUSION.md')),indent=2))


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();finalize(args.run.resolve())
