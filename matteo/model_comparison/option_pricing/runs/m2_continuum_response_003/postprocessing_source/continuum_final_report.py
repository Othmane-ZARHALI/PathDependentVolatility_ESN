from pathlib import Path
import json,hashlib,shutil
import pandas as pd
from scipy.special import betainc
b=Path('/Users/matteogambara/Desktop/code/PathDependentVolatility_ESN/matteo/model_comparison/option_pricing');r=b/'runs/m2_continuum_response_003';s=json.loads((r/'fresh_check/scores.json').read_text());tr={z['name']:z for z in json.loads((r/'training/scores.json').read_text())};sel=json.loads((r/'selection.json').read_text());audit=json.loads((r/'audit/audit.json').read_text());ma=json.loads((r/'multilevel_audit/audit.json').read_text());assert audit['status']==ma['status']=='PASS';ledger=json.loads((r/'optimizer_ledger.json').read_text());best=next(z for z in s if z['name']=='continuous');side=next((z for z in s if z['name']=='interior'),None);boundary=next(z for z in s if z['name']=='boundary');pairs=json.loads((r/'fresh_check/paired_losses.json').read_text())
text='''# Continuous M2/stress measure — completed VIX-only test, 27 October 2026

The continuous construction is implemented and tested with four VIX forward equality constraints. In the bounded searches it has not improved the two-component reference: the best Beta measure concentrates near the endpoints, while a broader continuous measure fits less well and reaches the stress-strength bound. Keep the previous working VIX fit; this experiment is not adopted. SPX remains paused.

## Model

For u in [0,1],

    R_t(u)=G_t^(1-u)/E_Q[G_t^(1-u)]
           * exp(b*u*Z_t - .5*(b*u)^2*q(t)),
    v_t=f+(xi(t)-f)*integral R_t(u) mu(du),
    mu=Beta(m*kappa,(1-m)*kappa), q(t)=1-exp(-2*a*t).

G is the unchanged capped M2 readout and Z is one independent OU. Components share the original reservoir and Z; no bank per node. Every power has its own normalizer and conditional expectation. The integral is formed before the VIX square root and option payoff. mu is a deterministic probability measure under Q, whose mass describes expected non-floor variance contributions, not regime probabilities.

All217 native VIX quotes and all four forwards (27,55,118,153days). Seven independent active xi levels of the16-piece curve; nine inactive levels unchanged. Floor volatility5–12%; b0–5; m.01–.99; kappa.01–30; half-life1–365days. Original readout/cap, positive beta_J, architecture and premiums frozen. No clipping, quote removal, curve smoothing, jumps, SPX repricing or baseline replacement.

## Forward condition and estimator repair

The optimizer minimizes native bid/ask scaled option-price distance subject to four forward equalities, not a forward penalty. Equality is enforced on the calibration estimator to numerical tolerance; a finite simulation does not certify exact population means.

The initial experiment002 used8192 powered paths/expiry and had a55-day checking gap-.1656 versus.0248SE. The corresponding two-atom refit also failed. Its training forward estimates were noisy, partly because legacy frozen controls were poorly matched to the new response. That calibration is rejected as reliable fulfillment of the forward condition; its outputs and independent model audit are preserved. See CONTINUUM_INITIAL_RESULTS.md and run002/forward_estimator_diagnostics.json.

Run003 keeps the exact model and uses the multilevel identity

    E[P_cont]=E[P_boundary]+E[P_cont-P_boundary].

The reference boundary uses the same b,m,half-life,floor and curve. Calibration estimates its cheap first term on32768 existing paths, plus the paired difference on8192 powered paths. Matching eight-net blocks retain the covariance of the nested samples. All model prices still come from the continuous variance law; subtracting the numerical reference is variance reduction, not a second model component or a mixture of finished option laws.

Checking uses65536 NEW paths/expiry for the boundary term (eight8192-point scrambles, seed20261027117), plus the existing16384 paired-difference checking paths. The difference stream is deliberately reused after002; this is not a wholly independent fresh study. Models and selection are frozen before generating the new large reference stream. Checking uncertainty and the reuse limitation are explicit.

## Fit comparison

Smaller price loss is better. Total comparison scores also include the historical forward-error penalty, retained for comparability; this penalty is not the new optimization objective.

| Model | Training option loss | Checking option loss | Checking total score | Within bid/ask | Max checking forward gap |
|---|---:|---:|---:|---:|---:|
'''
for z in s:text+=f"| {z['name']} | {tr[z['name']]['option_loss']:.8f} | {z['option_loss']:.8f} | {z['loss']:.8f} | {z['inside_bidask']}/217 | {z['maximum_forward_error']:.8f} |\n"
text+='\nPrevious is the last fitted-floor model repriced on the same checking streams. Boundary is refitted using the same large calibration estimator. Continuous is selected by forward-feasible training loss; interior is the best completed feasible start with kappa>=.1, using training only. All starts and local minima are retained; no global optimum is claimed.\n\n'
text+='## Four forward checks\n\n| Expiry | Observed | Continuous checking estimate | Gap | Eight-net SE |\n|---|---:|---:|---:|---:|\n'
f=pd.read_csv(r/'fresh_check/forwards_continuous.csv')
for z in f.itertuples():text+=f'| {z.expiry_days} | {z.forward_vix_points:.8f} | {z.model_forward:.8f} | {z.forward_error:+.8f} | {z.mc_se:.8f} |\n'
text+='\nCalibration gaps are in optimizer_ledger.json; calibration SEs are reported in training/forwards_*.csv. Numerical equalities should be read with both calibration and checking uncertainty. No forward shifting is applied to plotted prices. Canonical call IVs use the market forward for every curve; native IVs use each model forward and preserve unavailable inversions.\n\n'
text+='## Fitted continuous measures\n\n'
for z in [best]+([side] if side else []):
 m=z['model'];mass=betainc(m['alpha'],m['beta'],.9)-betainc(m['alpha'],m['beta'],.1);text+=f"{z['name']}: m={m['measure_mean']:.8f}, kappa={m['kappa']:.8f}, b={m['b']:.8f}, half-life={m['half_life_days']:.6f}days; mass in u=[.1,.9] is {100*mass:.3f}%. Floor={100*m['floor_volatility']:.6f}%.\n\n"
text+='The main kappa lower-bound hit approximates the previous two-atom construction. The broader measure remains genuinely spread, but its b=5 ceiling and12% floor hit limit conclusions about its best possible fit. Forward-variance bounds are also visible in curve_levels.csv. These results do not exclude richer continuous measures, a different interpolation path or wider stress-strength bounds.\n\n'
text+='## Paired checking differences\n\n```json\n'+json.dumps(pairs,indent=2)+'\n```\n\nEight-net paired jackknife errors measure conditional numerical checking variability, not fitting variability, global optimality or market out-of-sample performance.\n\n'
text+='## Verification and unresolved limits\n\n- Original M2 and pure-stress endpoints, Beta moments, independent analytic Gaussian power oracle,12-direction continuum and12-direction multilevel price/forward Jacobians passed.\n- Independent NumPy smooth-cap/power oracle, off-grid power checks and combined u33/measure24/GH32/time8 versus u65/measure48/GH48/time16 probes passed.\n- Independent direct model-VIX replay is scale-aware (absolute1e-9 OR relative1e-11): rare tiny-weight proposal paths have very large VIX values; absolute, relative and weighted-mean errors are retained. No unweighted maximum is used as tail evidence.\n- Independent reconstruction of base-plus-paired-difference prices/forwards, all eight blocks and both candidate models passed.\n- All32 new large reference scrambles and32 reused difference scrambles regenerated exactly. Source/input hashes, floor support and all nine inactive xi levels checked.\n- Signed native estimates and unavailable IVs retained. Separate raw positive-payoff put diagnostics and negative-price/SE diagnostics are saved, without substituting them into calibration.\n- Second-moment targets are exact model identities; numerical moment estimates and SEs appear in multilevel_audit/second_moments.json. Full bid/ask feasibility is not established.\n- Sampled deterministic quadrature checks are not uniform error bounds. No new Lean theorem or exact population-forward certificate.\n\n'
text+='Combined sampled refinement:\n\n```json\n'+json.dumps(audit['refinements'],indent=2)+'\n```\n\n'
text+='## Reproducible artifacts\n\n- Main run: runs/m2_continuum_response_003. Initial rejected calibration: run002; initial overly strict Beta-moment assertion: run001, preserved.\n- Mathematical and numerical core: continuum_response_core.py, continuum_power_kernel.cpp, continuum_multilevel_core.py, continuum_multilevel_study.py.\n- Independent audits: audit_continuum_scaled_fast.py, audit_continuum_response_v2.py, audit_multilevel_continuum.py.\n- Experimental models: experimental_model_continuous.json and experimental_model_interior.json.\n- Plots (PNG/SVG): visualizations/continuum-response023, including all strikes, right wing, native puts, forward errors/curve, and measure comparisons.\n- Frozen calibration sources and postprocessing sources are saved; artifact_manifest.json records hashes.\n\nThe original selected production baseline and previous working VIX candidate are preserved.\n'
(b/'CONTINUUM_RESPONSE_RESULTS.md').write_text(text)
checkpoint='''## 27 October 2026 — continuous measure, forward constraints and multilevel repair

Completed VIX-only continuous M2/stress Beta-measure experiment003, after rejecting002's noisy55-day forward calibration. The model is unchanged: normalized G^(1-u) times the normalized OU stress exponential of strength b*u, integrated against one Beta measure. Four forward equalities use32768-path boundary plus8192-path paired difference; checking uses65536 NEW boundary paths plus16384 reused difference paths. Do not claim wholly fresh validation or exact population equality.

Best bounded fit remains near kappa=.01; a broader measure is retained as a training-selected diagnostic, with b=5 at its bound. No clear improvement established over the two-atom reference. All217 quotes, four maturities, seven free active xi levels and fitted floor preserved. SPX paused; original baseline and previous working VIX candidate unchanged. See CONTINUUM_RESPONSE_RESULTS.md, runs/m2_continuum_response_003 and visualizations/continuum-response023.

'''
for name in ['NEXT_STEPS.md','CONTINUE_FITTING_PROMPT.md','RESEARCH_LEDGER.md','STOCHASTIC_VOL_OF_VOL_RESEARCH.md','AUDIT.md']:
 p=b/name
 if p.exists():p.write_text(checkpoint+p.read_text())
p=b/'active_scope.json';a=json.loads(p.read_text());a['decision']='Continuous M2/stress measure tested with four multilevel-estimator forward equalities; no clear fit improvement established. Previous working VIX fit retained, SPX paused, original baseline unchanged.';a['continuous_measure_experiment']=dict(run='runs/m2_continuum_response_003',primary_model='runs/m2_continuum_response_003/experimental_model_continuous.json',secondary_model='runs/m2_continuum_response_003/experimental_model_interior.json',status='completed_exploratory_not_adopted',report='CONTINUUM_RESPONSE_RESULTS.md',forward_treatment='Four exact calibration-estimator equalities; checking gaps and SE reported; reused difference stream disclosed',figure_folder='visualizations/continuum-response023');p.write_text(json.dumps(a,indent=2)+'\n')
post=r/'postprocessing_source';post.mkdir(exist_ok=True)
for path in [*b.glob('*continuum*.py'),b/'plot_multilevel_response.py',b/'plot_multilevel_comparison.py',b/'continuum_power_kernel.cpp',b/'CONTINUUM_RESPONSE_PLAN.md',b/'CONTINUUM_MULTILEVEL_PLAN.md',b/'CONTINUUM_RESPONSE_MATH.md',b/'CONTINUUM_RESPONSE_RESULTS.md',b/'CONTINUUM_INITIAL_RESULTS.md',Path('/private/tmp/continuum_final_report.py')]:shutil.copy2(path,post/path.name)
d=json.loads((r/'run.json').read_text());d['status']='completed with model-VIX and multilevel independent audits PASS';d['audit']='audit/audit.json';d['multilevel_audit']='multilevel_audit/audit.json';d['adopted']=False;(r/'run.json').write_text(json.dumps(d,indent=2)+'\n')
manifest={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(r.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'};figs={str(p.relative_to(b)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((b/'visualizations/continuum-response023').glob('*')) if p.is_file()};(r/'artifact_manifest.json').write_text(json.dumps(dict(sha256=manifest,figure_sha256=figs,count=len(manifest),figure_count=len(figs)),indent=2)+'\n');print('Finalized',len(manifest),'run artifacts and',len(figs),'figures')
