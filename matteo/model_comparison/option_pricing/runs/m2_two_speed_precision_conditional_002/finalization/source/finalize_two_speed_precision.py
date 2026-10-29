"""Archive the completed protected pricing phase and update continuation notes."""
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import HERE,dump


def run():
    runs=HERE/'runs';main=runs/'m2_two_speed_precision_002';conditional=runs/'m2_two_speed_precision_conditional_002';pilot=runs/'m2_two_speed_precision_conditional_001';plots=HERE/'visualizations/two-speed-precision027';review=HERE/'visualizations/two-speed-precision027-review'
    if review.exists():
        for f in review.iterdir():shutil.copy2(f,plots/f.name)
        shutil.rmtree(review)
    for d in [main,conditional,pilot]:assert json.loads((d/'audit.json').read_text())['status'].startswith('implementation_checks_PASS')
    frame=pd.read_csv(plots/'precision.csv');q=pd.read_csv(conditional/'check/prices_conditional.csv');primary=q[q.priority_selected];assert len(q)==1751 and len(primary)==1294
    ratios=frame.old/frame.conditional
    table='| Days | Previous median SE / scale | New median SE / scale | Reduction | Expensive component | Cheap component |\n|---:|---:|---:|---:|---:|---:|\n'
    for r in frame.itertuples():table+=f'| {r.day} | {r.old:.3f} | {r.conditional:.3f} | {r.old/r.conditional:.1f}x | {r.residual:.3f} | {r.cheap:.3f} |\n'
    budget=json.loads((plots/'sampling-budget-projection.json').read_text());bt='| Days | Projected median expensive paths | Projected worst-quote expensive paths |\n|---:|---:|---:|\n'
    for r in budget['per_maturity']:bt+=f"| {r['day']} | {r['median_expensive_paths']:,.0f} | {r['max_expensive_paths']:,} |\n"
    text=table+f'\nAll figures are checking-role results for the same frozen VIX-only model. The estimated median SE reduction ranges from{ratios.min():.1f}x to{ratios.max():.1f}x across maturities. That is a reduction in SE, not in variance or model error. Only{int(primary.SE_gate_pass.sum())}/1294primary quotes passes the0.125SE/scale screen. Negative native checking estimates:{int((q.model_price_U<0).sum())}/1751. No option optimizer was launched.\n\n'
    text+='[Precision comparison](visualizations/two-speed-precision027/precision-comparison.png) and [all seven native SPX smiles](visualizations/two-speed-precision027/frozen-model-spx-smiles.png). The corrected IV plot uses percentage points directly from the existing inversion routine; it does not multiply percentages by100 again. Its CSV retains every primary quote and flags invalid model IVs.\n\n'
    text+='The clearer SPX plots still show a substantially more U-shaped model smile than the observed downward skew. This is a fixed-parameter discrepancy; it is not a proof that refitting the nine allowed curve entries cannot help, or that this two-speed model family is infeasible. The small plotted sampling bars do not explain the broad current smile mismatch.\n\n'
    text+='The expensive correction is the main remaining sampling bottleneck. Taking only the cheap reference to infinite precision would leave the median expensive-component SE above the screen at every maturity. A simple equal-variance-budget projection gives:\n\n'+bt
    text+='\nThese are pilot-based sqrt(N) projections using the same expensive estimator/proposal, not guaranteed requirements or achieved precision. Each component is assigned SE0.125/sqrt(2). For the cheap reference, the companion CSV projects more independently scrambled4096-point nets; it does not assume a convergence rate when enlarging one Sobol net. The worst-quote expensive projection is about1.2million paths, versus512 currently stored. This makes blind expansion of the full fine-path cache unattractive: use streaming, adaptive allocation and a coupled multilevel estimator before committing to that cost. The projections exclude timestep and normalizer bias and can themselves be unstable for concentrated tails.\n\n'
    text+='Primary native interval loss (equal maturity weight) of the frozen model with the new estimator: '+f"{np.mean([np.mean(g.interval_residual**2) for _,g in primary.groupby('expiry_days')]):.3f}. This is a diagnostic of the unchanged model, not a new calibration result.\n\n"
    report=HERE/'TWO_SPEED_PRECISION_RESULTS.md';s=report.read_text();assert '<!-- GENERATED_RESULTS -->' in s;report.write_text(s.replace('<!-- GENERATED_RESULTS -->',text))
    checkpoint='''## 28 October 2026 — protected two-speed SPX precision phase

User authorized improving SPX precision with the successful two-speed VIX-only model fixed, then a protected fit of only nine off-window xi levels, with stress/equity correlation a later possible experiment. All model parameters remain unchanged in this phase: no joint candidate adoption, correlation, factor, readout or floor change. All1751eligible SPX quotes repriced,1294primary at all seven retained maturities; all217VIX quotes/four forwards remain protected. See [TWO_SPEED_PRECISION_RESULTS.md](TWO_SPEED_PRECISION_RESULTS.md), [TWO_SPEED_PRECISION_PLAN.md](TWO_SPEED_PRECISION_PLAN.md), `runs/m2_two_speed_precision_002`, `runs/m2_two_speed_precision_conditional_002`, and `visualizations/two-speed-precision027`.

Implemented an independent-stress auxiliary stock, likelihood-weighted differences, terminal/dynamic hedge controls, then an exact discrete one-step conditional Black-value martingale correction and independently scrambled OU-bridge reference integration. Checking median SE/scale improves about8–29x, from4.68/15.03/12.17/12.76/13.08/5.74/5.01 to0.579/0.550/0.903/1.153/1.239/0.435/0.174. Only1/1294primary quotes passes the existing sampling screen. The nine-level price optimizer was therefore NOT launched. Broad frozen-model SPX smile-shape mismatch remains visible; it is not a model-family infeasibility theorem.

Remaining bottleneck is the512-path full-M2 correction. More cheap stress paths alone cannot clear the gate. Pilot equal-component budget projections give median expensive counts about1.8k–93k across maturities and worst quote about1.2million, excluding timestep/normalizer errors; these are not guaranteed budgets. Next implement streaming fresh full-M2 corrections with adaptive/multilevel allocation and coupled stock-grid checks, rather than expanding the full cache blindly. Then free ONLY indices0,5,6,9,11,12,13,14,15; keep floor, both stresses and VIX-sensitive indices1,2,3,4,7,8,10 fixed. The exact VIX parameter map is invariant to those nine changes. Existing forward accuracy is preserved, not upgraded to population equality. Independent target/reference/conditional/RNG/covariance/full-horizon derivative and finite-law control checks PASS. Checking streams reused; not market OOS, full convergence clearance or a Lean proof. Current chat can continue; these notes are the restart handoff if a new context is chosen.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        f=HERE/name;f.write_text(checkpoint+f.read_text())
    dump(runs/'m2_two_speed_precision_001/run.json',dict(status='preflight stopped before price experiment',reason='Identical-model test fixture combined fine target stress and held reference stress; corrected by using identical fine reference grid in that limit.',successor='m2_two_speed_precision_002',model_changed=False))
    final=conditional/'finalization';final.mkdir(exist_ok=False);(final/'source').mkdir()
    for name in ['two_speed_precision.py','two_speed_precision_kernel.cpp','audit_two_speed_precision.py','two_speed_precision_conditional.py','two_speed_precision_conditional.cpp','audit_two_speed_precision_conditional.py','diagnose_two_speed_precision_mean_control.py','plot_two_speed_precision.py','finalize_two_speed_precision.py','TWO_SPEED_PRECISION_RESULTS.md','TWO_SPEED_PRECISION_PLAN.md','NEXT_STEPS.md','CONTINUE_FITTING_PROMPT.md','RESEARCH_LEDGER.md']:shutil.copy2(HERE/name,final/'source'/name)
    (final/'logs').mkdir()
    for f in Path('/private/tmp').glob('m2-two-speed-precision-*.log'):shutil.copy2(f,final/'logs'/f.name)
    artifacts={}
    for folder in [main,pilot,conditional,plots]:
        for f in sorted(folder.rglob('*')):
            if f.is_file() and f!=final/'artifact_manifest.json':artifacts[str(f.relative_to(HERE))]=hashlib.sha256(f.read_bytes()).hexdigest()
    dump(final/'artifact_manifest.json',artifacts)
    print(json.dumps(dict(archived_artifacts=len(artifacts),SE_reduction_range=[float(ratios.min()),float(ratios.max())],primary_sampling_passes=int(primary.SE_gate_pass.sum()),quotes=len(primary),VIX_model_changed=False,report=str(report)),indent=2))


if __name__=='__main__':run()
