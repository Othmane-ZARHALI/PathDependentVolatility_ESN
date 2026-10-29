"""Preserve numerical evidence and durable conclusions without adopting the pilot."""
import hashlib,json,shutil
from pathlib import Path
import pandas as pd
from .slow_readout import HERE,OUT,dump
WEAK=HERE/'runs/m2_slow_readout_002'

def main():
    plot=HERE/'visualizations/slow-readout028';summary=json.loads((plot/'summary.json').read_text());s=summary['scores'];v=summary['VIX_scores'];m=pd.read_csv(plot/'maturity_metrics.csv').query('day==29').set_index('candidate');note=f'''## 28 October 2026 — slower existing-bank downside signals tested

Read [SLOW_READOUT_RESULTS.md](SLOW_READOUT_RESULTS.md); runs `runs/m2_slow_readout_001` and `runs/m2_slow_readout_002`, figures `visualizations/slow-readout028`. The previously proposed slow-readout hypothesis is now tested. Existing feedback modes nearest log-half-life targets2/8/32days have actual half-lives3.329/9.481/27.004days. Original fast M2 function and positive spike shift retained; no new Gaussian bank or endogenous return/volatility feedback. Stock theta=.8 fixed. One slow centered unit-stationary-SD mode at a time enters D=exp(alpha softplus(-X)). Conditional-stress normalization H=(aR+wA)D(a+wA)/(a E[RD]+wA E[D]) preserves E[H|stress]=a+wA, not the VIX law. R and D are correlated; product moments use their full joint Gaussian law.

Exploratory alpha.5/1 then reactive.125/.25 grids, fixed curve and same seven active curve levels refitted under four training-forward equalities; nine off-window levels frozen.512paired pilot stock paths and512training VIX states/expiry; previously examined4096checking states. Frozen weak selection9.481days/alpha.125 before1024new stock paths and2048new VIX states/expiry. Five fresh controls include weak fixed curve, weak refit, curve-only and strong3.329day/alpha.5refit. All1751eligible/1294primarySPX/sevenexpiries and217VIX/fourforwards retained, signed prices/no smoothing/clipping.

Fine256step results: baseline SPX loss{s[0]['loss']:.2f}/IVMAE{s[0]['IV_MAE_pp']:.3f}pp/VIXloss{v[0]['option_loss']:.3f}; weak fixed{s[1]['loss']:.2f}/{s[1]['IV_MAE_pp']:.3f}/{v[1]['option_loss']:.3f}; weak refit{s[2]['loss']:.2f}/{s[2]['IV_MAE_pp']:.3f}/{v[2]['option_loss']:.3f}; curve alone{s[3]['loss']:.2f}/{s[3]['IV_MAE_pp']:.3f}/{v[3]['option_loss']:.3f}; strong refit{s[4]['loss']:.2f}/{s[4]['IV_MAE_pp']:.3f}/{v[4]['option_loss']:.3f}. Weak refit improves SPX price loss{summary['price_loss_reduction_percent']:.2f}%; approximate paired eight-batch bootstrap loss-change interval[-95.71,-7.82], excluding deterministic error. Strong contrast interval includes0.29day right-minus-left IV gap baseline{m.loc[0,'model_right_minus_left_IV_pp']:.2f},weak{m.loc[2,'model_right_minus_left_IV_pp']:.2f},strong{m.loc[4,'model_right_minus_left_IV_pp']:.2f}pp versusmarket{m.loc[0,'market_right_minus_left_IV_pp']:.2f}. Missing skew persists. Weak refit checking forward gap max{v[2]['max_forward_gap']:.4f}points: training equality is not exact population matching.

Alpha-zero old-kernel replay2.22e-16; independent conditional tensor product-moment checks<6e-7relative; GH48/64normalizers<6.66e-7; independent97step scalar five-candidate stock recurrence0error. Coupled128/256stock median changes.117/.135/.128/.118/.183quote scales; matched64state/expiry VIX GH32/time8 versus48/time12 max7.32e-6points. Numerical probes only, no Lean certification/uniform tail bound. Stock/VIX likelihood-weighted payoff concentration, SEs, stock means and invalid IVs retained.

Durable conclusion: this simple slow one-sided multiplier provides a small SPX response but no satisfactory double fit. Stronger changes damage VIX despite conditional-mean normalization. Do not repeat this same strength scan or describe the slow-readout hypothesis as untested. This is not a no-go theorem for other slower projections/readout structures. Production model unchanged; no candidate adopted. Keep fast Zumbach-motivated component in any subsequent design and rebuild downstream prices whenever a readout changes.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        path=HERE/name;path.write_text(note+path.read_text())
    path=HERE/'formalization/CORRESPONDENCE.md';old=path.read_text();addition='''## 28 October 2026: slower existing-bank downside multiplier

`slow_readout.cpp` and `slow_readout.py` retain the frozen exogenous Gaussian bank and fast M2 function while adding a selected existing slow OU mode to a one-sided readout. The conditional identity E[Hnew|stress]=a+wA follows by integrating correlated R,D jointly over the independent old bank. E[RD] must not be factored. Unit first moment is not invariance of the conditional VIX distribution. Full-state four-projection Gaussian means/covariances, conditional-X-on-Y/F quadrature, scalar tensor comparisons, alpha-zero replay and scalar stock recurrence are numerically audited in `../runs/m2_slow_readout_001` and `002`. These are finite numerical checks, not a new Lean certificate, continuous-stock martingale theorem, population-forward certificate or model-family impossibility result.

''';path.write_text(old.replace('# Lean checks for the joint curve fit\n\n','# Lean checks for the joint curve fit\n\n'+addition,1))
    # Repair the run.json evaluator labels; CSVs and score files already carry actual128/256 labels.
    for root in [OUT,WEAK]:
        r=json.loads((root/'run.json').read_text());r.update(status='completed exploratory slow-readout trial; independent checking and reports finalized',SPX_scores=json.loads((root/'SPX/SPX_scores.json').read_text()),production_baseline_unchanged=True,adopted=False);dump(root/'run.json',r)
        for name in ['SLOW_READOUT_PLAN.md','SLOW_READOUT_RESULTS.md','slow_readout_audit.py','audit_slow_recurrence.py','slow_readout_check.py','report_slow_readout.py','finalize_slow_readout.py']:
            shutil.copy2(HERE/name,root/'source'/name)
        # Preserve the exact original script snapshot as well as its final generalized version.
        shutil.copy2(HERE/'slow_readout.py',root/'source/slow_readout_final.py')
        if root==OUT:shutil.copy2('/private/tmp/m2_slow_readout_001.log',root/'run.log')
        else:shutil.copy2('/private/tmp/m2_slow_readout_002.log',root/'run.log');shutil.copy2('/private/tmp/m2_slow_readout_check.log',root/'fresh_validation/run.log')
        des=json.loads((root/'design.json').read_text())
        for file,h in des['source_hashes'].items():
            snap=root/'source'/Path(file).name;assert hashlib.sha256(snap.read_bytes()).hexdigest()==h,(root.name,file)
    counts={}
    for name in ['m2_leverage_scan_001','m2_asymmetric_readout_001','m2_targeted_skew_001','m2_forward_leverage_001']:
        root=OUT.parent/name;manifest=json.loads((root/'artifact_manifest.json').read_text())
        for file,h in manifest.items():assert hashlib.sha256((root/file).read_bytes()).hexdigest()==h,(name,file)
        counts[name]=len(manifest)
    for root in [OUT,WEAK]:
        dump(root/'finalization.json',dict(status='completed and reviewed',adopted=False,report='SLOW_READOUT_RESULTS.md',prior_artifact_hash_counts_verified=counts,initial_source_snapshots_hash_verified=True,final_generalized_script_separately_saved=True,visual_inspection='all seven SPX and four VIX panels reviewed',production_model_sha256=hashlib.sha256((HERE/'runs/m2_maturity_curve_fit_001/selected_model.json').read_bytes()).hexdigest(),visualization_hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in plot.iterdir() if p.is_file()}))
        manifest={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'};dump(root/'artifact_manifest.json',manifest)
        for file,h in manifest.items():assert hashlib.sha256((root/file).read_bytes()).hexdigest()==h
        print('SLOW FINALIZATION PASS',root.name,len(manifest),'artifacts; preserved prior',counts,flush=True)
if __name__=='__main__':main()
