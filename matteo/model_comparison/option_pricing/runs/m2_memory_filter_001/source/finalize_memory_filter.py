"""Preserve the corrected linear-memory experiment and its bounded conclusion."""
import hashlib,json,shutil
from pathlib import Path
import pandas as pd
from .memory_filter import HERE,OUT,dump

def main():
    plot=HERE/'visualizations/memory-filter029';summary=json.loads((plot/'summary.json').read_text());s=summary['scores'];v=summary['VIX_scores'];m=pd.read_csv(plot/'maturity_metrics.csv').query('day==29').set_index('candidate');note=f'''## 29 October 2026 — normalized linear-memory filters tested

Read [MEMORY_FILTER_RESULTS.md](MEMORY_FILTER_RESULTS.md); corrected run `runs/m2_memory_filter_001`, figures `visualizations/memory-filter029`. User authorized the Guyon-slide-inspired linear memory experiment. No nonlinear GLE or actual-return feedback was implemented. Gaussian OU dynamics, fast capped/Q-normalized M2 with positive spike shift, two stress factors, stock theta=.8 and floor stay fixed. For J=R,As,Af, dMJ=a(J-MJ)dt with MJ0=1. Test M2-only filtering and whole-multiplier filtering with half-lives3/10/30days, shares.25/.5, uniform across maturities. Variance=f+(xi-f)H; unit means preserve E[v]=xi algebraically, not exact population VIX forwards. This filters nonlinear signal history rather than slower Gaussian states. All downstream prices are rebuilt with shared factors.

Pilot26configurations: baseline and12nonzero fixed-curve controls plus13refits of the same seven active curve levels under four training-forward equalities; nine off-window levels frozen.512pilotstock paths,512weighted terminal-plus-history VIX samples/expiry. Selection M2-only3days/share.25, frozen before1024new stock paths and1024new joint VIX terminal states/conditional histories per expiry. Conditional VIX includes realized memory states, not just current OU projections. Keep1751eligible/1294primarySPX/sevenmaturities and217VIX/fourforwards, signed estimates/no clipping.

Fresh256step results: baseline SPX loss{s[0]['loss']:.2f}/IVMAE{s[0]['IV_MAE_pp']:.3f}pp/VIXloss{v[0]['option_loss']:.5f}; M2 fixed{s[1]['loss']:.2f}/{s[1]['IV_MAE_pp']:.3f}/{v[1]['option_loss']:.5f}; M2 refit{s[2]['loss']:.2f}/{s[2]['IV_MAE_pp']:.3f}/{v[2]['option_loss']:.5f}; curve alone{s[3]['loss']:.2f}/{s[3]['IV_MAE_pp']:.3f}/{v[3]['option_loss']:.5f}; whole3days/share.25refit{s[4]['loss']:.2f}/{s[4]['IV_MAE_pp']:.3f}/{v[4]['option_loss']:.5f}. M2refit price-loss reduction{summary['price_loss_reduction_percent']:.3f}%, approximate paired eight-batch bootstrap loss-change interval[-18.28,15.92], includes0.29day right-minus-left IV gap baseline{m.loc[0,'model_right_minus_left_IV_pp']:.2f},M2refit{m.loc[2,'model_right_minus_left_IV_pp']:.2f},whole{m.loc[4,'model_right_minus_left_IV_pp']:.2f}pp versusmarket{m.loc[0,'market_right_minus_left_IV_pp']:.2f}. No meaningful SPX improvement; missing skew persists. Fresh M2refit max forward gap{v[2]['max_forward_gap']:.5f}points; all four training equalities remain sample constraints.

History conditioning numerical issue was caught and corrected before final conclusions: direct inversion was unreliable; a whitened cutoff1e-14 pilot was rejected for projected history covariance error.02146 and preserved as `m2_memory_filter_001.conditioning_diagnostic001`. Corrected cutoff1e-11 retains21/21/22/22terminal directions, covariance-gated below2e-6; measured maximum5.53e-7. Actual terminal projection error<9.5e-9. Stable cutoff1e-10vs1e-11 gives identical matched-history probes. Scalar stock/filter recurrence and zero-share full-horizon replay errors0; independent scalar bridge-history recurrence2.22e-15. No Lean certification/uniform pricing bound.

Historical64/128 and futureGH32/time8/32lag versus48/time12/64lag refinements have large unweighted tail maxima. Do not interpret these without likelihood weights: maximum weighted absolute VIX change.000956points and native option-price change.001277quote scales in the declared probes. All extreme states retained. Stock grids128/256, paired SE/bootstrap intervals, raw/weighted payoff concentration, memory means/SEs and invalid IVs are saved. Reused65536-path conditional-lognormal auxiliary reference and same market snapshot limit independence claims.

Durable conclusion: positive linear averaging preserves mean variance and adds persistence, but it weakens the filtered signal's direct Brownian response. These tested filters do not supply the missing SPX skew. Longer whole-multiplier filtering also damages the VIX training fit. Do not repeat this same filter grid as untested; this is not a no-go theorem for nonlinear hysteresis or other exogenous Gaussian readouts. No candidate adopted; production model and all previous experiment artifacts remain unchanged.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        p=HERE/name;p.write_text(note+p.read_text())
    p=HERE/'formalization/CORRESPONDENCE.md';old=p.read_text();addition='''## 29 October 2026: normalized linear memory

`memory_filter.cpp/.py` filter the realized nonlinear normalized M2 and stress signals while retaining the exogenous Gaussian factor dynamics. Constant-input and mean-one filter identities are algebraic; they do not imply VIX-law invariance or certified population-forward matching. Conditional pricing uses realized filter histories plus the full-state conditional future signal law. Numerical scalar stock/history recurrences, zero-share replay, whitened Gaussian bridge covariance/terminal probes and likelihood-weighted payoff refinement checks are in `../runs/m2_memory_filter_001`. The corrected bridge uses a rank-truncated Gaussian conditioning approximation; rejected high-rank conditioning remains separate. None of these numerical checks is a new Lean certificate, continuous-stock martingale theorem or model-family impossibility result.

''';p.write_text(old.replace('# Lean checks for the joint curve fit\n\n','# Lean checks for the joint curve fit\n\n'+addition,1))
    for name in ['MEMORY_FILTER_RESULTS.md','audit_memory_filter.py','audit_memory_bridge.py','audit_memory_payoffs.py','check_memory_filter.py','report_memory_filter.py','finalize_memory_filter.py']:
        shutil.copy2(HERE/name,OUT/'source'/name)
    for source,h in json.loads((OUT/'design.json').read_text())['source_hashes'].items():assert hashlib.sha256((OUT/'source'/Path(source).name).read_bytes()).hexdigest()==h,source
    for name in ['asymmetric_readout.py','protected_curve_check.py','maturity_curve.py','forward_leverage.py','leverage_scan.py']:
        shutil.copy2(HERE/name,OUT/'source'/name)
    shutil.copy2('/private/tmp/m2_memory_filter.log',OUT/'run.log');shutil.copy2('/private/tmp/m2_memory_filter_check.log',OUT/'fresh_validation/run.log');shutil.copy2('/private/tmp/m2_memory_filter_audit.log',OUT/'implementation_audit.log');shutil.copy2('/private/tmp/m2_memory_bridge_audit.log',OUT/'bridge_precision_audit.log')
    counts={}
    for name in ['m2_leverage_scan_001','m2_asymmetric_readout_001','m2_targeted_skew_001','m2_forward_leverage_001','m2_slow_readout_001','m2_slow_readout_002']:
        root=OUT.parent/name;manifest=json.loads((root/'artifact_manifest.json').read_text())
        for file,h in manifest.items():assert hashlib.sha256((root/file).read_bytes()).hexdigest()==h,(name,file)
        counts[name]=len(manifest)
    r=json.loads((OUT/'run.json').read_text());r.update(status='completed corrected memory-filter experiment, fresh checking, numerical audits, reviewed figures and durable conclusion',adopted=False,report='MEMORY_FILTER_RESULTS.md');dump(OUT/'run.json',r)
    dump(OUT/'finalization.json',dict(status='completed and reviewed',adopted=False,prior_artifact_hash_counts_verified=counts,initial_source_snapshots_verified=True,visual_inspection='all seven SPX and four VIX panels reviewed',production_model_sha256=hashlib.sha256((HERE/'runs/m2_maturity_curve_fit_001/selected_model.json').read_bytes()).hexdigest(),excluded_conditioning_run='m2_memory_filter_001.conditioning_diagnostic001',visualization_hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in plot.iterdir() if p.is_file()}))
    manifest={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'};dump(OUT/'artifact_manifest.json',manifest)
    for file,h in manifest.items():assert hashlib.sha256((OUT/file).read_bytes()).hexdigest()==h
    print('MEMORY FILTER FINALIZATION PASS',len(manifest),'artifacts; preserved',counts,flush=True)
if __name__=='__main__':main()
