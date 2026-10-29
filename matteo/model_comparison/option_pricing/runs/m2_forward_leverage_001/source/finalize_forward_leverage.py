"""Freeze the future-leverage diagnostic and preserve the negative result."""
import hashlib,json,shutil
from pathlib import Path
import pandas as pd
from .forward_leverage import OUT,TARGET,HERE,dump

def main():
    plot=HERE/'visualizations/forward-leverage028';summary=json.loads((plot/'summary.json').read_text());memory=json.loads((OUT/'feedback_memory_diagnostic.json').read_text());scores=summary['scores'];m=pd.read_csv(plot/'maturity_metrics.csv');day29=m[m.day==29].set_index('candidate');base=next(s for s in scores if s['candidate']==0);forward=next(s for s in scores if s['candidate']==3);slope0=day29.loc[0,'model_right_minus_left_IV_pp'];slope3=day29.loc[3,'model_right_minus_left_IV_pp'];market=day29.loc[0,'market_right_minus_left_IV_pp'];losschange=100*(forward['loss']/base['loss']-1);report=HERE/'FORWARD_LEVERAGE_RESULTS.md';s=report.read_text();extra=f'''## Assessment and memory diagnostic

This particular30day exposure rule does not materially improve the SPX fit. On these paired draws, the8update price loss changes {base['loss']:.2f}to{forward['loss']:.2f}({losschange:+.2f}%), while native-IV MAE changes only{base['IV_MAE_pp']:.3f}to{forward['IV_MAE_pp']:.3f}pp. Its29day quoted IV gap becomes less negative, {slope0:.2f}to{slope3:.2f}pp versus market{market:.2f}pp. All paired price-loss sampling intervals include zero. There is no demonstrated improvement to adopt, and the missing left-wing slope remains visible.

Together with the earlier constant-rotation and instantaneous-gradient tests, this gives no empirical reason to keep adjusting stock-shock direction as the next priority. It does not prove that every possible state-dependent correlation or lookahead kernel fails. The tests have held the variance readout, bank projections and stress distributions fixed.

A separate, exact finite-dimensional Gaussian diagnostic exposes a plausible readout-memory limitation. For the existing feedback projection Y, stationary covariance gives autocorrelation{memory['correlation'][1]:.3f}after1day,{memory['correlation'][2]:.3f}after7days and{memory['correlation'][3]:.3f}after30days. This is a model stationary-memory calculation, not a dated option calibration or historical correlation estimate. Its projection weights are fixed and strongly favour fast modes; the slow OU modes already exist in the bank.

A concrete next hypothesis is to give a slower projection of those existing modes an independent role in the asymmetric readout, while retaining the current fast component and the Zumbach-motivated design. Gaussian state paths remain reusable; new readout normalization, stock paths and VIX conditional expectations would need rebuilding. This readout-memory hypothesis is untested and must not be presented as the diagnosed sole cause of the skew gap. It is more distinct from the failed exposure scans than another constant-theta/global-multiplier variation. No new bank, readout or model is adopted in this experiment.

''';s=s.replace('## What this establishes',extra+'## What this establishes');report.write_text(s)
    note='''## 28 October 2026 — future-integrated-variance stock exposure tested

Read [FORWARD_LEVERAGE_RESULTS.md](FORWARD_LEVERAGE_RESULTS.md); run `runs/m2_forward_leverage_001`, figures `visualizations/forward-leverage028`. User authorized the proposed persistent-leverage test. Stock follows negative Brownian sensitivity of conditional next30day integrated variance; same horizon for every maturity. Existing full24mode Gaussian bank, variance readout/positive spike shift, stress parameters, floor and16-level curve stay fixed. VIX distribution/options/forwards unchanged, including inherited numerical error. No endogenous factor feedback, new bank, readout fit or constanttheta sweep.

Full-state conditional Gaussian gradient: independent C convolution GH48, tensor Y/F GH12,0plus22geometric lag knots; hat weights integrate xi-floor exactly across curve breaks. Fixed last xi continued only for signal through1060days. Direction updated4or8times/tradingday on shared factors; paired128/256stock grids. GH12vs32 max probe direction difference.000514;22vs44lagknots.01269; all13primitive finite differences9.3e-12; scalar97step stock recursion/unitnorm2.2e-16; frozen baseline/instant full replay0. No new Lean or continuous-martingale theorem.

1024freshstock paths/eightantithetic batches, prior65536auxiliary reference reused. At256steps: fixed price loss1909.19/IVMAE3.201pp; instantaneous1978.79/3.163; future4updates1969.78/3.175; future8updates1941.31/3.157. Invalid primary IVs1/1/2/3 retained. Future8 loss is1.68%higher; all approximate paired loss-difference intervals include0.29day quoted right-minus-left IV gap nearK/F=.9,1.1 fixed-6.39pp,future8-4.51,market-14.05. No demonstrated fit improvement; missing skew persists. VIX217quotes/four forwards inherited unchanged; no new exact population-forward claim. All1751eligible/1294primarySPX/sevenmaturities retained, signed prices/no clipping. Update4to8 median.253quote scales versus paired SE.499; future8 stock128/256 median.144scales. Final sampling/tail accuracy unresolved.

Durable conclusion: do not repeat this same30day stock-direction rule, constanttheta grids or the global multiplier scan. This is not a no-go theorem for all correlations/horizons. Next untested hypothesis: change readout memory using slower projections of the existing OU modes, preserving the fast Zumbach-motivated component. Exact stationary Gaussian memory diagnostic of current feedbackY gives correlation.170after1day,.0479after7days,.0173after30days; this suggests a target but does not prove the cause of SPX misfit. No new projection/readout/bank adopted; baseline artifacts preserved.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        path=HERE/name;path.write_text(note+path.read_text())
    path=HERE/'formalization/CORRESPONDENCE.md';s=path.read_text();addition='''## 28 October 2026: future-integrated-variance sensitivity

`forward_leverage.cpp` computes the Brownian loading of conditional next30day integrated variance using full24mode decays, conditional Gaussian readout derivatives, analytic stresses and exact curve-break hat weights. Integration/recurrence numerical audits are in `../runs/m2_forward_leverage_001`. Primitive finite differences independently reprice the finite-rule conditional variance target; this checks the derivative correspondence at declared states, not an exact infinite-dimensional quadrature theorem. Unit predictable stock loading preserves the unchanged exogenous factor/variance law, so VIX is inherited unchanged. No new Lean proof, option-skew optimum or continuous-stock martingale result is claimed. Fixed projection memory is separately diagnosed from stationary Gaussian covariance; it is not proof of calibration infeasibility or of the proposed slower-readout remedy.

''';path.write_text(s.replace('# Lean checks for the joint curve fit\n\n','# Lean checks for the joint curve fit\n\n'+addition,1))
    sources=['audit_forward_leverage.py','audit_forward_recurrence.py','run_forward_leverage.py','report_forward_leverage.py','finalize_forward_leverage.py','FORWARD_LEVERAGE_RESULTS.md']
    for name in sources:shutil.copy2(HERE/name,OUT/'source'/name)
    for file,h in json.loads((OUT/'design.json').read_text())['source_hashes'].items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==h,file
    counts={}
    for name in ['m2_leverage_scan_001','m2_asymmetric_readout_001','m2_targeted_skew_001']:
        folder=OUT.parent/name;manifest=json.loads((folder/'artifact_manifest.json').read_text())
        for file,h in manifest.items():assert hashlib.sha256((folder/file).read_bytes()).hexdigest()==h,(name,file)
        counts[name]=len(manifest)
    dump(OUT/'finalization.json',dict(status='completed frozen future-leverage comparison, integration and scalar audits, plots, memory diagnostic, durable conclusion',adopted=False,authoritative_report='FORWARD_LEVERAGE_RESULTS.md',prior_artifact_hash_counts_verified=counts,initial_source_fingerprints_unchanged=True,visualization_hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in plot.iterdir() if p.is_file()}))
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json');manifest={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files};dump(OUT/'artifact_manifest.json',manifest)
    for file,h in manifest.items():assert hashlib.sha256((OUT/file).read_bytes()).hexdigest()==h
    print('FORWARD FINALIZATION PASS',len(manifest),'artifacts; prior',counts,flush=True)

if __name__=='__main__':main()
