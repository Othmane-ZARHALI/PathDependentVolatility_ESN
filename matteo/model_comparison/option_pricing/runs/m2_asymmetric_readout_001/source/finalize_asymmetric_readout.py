"""Freeze the readout trial, mechanism ablation and durable handoff."""
import hashlib,json,shutil
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .asymmetric_readout import HERE,OUT,dump

def main():
    plot=HERE/'visualizations/asymmetric-readout028';a=json.loads((OUT/'ablation/results.json').read_text());rows=a['rows'];frames={i:pd.read_csv(OUT/f'ablation/SPX/SPX_prices_c{i}_128.csv') for i in range(4)};fig,axes=plt.subplots(2,4,figsize=(16,8));labels=['Original','Multiplier only','Curve only','Both'];colors=['#999999','#d07a23','#9861ad','#187d91'];days=sorted(frames[0].expiry_days.unique())
    for ax,day in zip(axes.flat,days):
        g=frames[0][(frames[0].expiry_days==day)&frames[0].priority_selected];ax.scatter(g.K_over_F,g.market_IV,color='black',s=8,label='Market')
        for c,color,label in zip(range(4),colors,labels):
            g=frames[c][(frames[c].expiry_days==day)&frames[c].priority_selected];ax.scatter(g.K_over_F,g.model_IV,color=color,s=7,label=label)
        ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=7);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'α=1.5 mechanism ablation\nIdentical fresh factor streams\nθ=.8, no new bank\n\nSeparate multiplier / curve effects\nAll seven maturities retained\nInvalid IVs retained in CSV\n\nPost-selection diagnosis\nSampling / grid errors unresolved',va='top');fig.suptitle('The short-maturity improvement requires both readout and curve changes');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(plot/f'SPX-readout-curve-ablation.{ext}',dpi=170)
    plt.close(fig)
    table='\n'.join(f"| {r['name'].replace('_',' ')} | {r['IV_MAE_pp']:.3f} | {r['day29_right_minus_left_IV_pp']:.2f} | {r['VIX_loss']:.3f} | {r['invalid_IV']} |" for r in rows)
    extra=f'''## Separating readout and curve effects

The 29-day improvement in the combined fit must not be attributed to the multiplier alone. A four-way alpha1.5 ablation replays the exact same fresh stock Gaussian streams and uses the same fresh VIX states:

| Variant | SPX IV MAE, pp | 29-day IV(right)-IV(left), pp | VIX loss | Invalid primary SPX IVs |
|---|---:|---:|---:|---:|
{table}

The 29-day gap compares actual quotes nearest K/F=.9 and1.1; its market value is -14.05pp. At fixed curve, the multiplier changes this gap only from -5.86 to -6.05pp. Curve adjustment alone gives -8.51pp, while both give -13.44pp. The combined effect is nonlinear: neither isolated modification reproduces it. The first VIX-sensitive interval spans days27–29, so its refit also directly affects the 29-day SPX expiry; the interpretation must acknowledge this calendar interaction. No maturity-specific model parameter was introduced. The same parameters apply to all seven expiries, where shape problems remain.

The multiplier-only model has a maximum fresh VIX-forward gap of2.74points; combined readout/curve brings this down to.337points without restoring the VIX smiles. The curve-only ablation is a diagnostic control, not a forward-constrained candidate in its own right. Full-candidate stock replay in the ablation agrees exactly with the main fresh comparison. This is a post-selection mechanism check, not an independent model-selection sample.

All approximate paired price-loss confidence intervals for positive strengths versus the alpha0-refit control contain zero. Strong-readout stock sampling and timestep errors remain material: for alpha1.5 the median absolute64/128 price change is1.50quote scales, with maximum11.23. Thus the observed available-IV improvements and quoted skew gaps are finite-grid numerical diagnostics, not accuracy-certified market fits. The matched-state VIX GH32/time8 versus GH48/time12 probe changes VIX by at most8.76e-6points across the declared128-state probes per expiry; this does not resolve rare-tail sampling.

The tested readout/factor combination therefore provides no satisfactory joint calibration. Its useful short-maturity interaction is saved for future work; it does not justify adopting the model or concluding that asymmetric readouts in general fail. See [four-way ablation](visualizations/asymmetric-readout028/SPX-readout-curve-ablation.png).

'''
    f=HERE/'ASYMMETRIC_READOUT_RESULTS.md';s=f.read_text();s=s.replace('## Conclusion to retain',extra+'## Conclusion to retain');f.write_text(s)
    note='''## 28 October 2026 — asymmetric whole-response readout trial completed

User rejected further theta scans and authorized the asymmetric readout. Read [ASYMMETRIC_READOUT_RESULTS.md](ASYMMETRIC_READOUT_RESULTS.md); run `runs/m2_asymmetric_readout_001`, figures `visualizations/asymmetric-readout028`. New D_alpha=exp(alpha softplus(-X)), X=(Y-E_QY)/stationary_sd(Y), uses the existing equity-linked feedback bank. Normalize the whole H0*D product correctly. Theta=.8 and all Gaussian dynamics/readout/stress parameters/floor stay fixed; no endogenous feedback or new bank. Test alpha0,.25,.5,1,1.5,2; seven VIX-sensitive variance levels refit per strength, nine off-window levels frozen. All four training controlled-forward equalities hold within8.5e-10, not exact population-forward matching.

Fresh1024stock paths and fresh4096full-factor VIX draws/expiry: alpha0-refit SPX IVMAE3.278pp, VIX loss.0468; alpha1.5 combined IVMAE2.832pp, VIX loss25.725; six primary signed prices negative/invalid IVs and median stock SE7.88quote scales. Whole-product normalization changes mean stress excess-variance share from92.3% to53.2% at day153, despite unchanged nominal w. The VIX-preserving preselection retains alpha0. No model adopted or original baseline overwritten. All seven SPX maturities/all1751eligible/1294primary and all217VIX quotes/four forwards retained; no clipping/smoothing.

Important four-way paired ablation: at alpha1.5,29-day IV(right)-IV(left) near K/F=.9,1.1: original -5.86pp; multiplier only -6.05; curve only -8.51; both -13.44 versus market -14.05. Strong short-end improvement needs both changes and their interaction; the refitted27–29day interval also affects the29daySPX expiry. Do not attribute the combined result to the multiplier alone. This is no new maturity-specific model parameter. Multiplier-only VIX loss21.278 and maxforwardgap2.74points; combined freshgap.337 despite training equality. Other SPX maturity shapes remain mismatched.

Checks: independent direct-tensor/C-convolution3.7e-11scaled; scalarstock1.1e-16; zero-alpha fullreplay5.6e-14;7forwardderivatives1.4e-9; normalizerGH48vs64 max3.1e-8; VIX GH/time matched probe max8.8e-6points. Strong-stock grid error unresolved (alpha1.5 median1.50,max11.23quote scales for64/128); approximate paired price-loss intervals include zero. Stocks/auxiliary reference and VIX anchor uncertainty remain explicit. Fresh VIX factor draws use full38-dimensional defensive-mixture weights. No new Lean or continuous-martingale theorem.

Durable conclusion: this bounded whole-response multiplier on the existing feedback factor does not solve double calibration; do not repeat the same alpha/theta scan. A possible next hypothesis is component-specific joint normalization that preserves mean M2/stress allocations, or a different equity-linked memory factor. Neither has been tested/adopted here. Full stress/readout refitting was outside this isolated trial, so no model-family impossibility claim follows.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        f=HERE/name;f.write_text(note+f.read_text())
    f=HERE/'formalization/CORRESPONDENCE.md';s=f.read_text();heading='## 28 October 2026: asymmetric readout numerical correspondence\n\n';body='`asymmetric_readout.cpp` implements H0*D/E[H0*D] using the existing exogenous feedback projection. Because old-bank/stress factors are independent, numerator and normalizer are decomposed with E[R*D], not E[R]E[D]. Direct tensor quadrature, scalar stock recursion, zero-strength replay and forward-derivative numerical audits pass; see `../runs/m2_asymmetric_readout_001/implementation_audit.json`. No new Lean theorem was added. Existing normalization identities have their original scoped premises; they do not certify this compiled sampler, new conditional quadrature or continuous-time stock martingality. The new VIX law is repriced rather than claimed invariant.\n\n';f.write_text(s.replace('# Lean checks for the joint curve fit\n\n','# Lean checks for the joint curve fit\n\n'+heading+body,1))
    sources=['audit_asymmetric_readout.py','check_asymmetric_readout.py','plot_asymmetric_readout.py','ablate_asymmetric_readout.py','finalize_asymmetric_readout.py','ASYMMETRIC_READOUT_RESULTS.md']
    for name in sources:shutil.copy2(HERE/name,OUT/'source'/name)
    original=json.loads((OUT/'design.json').read_text())['source_hashes']
    for file,h in original.items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==h,file
    dump(OUT/'finalization.json',dict(status='completed readout scan, fresh stock/VIX checks, independent implementation audit and four-way mechanism ablation',adopted=False,selected_alpha=0.,authoritative_report='ASYMMETRIC_READOUT_RESULTS.md',fresh_VIX_followup='fresh_VIX/run.json',all_original_source_fingerprints_unchanged=True,visualization_hashes={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in plot.iterdir() if p.is_file()}))
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json');manifest={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files};dump(OUT/'artifact_manifest.json',manifest)
    for file,h in manifest.items():assert hashlib.sha256((OUT/file).read_bytes()).hexdigest()==h
    print('FINALIZED',len(manifest),'artifact hashes; original source fingerprints PASS',flush=True)

if __name__=='__main__':main()
