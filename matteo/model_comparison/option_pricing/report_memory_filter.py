"""Plot and explain linear-filter experiments, including paired numerical limits."""
import hashlib,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .memory_filter import HERE,OUT,dump,CONFIG,HALF
from .protected_curve_check import loss
from .plot_joint_curve_fit import implied_vol
CHECK=OUT/'fresh_validation'
PLOT=HERE/'visualizations/memory-filter029'
NAMES=['Baseline','M2 filter, fixed curve','M2 filter + curve refit','Curve refit alone','Whole filter + curve refit']
COLORS=['#888888','#1d8da4','#d17426','#846bb0','#b84b66']

def main():
    PLOT.mkdir(exist_ok=False,parents=True);sel=json.loads((CHECK/'selection.json').read_text());vix=json.loads((CHECK/'VIX/scores.json').read_text());scores=pd.read_csv(CHECK/'SPX/SPX_scores.csv');fine=scores[scores.steps==256].set_index('candidate');frames={c:pd.read_csv(CHECK/f'SPX/SPX_prices_c{c}_256.csv') for c in range(5)};days=sorted(frames[0].expiry_days.unique());metrics=[]
    def save(fig,name):
        fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(PLOT/f'{name}.{ext}',dpi=170)
        plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(17,8))
    for ax,day in zip(axes.flat,days):
        a=frames[0][(frames[0].expiry_days==day)&frames[0].priority_selected];ax.scatter(a.K_over_F,a.market_IV,s=9,c='black',alpha=.6,label='Market')
        for c in range(5):
            a=frames[c][(frames[c].expiry_days==day)&frames[c].priority_selected];ax.scatter(a.K_over_F,a.model_IV,s=8,c=COLORS[c],alpha=.65,label=NAMES[c]);lo=a.loc[(a.K_over_F-.9).abs().idxmin()];hi=a.loc[(a.K_over_F-1.1).abs().idxmin()];metrics.append(dict(candidate=c,day=int(day),IV_MAE_pp=float(np.nanmean(abs(a.model_IV-a.market_IV))),invalid_IV=int(a.model_IV.isna().sum()),stock_mean=float(a.stock_mean.iloc[0]),stock_mean_SE=float(a.stock_mean_SE.iloc[0]),left_K_over_F=float(lo.K_over_F),right_K_over_F=float(hi.K_over_F),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV)))
        ax.set(title=f'{day} days',xlabel='Strike / forward',ylabel='Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=6);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'1024 new stock paths\n256 steps / trading day\nAll seven SPX maturities\n\n'+sel['M2_name']+'\n'+sel['whole_name']+'\n\nLinear memory of nonlinear signals\nOriginal fast M2 and Gaussian dynamics retained\nNo actual-return feedback\nSigned prices and IV failures retained',va='top',fontsize=9);fig.suptitle('SPX: linear memory filters on frozen OU-driven signals');save(fig,'SPX-memory-filter');pd.DataFrame(metrics).to_csv(PLOT/'maturity_metrics.csv',index=False)
    vf={};fig,axes=plt.subplots(1,4,figsize=(17,4.5))
    for c in range(5):
        a=pd.read_csv(CHECK/f'VIX/VIX_prices_c{c}.csv');a['market_IV']=[implied_vol(r.price_mid_U,r.forward,r.strike,r.T_years,r.option_type) for r in a.itertuples()];a['model_IV']=[implied_vol(r.model_price_U,r.forward,r.strike,r.T_years,r.option_type) for r in a.itertuples()];a.to_csv(CHECK/f'VIX/native_IV_c{c}.csv',index=False);vf[c]=a
    for ax,(day,g) in zip(axes,vf[0].groupby('expiry_days',sort=True)):
        ax.scatter(g.strike,g.market_IV,c='black',s=12,label='Market')
        for c in range(5):
            a=vf[c][vf[c].expiry_days==day];ax.scatter(a.strike,a.model_IV,c=COLORS[c],s=12,label=NAMES[c])
        ax.set(title=f'{day} days',xlabel='VIX strike',ylabel='Implied volatility (%)');ax.grid(alpha=.15)
    axes[0].legend(fontsize=6);fig.suptitle('VIX: realized memory histories included in conditional pricing');save(fig,'VIX-memory-filter')
    pilot=pd.read_csv(OUT/'SPX/SPX_scores.csv');cfg=json.loads((OUT/'candidates.json').read_text());ledger=json.loads((OUT/'VIX_fit_ledger.json').read_text());records=[]
    for z in pilot[pilot.steps==256].itertuples():
        c=int(z.candidate);records.append(dict(candidate=c,name=cfg['names'][c],SPX_loss=float(z.loss),SPX_IV_MAE_pp=float(z.IV_MAE_pp),training_VIX_refit_loss=ledger[cfg['exposure_indices'][c]]['loss']))
    pd.DataFrame(records).to_csv(PLOT/'pilot_grid.csv',index=False)
    B=np.load(CHECK/'SPX/SPX_blocks.npz');base=frames[0];mask=base.priority_selected.to_numpy();pq=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();rng=np.random.default_rng(20261029832);intervals=[];grid=[]
    for c in range(5):
        coarse=pd.read_csv(CHECK/f'SPX/SPX_prices_c{c}_128.csv');dv=abs(frames[c].model_price_U-coarse.model_price_U)/base.price_scale_U;grid.append(dict(candidate=c,median_change_quote_scales=float(np.median(dv[mask])),max_change_quote_scales=float(np.max(dv[mask]))))
        if c==0:continue
        changes=[]
        for _ in range(1024):
            ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);C=B['cheap'][ic].mean(0);p=(B[f'c{c}_256'][ie].mean(0)+C-parity)[mask];p0=(B['c0_256'][ie].mean(0)+C-parity)[mask];changes.append(loss(pq,p)-loss(pq,p0))
        intervals.append(dict(candidate=c,loss_delta=float(fine.loc[c,'loss']-fine.loc[0,'loss']),approximate_batch_bootstrap_CI95=np.quantile(changes,[.025,.975]).tolist()))
    dump(PLOT/'paired_intervals.json',intervals);dump(PLOT/'grid_check.json',grid)
    stocks=np.load(CHECK/'stocks.npz')['stocks'];tails=[]
    for c in range(5):
        for j,day in enumerate(days):
            a=base[(base.expiry_days==day)&base.priority_selected];S=stocks[:,j,c,1]
            for side in ['put','call']:
                q=a[a.option_type==side];K=float(q.K_over_F.min() if side=='put' else q.K_over_F.max());pay=np.maximum(K-S,0) if side=='put' else np.maximum(S-K,0);tot=pay.sum();tails.append(dict(candidate=c,day=int(day),side=side,K_over_F=K,raw_payoff_mean=float(pay.mean()),payoff_concentration_ESS=float(tot*tot/(pay@pay)) if tot else 0.,top10_share=float(np.sort(pay)[-10:].sum()/tot) if tot else 0.))
    dump(CHECK/'stock_tail_diagnostics.json',tails);m=pd.DataFrame(metrics);m29=m[m.day==29].set_index('candidate');table='\n'.join(f"| {NAMES[c]} | {fine.loc[c,'loss']:.2f} | {fine.loc[c,'IV_MAE_pp']:.3f} | {vix[c]['option_loss']:.3f} | {vix[c]['max_forward_gap']:.4f} | {int(fine.loc[c,'invalid_IV'])} |" for c in range(5));gain=100*(1-fine.loc[2,'loss']/fine.loc[0,'loss']);audit=json.loads((OUT/'implementation_audit.json').read_text());bridge=json.loads((OUT/'bridge_precision_audit.json').read_text());history=json.loads((CHECK/'history_grid_audit.json').read_text());quad=json.loads((CHECK/'quadrature_audit.json').read_text());summary=dict(scores=fine.reset_index().to_dict('records'),VIX_scores=vix,price_loss_reduction_percent=gain,paired_intervals=intervals,grid=grid,selection=sel,adopted=False);dump(PLOT/'summary.json',summary)
    report=r'''# Linear memory filters of normalized variance signals — 29 October 2026

This experiment implements the linear-memory idea inspired by the two user-supplied Guyon slides. It tests M2-only filtering first, alongside filtering the whole multiplier. It does not implement the nonlinear GLE paper, feedback from actual returns, or a new Gaussian bank. The original fast M2 function and both stress factors remain. This report records a bounded numerical experiment, not a model-family impossibility or a new Lean certificate.

## Construction and interpretation

The frozen baseline is H0=(1-w)R+wA, where R is the original capped and Q-normalized M2 readout and A=(1-r)As+rAf is the existing two-speed lognormal stress mixture. Strictly positive spike shift and qJ are unchanged. Write b=1-w. For each normalized signal J in {R,As,Af}, add deterministic filtering of its realized path:

$$dM^J_t=a(J_t-M^J_t)dt,\qquad M^J_0=1,\qquad a=\log(2)/h.$$

$$M^J_t=e^{-at}+a\int_0^t e^{-a(t-s)}J_sds.$$

M2-only filtering uses

$$H^{M2}_t=b[(1-\chi)R_t+\chi M^R_t]+wA_t.$$

Whole-multiplier filtering uses

$$H^{all}_t=(1-\chi)H^0_t+\chi[bM^R_t+w((1-r)M^{As}_t+rM^{Af}_t)].$$

In both cases V=f+(xi(t)-f)H. Since the inputs and initial filters have mean one, E[H]=1 and E[V]=xi(t), assuming exact Q normalization. Positivity is preserved. The dated Monte Carlo implementation inherits the existing normalizer's approximation error. No new cap or sample normalization is introduced. Share0 recovers the baseline. Share chi controls how much current signal is replaced by its history; h is the effective filter half-life, not the original feedback-kernel decay constant in the slides.

This is a normalized-signal adaptation of linear hysteresis. Literal feedback on level variance with a time-varying xi would require a different deterministic input curve. Filtering a nonlinear signal also differs from passing a slower Gaussian state through a nonlinear readout: average(g(Z)) generally differs from g(average(Z)). Gaussian factors remain exogenous and reusable, but histories, stock paths and VIX exposures must be recomputed. Filters are ordinary path integrals; they add no Brownian motion or volatility feedback to the OU equations.

## Conditional VIX calculation

At expiry T the current OU state alone is insufficient: M_T is a realized random state. Retain that history. For u>=0,

$$E_T[M^J_{T+u}]=e^{-au}M^J_T+a\int_0^u e^{-a(u-s)}E_T[J_{T+s}]ds.$$

Future R means come from the existing full-state conditional Gaussian law, not from treating its three projections as Markov. Stress means are analytic. Integrate conditional filter means and the deterministic xi curve over the same 30-day VIX window. Four VIX forwards are constrained numerically in training; preserving E[V] alone is not a population-forward certificate or invariance of VIX smiles.

To obtain M_T jointly with each importance-sampled terminal state, simulate a Gaussian bridge history, then filter the actual nonlinear R and stress signals. Bridge residuals are sampled independently of the terminal state, with its exact target/proposal likelihood weight retained. There is no replacement of M_T by its conditional mean and no sample-weight normalization.

The terminal covariance has nearly redundant OU modes. Whitening and a relative eigenvalue cutoff1e-11 retain21,21,22,22 directions at the four expiries. This is a documented finite-precision Gaussian conditioning approximation. It is checked against the unconditional covariance of the implemented history, as well as cross-covariances and terminal projections. An earlier cutoff1e-14 run was rejected because its projected history covariance error reached .02146 despite tiny endpoint errors; it is preserved separately in `m2_memory_filter_001.conditioning_diagnostic001` and excluded from all tables and plots here. Initial preparation failures are retained separately. Calibration and checking were repeated after correction, using new checking seeds. The algebra did not change.

## Experiment

Generic half-lives 3, 10 and 30 days, and shares .25 and .5, applied uniformly across maturities; M2-only and whole filtering give12 nonzero settings. Each has a fixed-curve control and a fit of the same seven active variance levels under four training-forward equalities. Nine off-window levels, M2/stress parameters, stock theta=.8 and all Gaussian dynamics remain frozen. There are26 pilot stock configurations including baseline and baseline-refit controls. All 1,751 eligible SPX quotes (1,294 in the primary set), seven SPX maturities, 217 VIX quotes and four VIX forwards are retained.

The pilot uses 512 stock paths in eight independent antithetic batches, and 512 training VIX terminal states per expiry, augmented with independently simulated conditional histories. The terminal training stream is reused research data. Coupled historical grids with 64/128 steps per trading day and stock grids with 128/256 steps are saved. Every nonzero fit is attempted; unsuccessful or uncompetitive settings are retained.

Selection is frozen before new checking draws: training-forward gap<1e-5, training VIX loss<=alpha-zero-refit+.1, then smallest pilot fine-grid native-IV error among refits. The global selection is GLOBAL. Additional checking compares baseline, M2 fixed, M2 refit, the same refitted curve with no filter, and a whole-filter refit. M2 candidate M2NAME; whole candidate WHOLENAME. If a scope has no feasible candidate, its lowest-VIX-loss diagnostic is retained with its feasibility status. The tolerance is an exploratory gate, not a satisfactory-fit certificate.

New checking uses 1,024 stock paths and 1,024 joint VIX terminal states plus independent conditional histories per expiry. All four training equalities are numerical sample statements; checking gaps below remain visible.

| Configuration | SPX interval loss | SPX IV error, pp | VIX interval loss | Max VIX forward gap, points | Invalid primary SPX IVs |
|---|---:|---:|---:|---:|---:|
TABLE

M2 refit SPX price-loss reduction versus baseline is GAIN percent. IV error equally averages available native IVs across seven maturities; unavailable IVs and signed negative prices remain in the data. At 29 days, quoted right-minus-left IV gaps near strike/forward .9 and 1.1 are SLOPES percentage points for the five rows, versus MARKET in the market. Both wing levels and skew must be assessed; aggregate improvement alone does not establish a satisfactory fit.

Approximate paired 1,024-resample eight-batch bootstrap intervals for candidate-minus-baseline price loss are INTERVALS. Candidate and baseline share expensive draws; cheap reference batches are resampled independently and shared across candidates. These intervals exclude deterministic errors, fitting variability and market uncertainty. The existing 65,536-path stress-only auxiliary reference is reused, with valid conditional-lognormal call prices and fixed deterministic delta controls. Full checking is therefore not wholly independent of previous research even though the target draws are new; no market out-of-sample claim is made.

## Numerical audits and limits

Independent 97-step scalar stock recurrence error STOCK; full-horizon zero-share replay error REPLAY; independent scalar bridge-history/filter recurrence error HISTORY; constant-input identity error CONSTANT. These are implementation checks, not formal simulator certification. No new Lean theorem is claimed.

Corrected bridge probe maximum unconditional projected covariance error COV, maximum actual-terminal projection error ENDPOINT, and stable cutoff1e-10vs1e-11 matched-history change RANK. The bridge's dropped covariance directions and residual cross-covariance errors remain in `bridge_audit.json`; this is not exact full-rank conditioning or a uniform tail bound.

Historical 64/128 grid comparison: max VIX change HMAX points, max absolute weighted-mean change HMEAN. Future conditional GH32/time8/32lag knots versus GH48/time12/64lag knots on 64 matched states per expiry: max VIX change QMAX, max absolute weighted-mean change QMEAN. Coupled stock 128/256 median price changes in quote-scale units are GRID. These are probes, not rigorous error bounds.


Large unweighted maxima above occur in extreme importance-sampled states and must be interpreted with their likelihood weights. Across these declared refinement probes, the maximum likelihood-weighted absolute VIX change is 0.000955629 points. The maximum absolute native option-price change is 0.0012763 quote scales. Future quadrature changes include the paired baseline control; history changes use all 1,024 checking histories, while future refinement uses 64 matched states per expiry. No extreme states are clipped. Detailed paired payoff changes and SEs are in `fresh_validation/weighted_numerical_payoff_audit.json`; these are finite probes, not uniform tail bounds.

Raw stock payoff concentration, importance-weighted VIX tail concentration, likelihood-weight sample mass, memory means/SEs, stock means/SEs, signed prices, primary membership, invalid IVs, and separate anchor/difference forward SEs are saved. No clipping, smoothing, sample-forward renormalization or hidden quote deletion is used. Rare-history sampling can remain uncertain even when the terminal covariance and deterministic quadrature agree.

## Assessment

ASSESSMENT

The mechanism is persistent averaging. In particular dM has no direct Brownian shock, so replacing a fraction of a signal by its filtered history reduces that signal's instantaneous Brownian loading. Longer persistence does not automatically strengthen short-maturity downside leverage. The numerical result is specific to these positive linear filters, shares, memories and frozen base parameters; it does not exclude nonlinear hysteresis, other readouts, or exogenous Gaussian models in general.

No candidate is adopted. The production model and previous experiments remain unchanged. Do not repeat this same filter grid as an untested hypothesis. Figures: [seven SPX smiles](visualizations/memory-filter029/SPX-memory-filter.png), [four VIX smiles](visualizations/memory-filter029/VIX-memory-filter.png). Raw run: `runs/m2_memory_filter_001`.
'''
    if abs(gain)<2:assessment='The M2 filter produces no material improvement in this checking sample. The missing SPX skew remains. Whole filtering and curve adjustments must be read together with the VIX costs and paired uncertainty above; no satisfactory double calibration has been demonstrated.'
    else:assessment='The aggregate SPX response must be read alongside the wing gaps, curve-only control, VIX costs and paired uncertainty above. The tested filter has not demonstrated a satisfactory double calibration.'
    replacements={'GLOBAL':sel['global_pilot_selection']['name'],'M2NAME':sel['M2_name'],'WHOLENAME':sel['whole_name'],'TABLE':table,'GAIN':f'{gain:.2f}','SLOPES':', '.join(f"{m29.loc[c,'model_right_minus_left_IV_pp']:.2f}" for c in range(5)),'MARKET':f"{m29.loc[0,'market_right_minus_left_IV_pp']:.2f}",'INTERVALS':json.dumps(intervals),'STOCK':str(audit['scalar_stock_error']),'REPLAY':str(audit['zero_share_full_stock_replay_error']),'HISTORY':str(audit['scalar_conditional_history_error']),'CONSTANT':str(audit['constant_input_filter_error']),'COV':f"{max(a['max_unconditional_projected_bridge_covariance_error'] for a in bridge):.3g}",'ENDPOINT':f"{max(a['actual_terminal_projection_max_error'] for a in bridge):.3g}",'RANK':f"{max(a['rank_cutoff_1e10_vs1e11_max_memory_change'] for a in bridge):.3g}",'HMAX':f"{max(a['max_VIX_change'] for a in history):.4g}",'HMEAN':f"{max(abs(a['weighted_mean_change']) for a in history):.4g}",'QMAX':f"{max(a['max_VIX_change'] for a in quad):.4g}",'QMEAN':f"{max(abs(a['weighted_mean_change']) for a in quad):.4g}",'GRID':', '.join(f"{a['median_change_quote_scales']:.3f}" for a in grid),'ASSESSMENT':assessment}
    for a,b in replacements.items():report=report.replace(a,b)
    (HERE/'MEMORY_FILTER_RESULTS.md').write_text(report);print('MEMORY REPORT COMPLETE',gain,flush=True)
if __name__=='__main__':main()
