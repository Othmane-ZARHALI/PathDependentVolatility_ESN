"""Native smiles, paired uncertainty and durable future-leverage conclusion."""
import hashlib,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .forward_leverage import OUT,TARGET,HERE,arr,dump
from .asymmetric_readout import VBASE
from .protected_curve_check import loss
from .vix_only_curve_fit import DAYS

def main():
    plot=HERE/'visualizations/forward-leverage028';plot.mkdir(exist_ok=False,parents=True);names=['Fixed exposure','Instantaneous variance','Future30days,4updates','Future30days,8updates'];colors=['#999999','#b4772c','#8659a1','#167e96'];frames={c:pd.read_csv(OUT/f'prices/SPX_prices_c{c}_256.csv') for c in range(4)};scores=pd.read_csv(OUT/'prices/SPX_scores.csv');scores=scores[scores.steps==256];days=sorted(frames[0].expiry_days.unique());metrics=[]
    def save(fig,name):
        fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(plot/f'{name}.{ext}',dpi=170)
        plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,day in zip(axes.flat,days):
        g=frames[0][(frames[0].expiry_days==day)&frames[0].priority_selected];ax.scatter(g.K_over_F,g.market_IV,color='black',s=8,alpha=.55,label='Market')
        for c,color,label in zip(range(4),colors,names):
            a=frames[c][(frames[c].expiry_days==day)&frames[c].priority_selected];ax.scatter(a.K_over_F,a.model_IV,color=color,s=7,alpha=.7,label=label);lo=a.loc[(a.K_over_F-.9).abs().idxmin()];hi=a.loc[(a.K_over_F-1.1).abs().idxmin()];metrics.append(dict(candidate=c,name=label,day=int(day),IV_MAE_pp=float(np.nanmean(abs(a.model_IV-a.market_IV))),invalid_IV=int(a.model_IV.isna().sum()),stock_mean=float(a.stock_mean.iloc[0]),stock_mean_SE=float(a.stock_mean_SE.iloc[0]),left_K_over_F=float(lo.K_over_F),right_K_over_F=float(hi.K_over_F),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV)))
        ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=6.4);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'1024 fresh stock paths\n256steps / trading day\nPaired128/256variance grids\nSame factors / readout / curve\n\nStock exposure targets conditional\nnext30day integrated variance\n4vs8updates per trading day\n\nVIX law unchanged for all models\nSigned prices / invalid IVs retained',va='top');fig.suptitle('SPX: exposure to shocks that raise future integrated variance');save(fig,'SPX-forward-leverage');pd.DataFrame(metrics).to_csv(plot/'maturity_metrics.csv',index=False)
    # Existing VIX prices are inherited, not newly estimated by this experiment.
    q=pd.read_csv(VBASE/'fresh_check/native_iv_two_speed.csv');fig,axes=plt.subplots(1,4,figsize=(15,4))
    for ax,(day,g) in zip(axes,q.groupby('expiry_days',sort=True)):
        ax.scatter(g.strike,g.market_native_IV,color='black',s=10,label='Market');ax.plot(g.strike,g.model_native_IV,'.-',color='#167e96',label='All four couplings');ax.set_title(f'{day} days');ax.set_xlabel('VIX strike');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes[0].legend(fontsize=7);fig.suptitle('VIX is inherited unchanged: only the stock shock exposure changes');save(fig,'VIX-unchanged')
    B=np.load(OUT/'prices/SPX_blocks.npz');base=frames[0];mask=base.priority_selected.to_numpy();pq=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();rng=np.random.default_rng(20261028707);intervals=[];grids=[]
    for c in range(4):
        coarse=pd.read_csv(OUT/f'prices/SPX_prices_c{c}_128.csv');fine=frames[c];delta=abs(fine.model_price_U-coarse.model_price_U)/fine.price_scale_U;grids.append(dict(candidate=c,median_grid_change_scales=float(np.median(delta[mask])),max_grid_change_scales=float(np.max(delta[mask]))))
        if c==0:continue
        E=B[f'c{c}_256'];E0=B['c0_256'];differences=[]
        for _ in range(1024):
            ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);C=B['cheap'][ic].mean(0);differences.append(loss(pq,(E[ie].mean(0)+C-parity)[mask])-loss(pq,(E0[ie].mean(0)+C-parity)[mask]))
        intervals.append(dict(candidate=c,name=names[c],price_loss_delta_CI95=np.quantile(differences,[.025,.975]).tolist()))
    a=frames[2];b=frames[3];difference=abs(b.model_price_U-a.model_price_U)/b.price_scale_U;E=B['c3_256']-B['c2_256'];update=dict(median_price_change_scales=float(np.median(difference[mask])),maximum_price_change_scales=float(np.max(difference[mask])),median_paired_SE_scales=float(np.median(E.std(0,ddof=1)[mask]/np.sqrt(8)/b.price_scale_U.to_numpy()[mask])));dump(plot/'update_check.json',update);dump(plot/'grid_check.json',grids);dump(plot/'paired_intervals.json',intervals)
    # Matched paired update differences, without changing signs or forwards.
    paired=b[['expiry_days','strike','option_type','priority_selected','price_scale_U']].copy();paired['eight_minus_four_price_U']=b.model_price_U-a.model_price_U;paired['paired_SE_U']=E.std(0,ddof=1)/np.sqrt(8);paired.to_csv(plot/'paired_update_prices.csv',index=False)
    stocks=np.load(OUT/'stocks.npz')['stocks'];tail=[]
    for c in range(4):
        for j,day in enumerate(days):
            g=base[(base.expiry_days==day)&base.priority_selected];S=stocks[:,j,c,1]
            for side in ['put','call']:
                gg=g[g.option_type==side]
                if len(gg)==0:continue
                K=float(gg.K_over_F.min() if side=='put' else gg.K_over_F.max());pay=np.maximum(K-S,0) if side=='put' else np.maximum(S-K,0);total=pay.sum();tail.append(dict(candidate=c,day=int(day),side=side,K_over_F=K,raw_payoff_mean=float(pay.mean()),payoff_concentration_ESS=float(total**2/(pay@pay)) if total else 0.,top10_share=float(np.sort(pay)[-10:].sum()/total) if total else 0.,note='raw target-Q payoffs; concentration statistic not independent sample count; antithetic pairs determine uncertainty'))
    dump(OUT/'tail_diagnostics.json',tail);audit=json.loads((OUT/'integration_audit.json').read_text());rec=json.loads((OUT/'recurrence_audit.json').read_text());replay=json.loads((OUT/'baseline_replay.json').read_text());run=json.loads((OUT/'run.json').read_text());m=pd.DataFrame(metrics);m29=m[m.day==29].set_index('candidate');table='\n'.join(f"| {names[int(r.candidate)]} | {r.loss:.2f} | {r.IV_MAE_pp:.3f} | {int(r.invalid_IV)} | {r.median_SE_scales:.2f} |" for r in scores.itertuples());s0=scores[scores.candidate==0].iloc[0];s3=scores[scores.candidate==3].iloc[0];gain=100*(1-s3.loss/s0.loss)
    report=f'''# Stock exposure to future integrated variance — 28 October 2026

Implemented the user-authorized follow-up to TARGETED_SKEW_RESULTS.md. All variance/factor parameters stay fixed; only the stock shock direction changes. No new bank, endogenous factor feedback, readout fitting or constant-theta search. All seven SPX maturities/1751eligible/1294primary quotes and217VIX quotes/four forwards remain. This is an exploratory fixed-model comparison, not an adopted double calibration.

## What the stock now follows

For the same30day horizon L=30/365 at every time/maturity, define

$$F(t,X_t,Z_t)=E_Q[\\int_t^{{t+L}}v_u du\\mid\\mathcal F_t],\\qquad dF_t=\\text{{drift}}dt+\\mathbf K_t\\cdot d\\mathbf B_t.$$

The new stock direction is

$$dS_t/S_t=\\sqrt{{v_t}}\\mathbf a_t\\cdot d\\mathbf B_t,\\qquad \\mathbf a_t=-\\mathbf K_t/\\|\\mathbf K_t\\|.$$

At zero K, use the old mixed theta=.8 loading. The state derivative of F uses all24retained Gaussian modes and their individual decay speeds, plus both stresses, in13independent primitive Brownian coordinates. It is not a three-projection Markov approximation and does not use realized future shocks. Compare the fixed mixed stock loading, previous instantaneous-variance loading, future-variance loading updated4and8times/trading-day. Directions are predictable and unit length between updates. The purpose is to test persistent leverage, not to prove that this choice maximizes finite-wing implied skew.

The original M2 cap/readout, positive spike shift, stress rates/amplitudes/weights, floor,16-level curve, initial state and Q law remain fixed. The entire conditional variance-factor law is unchanged, so VIX options/forwards are identical to the original model. Existing numerical forward gaps (up to.0552points) are inherited; this experiment does not newly establish exact population forward matching. Stock martingality in continuous time is separate from unit Brownian loading.

The last variance level is continued flat beyond day1030 solely to define the common30day signal through day1060. No observed-horizon variance value changes and no new level is fitted. The existing normalizer covers this lookahead.

## Calculation and validation

For each future lag, project the complete centered24mode state with its exponential decays and add the deterministic Q mean. Integrate the independent cluster C dimension of the capped readout derivative using48point Gaussian convolution and an8193knot cubic-Hermite table; integrate feedback/spike Y/F with12point tensor Gaussian quadrature. Stress sensitivities are analytic. Sum the sensitivities on0plus22geometric lag knots from1e-5years to30/365. Piecewise-linear hat weights integrate xi(t+lag)-floor exactly across every deterministic curve break. Time means, normalizers and decayed diffusion matrices are precomputed, not fitted.

Independent probes compare12/24/32Gaussian nodes,22/44lag knots and direct finite differences of independently repriced F in every primitive direction. At60declared Q/tail-amplified probe states, max unit-direction distance GH12vs32 is {audit['GH12_vs32_max_direction_distance']:.3g};22vs44lag knots gives {audit['lag22_vs44_max_direction_distance']:.3g}. These meet the declared exploratory.03direction-distance gate, not a uniform error bound. Full primitive finite-difference sensitivity error is {audit['finite_difference_primitive_sensitivity_max_scaled_error']:.3g}; exact curve hat-weight integral error is {audit['hat_curve_integral_max_error']:.3g}. The signal's Gaussian old-bank gradient is obtained by differentiating the conditional expectation, not by evaluating the raw readout derivative at the conditional mean.

Missing primitive Brownian increments use the previously audited conditional OU bridge. Frozen baseline/instantaneous full-horizon replay matches earlier stocks to {replay['maximum_baseline_and_instant_stock_error']:.3g} on four declared paths, with original counter stride128. An independent97step scalar recurrence checks both update schedules, both stock grids and antithetic partners to {rec['scalar_stock_recursion_max_absolute_error']:.3g}; unit forward direction error is {rec['unit_forward_stock_loading_max_error']:.3g}. No new Lean or compiled-simulator theorem is claimed.

## Fresh paired results

1024new full-model stock paths,512antithetic pairs/eightindependent batches, coupled128/256stock grids. The4/8update schedules share the identical factors and variance. All choices were fixed before pricing; there is no market-selected winner. The table uses256steps/trading-day; IV MAE equally averages available native IV errors across seven maturities.

| Stock exposure | SPX interval loss | IV MAE,pp | Invalid primary IVs | Median MC SE / quote scale |
|---|---:|---:|---:|---:|
{table}

Future30day/8update point-estimate price loss relative change versus fixed exposure: {-gain:+.2f}%(positive means worse). Approximate paired1024resample batch-bootstrap95%loss-change intervals: {json.dumps(intervals)}. Eight batches underpin these sampling intervals; deterministic integration, stock-grid, normalizer and market errors are excluded. Invalid IVs/negative signed price estimates are retained. A lower available-IV error with omitted failed IVs is not a final calibration certificate.

At29days, quoted K/F near.9,1.1 gives IV(right)-IV(left) {m29.loc[0,'model_right_minus_left_IV_pp']:.2f}pp for fixed exposure, {m29.loc[1,'model_right_minus_left_IV_pp']:.2f}pp for instantaneous exposure, {m29.loc[3,'model_right_minus_left_IV_pp']:.2f}pp for future30day/8updates, versus market {m29.loc[0,'market_right_minus_left_IV_pp']:.2f}pp. All seven maturity panels and per-maturity metrics are saved; assess the slope and both wings, not only the aggregate price loss.

Update-frequency diagnostic4to8: median absolute price change {update['median_price_change_scales']:.3f}quote scales, maximum {update['maximum_price_change_scales']:.3f}; median paired SE of the change {update['median_paired_SE_scales']:.3f}scales. Stock128/256grid diagnostics: {json.dumps(grids)}. Nonzero changes remain numerical approximation error, not a separate calibrated model benefit.

## Estimator and limits

Use the same65536path stress-only auxiliary reference and its correct conditional lognormal call prices; that auxiliary stock model has not acquired the new coupling. The expensive correction is target/reference payoff difference minus a deterministic delta times stock difference. Include independent reference-batch uncertainty. No future-stress-dependent delta, sample-forward renormalization, clipping or price smoothing. The reference is reused from earlier research, so the complete estimator is not newly independent even though the1024factor paths are new. This is the same market snapshot, not market out-of-sample validation. Raw downside/upside payoff concentration, stock means/SEs, quote membership, prices, invalid IVs and paired differences remain in the run.

VIX prices/forwards are inherited unchanged, not newly re-estimated here. The unchanged factor/readout law is the invariance argument. No new claim about precise population-forward matching, continuous-stock martingality, optimum leverage direction for option prices or impossibility of fitting the model follows.

## What this establishes

The implemented direction targets expected future integrated variance with verified full-state sensitivities and retains the successful VIX law. The response and limits must be read from the paired prices, all-maturity smiles and uncertainty above. This is a diagnostic of one common30day horizon and the existing variance model, not a broad search over horizon functions or readouts. Do not revert to another constant-theta/global-multiplier scan or treat a favorable aggregate metric as a satisfactory joint calibration. Any extension should state which missing maturity/wing mechanism it addresses and keep factor reuse and VIX protection explicit.

Figures: [all seven SPX smiles](visualizations/forward-leverage028/SPX-forward-leverage.png), [unchanged VIX smiles](visualizations/forward-leverage028/VIX-unchanged.png). Detailed run: `runs/m2_forward_leverage_001`. Stock simulation/pricing wall time {run['wall_seconds']:.1f}seconds, excluding preparation and audits.
'''
    (HERE/'FORWARD_LEVERAGE_RESULTS.md').write_text(report);dump(plot/'summary.json',dict(scores=scores.to_dict('records'),price_loss_reduction_percent=gain,paired_intervals=intervals,grid=grids,update=update,adopted=False));print('FORWARD REPORT COMPLETE',gain,flush=True)

if __name__=='__main__':main()
