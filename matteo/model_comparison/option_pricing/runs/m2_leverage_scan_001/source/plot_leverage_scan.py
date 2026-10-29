"""Native smiles, grid diagnostics and transparent experiment report."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .leverage_scan import HERE,OUT,THETAS,stock,blocks,dump

def main():
    out=HERE/'visualizations/leverage028';out.mkdir(exist_ok=False,parents=True);scores=pd.read_csv(OUT/'fresh/scores.csv');sel=json.loads((OUT/'selection.json').read_text())['selected'];theta=sel['theta'];direction=sel['direction'];base=pd.read_csv(OUT/'fresh/prices_mixed_0_128.csv');chosen=pd.read_csv(OUT/f'fresh/prices_{direction}_{theta:g}_128.csv');days=sorted(base.expiry_days.unique());primary=base.priority_selected.to_numpy();rows=[]
    for name,q in [('theta0',base),(f'{direction}_{theta:g}',chosen)]:
        for day,g in q[q.priority_selected].groupby('expiry_days'):
            lo=g.loc[(g.K_over_F-.9).abs().idxmin()];hi=g.loc[(g.K_over_F-1.1).abs().idxmin()];res=np.where(g.model_price_U<g.price_bid_U,g.model_price_U-g.price_bid_U,np.where(g.model_price_U>g.price_ask_U,g.model_price_U-g.price_ask_U,0))/g.price_scale_U;rows.append(dict(model=name,day=int(day),loss=float(np.mean(res**2)),IV_MAE_pp=float(np.nanmean(abs(g.model_IV-g.market_IV))),invalid_IV=int(g.model_IV.isna().sum()),stock_mean=float(g.stock_mean.iloc[0]),stock_mean_SE=float(g.stock_mean_SE.iloc[0]),Kleft=float(lo.K_over_F),Kright=float(hi.K_over_F),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV)))
    pd.DataFrame(rows).to_csv(out/'maturity_metrics.csv',index=False)
    def save(fig,name):
        fig.tight_layout()
        for ext in ['png','svg']:fig.savefig(out/f'{name}.{ext}',dpi=160)
        plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(16,8));colors=plt.cm.viridis(np.linspace(.05,.95,6))
    for ax,day in zip(axes.flat,days):
        g=base[(base.expiry_days==day)&base.priority_selected];ax.scatter(g.K_over_F,g.market_IV,c='black',s=8,alpha=.5,label='Market')
        for t,color in zip(THETAS,colors):
            a=pd.read_csv(OUT/f'fresh/prices_mixed_{t:g}_128.csv');a=a[(a.expiry_days==day)&a.priority_selected];ax.scatter(a.K_over_F,a.model_IV,color=color,s=6,alpha=.7,label=f'θ = {t:g}')
        ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=7);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'Mixed stress leverage\nFixed variance model and curve\nFresh 1024 full-model paths\n65536 auxiliary reference paths\n128 steps / trading day\n\nAll primary quotes displayed\nInvalid IVs retained in CSV\nSame market snapshot',va='top');fig.suptitle('Increasing negative stock / stress dependence changes the SPX smile');save(fig,'SPX-theta-fan')
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,day in zip(axes.flat,days):
        g=base[(base.expiry_days==day)&base.priority_selected];ax.scatter(g.K_over_F,g.market_IV,c='black',s=9,alpha=.6,label='Market');ax.scatter(g.K_over_F,g.model_IV,c='#999999',s=7,label='θ = 0')
        for name,color in [('mixed','#147d92'),('slow','#bf642d'),('fast','#8c4ba4')]:
            a=pd.read_csv(OUT/f'fresh/prices_{name}_{theta:g}_128.csv');a=a[(a.expiry_days==day)&a.priority_selected];ax.scatter(a.K_over_F,a.model_IV,c=color,s=7,label=name)
        ax.set_title(f'{day} days');ax.set_xlabel('Strike / forward');ax.set_ylabel('Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=7);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,f'Common θ = {theta:g}\nChosen on replay price loss\nNo curve/readout fitting\nVIX factor law unchanged\n\nSlow/fast panels isolate direction\nSampling error remains material',va='top');fig.suptitle('Which stress shocks should also move the stock?');save(fig,'SPX-stress-directions')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for name,color in [('mixed','#147d92'),('slow','#bf642d'),('fast','#8c4ba4')]:
        a=scores[(scores.direction==name)&(scores.steps==128)]
        axes[0].plot(a.theta,a.loss,'.-',color=color,label=name);axes[1].plot(a.theta,a.IV_MAE_pp,'.-',color=color,label=name)
    for ax in axes:ax.set_xlabel('θ');ax.legend();ax.grid(alpha=.15)
    axes[0].set_ylabel('Equal-maturity price interval loss');axes[1].set_ylabel('Equal-maturity IV MAE (percentage points)');fig.suptitle('Fresh expensive paths; the cheap reference sample is shared with screening');save(fig,'SPX-response')
    # Paired grid differences, with the same auxiliary blocks cancelling.
    grid=[];a=pd.read_csv(OUT/f'fresh/prices_{direction}_{theta:g}_64.csv');b=chosen;d=abs(a.model_price_U-b.model_price_U)/b.price_scale_U
    for day in days:
        mask=primary & b.expiry_days.eq(day).to_numpy();grid.append(dict(day=int(day),median_abs_grid_difference_scales=float(np.median(d[mask])),max_abs_grid_difference_scales=float(np.max(d[mask]))))
    pd.DataFrame(grid).to_csv(out/'grid_metrics.csv',index=False)
    # Reference forward diagnostic: analytically integrate independent W0.
    ref=np.load(OUT/'reference.npz')['stats'];x=np.array(json.loads((OUT/'design.json').read_text())['parameters']);cs=np.sqrt(1-x[22]);cf=np.sqrt(x[22]);mu=np.exp(-.5*theta**2*ref[:,:,0]-theta*(cs*ref[:,:,2]+cf*ref[:,:,3]));fwd=[]
    for j,day in enumerate(days):
        v=mu[:,j];block=blocks(v[:,None]).ravel();fwd.append(dict(day=int(day),conditional_reference_stock_mean=float(v.mean()),batch_SE=float(block.std(ddof=1)/np.sqrt(8)),largest_mean_contribution_share=float(v.max()/v.sum())))
    pd.DataFrame(fwd).to_csv(out/'reference_forward_diagnostics.csv',index=False)
    s0=scores[(scores.direction=='mixed')&(scores.theta==0)&(scores.steps==128)].iloc[0];s1=scores[(scores.direction==direction)&(scores.theta==theta)&(scores.steps==128)].iloc[0];gain=100*(1-s1.loss/s0.loss);audit=json.loads((OUT/'implementation_audit.json').read_text());run=json.loads((OUT/'run.json').read_text())
    report=f'''# Stock/stress leverage scan — 28 October 2026

The experiment introduces negative stock/stress dependence while keeping the full variance-factor law fixed. It materially improves both price loss and native implied-volatility error. It does not yet give a satisfactory double calibration.

## What changed

For fixed stress direction (cs,cf), use

$$dW^{{stock}}_\\theta=\\sqrt{{1-\\theta^2}}dW^0-\\theta(c_s dB^s+c_f dB^f).$$

A positive stress shock raises variance and now tends to push the stock down. We tested theta = 0, .2, .4, .6, .8, .95 with mixed, slow-only and fast-only directions. Mixed uses (sqrt(1-r),sqrt(r)); the other directions isolate each factor. All OU dynamics, their joint law, M2 readout, floor and protected-candidate variance curve remain fixed. Original positive spike shift remains {audit['positive_spike_shift']:.9f}. There is no endogenous volatility/return feedback. The old stock/M2 correlations shrink by sqrt(1-theta²).

The variance-factor process and conditional future variance law are unchanged. Thus the VIX distribution, its options and forwards are protected algebraically. Existing finite-sample VIX forward error is inherited; this is not a new population-forward equality certificate. Baselines and selected models from earlier runs are not overwritten.

## Fresh numerical comparison

Selection was frozen before fresh expensive checking: **{direction}, theta={theta:g}**, based on the replay sample's fine-grid equal-maturity price loss. All 1751 eligible SPX quotes / all seven maturities remain in CSV; the original 1294 primary quotes determine scores.

| Metric, 128 steps/trading day | Theta 0 | Selected theta {theta:g} |
|---|---:|---:|
| Equal-maturity price interval loss | {s0.loss:.2f} | {s1.loss:.2f} |
| Equal-maturity IV MAE, percentage points | {s0.IV_MAE_pp:.3f} | {s1.IV_MAE_pp:.3f} |
| Primary prices inside bid–ask | {int(s0.inside)} / 1294 | {int(s1.inside)} / 1294 |
| Invalid primary IVs | {int(s0.invalid_IV)} | {int(s1.invalid_IV)} |
| Negative primary price estimates | {int(s0.negative)} | {int(s1.negative)} |

Price loss falls **{gain:.1f}%**. Approximate paired 95% batch-bootstrap loss-change interval: {s1.paired_loss_delta_CI95}. This interval covers sampling of expensive and cheap blocks; it excludes timestep, normalizer and market errors. Only eight independent batches underpin this approximation. IV MAE omits unavailable IVs, whose counts remain explicit. No signed estimate is clipped or smoothed.

Fresh full-grid scores are in [scores.csv](runs/m2_leverage_scan_001/fresh/scores.csv). Looking across that grid after checking is exploratory; it is not an independent validation of a newly selected winner. Some maturities remain visibly mismatched; inspect the seven maturity panels and maturity_metrics.csv rather than treating the aggregate gain as a good calibration.

## Computation and checks

Each factor stream is built once, then caches I=sum(v dt), U0=sum(sqrt(v)dW0), Us=sum(sqrt(v)dBs), Uf=sum(sqrt(v)dBf). For every theta/direction,

$$\\log S_T=-I/2+\\sqrt{{1-\\theta^2}}U_0-\\theta(c_sU_s+c_fU_f).$$

Replay uses 512 earlier paths; checking uses 1024 fresh full-model paths (512 antithetic pairs). Both reuse the same 65536-path independent auxiliary stress reference, with its separate batch uncertainty included. Thus the expensive checking streams are fresh, while the complete estimator is not wholly independent of screening. These are the same market quotes, not market out-of-sample data. Coupled 64/128 grids diagnose timestep sensitivity. Total scan wall time: {run['wall_seconds']:.1f} seconds, including compile, both factor samples, reference simulation, pricing all configurations, IV inversion and paired bootstrap.

Stress Brownian increments are sampled conditionally on OU innovations using their exact joint covariance, including the independent residual. Full replay at theta=0 matches stored stocks to {run['replay_theta0_max_absolute_stock_error']:.3g}. Independent scalar target/reference sufficient-statistic oracles agree to {audit['scalar_target_statistics_max_error']:.3g} and {audit['scalar_reference_statistics_max_error']:.3g}; standardized stress covariance error is {audit['maximum_standardized_stress_increment_covariance_error']:.3g}.

The old future-stress-dependent delta control is not reused. The new auxiliary call is integrated conditionally with log mean -I/2-theta(cs Us+cf Uf) and log variance (1-theta²)I. Add the paired full/reference payoff difference and subtract a deterministic strike-dependent Black delta times their stock difference. Each finite predictable log-Euler stock has mathematical mean one; no sample-forward renormalization is applied. Continuous-time martingality and tail convergence are separate questions, not established here. Stock/reference forward diagnostics are saved, not used to repair prices.

Existing Lean results prove coefficient normalization, a covariance-sign identity and an abstract discrete reconstruction identity. These numerical implementation checks do not turn them into a full simulator, VIX-invariance or continuous-martingale theorem.

## What to do next

Keep the mixed-leverage construction as a promising experimental extension. Before fitting additional parameters, improve independent reference/tail precision for the promising theta range and inspect remaining maturity-specific shape errors. Then a joint SPX/VIX fit can vary the previously authorized curve/readout parameters within explicit VIX protection. Do not introduce endogenous OU feedback or adopt a new production model from this exploratory scan alone.

Figures: [theta response](visualizations/leverage028/SPX-response.png), [all seven mixed-direction smiles](visualizations/leverage028/SPX-theta-fan.png), [stress-direction comparison](visualizations/leverage028/SPX-stress-directions.png). Detailed maturity, grid and reference-forward diagnostics accompany them.
'''
    (HERE/'LEVERAGE_SCAN_RESULTS.md').write_text(report)
    dump(out/'summary.json',dict(selected=sel,fresh_selected=s1.to_dict(),fresh_theta0=s0.to_dict(),price_loss_improvement_percent=gain,grid=grid,reference_forward=fwd))
    print('REPORT',gain,flush=True)

if __name__=='__main__':main()
