"""All-maturity diagnostics and durable record for the protected correction."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .m2_protected_variance import *
from .protected_curve_check import loss

PLOT=HERE/'visualizations/protected-variance029'

def main():
    PLOT.mkdir(exist_ok=False,parents=True);run=json.loads((OUT/'run.json').read_text());design=json.loads((OUT/'design.json').read_text());selection=json.loads((OUT/'selection.json').read_text());audit=json.loads((OUT/'implementation_audit.json').read_text());recurrence=json.loads((OUT/'recurrence_audit.json').read_text())
    root=OUT/'fresh/SPX';scores=pd.read_csv(root/'SPX_scores.csv');fine=scores[scores.steps==256].set_index('candidate');frames=[pd.read_csv(root/f'SPX_prices_c{c}_256.csv') for c in range(2)];base=frames[0];days=sorted(base.expiry_days.unique());metrics=[]
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,day in zip(axes.flat,days):
        market=base[(base.expiry_days==day)&base.priority_selected];ax.scatter(market.K_over_F,market.market_IV,s=10,color='black',label='Market')
        for c,(label,color) in enumerate([('M2 baseline','#929292'),('Protected correction','#c26925')]):
            q=frames[c].query('expiry_days==@day and priority_selected');ax.scatter(q.K_over_F,q.model_IV,s=9,color=color,label=label)
            lo=q.loc[(q.K_over_F-.9).abs().idxmin()];hi=q.loc[(q.K_over_F-1.1).abs().idxmin()]
            metrics.append(dict(candidate=c,day=int(day),loss=loss(q.reset_index(drop=True),q.model_price_U.to_numpy()),IV_MAE_pp=float(np.nanmean(abs(q.model_IV-q.market_IV))),left_K_over_F=float(lo.K_over_F),right_K_over_F=float(hi.K_over_F),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV),stock_mean=float(q.stock_mean.iloc[0]),stock_mean_SE=float(q.stock_mean_SE.iloc[0]),invalid_IV=int(q.model_IV.isna().sum())))
        ax.set(title=f'{day} days',xlabel='Strike / forward',ylabel='Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=8);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'1024 fresh stock paths\n256 steps / trading day\n\nThree existing OU modes tested\nFrozen selection: 27-day half-life\nFull positive variance, no clipping\nAll four quoted VIX values protected\nBaseline VIX fitting errors inherited',va='top');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(PLOT/f'SPX-protected-variance.{ext}',dpi=170)
    plt.close(fig);pd.DataFrame(metrics).to_csv(PLOT/'maturity_metrics.csv',index=False)
    fig,axes=plt.subplots(2,1,figsize=(12,6));t=np.linspace(0,1030/365,4000);a,ap=envelope(t);axes[0].plot(t*365,a);axes[0].set(ylabel='Common envelope a(t)',xlim=(0,220));rates=np.log(2)*365/np.array(design['half_life_days'])
    for rate,half in zip(rates,design['half_life_days']):axes[1].plot(t*365,correction(t,-1,rate),label=f'{half:.1f}-day mode, U=-1')
    axes[1].set(xlabel='Calendar days',ylabel='Correction h(t, U)',xlim=(0,220));axes[1].legend()
    for ax in axes:
        for day in ENDS*365:ax.axvline(day,color='grey',alpha=.2)
        ax.grid(alpha=.15)
    fig.tight_layout();fig.savefig(PLOT/'correction-envelope.png',dpi=170);plt.close(fig)
    B=np.load(root/'SPX_blocks.npz');mask=base.priority_selected.to_numpy();pq=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();rng=np.random.default_rng(20261029044);changes=[]
    for _ in range(1024):
        ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);cheap=B['cheap'][ic].mean(0);p0=(B['c0_256'][ie].mean(0)+cheap-parity)[mask];p1=(B['c1_256'][ie].mean(0)+cheap-parity)[mask];changes.append(loss(pq,p1)-loss(pq,p0))
    interval=np.quantile(changes,[.025,.975]).tolist();grid=[]
    for c in range(2):
        coarse=pd.read_csv(root/f'SPX_prices_c{c}_128.csv');delta=abs(frames[c].model_price_U-coarse.model_price_U)/base.price_scale_U
        grid.append(dict(candidate=c,coarse_loss=float(scores.query('candidate==@c and steps==128').loss.iloc[0]),fine_loss=float(fine.loc[c,'loss']),median_price_change_quote_scales=float(np.median(delta[mask])),max_price_change_quote_scales=float(np.max(delta[mask]))))
    # Discretization of the correction itself on the two stock grids. These
    # conditional sums audit the sampled h term, not the old baseline VIX error.
    d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));k=Kernel(OUT/'protected.dylib');z,w=hermgauss(96);z=z*np.sqrt(2);w=w/np.sqrt(np.pi);discrete=[]
    rate=rates[2]
    for steps in [128,256]:
        idx=np.flatnonzero(d['coarse']) if steps==128 else np.arange(len(d['time']));tt=d['time'][idx];right=np.r_[tt[1:],d['time'][-1]+d['dt'][-1]]
        for day in [27,55,118,153]:
            T=day/365;end=T+30/365;active=(tt>=T-1e-13)&(tt<end-1e-13);times=tt[active];weights=np.minimum(right[active],end)-times;assert abs(times[0]-T)<1e-11
            for state in [-4.,-1.,0.,1.,4.]:
                lag=np.maximum(times-T,0);u=np.exp(-rate*lag[:,None])*state+np.sqrt(-np.expm1(-2*rate*lag[:,None]))*z
                eh=k.probe(times[:,None],u,rate)[...,2]@w;integral=float(weights@eh);discrete.append(dict(steps=steps,day=day,state=state,conditional_h_integral=integral,conditional_integrated_variance_change=float(.0072*integral),VIX_squared_change_points_squared=float(10000*.0072*integral/(30/365))))
    dump(OUT/'discrete_window_audit.json',discrete)
    stocks=np.load(OUT/'fresh/stocks.npz')['stocks'];tails=[]
    for c in range(2):
        for j,day in enumerate(days):
            q=base[(base.expiry_days==day)&base.priority_selected];S=stocks[:,j,c,1]
            for side in ['put','call']:
                g=q[q.option_type==side];K=float(g.K_over_F.min() if side=='put' else g.K_over_F.max());pay=np.maximum(K-S,0) if side=='put' else np.maximum(S-K,0);total=pay.sum()
                tails.append(dict(candidate=c,day=int(day),side=side,K_over_F=K,raw_payoff_mean=float(pay.mean()),payoff_ESS=float(total**2/(pay@pay)) if total else 0.,top10_share=float(np.sort(pay)[-10:].sum()/total) if total else 0.))
    dump(OUT/'stock_tail_diagnostics.json',tails)
    gain=float(100*(1-fine.loc[1,'loss']/fine.loc[0,'loss']));summary=dict(fresh_scores=fine.reset_index().to_dict('records'),price_loss_reduction_percent=gain,paired_loss_delta_CI95=interval,grid=grid,selection=selection,adopted=False,VIX_prices='inherited unchanged at four quoted expiries by conditional Dynkin identity; independently audited compiled correction, not recomputed unchanged baseline estimates',max_discrete_VIX_squared_change=float(max(abs(r['VIX_squared_change_points_squared']) for r in discrete)))
    dump(PLOT/'summary.json',summary);m=pd.DataFrame(metrics).query('day==29').set_index('candidate');table='\n'.join(f"| {label} | {fine.loc[c,'loss']:.2f} | {fine.loc[c,'IV_MAE_pp']:.3f} | {m.loc[c,'model_right_minus_left_IV_pp']:.2f} |" for c,label in enumerate(['Baseline','Protected correction']))
    report=fr'''# VIX-protected additive M2 variance correction — 29 October 2026

The first test produces a modest SPX improvement while analytically protecting the four quoted VIX expiry values. It does not solve the missing SPX skew. No candidate is adopted; the production model and all previous runs are unchanged. Figures: [all seven SPX smiles](visualizations/protected-variance029/SPX-protected-variance.png), [correction envelope](visualizations/protected-variance029/correction-envelope.png). Full artifacts: `runs/m2_protected_variance_001`.

## What changed, in plain language

Keep the existing M2 component, both stress factors, their mixture, cap, normalizers, positive spike shift, forward-variance curve and stock driver. Add a small signed correction directly to variance. It can redistribute volatility along stock paths, but its conditional expected total over each quoted VIX window is zero. This resembles moving water within a container: its distribution can change while its expected total remains fixed. That cancellation is conditional on the entire current model state, stronger than cancelling only in the unconditional average.

Let the baseline variance be $v^0_t$. Existing centered OU modes satisfy

$$dU^j_t=-\kappa_j U^j_tdt+\sqrt{{2\kappa_j}}dB^j_t,\quad U^j_0=0.$$

Use the same geometric feedback bank, modes 7,6,5 with half-lives {design['half_life_days']}. There is no new bank or feedback of returns or variance into these dynamics. Stock paths are recomputed; Gaussian factors/random draws remain reusable.

$$v_t=v^0_t+\frac f2\sum_j c_jh_j(t,U^j_t),\qquad \sum_j|c_j|\le1,\quad f=0.0144.$$

$$a(t)=\prod_{{b\in B}}\tanh^2((t-b)/\tau),\quad B=\{{27,55,57,85,118,148,153,183\}}/365,\quad \tau=10/365.$$

$$q_j(t,u)=a(t)\tanh(u)/K_j,\quad K_j=1/\tau+\kappa_j[1/2+4/(3\sqrt3)],$$

$$h_j=(\partial_t+\mathcal L_j)q_j=\frac{{a'\tanh u-\kappa_j a(u+2\tanh u)\operatorname{{sech}}^2u}}{{K_j}}.$$

The common calendar envelope treats every quoted expiry by the same rule. This is deliberately an initial finite-expiry protection test, not a universal calendar-independent VIX model. Smoothness here does not mean a small effect is guaranteed or that smile fit follows automatically.

## Why VIX is protected

For each quoted $T_i=27,55,118,153$ days and $\Delta=30/365$, both endpoint functions vanish. Conditional Dynkin's identity gives

$$E_{{T_i}}\int_{{T_i}}^{{T_i+\Delta}}h_j(s,U^j_s)ds
=E_{{T_i}}q_j(T_i+\Delta,U^j_{{T_i+\Delta}})-q_j(T_i,U^j_{{T_i}})=0.$$

Thus $\mathrm{{VIX}}^2_{{T_i}}=10000\Delta^{{-1}}E_{{T_i}}\int_{{T_i}}^{{T_i+\Delta}}v_sds$ is unchanged state by state for the unchanged full Gaussian filtration. All VIX strikes and forwards at these four dates inherit their baseline values. Baseline VIX loss is approximately0.05008, with maximum checking-forward gap0.05517points: protection does not newly make those baseline market errors zero. Baseline VIX estimates are reused under this identity, not presented as a new independent calibration. Unquoted dates and the full VIX process law are not protected. At the protected dates even their joint law is unchanged, so this is more restrictive than the paper's approach of fixing only expiry marginal distributions.

Oddness of $h_j$ in centered Gaussian $U^j$ also preserves $E[v_t]=\xi(t)$ in the ideal normalized baseline. A global bound $|h_j|\le1$ and the amplitude budget imply $v_t\ge f/2=0.0072$: there is no clipping. This retains strict positivity but lowers the guaranteed floor from the original $f$; it does not retain the original exact hard floor. Since the correction is additive, it is not multiplied by $\xi-f$ and is active even on baseline floor intervals.

## Actual experiment and results

The pilot evaluated baseline plus twelve axis configurations: each of three modes, both signs, and strengths0.5/1 under the budget. There are512pilot stock paths. The minimum all-primary SPX interval-loss configuration was frozen before1024fresh stock paths: $c=(0,0,1)$, the27-day mode at its allowed amplitude. No curve refit or VIX-forward fit was performed. This is a small discrete parameter search, not optimization over all combinations or a proof of optimality.

Keep all1751eligible SPX quotes,1294primary quotes and seven expiries29/57/119/211/302/666/1030days. Four VIX expiries and217quotes remain protected. Signed prices and all native-IV failures remain in artifacts. The table uses fresh paths and256steps per trading day. The last column is IV near strike/forward1.1 minus IV near0.9 at29days, in percentage points: more negative means a stronger downward skew.

| Model | SPX price loss | IV mean absolute error (pp) | 29-day right minus left IV (pp) |
|---|---:|---:|---:|
{table}

The market29-day gap is{m.loc[0,'market_right_minus_left_IV_pp']:.2f}pp. Price loss improves{gain:.2f}%. The approximate paired1024-resample eight-batch95% bootstrap interval for candidate-minus-baseline loss is[{interval[0]:.2f},{interval[1]:.2f}]. It resamples shared full-model batches and independently resamples the shared65536-path auxiliary-reference batches. The new target simulations are fresh; the auxiliary reference and market snapshot are reused. These are sampling diagnostics, not a market out-of-sample certificate, and exclude deterministic error and fitting uncertainty.

The result separates two questions: this construction can change SPX prices while protecting quoted-date VIX; this first conservative correction is too weak to recover the observed SPX skew. A small aggregate loss improvement is not a satisfactory double calibration. The conservative bound constrains amplitude globally even where the baseline variance has more room. Also, the common envelope is suppressed near the densely spaced early VIX endpoints. Those are plausible limitations of this construction, not numerically established explanations or a no-go theorem for exogenous M2.

## Verification and limits

Compiled/Python formula maximum discrepancy{audit['compiled_formula_max_error']:.3g}. Conditional cancellation probes cover three modes, four dates, five conditioning states and GH64/96; maximum integral residual{audit['max_conditional_integral_error']:.3g}. Zero-amplitude compiled stock replay error{audit['zero_amplitude_stock_replay_max_error']:.3g}. Independent scalar recurrence on a97-step synthetic clock spanning envelope features, five configurations, both antithetic partners and both grids gives maximum stock discrepancy{recurrence['maximum_scalar_stock_error']:.3g}. These are numerical correspondence checks, alongside the stated analytic identity, with no new Lean proof or continuous-stock martingale theorem.

Both128/256stock grids are retained. Grid diagnostics are in `visualizations/protected-variance029/summary.json`. A separate GH96 conditional left-sum audit on both actual grids, selected mode and five states at all four dates is in `discrete_window_audit.json`; its largest correction-induced VIX-squared change is{summary['max_discrete_VIX_squared_change']:.6g}VIX-points squared. Exact cancellation is a continuous-time claim; the grid sum has small nonzero error. This audit concerns the new correction, not the baseline VIX quadrature. No finite probe establishes a uniform discretization bound for every state.

Per-maturity native-IV/skew/stock-mean diagnostics, paired price errors and raw tail payoff concentration are saved. No stock-mean renormalization or post-pricing smile smoothing was used. The old valid conditional-lognormal auxiliary reference and fixed deterministic delta controls remain unchanged.

## Research conclusion

Retain the protected-correction idea as a demonstrated way to separate part of SPX path dependence from four quoted-date VIX constraints, but do not adopt this first bounded-tanh specification as a solution. Do not repeat this exact mode/strength grid as untested. Before adding banks, a useful next design would investigate better use of the available positivity budget and less restrictive time/state bases while preserving the explicit conditional cancellation. Any enlarged correction requires a new positivity argument and implementation checks; simply increasing the coefficient beyond the tested bound is not justified by this experiment.
'''
    (HERE/'PROTECTED_VARIANCE_RESULTS.md').write_text(report)
    note=f'''## 29 October 2026 — VIX-protected additive variance correction tested

Read [PROTECTED_VARIANCE_RESULTS.md](PROTECTED_VARIANCE_RESULTS.md), run `runs/m2_protected_variance_001`, figures `visualizations/protected-variance029`. Retain full M2, stress mixture, curve, cap, OU dynamics and theta=.8. Add v=v0+(f/2)sum(c_j h_j), h_j=(partial_t+L_j)[a(t)tanh(U_j)/K_j], common envelope zero at all eight quoted VIX-window endpoints. Three existing modes7/6/5 have half-lives3.329/9.481/27.004days. L1(c)<=1 proves v>=f/2=.0072 without clipping; original exact floor f=.0144 is not retained. Centered odd correction preserves ideal mean variance. Conditional Dynkin identity protects all quoted-date VIX values/strikes/forwards state by state; inherits baseline VIXloss.05008/maxforwardgap.05517, does not newly match market forwards exactly, and does not protect unquoted VIX dates/full process law.

512pilotpaths, baseline+12axis configurations, selection frozen at27day mode/c=+1 before1024fresh paths. All1751eligible/1294primary SPX/sevenmaturities,217VIX/fourdates, signed estimates retained. Fresh256step price loss{fine.loc[0,'loss']:.2f}->{fine.loc[1,'loss']:.2f} ({gain:.2f}% improvement), IVMAE{fine.loc[0,'IV_MAE_pp']:.3f}->{fine.loc[1,'IV_MAE_pp']:.3f}pp;29day skew{m.loc[0,'model_right_minus_left_IV_pp']:.2f}->{m.loc[1,'model_right_minus_left_IV_pp']:.2f}pp versusmarket{m.loc[0,'market_right_minus_left_IV_pp']:.2f}. Paired bootstrap loss-changeCI[{interval[0]:.2f},{interval[1]:.2f}], excludes deterministic/model error. Compiled formula, all-mode conditional cancellation, zero-correction replay and independent nonzero scalar recurrence audited; stock128/256refinement and discrete-window audit saved. No new Lean theorem. Fresh targets use reused65536path auxiliary, same market snapshot.

Conclusion: modest SPX improvement with explicit quoted VIX protection, missing skew remains. No candidate adopted or production/old run changes. Do not repeat this exact test as untested; investigate positivity-budget efficiency and richer protected bases before new banks. New amplitudes require a fresh global positivity argument. The following proposal entry is historical and has now been tested.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        path=HERE/name;path.write_text(note+path.read_text())
    path=HERE/'formalization/CORRESPONDENCE.md';path.write_text(path.read_text()+'''\n## 29 October 2026: quoted-VIX-protected variance correction\n\n`m2_protected_variance.cpp/.py` use a bounded additive generator correction on existing centered OU modes. Conditional Dynkin cancellation over four quoted windows is an analytic argument under the unchanged filtration and integrability assumptions. It protects model values, not exact market forwards. `runs/m2_protected_variance_001` contains compiled/Python formula checks, conditional quadrature probes, zero-amplitude replay, independent signed-correction stock recurrence and two-grid conditional left-sum audits. These finite diagnostics are not new Lean proofs, a uniform numerical error bound, a stock martingale theorem or an impossibility claim. Positivity uses a lower bound f/2 rather than the original exact f.\n''')
    sources=['audit_protected_variance.py','report_protected_variance.py','PROTECTED_VARIANCE_RESULTS.md','asymmetric_readout.py','forward_leverage.py','leverage_scan.py','protected_curve_check.py','audit_leverage_scan.py']
    for name in sources:shutil.copy2(HERE/name,OUT/'source'/name)
    shutil.copy2('/private/tmp/m2-protected-variance-001.log',OUT/'run.log')
    for source,h in json.loads((OUT/'source_hashes.json').read_text()).items():assert hashlib.sha256(Path(source).read_bytes()).hexdigest()==h,source
    counts={}
    for name in ['m2_leverage_scan_001','m2_asymmetric_readout_001','m2_targeted_skew_001','m2_forward_leverage_001','m2_slow_readout_001','m2_slow_readout_002','m2_memory_filter_001','m2_protected_variance_design_001']:
        prior=OUT.parent/name;manifest=json.loads((prior/'artifact_manifest.json').read_text())
        for file,h in manifest.items():assert hashlib.sha256((prior/file).read_bytes()).hexdigest()==h,(name,file)
        counts[name]=len(manifest)
    run.update(status='completed pricing, implementation audits and documented conclusion; figures pending visual review',adopted=False,summary=summary,prior_artifact_hash_counts_verified=counts);dump(OUT/'run.json',run)
    dump(OUT/'artifact_manifest.json',{str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'})
    print('PROTECTED VARIANCE REPORT COMPLETE',json.dumps(summary),flush=True)

if __name__=='__main__':main()
