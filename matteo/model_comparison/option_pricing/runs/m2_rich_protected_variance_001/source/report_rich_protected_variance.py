"""Compare richer protected shapes against baseline and the previous correction."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .m2_rich_protected_variance import *
from .protected_curve_check import loss

PLOT=HERE/'visualizations/rich-protected-variance029'
NAMES=['Baseline','Previous protected correction','Richer protected correction']
COLORS=['#999999','#318b9c','#c06b2b']

def readable_report(report,selection,intervals,grid):
    selected=selection['new_shape'];state='tanh' if selected['kind']==0 else 'error-function'
    report=report.replace(json.dumps(selected),f"the {state} response on the existing approximately {[3.3,9.5,27.][selected['mode']]}-day OU mode, with width {selected['width']:g}")
    rows='\n'.join(f"| {NAMES[r['candidate']]} minus {NAMES[r['reference']]} | {r['loss_delta']:.2f} | [{r['approximate_CI95'][0]:.2f}, {r['approximate_CI95'][1]:.2f}] |" for r in intervals)
    report=report.replace(json.dumps(intervals),f"\n\n| Comparison | Loss change | Approximate 95% interval |\n|---|---:|---:|\n{rows}\n\n")
    rows='\n'.join(f"| {NAMES[r['candidate']]} | {r['coarse_loss']:.2f} | {r['fine_loss']:.2f} | {r['median_price_change_quote_scales']:.3f} | {r['max_price_change_quote_scales']:.3f} |" for r in grid)
    report=report.replace(json.dumps(grid),f"\n\n| Model | 128-step loss | 256-step loss | Median price change / quote scale | Maximum |\n|---|---:|---:|---:|---:|\n{rows}\n\n")
    for a,b in {'half-lives are3.329/9.481/27.004days':'half-lives are 3.329, 9.481 and 27.004 days','a10-day':'a 10-day','from0to1, then keep1':'from 0 to 1, then keep it at 1','approximately0.05008':'approximately 0.05008','gap0.05517points':'gap 0.05517 points','has512stock paths and26configurations':'uses 512 stock paths and 26 configurations','and24new choices':'and 24 new choices','before1024fresh':'before 1,024 fresh','existing65536-path':'existing 65,536-path','All1751eligible/1294primary':'All 1,751 eligible and 1,294 primary','seven expiries29/57/119/211/302/666/1030days':'seven expiries (29, 57, 119, 211, 302, 666 and 1,030 days)','All217VIX':'All 217 VIX','the256steps/trading-day':'the 256 steps per trading day','The29-day':'The 29-day','percent relative':'% relative','Each uses1024resamples':'Each uses 1,024 resamples','all13basis':'all 13 basis','Both128/256stock grids':'Both 128- and 256-step stock grids','states-4/-1/+1/+4':'states -4, -1, +1, +4','states-4/-1/0/+1/+4':'states -4, -1, 0, +1, +4','top10concentration':'top-10 concentration','factor or coefficient':'factor or coefficient'}.items():report=report.replace(a,b)
    return report

def main():
    PLOT.mkdir(exist_ok=False,parents=True);run=json.loads((OUT/'run.json').read_text());design=json.loads((OUT/'design.json').read_text());selection=json.loads((OUT/'selection.json').read_text());audit=json.loads((OUT/'implementation_audit.json').read_text());recurrence=json.loads((OUT/'recurrence_audit.json').read_text());rates=np.log(2)*365/np.array(design['half_life_days']);root=OUT/'fresh/SPX'
    scores=pd.read_csv(root/'SPX_scores.csv');fine=scores[scores.steps==256].set_index('candidate');frames=[pd.read_csv(root/f'SPX_prices_c{c}_256.csv') for c in range(3)];base=frames[0];days=sorted(base.expiry_days.unique());metrics=[]
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,day in zip(axes.flat,days):
        market=base.query('expiry_days==@day and priority_selected');ax.scatter(market.K_over_F,market.market_IV,s=10,color='black',label='Market')
        for c in range(3):
            q=frames[c].query('expiry_days==@day and priority_selected');ax.scatter(q.K_over_F,q.model_IV,s=9,color=COLORS[c],label=NAMES[c]);lo=q.loc[(q.K_over_F-.9).abs().idxmin()];hi=q.loc[(q.K_over_F-1.1).abs().idxmin()]
            metrics.append(dict(candidate=c,day=int(day),loss=loss(q.reset_index(drop=True),q.model_price_U.to_numpy()),IV_MAE_pp=float(np.nanmean(abs(q.model_IV-q.market_IV))),left_K_over_F=float(lo.K_over_F),right_K_over_F=float(hi.K_over_F),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV),stock_mean=float(q.stock_mean.iloc[0]),stock_mean_SE=float(q.stock_mean_SE.iloc[0]),invalid_IV=int(q.model_IV.isna().sum())))
        ax.set(title=f'{day} days',xlabel='Strike / forward',ylabel='Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=7);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,'1024 fresh stock paths\n256 steps / trading day\n\n'+selection['new_shape']['name']+'\n\nSame OU bank and stock driver\nSame positive lower bound f/2\nFour quoted VIX expiries protected\nBaseline VIX fit errors inherited',va='top',fontsize=9);fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(PLOT/f'SPX-rich-protected.{ext}',dpi=170)
    plt.close(fig);pd.DataFrame(metrics).to_csv(PLOT/'maturity_metrics.csv',index=False)
    cfg=selection['new_shape'];t=np.linspace(0,1030/365,5000);fig,axes=plt.subplots(2,1,figsize=(12,6))
    for c,index in enumerate([0,int(np.flatnonzero(selection['coefficients'][2])[0])]):
        cf=CONFIG[index];rate=rates[cf['mode']];q,qd=time_basis(t,rate,cf['width'],cf['kind']);label=NAMES[c+1];axes[0].plot(t*365,q,label=label);axes[1].plot(t*365,hvalue(t,-1,rate,cf['width'],cf['kind']),label=label)
    axes[0].set(ylabel='Potential q(t, -1) / phi(-1)',xlim=(0,220));axes[1].set(xlabel='Calendar days',ylabel='Correction h(t, U=-1)',xlim=(0,220))
    for ax in axes:
        for day in ENDS*365:ax.axvline(day,color='grey',alpha=.2)
        ax.grid(alpha=.15);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(PLOT/'local-scaling.png',dpi=170);plt.close(fig)
    # Paired full-model and independently resampled cheap reference batches.
    B=np.load(root/'SPX_blocks.npz');mask=base.priority_selected.to_numpy();pq=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();rng=np.random.default_rng(20261029133);intervals=[];grid=[]
    for candidate,reference in [(1,0),(2,0),(2,1)]:
        changes=[]
        for _ in range(1024):
            ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);cheap=B['cheap'][ic].mean(0);p=(B[f'c{candidate}_256'][ie].mean(0)+cheap-parity)[mask];p0=(B[f'c{reference}_256'][ie].mean(0)+cheap-parity)[mask];changes.append(loss(pq,p)-loss(pq,p0))
        intervals.append(dict(candidate=candidate,reference=reference,loss_delta=float(fine.loc[candidate,'loss']-fine.loc[reference,'loss']),approximate_CI95=np.quantile(changes,[.025,.975]).tolist()))
    for c in range(3):
        coarse=pd.read_csv(root/f'SPX_prices_c{c}_128.csv');delta=abs(frames[c].model_price_U-coarse.model_price_U)/base.price_scale_U
        grid.append(dict(candidate=c,coarse_loss=float(scores.query('candidate==@c and steps==128').loss.iloc[0]),fine_loss=float(fine.loc[c,'loss']),median_price_change_quote_scales=float(np.median(delta[mask])),max_price_change_quote_scales=float(np.max(delta[mask]))))
    # Selected correction conditional left sums on both actual stock grids.
    d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));k=Kernel(OUT/'rich.dylib');z,w=hermgauss(96);z=z*np.sqrt(2);w=w/np.sqrt(np.pi);discrete=[];rate=rates[cfg['mode']]
    for steps in [128,256]:
        idx=np.flatnonzero(d['coarse']) if steps==128 else np.arange(len(d['time']));tt=d['time'][idx];right=np.r_[tt[1:],d['time'][-1]+d['dt'][-1]]
        for day in [27,55,118,153]:
            T=day/365;end=T+30/365;active=(tt>=T-1e-13)&(tt<end-1e-13);times=tt[active];weights=np.minimum(right[active],end)-times;assert abs(times[0]-T)<1e-11
            for state in [-4.,-1.,0.,1.,4.]:
                lag=np.maximum(times-T,0);u=np.exp(-rate*lag[:,None])*state+np.sqrt(-np.expm1(-2*rate*lag[:,None]))*z;eh=k.probe(times[:,None],u,rate,cfg['width'],cfg['kind'])[...,2]@w;integral=float(weights@eh);discrete.append(dict(steps=steps,day=day,state=state,conditional_h_integral=integral,VIX_squared_change_points_squared=float(10000*.0072*integral/(30/365))))
    dump(OUT/'discrete_window_audit.json',discrete)
    # Independent finite differences of q in time and state verify generator sign.
    rng=np.random.default_rng(20261029134);fd=[]
    for config in CONFIG:
        rate=rates[config['mode']];width=config['width'];kind=config['kind'];t=rng.uniform(.001,1030/365,100);u=rng.uniform(-3,3,100);keep=np.min(abs(t[:,None]-np.r_[ENDS,ENDS[-1]+10/365]),axis=1)>1e-4;t=t[keep];u=u[keep];dt=1e-7;du=1e-4
        def potential(t,u):
            q,_=time_basis(t,rate,width,kind);return q*(np.tanh(u/width) if kind<=0 else erf(u/(np.sqrt(2)*width)))
        q=potential(t,u);dtime=(potential(t+dt,u)-potential(t-dt,u))/(2*dt);du1=(potential(t,u+du)-potential(t,u-du))/(2*du);du2=(potential(t,u+du)-2*q+potential(t,u-du))/du**2;expected=dtime+rate*(du2-u*du1);actual=k.probe(t,u,rate,width,kind)[:,2];err=float(abs(actual-expected).max());assert err<2e-6;fd.append(dict(config=config,max_generator_error=err))
    dump(OUT/'finite_difference_generator_audit.json',fd)
    stocks=np.load(OUT/'fresh/stocks.npz')['stocks'];tails=[]
    for c in range(3):
        for j,day in enumerate(days):
            q=base.query('expiry_days==@day and priority_selected');S=stocks[:,j,c,1]
            for side in ['put','call']:
                g=q[q.option_type==side];K=float(g.K_over_F.min() if side=='put' else g.K_over_F.max());pay=np.maximum(K-S,0) if side=='put' else np.maximum(S-K,0);total=pay.sum();tails.append(dict(candidate=c,day=int(day),side=side,K_over_F=K,raw_payoff_mean=float(pay.mean()),payoff_ESS=float(total**2/(pay@pay)) if total else 0.,top10_share=float(np.sort(pay)[-10:].sum()/total) if total else 0.))
    dump(OUT/'stock_tail_diagnostics.json',tails)
    gain=float(100*(1-fine.loc[2,'loss']/fine.loc[0,'loss']));extra=float(100*(1-fine.loc[2,'loss']/fine.loc[1,'loss']));summary=dict(fresh_scores=fine.reset_index().to_dict('records'),price_loss_reduction_percent=gain,extra_reduction_over_previous_percent=extra,paired_intervals=intervals,grid=grid,selection=selection,adopted=False,max_discrete_VIX_squared_change=float(max(abs(r['VIX_squared_change_points_squared']) for r in discrete)),max_finite_difference_generator_error=max(r['max_generator_error'] for r in fd));dump(PLOT/'summary.json',summary)
    m=pd.DataFrame(metrics).query('day==29').set_index('candidate');table='\n'.join(f"| {name} | {fine.loc[c,'loss']:.2f} | {fine.loc[c,'IV_MAE_pp']:.3f} | {m.loc[c,'model_right_minus_left_IV_pp']:.2f} | {int(fine.loc[c,'negative'])} |" for c,name in enumerate(NAMES));assessment='The richer shape improves aggregate SPX fit in this checking sample, but the observed downward skew remains far from market.' if gain>0 else 'The richer shape does not improve aggregate SPX fit in this checking sample. The missing skew remains.'
    report=r'''# Richer VIX-protected variance shapes — 29 October 2026

ASSESSMENT No candidate is adopted. The production model and previous runs are unchanged. [Seven SPX smiles](visualizations/rich-protected-variance029/SPX-rich-protected.png), [time scaling](visualizations/rich-protected-variance029/local-scaling.png), full run `runs/m2_rich_protected_variance_001`.

## What this tests

The previous protected correction used a single global normalization and a product envelope. This experiment tests a more efficient local time normalization and two smooth state shapes. It keeps the same M2 readout, cap, Q normalizer, strictly positive spike shift, both stresses, mixture, stock theta=.8 and forward-variance curve. No new bank, curve refit or feedback of returns/variance into the Gaussian OU dynamics is added. Stock paths are recomputed from fresh primitive draws; the underlying exogenous factor construction is reusable.

Let the existing centered modes obey $dU=-\kappa Udt+\sqrt{2\kappa}dB$, $U_0=0$. Their half-lives are3.329/9.481/27.004days. As before,

$$v_t=v^0_t+\frac f2\sum_j c_jh_j(t,U^j_t),\quad \sum_j|c_j|\le1,\quad f=0.0144,$$

$$q_j(t,u)=b_j(t)\phi_j(u),\qquad h_j=(\partial_t+\mathcal L_j)q_j=b'_j\phi_j+\kappa_jb_j(\phi_j''-u\phi_j').$$

Here $\phi$ is either $\tanh(u/s)$ or $\operatorname{erf}(u/(\sqrt2s))$, with widths $s=1,2$. The width changes how rapidly the response saturates as the Gaussian state moves away from zero; it adds no randomness.

## Local time scaling and positivity

Every quoted VIX-window endpoint is in $B=\{27,55,57,85,118,148,153,183\}/365$. On each successive interval $[l,r]$ take

$$a(t)=\sin^2(\pi(t-l)/(r-l)),\qquad A=\pi/(r-l).$$

Before the first endpoint, $a(t)=\cos^2(\pi t/(2T_1))$ with $A=\pi/(2T_1)$. After the last endpoint, use a10-day $\sin^2$ ramp from0to1, then keep1, with $A=\pi/(2\tau)$ and $\tau=10/365$. These are boundary extensions of a single calendar construction; no expiry is given its own fitted factor or coefficient. It remains dependent on the finite quoted calendar, as the first test was.

For either shape $|\phi|\le1$ and $|\phi''-u\phi'|\le C_s$, with the analytic bounds

$$C_s^{\tanh}=\frac12+\frac{4}{3\sqrt3s^2},\qquad C_s^{\operatorname{erf}}=\sqrt{\frac2{\pi e}}(1+s^{-2}).$$

On each interval use the constant

$$K=\frac{\kappa C_s}{2}+\sqrt{A^2+\frac{\kappa^2C_s^2}{4}},\qquad b(t)=a(t)/K.$$

To see the bound, write $a=(1-\cos2\theta)/2$, $a'=A\sin2\theta$. For fixed $u$, the numerator $a'\phi+a\mathcal L\phi$ is a constant plus a sine and cosine. Its absolute maximum is at most $|\mathcal L\phi|/2+\sqrt{A^2\phi^2+(\mathcal L\phi)^2/4}$, which is at most $K$. The cosine boundary segment has the same bound, and the flat extension satisfies it too. Thus $|h_j|\le1$ globally and $v\ge f/2=0.0072$. This is a proved analytic bound under the stated OU generator, not a bound inferred from observed paths. No variance clipping is used. It preserves the same conservative positive floor as the first correction, rather than lowering it again.

Although $K$ changes between intervals, both $a$ and $a'$ vanish at their joins. Therefore $b$ is continuously differentiable, and the potential has no jump term. A large local scale cannot be inserted without this endpoint check. The two possible state shapes and their generators are implemented directly; there is no post-pricing smile smoothing.

## VIX protection and its limits

For every quoted expiry $T_i=27,55,118,153$days and $\Delta=30/365$, $q(T_i,u)=q(T_i+\Delta,u)=0$ for all states. Conditional Dynkin gives

$$E_{T_i}\int_{T_i}^{T_i+\Delta}h(s,U_s)ds=E_{T_i}q(T_i+\Delta,U_{T_i+\Delta})-q(T_i,U_{T_i})=0.$$

This protects model VIX values state by state, hence all quoted-date VIX options and forwards, for the unchanged full Gaussian filtration. It inherits the baseline VIX fit: loss approximately0.05008 and maximum checking-forward gap0.05517points. No new claim of exact market/population forwards is made and no unchanged baseline exposure estimate is mislabelled as fresh simulation. Unquoted dates and the full VIX process law are not protected. Quoted-date joint VIX distributions are protected too, a stronger restriction than preserving only marginal distributions.

Both shapes are odd in the centered Gaussian state, so $E[h_t]=0$ and ideal $E[v_t]=\xi(t)$ is unchanged. Positivity and conditional cancellation are separate requirements; neither by itself guarantees a satisfactory SPX smile. No new continuous-stock martingale theorem or Lean proof is asserted.

## Selection and independent checking

The pilot has512stock paths and26configurations: baseline, the previous selected correction, and24new choices from three OU modes, two state shapes, two widths and both signs, at the same maximum L1 budget. This is a discrete shape search, not full coefficient optimization. The best new choice by all-primary fine-grid SPX interval loss was frozen before1024fresh stock paths: SELECTED. These target draws differ from the first experiment; all three checking configurations share draws for precise differences. The same existing65536-path auxiliary reference and market snapshot are reused, so this is not market out-of-sample validation.

All1751eligible/1294primary SPX quotes and seven expiries29/57/119/211/302/666/1030days remain. All217VIX quotes at four dates are protected. The table is on the256steps/trading-day grid. The29-day wing gap is implied volatility near strike/forward1.1 minus that near0.9, in percentage points; market gap is MARKETpp.

| Model | SPX price loss | IV mean absolute error (pp) | 29-day wing gap (pp) | Signed negative primary prices |
|---|---:|---:|---:|---:|
TABLE

The richer correction changes price loss by GAINpercent relative to baseline and EXTRApercent relative to the previous correction. Approximate paired eight-batch bootstrap intervals for candidate-minus-reference loss are INTERVALS. Each uses1024resamples of paired full-model batches with independent resampling of the common cheap-reference batches. They exclude deterministic grid error, fitting variability and market uncertainty. Signed negative prices are retained in the price score; unavailable IVs are reported and omitted only from the available-IV mean. IV means equally weight the seven maturities and are not a replacement for wing-by-wing inspection.

## Implementation and numerical evidence

Compiled/Python formula discrepancy FORMULA; conditional-window maximum residual CONDITIONAL across all13basis functions, four quoted dates, states-4/-1/+1/+4 and GH64/96; endpoint potential/derivative discrepancy ENDPOINT; baseline and previous-correction replay discrepancy REPLAY. Independent nonzero scalar stock recurrence error RECURRENCE. An additional independent finite-difference generator check has maximum error FD. These finite numerical checks support implementation correspondence; they do not establish a uniform error bound over all states.

Both128/256stock grids are retained; summary GRID. The selected correction's GH96 conditional left-sum audit on actual grids and states-4/-1/0/+1/+4 at all four dates has maximum correction-induced VIX-squared change DISCRETEpoints squared. Exact cancellation belongs to the continuous-time model; the discretized correction has small nonzero error. Baseline VIX quadrature is a separate inherited error. Paired option-price errors, native-IV/stock-mean diagnostics by maturity, and raw tail payoff ESS/top10concentration are saved. No sample stock-mean renormalization was used.

## Research conclusion

ASSESSMENT This test distinguishes richer state/time shape from adding another stochastic bank, with the same variance budget and positivity floor. The allowed corrections remain small compared with stressed baseline variance, and finite-date protection restricts their calendar shape. These are limitations of this tested family, not an impossibility theorem for M2 or all exogenous readouts. Do not repeat this exact shape grid as untested, and do not increase amplitudes past the established bound without a new positivity argument. Any next step should address the remaining wing-level mismatch explicitly rather than treating aggregate price improvement as a solved double calibration.
'''
    replacements={'ASSESSMENT':assessment,'SELECTED':json.dumps(selection['new_shape']),'MARKET':f"{m.loc[0,'market_right_minus_left_IV_pp']:.2f}",'TABLE':table,'GAIN':f'{gain:.2f}','EXTRA':f'{extra:.2f}','INTERVALS':json.dumps(intervals),'FORMULA':f"{audit['compiled_formula_max_error']:.3g}",'CONDITIONAL':f"{audit['max_conditional_integral_error']:.3g}",'ENDPOINT':f"{audit['endpoint_q_and_dtq_max_error']:.3g}",'REPLAY':f"{audit['baseline_and_previous_correction_replay_error']:.3g}",'RECURRENCE':f"{recurrence['maximum_scalar_stock_error']:.3g}",'FD':f"{summary['max_finite_difference_generator_error']:.3g}",'GRID':json.dumps(grid),'DISCRETE':f"{summary['max_discrete_VIX_squared_change']:.4g}"}
    for a,b in replacements.items():report=report.replace(a,b)
    report=readable_report(report,selection,intervals,grid)
    (HERE/'RICH_PROTECTED_VARIANCE_RESULTS.md').write_text(report)
    note=f'''## 29 October 2026 — richer protected variance shapes tested

Read [RICH_PROTECTED_VARIANCE_RESULTS.md](RICH_PROTECTED_VARIANCE_RESULTS.md), run `runs/m2_rich_protected_variance_001`, figures `visualizations/rich-protected-variance029`. Same frozen M2/OU bank/stresses/stock theta=.8/curve. Replace first product-envelope/global scale by interval sine-squared potentials b(t)phi(U), normalized locally by K=kappa*C/2+sqrt(A^2+(kappa*C)^2/4). At all quoted window endpoints b=b'=0, so jumps of the interval scale add no jump term. Tanh/erf shapes, widths1/2, existing modes7/6/5, both signs. Global |h|<=1 and L1 coefficient budget retain the same v>=f/2=.0072; no clipping or new bank. Conditional Dynkin still protects all217VIX quotes/four quoted dates state by state, inherits baseline errors, not exact market forwards. No unquoted/full-process invariance.

512pilotpaths/26choices; frozen best new choice {selection['new_shape']};1024fresh checking paths and reused65536path auxiliary. Fresh256step SPX losses baseline{fine.loc[0,'loss']:.2f},previous{fine.loc[1,'loss']:.2f},richer{fine.loc[2,'loss']:.2f}; richer gain{gain:.2f}%vsbaseline/{extra:.2f}%vsprevious. IVMAE baseline{fine.loc[0,'IV_MAE_pp']:.3f}/previous{fine.loc[1,'IV_MAE_pp']:.3f}/richer{fine.loc[2,'IV_MAE_pp']:.3f}pp.29day gaps{m.loc[0,'model_right_minus_left_IV_pp']:.2f}/{m.loc[1,'model_right_minus_left_IV_pp']:.2f}/{m.loc[2,'model_right_minus_left_IV_pp']:.2f}pp vsmarket{m.loc[0,'market_right_minus_left_IV_pp']:.2f}. {assessment} All1751eligible/1294primary/sevenSPXmaturities, signed estimates/invalid IVs retained. Paired intervals, grid sensitivity and tails saved. Compiled formula/conditional cancellation, old-model replay, independent scalar recurrence and finite-difference generator checked. No new Lean proof or full positivity certificate from numerical probes alone; positivity rests on the analytic bound in the report.

No candidate adopted; production model and old artifacts unchanged. Preserve this failed-or-partial shape search rather than repeating it. The following entries are historical.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        path=HERE/name;path.write_text(note+path.read_text())
    path=HERE/'formalization/CORRESPONDENCE.md';path.write_text(path.read_text()+'''\n## 29 October 2026: interval-scaled protected corrections\n\n`m2_rich_protected_variance.cpp/.py` preserve four conditional VIX-window integrals using piecewise C1 time potentials and bounded tanh/erf state potentials. Local analytic normalization bounds each generator correction by1 and the same L1 budget proves variance>=f/2. Endpoint zeros of both potential and time derivative justify joining interval scales without distributional jump terms. Numerical formula, endpoint, conditional quadrature, replay, independent stock recurrence, finite-difference generator and two-grid left-sum audits are in `../runs/m2_rich_protected_variance_001`. These checks are not new Lean proofs, a uniform pricing bound, exact market-forward calibration or a continuous-stock martingale theorem.\n''')
    for name in ['audit_rich_protected_variance.py','report_rich_protected_variance.py','RICH_PROTECTED_VARIANCE_RESULTS.md','asymmetric_readout.py','forward_leverage.py','leverage_scan.py','protected_curve_check.py','audit_leverage_scan.py']:shutil.copy2(HERE/name,OUT/'source'/name)
    shutil.copy2('/private/tmp/m2-rich-protected-001.log',OUT/'run.log');shutil.copy2('/private/tmp/m2-rich-recurrence.log',OUT/'recurrence.log')
    for source,h in json.loads((OUT/'source_hashes.json').read_text()).items():assert hashlib.sha256(Path(source).read_bytes()).hexdigest()==h,source
    counts={}
    for name in ['m2_leverage_scan_001','m2_asymmetric_readout_001','m2_targeted_skew_001','m2_forward_leverage_001','m2_slow_readout_001','m2_slow_readout_002','m2_memory_filter_001','m2_protected_variance_design_001','m2_protected_variance_001']:
        prior=OUT.parent/name;manifest=json.loads((prior/'artifact_manifest.json').read_text())
        for file,h in manifest.items():assert hashlib.sha256((prior/file).read_bytes()).hexdigest()==h,(name,file)
        counts[name]=len(manifest)
    run.update(status='completed pricing, audits and report; figures pending visual review',summary=summary,prior_artifact_hash_counts_verified=counts);dump(OUT/'run.json',run)
    dump(OUT/'artifact_manifest.json',{str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'})
    print('RICH REPORT COMPLETE',json.dumps(summary),flush=True)

if __name__=='__main__':main()
