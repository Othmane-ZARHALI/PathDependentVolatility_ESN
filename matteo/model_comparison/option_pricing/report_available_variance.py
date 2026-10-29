"""Report the stress-variance-budget experiment and preserve its finite scope."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .m2_available_variance import *
from .protected_curve_check import loss

PLOT=HERE/'visualizations/available-variance029'
NAMES=['Baseline','Previous rich correction','Stress-budget correction']
COLORS=['#999999','#318b9c','#c06b2b']

def main():
    PLOT.mkdir(exist_ok=False,parents=True);run=json.loads((OUT/'run.json').read_text());design=json.loads((OUT/'design.json').read_text());sel=json.loads((OUT/'selection.json').read_text());audit=json.loads((OUT/'implementation_audit.json').read_text());rec=json.loads((OUT/'recurrence_audit.json').read_text());rates=np.log(2)*365/np.array(design['half_life_days']);x=arr(design['baseline']);sr,bb=stress_parameters(x);cfg=sel['configs'][2];mode,li,eta,_=cfg;mode=int(mode);li=int(li);lam=LAM[li];rate=rates[mode];f=(.05+.07*x[19])**2;root=OUT/'fresh/SPX';scores=pd.read_csv(root/'SPX_scores.csv');fine=scores.query('steps==256').set_index('candidate');frames=[pd.read_csv(root/f'SPX_prices_c{c}_256.csv') for c in range(3)];base=frames[0];days=sorted(base.expiry_days.unique());metrics=[]
    fig,axes=plt.subplots(2,4,figsize=(16,8))
    for ax,day in zip(axes.flat,days):
        market=base.query('expiry_days==@day and priority_selected');ax.scatter(market.K_over_F,market.market_IV,s=10,color='black',label='Market')
        for c in range(3):
            q=frames[c].query('expiry_days==@day and priority_selected');ax.scatter(q.K_over_F,q.model_IV,s=9,color=COLORS[c],label=NAMES[c]);lo=q.loc[(q.K_over_F-.9).abs().idxmin()];hi=q.loc[(q.K_over_F-1.1).abs().idxmin()];atm=q.loc[(q.K_over_F-1).abs().idxmin()]
            metrics.append(dict(candidate=c,day=int(day),loss=loss(q.reset_index(drop=True),q.model_price_U.to_numpy()),IV_MAE_pp=float(np.nanmean(abs(q.model_IV-q.market_IV))),left_K_over_F=float(lo.K_over_F),right_K_over_F=float(hi.K_over_F),left_IV_pp=float(lo.model_IV),market_left_IV_pp=float(lo.market_IV),right_IV_pp=float(hi.model_IV),market_right_IV_pp=float(hi.market_IV),ATM_IV_pp=float(atm.model_IV),model_right_minus_left_IV_pp=float(hi.model_IV-lo.model_IV),market_right_minus_left_IV_pp=float(hi.market_IV-lo.market_IV),stock_mean=float(q.stock_mean.iloc[0]),stock_mean_SE=float(q.stock_mean_SE.iloc[0]),invalid_IV=int(q.model_IV.isna().sum())))
        ax.set(title=f'{day} days',xlabel='Strike / forward',ylabel='Implied volatility (%)');ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=7);axes.flat[-1].axis('off');axes.flat[-1].text(0,1,f'1024 fresh stock paths\n256 steps / trading day\n\nExisting {design["half_life_days"][mode]:.2f}-day OU mode\nDamping lambda={lam:g}, eta={eta:g}\n\nStress-scaled generator correction\nOriginal floor f retained\nFour quoted VIX expiries protected\nBaseline VIX fit errors inherited',va='top',fontsize=9);fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(PLOT/f'SPX-available-variance.{ext}',dpi=170)
    plt.close(fig);pd.DataFrame(metrics).to_csv(PLOT/'maturity_metrics.csv',index=False)
    B=np.load(root/'SPX_blocks.npz');mask=base.priority_selected.to_numpy();pq=base[mask].reset_index(drop=True);parity=base.option_type.eq('put').to_numpy()*(base.forward-base.strike).to_numpy();rng=np.random.default_rng(20261029277);intervals=[];grid=[]
    for c,r in [(1,0),(2,0),(2,1)]:
        delta=[]
        for _ in range(1024):
            ie=rng.integers(8,size=8);ic=rng.integers(8,size=8);cheap=B['cheap'][ic].mean(0);p=(B[f'c{c}_256'][ie].mean(0)+cheap-parity)[mask];p0=(B[f'c{r}_256'][ie].mean(0)+cheap-parity)[mask];delta.append(loss(pq,p)-loss(pq,p0))
        intervals.append(dict(candidate=c,reference=r,loss_delta=float(fine.loc[c,'loss']-fine.loc[r,'loss']),approximate_CI95=np.quantile(delta,[.025,.975]).tolist()))
    for c in range(3):
        coarse=pd.read_csv(root/f'SPX_prices_c{c}_128.csv');delta=abs(frames[c].model_price_U-coarse.model_price_U)/base.price_scale_U;grid.append(dict(candidate=c,coarse_loss=float(scores.query('candidate==@c and steps==128').loss.iloc[0]),fine_loss=float(fine.loc[c,'loss']),median_price_change_quote_scales=float(np.median(delta[mask])),max_price_change_quote_scales=float(np.max(delta[mask]))))
    # Conditional discretization check uses the independently audited analytic
    # Gaussian oracle, including the actual fitted curve coefficient.
    d=dict(np.load(HERE/'runs/m2_targeted_skew_001/grid_refinement/grid.npz'));discrete=[]
    for steps in [128,256]:
        idx=np.flatnonzero(d['coarse']) if steps==128 else np.arange(len(d['time']));tt=d['time'][idx];right=np.r_[tt[1:],d['time'][-1]+d['dt'][-1]]
        for day in [27,55,118,153]:
            T=day/365;end=T+30/365;active=(tt>=T-1e-13)&(tt<end-1e-13);times=tt[active];weights=np.minimum(right[active],end)-times;assert abs(times[0]-T)<1e-11;segments=np.minimum(np.searchsorted(ENDSAV,times,side='right'),15)
            for U,Z in [(-3.,-2.),(-1.,0.),(1.,2.),(3.,5.)]:
                for a,b,share in zip(sr,bb,[1-x[22],x[22]]):
                    vals=np.array([conditional(t,T,U,Z,rate,a,b,lam) for t in times]);I=float(weights@((.2-f)*x[segments]*x[17]*share*eta*vals));initial=np.exp(b*Z-.5*b*b*(-np.expm1(-2*a*T)));discrete.append(dict(steps=steps,day=day,U=U,Z=Z,stress_rate=float(a),conditional_integrated_variance_change=I,relative_to_initial_stress_and_max_curve=I/(initial*(.2-f)*max(x[:16])*x[17]*share),VIX_squared_change_points_squared=10000*I/(30/365)))
    dump(OUT/'discrete_window_audit.json',discrete)
    # Typical Q-state probes measure how much variance the new budget uses.
    budget=[]
    for day in [7,20,29,57,119,211,666]:
        t=day/365;seg=min(np.searchsorted(ENDSAV,t,side='right'),15);n=32768;U=np.sqrt(-np.expm1(-2*rate*t))*rng.normal(size=n);Z=np.sqrt(-np.expm1(-2*sr*t))[:,None]*rng.normal(size=(2,n));parts=[scalar(t,U,Z[j],rate,sr[j],bb[j],lam) for j in range(2)];h=(.2-f)*x[seg]*x[17]*eta*((1-x[22])*parts[0][:,1]+x[22]*parts[1][:,1]);avail=(.2-f)*x[seg]*x[17]*((1-x[22])*parts[0][:,2]+x[22]*parts[1][:,2]);assert np.all(abs(h)<=avail+1e-12);budget.append(dict(day=day,mean_absolute_correction=float(abs(h).mean()),p99_absolute_correction=float(np.quantile(abs(h),.99)),fraction_exceeding_previous_absolute_budget=float(np.mean(abs(h)>f/2)),max_observed_correction_over_stress_budget=float(np.max(np.divide(abs(h),avail,out=np.zeros_like(h),where=avail>0))),scope='finite direct-Q state probes, not rare-tail convergence or a uniform bound'))
    dump(OUT/'variance_budget_probes.json',budget)
    stocks=np.load(OUT/'fresh/stocks.npz')['stocks'];tails=[]
    for c in range(3):
        for j,day in enumerate(days):
            q=base.query('expiry_days==@day and priority_selected');S=stocks[:,j,c,1]
            for side in ['put','call']:
                g=q[q.option_type==side];K=float(g.K_over_F.min() if side=='put' else g.K_over_F.max());pay=np.maximum(K-S,0) if side=='put' else np.maximum(S-K,0);total=pay.sum();tails.append(dict(candidate=c,day=int(day),side=side,K_over_F=K,raw_payoff_mean=float(pay.mean()),payoff_ESS=float(total**2/(pay@pay)) if total else 0.,top10_share=float(np.sort(pay)[-10:].sum()/total) if total else 0.))
    dump(OUT/'stock_tail_diagnostics.json',tails)
    gain=float(100*(1-fine.loc[2,'loss']/fine.loc[0,'loss']));extra=float(100*(1-fine.loc[2,'loss']/fine.loc[1,'loss']));summary=dict(fresh_scores=fine.reset_index().to_dict('records'),price_loss_reduction_percent=gain,extra_reduction_over_previous_percent=extra,paired_intervals=intervals,grid=grid,selection=sel,adopted=False,max_discrete_window_error_scaled=float(max(abs(r['relative_to_initial_stress_and_max_curve']) for r in discrete)));dump(PLOT/'summary.json',summary)
    m=pd.DataFrame(metrics).query('day==29').set_index('candidate');table='\n'.join(f"| {name} | {fine.loc[c,'loss']:.2f} | {fine.loc[c,'IV_MAE_pp']:.3f} | {m.loc[c,'left_IV_pp']:.2f} | {m.loc[c,'right_IV_pp']:.2f} | {m.loc[c,'model_right_minus_left_IV_pp']:.2f} |" for c,name in enumerate(NAMES));cis='\n'.join(f"| {NAMES[r['candidate']]} minus {NAMES[r['reference']]} | {r['loss_delta']:.2f} | [{r['approximate_CI95'][0]:.2f}, {r['approximate_CI95'][1]:.2f}] |" for r in intervals)
    theory=r'''# Using existing stress variance as the protected-correction budget — 29 October 2026

## Question and construction

The earlier additive corrections were limited by a fixed variance budget f/2, regardless of how high current variance was. Here the budget grows with the existing stress contribution. This retains the full exogenous OU bank, old capped/Q-normalized M2 readout, positive spike shift, both stress factors, their mixture, stock driver and fitted forward-variance curve. It does not add actual-return or volatility feedback to the OU dynamics. There is no new bank or market-data refit.

Write the baseline as

$$v^0_t=f+c_M(t)R_t+c_s(t)A^s_t+c_f(t)A^f_t,$$

where $c_M=(\xi-f)(1-w)$, $c_s=(\xi-f)w(1-r)$, $c_f=(\xi-f)wr$, all nonnegative. Each stress follows $dZ=-aZdt+\sqrt{2a}dB$, $Z_0=0$, with $A(t,Z)=\exp(bZ-b^2Q(t)/2)$ and $Q(t)=1-e^{-2at}$. The selected pre-existing mode $U$ is driven independently of these two stress factors. Their stock correlations remain unchanged.

Simply setting h=current variance times the old correction would not generally retain cancellation. Instead scale the potential before applying the full generator. For one stress component define

$$y=Z-bQ(t),\quad g(t,Z)=e^{-\lambda y^2},\quad\phi(U)=\operatorname{erf}(U/\sqrt2),$$

$$P_j(t,U,Z)=c_j(t)\beta_j(t)\phi(U)A(t,Z)g(t,Z),\quad h_j=(\partial_t+\mathcal L)P_j.$$

The Gaussian factor g keeps the derivative terms bounded relative to available stress variance. Its center bQ(t) is the mean of Z under weighting by A, so it follows the states that contribute most to mean stress variance. It is a smooth readout weight, not a new state or clipping rule.

Put $D_\lambda(t,y)=a[(4\lambda^2+2\lambda)y^2-b(1+2\lambda Q(t))y-2\lambda]$. Direct differentiation gives

$$\frac{(\partial_t+\mathcal L_Z)(Ag)}{A}=gD_\lambda,$$

$$h_j=c_j A g\,[\beta'_j\phi+\beta_j\mathcal L_U\phi+\beta_j\phi D_\lambda],\qquad\mathcal L_U\phi=-2\kappa U\sqrt{2/\pi}e^{-U^2/2}.$$

There is no mixed U/Z derivative term because the existing OU mode and stress driver are independent. The date-dependent coefficients c_j are constant between fitted curve breakpoints. The potential and its first time derivative are set to zero at every curve breakpoint as well as every quoted VIX-window endpoint. Thus coefficient jumps introduce no omitted jump term. Where xi=f, c_s=c_f=0 and this correction is zero; this differs from the earlier f/2 correction, which was active even on floor intervals.

## Global positivity bound

We have $|\phi|\le1$, $|\mathcal L_U\phi|\le\kappa C$ with $C=2\sqrt{2/(\pi e)}$. Since $Q(t)\in[0,1]$ and b>0,

$$|gD_\lambda|\le D_*:=a\left[\frac{4\lambda+2}{e}+\frac{b(1+2\lambda)}{\sqrt{2\lambda e}}+2\lambda\right].$$

This uses $|y|e^{-\lambda y^2}\le1/\sqrt{2\lambda e}$ and $y^2e^{-\lambda y^2}\le1/(\lambda e)$. On each interval between endpoints use the same sine-squared clock as in the preceding experiment, with $a_0=\sin^2(\pi(t-l)/(r-l))$ and $A_0=\pi/(r-l)$. Before the first endpoint use a cosine-squared boundary segment. After the last use a10-day sine-squared ramp then a flat extension. Set

$$B_j=\kappa C+D_{*,j},\quad K_j=B_j/2+\sqrt{A_0^2+B_j^2/4},\quad\beta_j=a_0/K_j.$$

The same sine/cosine amplitude bound used previously now gives $|h_j|\le c_jA^j$. Consequently, with

$$v_t=v^0_t+\eta(h_s+h_f),\qquad |\eta|\le1,$$

$$v_t\ge f+c_MR+(1-|\eta|)(c_sA^s+c_fA^f)\ge f.$$

This retains the ORIGINAL floor f=0.0144, stronger than the f/2 guarantee in the previous experiment. It can exceed the old absolute correction budget f/2 in high-variance states. The inequality is analytic, not inferred from simulated minima. It covers this stress-scaled construction, not arbitrary multiplication of corrections by variance or the entire M2 component. No variance clipping is used. Gaussian damping makes the relevant potentials and derivatives integrable on finite horizons.

## Conditional VIX protection and mean variance

All potentials vanish at both ends of all four quoted VIX windows. For the unchanged full Gaussian filtration,

$$E_{T_i}\int_{T_i}^{T_i+\Delta}(h_s+h_f)dt=0,\quad T_i=27,55,118,153\text{ days},\quad\Delta=30/365.$$

Conditional Dynkin's identity therefore protects model VIX values state by state, hence all217quoted VIX options and four forwards at these dates. Baseline VIX loss approximately0.05008 and max checking-forward gap0.05517points are inherited, not newly corrected or certified zero. Unquoted dates and the full VIX process law are not protected. Quoted-date joint VIX distributions are unchanged too. Oddness in the centered Gaussian U and independence from stresses give E[h]=0, retaining ideal E[v_t]=xi(t). Finite numerical normalization errors of the original M2 remain.

## Numerical translation and evidence scope

The compiled formula is compared against independent Python algebra, finite differences of the full time/U/Z potential, and a Gaussian conditional-expectation oracle. For the latter, exponential tilting and Gaussian damping give E[Ag] and E[AgD] analytically; E[erf(U/sqrt2)] and E[L_U phi] are analytic too. The product of these conditional moments is checked against compiled two-dimensional GH64 quadrature. Time integration splits at every curve endpoint and tests all three existing modes, both damping values, both stress factors, four VIX dates and four pairs of initial U/Z states. The actual fitted curve coefficient is included. Conditional residuals are scaled by initial stress and maximum curve level, since raw errors at extreme stress states scale with A. The analytic identity supplies the general cancellation claim; finite probes do not certify a uniform numerical bound.

The new construction is compared with the previously selected rich correction. This is not an ablation changing only amplitude: the calendar now includes every curve breakpoint, the stress readout is damped, and the new family retains the stricter original floor. Keep these differences when interpreting the outcome. No new Lean theorem, continuous-stock martingale theorem or universal model-family impossibility is asserted.
'''
    result=f'''
## Experiment and outcome

512 pilot paths test38configurations: baseline, the previous rich correction, and36new choices from three existing modes, damping lambda0.125/0.5, strengths0.25/0.5/1 and both signs. The best new primary-SPX interval-loss choice was frozen before1,024fresh stock paths: half-life{design['half_life_days'][mode]:.3f}days, lambda{lam}, eta{eta}. All1,751eligible/1,294primary SPX quotes and seven expiries are retained. Signed prices and invalid IVs stay in the artifacts. The auxiliary65,536-path reference and market snapshot are reused, so this is not market out-of-sample validation. New target draws are independent of previous experiments; all checking candidates share these draws.

Results on256steps per trading day; left/right IVs are near strike/forward0.9/1.1 at29days. The market levels are{m.loc[0,'market_left_IV_pp']:.2f}% and{m.loc[0,'market_right_IV_pp']:.2f}%, giving a gap{m.loc[0,'market_right_minus_left_IV_pp']:.2f}pp.

| Model | SPX price loss | IV MAE (pp) | 29d left IV (%) | 29d right IV (%) | Gap (pp) |
|---|---:|---:|---:|---:|---:|
{table}

Price-loss reduction versus baseline:{gain:.2f}%; versus previous rich correction:{extra:.2f}%. Approximate1,024-resample eight-batch paired bootstrap intervals are below; they resample full-model batches together across candidates and resample cheap-reference batches independently. They exclude deterministic error, fitting variability and market uncertainty.

| Comparison | Loss change | Approximate95% interval |
|---|---:|---:|
{cis}

Compiled/Python relative formula error{audit['compiled_formula_relative_error']:.3g}; full-potential finite-difference generator error scaled by1+A:{audit['finite_difference_full_generator_error_scaled']:.3g}; conditional scaled-integral error{audit['conditional_scaled_integral_error']:.3g}; analytic conditional oracle versus compiled GH64 error{audit['conditional_oracle_GH64_error']:.3g}; baseline/previous-rich replay error{audit['baseline_and_rich_replay_error']:.3g}; independent signed stock recurrence error{rec['maximum_scalar_stock_error']:.3g}. The scalar recurrence uses a synthetic clock to cover all time features; it is a finite correspondence audit, not a separate pricing model.

Both128/256stock grids and all per-maturity native-IV/stock-mean diagnostics are saved in `visualizations/available-variance029/summary.json` and `maturity_metrics.csv`. The selected correction's conditional left-sum audit on actual grids is `runs/m2_available_variance_001/discrete_window_audit.json`; maximum error relative to initial stress and maximal curve budget is{summary['max_discrete_window_error_scaled']:.3g}years. Exact cancellation belongs to the continuous-time model. Discrete extreme-state errors must be assessed relative to stress size rather than presented as a universal absolute bound. Direct-Q variance-budget probes and raw option-tail payoff concentration/ESS are saved; their finite samples do not certify rare-tail convergence. No sample stock-mean renormalization or smile smoothing is used.

## Conclusion

The main mathematical answer is yes: this explicit correction can use the current stress variance while preserving quoted-date VIX and the original positive floor. The numerical table assesses whether this particular construction supplies the missing SPX shape; do not identify feasibility with calibration success. Both the remaining downside IV shortfall and any wing deterioration must be read from the plots, not just the aggregate price score. No candidate is adopted. Preserve this search rather than rerunning it as untested. Future changes must retain the full-generator terms, the curve-jump endpoint conditions and the positivity argument.

[Seven SPX smiles](visualizations/available-variance029/SPX-available-variance.png). Full artifacts: `runs/m2_available_variance_001`.
'''
    outcome=f'''
## Decision from this test

The new correction has a {gain:.2f}% lower sampled price loss than baseline, but its paired loss-change interval includes zero. Its sampled loss is higher than the previous rich correction, and that paired difference also has an interval including zero. This sample therefore does not establish an improvement over either comparator. The 29-day left-wing IV rises from {m.loc[0,'left_IV_pp']:.2f}% to {m.loc[2,'left_IV_pp']:.2f}%, while the market is {m.loc[0,'market_left_IV_pp']:.2f}%. The much larger downside-volatility shortfall remains.

The direct-Q probes confirm that the new correction does use more than the old absolute f/2 budget in some states. However, the probes at fitted curve endpoints are zero by construction and say nothing about typical interior amplitudes. At day7 the observed 99th percentile of absolute correction is {budget[0]['p99_absolute_correction']:.5f}, versus the old fixed limit {f/2:.5f}. These finite samples are descriptive rather than a tail-convergence certificate.

The feasibility question has a constructive positive answer; this numerical specification is not a successful calibration. Do not increase simulation count merely to establish a small aggregate gain while the wing-level error remains this large. Preserve the model as an exploratory family, with no adoption.
'''
    (HERE/'AVAILABLE_VARIANCE_RESULTS.md').write_text(theory+result+outcome)
    note=f'''## 29 October 2026 — existing stress variance as correction budget tested

Read [AVAILABLE_VARIANCE_RESULTS.md](AVAILABLE_VARIANCE_RESULTS.md), run `runs/m2_available_variance_001`, plots `visualizations/available-variance029`. Potential P_j=c_j(t)beta_j(t)erf(U/sqrt2)A_j exp[-lambda(Z_j-b_j Q_j(t))^2], j=slow/fast stress; h_j=(partial_t+L)P_j, v=v0+eta sumh. Include full time/stress generator terms; simple multiplication of old h by variance is not protected. Clock zeros at all16curve breakpoints (includes quoted VIX endpoints) remove curve-coefficient jump terms. Analytic |h_j|<=c_j A_j implies v>=ORIGINALf=.0144 for|eta|<=1; can exceed old fixed f/2 budget. OU law, oldM2/cap/Qnormalizer/positive spike shift, stresses, curve and stock driver fixed. Odd independent U preserves ideal mean variance. Four quoted VIX values protected state by state, baseline errors inherited; no exact market-forward claim, no unquoted/full-process protection. This is a whole stress-scaled construction, not a pure amplitude ablation versus prior rich correction; calendar/damping/floor differ too.

512pilotpaths/38configs, selected half-life{design['half_life_days'][mode]:.3f}d/lambda{lam}/eta{eta} frozen before1024fresh targets. All1751eligible/1294primary/sevenSPXmaturities,217VIX/fourdates retained. Fresh256step losses baseline{fine.loc[0,'loss']:.2f}/previousrich{fine.loc[1,'loss']:.2f}/new{fine.loc[2,'loss']:.2f}, gain{gain:.2f}%vsbaseline/{extra:.2f}%vsprevious.29d gaps{m.loc[0,'model_right_minus_left_IV_pp']:.2f}/{m.loc[1,'model_right_minus_left_IV_pp']:.2f}/{m.loc[2,'model_right_minus_left_IV_pp']:.2f}pp vsmarket{m.loc[0,'market_right_minus_left_IV_pp']:.2f}; new leftIV{m.loc[2,'left_IV_pp']:.2f}%vsmarket{m.loc[0,'market_left_IV_pp']:.2f}%. Paired intervals, stock128/256, conditional discrete sums, tail/budget probes saved. Compiled full formula, independent full-potential generator differences, conditional analytic/GH oracle, old replay and scalar stock recurrence audited. No new Lean proof, uniform numerical bound or martingale theorem. Auxiliary65536paths and market snapshot reused. No candidate adopted or production/old artifacts changed.

'''
    for name in ['NEXT_STEPS.md','RESEARCH_LEDGER.md','CONTINUE_FITTING_PROMPT.md']:
        p=HERE/name;p.write_text(note+p.read_text())
    p=HERE/'formalization/CORRESPONDENCE.md';p.write_text(p.read_text()+'''\n## 29 October 2026: existing-stress variance budget\n\n`m2_available_variance.cpp/.py` apply the full generator to stress-scaled, Gaussian-damped potentials, including all time/stress terms. Analytic bounds give |h_j|<=c_j A_j and retain original floor f for |eta|<=1. Zero potentials/time derivatives at all fitted curve breaks avoid jump terms and retain four conditional VIX-window cancellations. Numerical formula, full-generator differences, analytic/GH conditional oracle, old replay, scalar stock recurrence and conditional discrete sums are in `../runs/m2_available_variance_001`. These are scoped numerical correspondence checks, not new Lean proofs, a uniform error bound, a continuous-stock martingale theorem or exact market-forward calibration.\n''')
    for name in ['audit_available_variance.py','report_available_variance.py','AVAILABLE_VARIANCE_RESULTS.md','asymmetric_readout.py','forward_leverage.py','leverage_scan.py','protected_curve_check.py','audit_leverage_scan.py']:shutil.copy2(HERE/name,OUT/'source'/name)
    shutil.copy2('/private/tmp/m2-available-001.log',OUT/'run.log');shutil.copy2('/private/tmp/m2-available-recurrence.log',OUT/'recurrence.log')
    for source,h in json.loads((OUT/'source_hashes.json').read_text()).items():assert hashlib.sha256(Path(source).read_bytes()).hexdigest()==h,source
    counts={}
    for name in ['m2_leverage_scan_001','m2_asymmetric_readout_001','m2_targeted_skew_001','m2_forward_leverage_001','m2_slow_readout_001','m2_slow_readout_002','m2_memory_filter_001','m2_protected_variance_design_001','m2_protected_variance_001','m2_rich_protected_variance_001']:
        prior=OUT.parent/name;manifest=json.loads((prior/'artifact_manifest.json').read_text())
        for file,h in manifest.items():assert hashlib.sha256((prior/file).read_bytes()).hexdigest()==h,(name,file)
        counts[name]=len(manifest)
    run.update(status='completed pricing, numerical audits and documented theory; figures pending visual review',summary=summary,prior_artifact_hash_counts_verified=counts);dump(OUT/'run.json',run);dump(OUT/'artifact_manifest.json',{str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json'});print('AVAILABLE REPORT COMPLETE',json.dumps(summary),flush=True)

if __name__=='__main__':main()
