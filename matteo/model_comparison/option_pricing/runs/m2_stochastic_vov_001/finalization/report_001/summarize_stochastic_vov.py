"""Evidence-linked research report for the bounded dedicated-layer grid."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .full_readout_fit import HERE
from .cap_controlled_core import shape_metrics


def vf(score):return score['original_otm_components']['VIX']+score['original_otm_components']['forwards']


def report(run):
    profiles=json.loads((run/'profiles.json').read_text());selection=json.loads((run/'selection.json').read_text());audit=json.loads((run/'audit.json').read_text());design=json.loads((run/'design.json').read_text());names=selection['full_candidates']
    full={n:json.loads((run/'full'/f'score_{n}.json').read_text()) for n in names};fresh={n:json.loads((run/'fresh_check'/f'score_{n}.json').read_text()) for n in names};lines=[];add=lines.append
    add('# Dedicated stochastic volatility response: results');add();add('25 October 2026. Completed bounded fixed-curve experiment. The selected baseline remains `runs/m2_maturity_curve_fit_001/selected_model.json`. New run: `runs/m2_stochastic_vov_001`. No new model is adopted.');add()
    add('## What we tried, in everyday language');add();add('The new state is a sensitivity dial. A high setting makes the existing spike factor move volatility more strongly; its half-life determines how long that episode tends to persist. This changes the distribution of future volatility while keeping its expected variance curve fixed through a new model-specific normalizer. It is an extra continuous independent OU state, with no jumps or changes to the existing reservoir.');add()
    add('We crossed 1.5×/2× sensitivity ranges with 7/30/90-day half-lives and two outer transformations: the original cap and the previously tested logarithmic continuation of strength .5. Both transformations also have a zero-strength reference. All16 variance levels, floor .00264, original readout coefficients, OU rates/geometric grids/state/correlations and risk-price scenario remain EXACTLY fixed. There is no curve or floor refit in this experiment.');add()
    add('## Equations and mathematical scope');add();add('Under Q:');add();add('```text');add('dZ = -aZ Z dt + sqrt(2aZ) dBZ, Z0=0, BZ independent;');add('A(Z)=exp(s tanh Z), aZ=log(2)*365/half_life_days;');add('eta=cC C+qY(Y-betaY)^2+qJ A(Z)(F-betaJ)^2;');add('G=exp(2 h(eta)), M(t)=E_Q[G(t)];');add('v=f+(xi(t)-f) G/M(t).');add('```');add()
    add('The original cap has h=L-softplus(k(L-eta))/k, with L=5.502780056785304 and k=2. The continuation adds .5 log(1+softplus(k(eta-L))/k). A stays between exp(-s) and exp(s). At s=0 the layer disappears exactly. The extra Brownian variance response is sqrt(2aZ)*partial_Z(v), where partial_Z(eta)=qJ A s sech²(Z)(F-betaJ)². Thus the new source directly changes spike sensitivity; cap saturation can still damp it.');add()
    add('A is centered at one in log space; its arithmetic mean need not equal one. The G/M normalization preserves expected variance xi, not the average spike coefficient. The comparison tests this specified modulation as a whole, including both its random sensitivity and its nonlinear average-response effects.');add()
    add('The capped multiplier remains bounded. With the continuation, bounded A gives G at most quadratic growth in the enlarged Gaussian state. Gaussian multiplier moments remain finite, and the finite-horizon linear-growth/nonexplosion stock-martingale argument in CAP_TAIL_MATHEMATICS.md extends to this layer. These mathematical arguments are distinct from quadrature/sampling checks and are not new Lean proofs. They do not certify finite long-maturity stock-payoff or importance-estimator variance for the continued model. The new Q mean is zero; no P calibration or added risk-premium identification is claimed.');add()
    add('## Complete VIX grid, entire curve fixed');add();add('|Case|Cap continuation|Sensitivity range|Half-life days|Native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Unavailable wing IVs|');add('|---|---:|---:|---:|---:|---:|---:|---:|')
    shape_rows=[]
    for row in profiles:
        s=row['shape'];spec=row.get('spec');missing=sum(x['wing_quotes']-x['wing_available'] for x in s['rows'])
        add(f"|{row['name']}|{row['continuation']:g}|{spec[4] if spec else 1:g}×|{spec[3] if spec else '—'}|{vf(row['score']):.6f}|{s['wing_MAE']:.3f}|{s['ATM_MAE']:.3f}|{missing}|")
        shape_rows.extend(dict(case=row['name'],**v) for v in s['rows'])
    pd.DataFrame(shape_rows).to_csv(run/'grid_smile_metrics.csv',index=False)
    add();add('IV MAE is in percentage points, averaging maturities equally. Wing uses K/marketF>2; near ATM uses |log(K/marketF)|<=.2. Slopes and curvatures are descriptive regressions of raw implied volatilities, not a polynomial volatility model. Complete per-maturity values are saved in `grid_smile_metrics.csv` and `profiles.json`.');add()
    add('## Wing level is different from wing steepness');add()
    add('A lower wing IV error can come from lifting the whole smile, even when the right-wing slope becomes flatter. Selection by wing IV MAE does not select the best slope. Compare the actual slope and curvature below; they are measured against log strike over K/marketF>2.');add()
    market_table=pd.read_csv(run/'prices_canonical_cap0_ref.csv');market=shape_metrics(market_table,market_table.price_mid_U.to_numpy());market_byday={x['days']:x for x in market['rows']}
    add('|Case|Days|Mean wing IV|Wing slope|Wing curvature|');add('|---|---:|---:|---:|---:|')
    for x in market['rows']:add(f"|Market|{x['days']}|{x['wing_IV_level']:.3f}|{x['wing_slope']:.3f}|{x['wing_curvature']:.3f}|")
    for n in names:
        row=next(x for x in profiles if x['name']==n)
        for x in row['shape']['rows']:add(f"|{n}|{x['days']}|{x['wing_IV_level']:.3f}|{x['wing_slope']:.3f}|{x['wing_curvature']:.3f}|")
    add();add('The term “wing improvement” in MAE tables refers to closeness in IV level, not a claim of greater or better-matched steepness.');add()
    regions=pd.read_csv(run/'grid_region_diagnostics.csv');add('The native option-loss decomposition shows where the price-fit tradeoff occurs. These contributions sum to the saved VIX option loss; forward penalties are separate.');add();add('|Case|Near-ATM price loss contribution|Other strikes|Right-wing price loss contribution|');add('|---|---:|---:|---:|')
    for n in names:
        rr=regions[regions.candidate.eq(n)].groupby('region').native_price_loss_contribution.sum();add(f"|{n}|{rr['near_ATM']:.4f}|{rr['other']:.4f}|{rr['right_wing']:.4f}|")
    add()
    add('The normalizer matters here: a larger raw sensitivity does not increase every normalized volatility outcome. The new model divides by its own average multiplier. For the capped case, increasing M lowers the normalized height of already saturated states. This can help explain why smile levels and wing slopes respond differently; it is not a proof of a unique causal decomposition.');add()
    normalizers=pd.read_csv(run/'normalizer_comparison.csv');add('|Cap|Selected sensitivity/half-life|M_new / M_no_layer at153 days|');add('|---|---|---:|')
    for cap in ['cap0','cap05']:
        for name in selection['matched_specs']:
            v=normalizers[normalizers.cap.eq(cap)&normalizers.candidate.eq(name)&normalizers.days.eq(153)].iloc[0];add(f"|{cap}|{name}|{v.ratio:.4f}|")
    add();add('All six candidates and all four dates are saved in `normalizer_comparison.csv`. Expected variance still equals the same xi curve in the mathematical model; these ratios describe redistribution across states.');add()
    add('Training-only wing winners: '+', '.join('`'+x+'`' for x in selection['training_winners'])+'. Selection froze before SPX prices or fresh checking. The union of their sensitivity/half-life settings was crossed with BOTH caps, plus the two no-layer references. No aggregate, forward or near-ATM acceptance was imposed by this descriptive selection. The full-SPX subset does not establish the best joint point over the entire grid.');add()
    add('## Actual full SPX repricing');add();add('|Case|SPX loss|Native VIX + forwards|Native joint loss|Native prices inside bid/ask /1511|Negative native prices|Unavailable call IVs|');add('|---|---:|---:|---:|---:|---:|---:|')
    for n in names:
        s=full[n];add(f"|{n}|{s['original_otm_components']['SPX']:.6f}|{vf(s):.6f}|{s['original_otm_loss']:.6f}|{s['inside_original_otm_bidask']}|{s['native_otm_negative_prices']}|{s['iv_unavailable']}|")
    add();add('Every SPX entry above recomputes volatility and stock at all seven agreed maturities, using all4096 paths at every strike. Old exogenous OU states, equity increments, likelihood weights and frozen reference/hedge controls remain reusable; the terminal stock paths do not. The independent Z process is streamed at the SAME existing left-point mesh using exact OU transitions and parameter-independent counter-based Gaussian innovations. There is no additional multi-gigabyte all-time Z cache.');add()
    add('|Layer case|Comparison reference|SPX loss change ± paired block SE|Native joint change ± paired block SE|');add('|---|---|---:|---:|')
    for n in names:
        if n.endswith('_ref'):continue
        for suffix,reference in [('baseline','cap0_ref'),('same_cap',n.split('_')[0]+'_ref')]:
            z=json.loads((run/'full'/f'paired_{n}_{suffix}.json').read_text());add(f"|{n}|{reference}|{z['change'][1]:+.6f} ± {z['jackknife_SE'][1]:.6f}|{z['change'][0]:+.6f} ± {z['jackknife_SE'][0]:.6f}|")
    add();add('Paired eight-block sensitivities are descriptive numerical diagnostics, not a population confidence or convergence certificate. Loss scores measure fitting error; changes in scores are not proportional changes in option prices. Native/common prices and finite monotonicity/convexity flags are retained in score/CSV artifacts, with no clipping or smoothing.');add()
    add('Negative native estimates here occur in SPX. VIX can separately have unavailable call IVs because the graph deliberately uses the observed market forward while the model forward differs: a positive model call can fall below the intrinsic bound associated with that market-forward Black inversion. These two kinds of failure are recorded separately; neither is hidden by changing the plotting forward.');add()
    neg=[x for x in audit['negative_estimate_audit'] if x['role']=='SPX+VIX' and x['count']]
    if neg:add(f"The negative SPX estimates have magnitudes at most {max(x['largest_negative_magnitude_in_reported_SE'] for x in neg):.3f} of their reported sampling SE. They are compatible with finite-estimator noise around small positive prices; they are not evidence of negative theoretical option prices.");add()
    add('|Case|1030-day raw weighted stock mean (target1)|Raw iid SE|');add('|---|---:|---:|')
    for n in names:
        mm=next(x for x in full[n]['moments'] if x['instrument']=='SPX' and x['days']==1030);add(f"|{n}|{mm['raw_forward_ratio']:.6f}|{mm['raw_forward_ratio_se']:.6f}|")
    add();add('The long-maturity raw stock-mean discrepancy already exists in the baseline. Exact mean-one controls remain mathematically valid, but do not certify sufficient finite-sample accuracy. New negative signed estimates or missing IVs must be read with the reported sampling errors and loss sensitivities, not as model-level arbitrage or a stock-martingale counterexample.');add()
    add('## Genuinely new joint VIX checking draws');add();add('|Case|New native VIX + forwards|Wing IV MAE|Near-ATM IV MAE|Unavailable wing IVs|Native VIX loss change vs baseline ± paired net SE|');add('|---|---:|---:|---:|---:|---:|');fresh_shapes=[]
    for n in names:
        q=pd.read_csv(run/'fresh_check'/f'prices_canonical_{n}.csv');s=shape_metrics(q,q.model_price_U.to_numpy());fresh_shapes.append(dict(name=n,**s));z=json.loads((run/'fresh_check'/f'paired_{n}_baseline.json').read_text())
        missing=sum(x['wing_quotes']-x['wing_available'] for x in s['rows'])
        add(f"|{n}|{vf(fresh[n]):.6f}|{s['wing_MAE']:.3f}|{s['ATM_MAE']:.3f}|{missing}|{z['change'][2]:+.6f} ± {z['jackknife_SE'][2]:.6f}|")
    (run/'fresh_smile_metrics.json').write_text(json.dumps(fresh_shapes,indent=2)+'\n')
    add();add('After selection,32 new complete scrambled nets were generated, eight512-point nets per maturity, in the enlarged37-dimensional Gaussian space. Old defensive-mixture shifts were extended by a zero Z component, retaining exact likelihood weights without self-normalization. New joint draws improve the dimension pairing over the training independent-1D augmentation; neither construction guarantees adequate sampling of joint high-response tails. Old independently fitted control coefficients remain frozen. Checking cannot trigger reselection. The market quotes are the SAME in-sample observations; no market out-of-sample claim is made. No fresh SPX ensemble was generated.');add()
    add('## Weighted tails and variance targets');add();mom=pd.read_csv(run/'fresh_check/moments.csv');add('|Case|Days|E[VIX²] estimate|Fixed curve target|Eight-net SE|Difference / SE|');add('|---|---:|---:|---:|---:|---:|')
    for v in mom.itertuples():add(f'|{v.candidate}|{v.days}|{v.second:.4f}|{v.target:.4f}|{v.second_SE:.4f}|{(v.second-v.target)/v.second_SE:+.2f}|')
    add();add('Because the entire xi curve stays fixed, the target is identical for every case at each maturity. Historical 153-day second-moment discrepancies are preserved; a new estimate closer to target is not a convergence pass.');add()
    tail=pd.read_csv(run/'fresh_check/tails.csv');add('|Case|Days|Highest strike|Proposal hits|Weighted tail probability|Raw call ± SE|Payoff ESS|Top10 payoff share|');add('|---|---:|---:|---:|---:|---:|---:|---:|')
    for (n,day),g in tail.groupby(['candidate','days'],sort=False):
        v=g.loc[g.strike.idxmax()];add(f'|{n}|{day}|{v.strike:g}|{int(v.hits)}|{v.probability:.5g}|{v.raw_call:.6g} ± {v.call_SE:.2g}|{v.payoff_ESS:.1f}|{100*v.top10_share:.1f}%|')
    add();add('Raw hit counts are proposal events, not Q probabilities. Tail probabilities and calls use exact likelihood weights; raw calls exclude the signed control corrections used in actual option estimates. Low payoff ESS or concentrated contributions remain visible and limit interpretation.');add()
    add('## Verification and remaining numerical limits');add()
    checks=json.loads((run/'interpolation_checks.json').read_text());checks=[r for cap in checks.values() for r in cap]
    add(f"- Three independent tests check the Gaussian-response integral against elementary and dense-grid oracles, counter-normal replay/OU variance, and streamed stock paths against a separate Python implementation. Zero-strength stock prices nest the old kernel for both cap choices.")
    add(f"- Original baseline full price replay maximum absolute difference: {audit['baseline_full_maximum_price_difference']:.3g}. Independent saved-draw audit replays {len(audit['score_records'])} price/SE/parity/loss records, {audit['frozen_sources']} archived sources and all32 complete new joint scrambles. The selected baseline file is unchanged; selection timing is verified.")
    add(f"- Checked response interpolation uses degrees {sorted(set(x['degree'] for x in checks))}; maximum sampled off-node relative error {max(x['relative_error'] for x in checks):.3g}. Six off-node factors per table include selected extreme old-state means. This checks finite quadrature/interpolation, not a uniform error bound over all states.")
    add(f"- A separate numpy four-dimensional Gaussian implementation checks {len(audit['independent_4D_checks'])} conditional/unconditional expectations, maximum relative error {max(x['relative_error'] for x in audit['independent_4D_checks']):.3g}.")
    add(f"- Representative time order8-to16 checks at ordinary/extreme153-day checking states change VIX by at most {max(x['time_order8_to16_maximum_VIX_change'] for x in audit['time_order_checks']):.3g} points.")
    add(f"- Representative original-state normalizer GH48-to96 (also refining new-state192-to256) changes M by at most {max(abs(x['old_GH48_to96_and_new192_to256_relative_change']) for x in audit['normalizer_order_checks']):.3g} relatively; all six added-state cases and both caps are checked at27/153/1030 days.")
    for cap in ['cap0','cap05']:
        q=json.loads((run/f'quadrature_{cap}.json').read_text())
        for kind in sorted(set(x['kind'] for x in q)):
            add(f"- {cap}, {kind}: maximum representative relative change {max(x['relative_error'] for x in q if x['kind']==kind):.3g}.")
    add('- Conditional original-state GH32, normalizer original-state GH48, additional-state GH128/GH192, time order8 and cached stock512 mesh are approximations. Sampled refinements are not whole-surface price-convergence proofs. Standard errors omit integration/stock biases, pilot/parameter uncertainty and unsampled tails. New-state terminal moment diagnostics are saved in the audit.');add()
    add('## Decision');add()
    base=full['cap0_ref'];layer=[n for n in names if not n.endswith('_ref')]
    improving=[n for n in layer if full[n]['original_otm_components']['SPX']<base['original_otm_components']['SPX'] and vf(full[n])<vf(base)]
    if improving:add('These tested cases reduce both sampled native SPX and VIX-plus-forward losses versus the original baseline: '+', '.join(improving)+'. They remain experimental pending the numerical and cross-market checks above; no automatic adoption is made.')
    else:add('None of the matched layer representatives improves both sampled native SPX and VIX-plus-forward losses versus the selected baseline. Preserve the baseline. This bounded experiment does not establish that every continuous stochastic-response model fails.')
    for cap in ['cap0','cap05']:
        ref=next(x for x in profiles if x['name']==cap+'_ref');win=next(x for x in profiles if x['name'] in selection['training_winners'] and x['cap']==cap)
        direction='improves' if win['shape']['wing_MAE']<ref['shape']['wing_MAE'] else 'worsens'
        add(f"For {cap}, its wing-selected extra state {direction} training wing IV MAE from {ref['shape']['wing_MAE']:.3f} to {win['shape']['wing_MAE']:.3f}; near-ATM MAE changes from {ref['shape']['ATM_MAE']:.3f} to {win['shape']['ATM_MAE']:.3f}. Read this alongside actual forwards, price losses, SPX and new-stream evidence above.")
    add();add('This tests one spike-only bounded modulation and the predeclared two-cap grid, without fitting xi or the old coefficients. The earlier cluster experiment, a different modulation, and a possible future feedback modulation are distinct model questions. No full calibration convergence or structural impossibility follows.');add()
    add('Plan: [STOCHASTIC_VOV_PLAN.md](STOCHASTIC_VOV_PLAN.md). Figures and reading guide: [stochastic-vov016](visualizations/stochastic-vov016/README.md). Reproduce with the pinned pdv-options-venv and package PYTHONPATH using modules `stochastic_vov_study --output NEW_DIRECTORY`, `audit_stochastic_vov --run RUN`, `plot_stochastic_vov --run RUN --output PLOTS`, `summarize_stochastic_vov --run RUN`. Source snapshots and raw draws are in the run.');add()
    (HERE/'STOCHASTIC_VOV_RESULTS.md').write_text('\n'.join(lines)+'\n');print('Wrote STOCHASTIC_VOV_RESULTS.md',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);report(p.parse_args().run.resolve())
