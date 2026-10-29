# Tight VIX forwards and a fitted variance floor

Completed27October2026. User authorized points1and2: tighten forward fitting,
then let the minimum variance floor be calibrated, retaining free forward
variance. Run `runs/m2_forward_floor_001`; plan
[FORWARD_FLOOR_PLAN.md](FORWARD_FLOOR_PLAN.md). No new factor/readout was added.

Allowing the floor to move substantially improves the native option fit, with
much smaller fresh forward mismatches. The55-day lowest put is close to market;
the27-day left end improves but stays high. The lowest longer-dated puts become
too cheap. Right-wing IV error remains similar on average, with different
maturity trade-offs. This is an improved exploratory VIX candidate, not a complete
smile fit or accepted joint SPX/VIX model.

## Controlled comparison and scope

Same additive response v=f+(xi-f)[(1-w)R_old+w L_b] as the prior mixture winner.
The old readout/cap, banks, Q premiums, initial state and positive beta_J remain
fixed. Three candidates are replayed on common training and fresh checking draws:

1. Previous mixture, unchanged, including its forward mismatches.
2. Fixed old floor0.00264 (5.138% volatility), four forward equality constraints.
3. Same equality constraints with floor volatility free from5%to12%.

All16 curve intervals remain represented. Seven VIX-sensitive levels fit
independently through xi_j=f+(0.20-f)y_j,0<=y_j<=1; the other nine have no influence
on these VIX contracts and retain old values. All16 satisfy xi>=f. Curve shape,
b in[0,5], mixture weight in[0,1] and half-life in[1,365]days can change in both
fits. There is no curve smoothing or monotonicity penalty. All217 original native
OTM quotes/four forwards remain; SPX stays paused.

The calibration minimizes native bid–ask option-price distances with original
quote scales and equal total expiry weights, **subject to the four forward
estimates matching their archived targets**. Six fixed-floor and seven free-floor
starts use bounded SLSQP, at most120iterations/start, after a bounded curve-only
forward prefit. Only training forward gaps<=1e-5VIXpoint are selection-eligible.
Forward equalities hold on the training estimator; no exact population-forward
claim is made. Price comparisons below use the original option-plus-forward
score with0.5point forward scale, retaining comparability across all candidates.

Training explicitly reuses the previous mixture study's16384-point/expiry fresh
stream and Gaussian conditional nodes, now training data. Final selection is
frozen before new32768-point/expiry numerical checking (eight4096-point nets,
seed base20261026083, the same68-component joint37D importance proposal). Exact
importance weights and the original independent-pilot controls remain frozen;
there is no self-normalization. All market data remain calibration observations.

## Results

| Candidate | Training option-only loss | Fresh comparable total loss | Max fresh forward gap | ATM IV MAE | Wing IV MAE | Inside bid–ask |
|---|---:|---:|---:|---:|---:|---:|
| previous_mixture | 9.1787 | 10.2653 | 0.9001 | 9.37 | 8.61 | 49/217 |
| tight_forwards | 12.7156 | 12.7780 | 0.0322 | 8.81 | 9.20 | 42/217 |
| tight_forwards_free_floor | 1.7690 | 1.7830 | 0.0269 | 5.62 | 8.45 | 84/217 |

IV errors are volatility percentage points. The fixed-floor constrained fit
sacrifices option fit to match forwards. Releasing the floor lowers fresh total
loss10.2653->1.7830 relative to the previous mixture: paired change-8.4824,
eight-net jackknife SE0.1083. Relative to the fixed-floor constrained fit, change
is-10.9951, SE0.1142. These sensitivity estimates condition on fitted parameters
and the frozen pilot; they do not account for the full adaptive fitting history.

Training forward gaps are below1.9e-10point for the selected constrained fits.
The free-floor candidate's fresh gaps at27/55/118/153days are
-0.02445,-0.00845,+0.02694,-0.00884points, with scramble SEs
0.01331,0.01745,0.02117,0.01396 respectively. Sampling error in both training
and checking estimators remains. The prior largest mismatch was0.90012point.

The selected floor is **12% volatility**, f=0.0144, at the declared search ceiling.
It is the best tested value, not an identified optimum or a required fixed floor.
The selected sensitivity b=3.1710, stress-component weight0.66150 and half-life
104.215days are interior. Original M2 now contributes about33.85% of expected
non-floor variance, compared with0.43% in the previous mixture. Fitting the floor
permits a different response balance and persistence; this is a joint refit,
not the effect of changing f while holding all other parameters constant.

Six free-floor starts converge to a local option loss around4.4034; a start in
the slower-mixture region finds1.7690. Every outcome is retained in
`optimizer_ledger.json`; this finite search gives no global-optimum certificate.
One short curve level reaches the floor and the118–119day piece reaches0.20.
Full16-level values are saved in `curve_levels.csv` and experimental model JSONs.

## Left-end prices and right-wing trade-offs

These are original native put prices, without converting to calls. Signed
controlled estimates are retained, including small negative estimates.

| Expiry days | Lowest strike | Market put mid | Previous mixture | Tight forwards only | Tight forwards + free floor |
|---|---:|---:|---:|---:|---:|
| 27 | 15.5 | 0.091196 | 0.270502 | 0.320976 | 0.186498 |
| 55 | 14.5 | 0.049496 | 0.279971 | 0.314193 | 0.050529 |
| 118 | 12.5 | 0.043016 | 0.028503 | 0.041132 | -0.000011 |
| 153 | 12.5 | 0.045013 | 0.016896 | 0.032078 | -0.000078 |

At55days the lowest put moves close to the market mid; at27days its overpricing
shrinks. At118/153days the new lowest-strike puts are much too cheap. Five small
negative fresh controlled put estimates have magnitudes0.10–2.35reported SEs;
the largest magnitude is below0.000451price units. They remain saved and are not
clipped. Additional direct nonnegative payoff estimates are saved separately in
`raw_native_puts_*.csv`. They are zero at the two lowest long-dated strikes in
these finite samples; this does not prove population prices are zero. The
underpricing relative to market mids around0.04 is much larger than the small
signed control correction. This is a model fit limitation as well as a reason
to preserve estimator diagnostics.

Native model-forward IVs are unavailable at five fresh quotes; canonical
market-forward call IVs are unavailable at four. Missing values are marked in
the plots. We provide both conventions and native put-price plots, because even
a small forward mismatch can noticeably affect IV at low strikes. We never
replace model forwards with targets inside model prices or subtract their error
to manufacture a better fit. The earlier exact parity-error decomposition is in
`runs/m2_left_end_diagnostic_001`.

Wing slopes versus log(K/F), using all retained strikes above2F:

| Days | Market | Previous mixture | Tight forwards only | Tight forwards + free floor |
|---|---:|---:|---:|---:|
| 27 | 52.38 | 76.05 | 76.73 | 43.61 |
| 55 | 38.31 | 56.28 | 56.40 | 35.92 |
| 118 | 31.23 | 35.52 | 35.96 | 29.61 |
| 153 | 27.32 | 29.30 | 30.01 | 27.73 |

The averaged wing IV-level error stays near8.5points. At27days the far-right wing
level becomes too low, while the longer-dated wing levels improve. Hence the
left-end gain does not give a uniformly better IV fit at every strike/maturity.
Only84/217 native prices lie inside bid–ask; no complete-smile fit is claimed.

The free-floor exact E[VIX²] budgets are466.29,563.69,691.52,746.05. If exact
population forwards are required, the27-day budget is below the archived necessary
bound469.00 for all its quotes to fit within bid–ask. This conditional obstruction
is another reason not to declare a complete calibration; it does not turn
training-estimator forward equalities into population identities. Earlier
exact-rational payoff checks are in `runs/m2_vix_design_research_001`.

## Numerical checks and saved plots

Old additive prices replay within1e-9. All11analytic native-price and forward
Jacobian directions pass finite-difference checks: maximum relative discrepancies
4.55e-6 and1.19e-9 respectively. Independent physical-parameter nodal replay
checks six candidate/stream price sets, forwards, parity, floor support, all16
levels, unchanged off-window values and exact second moments; maximum native
price difference is below1.8e-15. All32fresh scrambles and proposal weights replay.
291source/fiveinput hashes pass. `audit.json` is PASS for these finite checks;
it does not assert option-price positivity or a successful complete calibration.

Sampled GH32/time8 to GH48/time16 refinement uses35–38states per expiry, including
large weighted payoffs. Max VIX difference is1.43e-6point for the free-floor
candidate,0.000856across all candidates. This is not a uniform integration-error
bound. Its likelihood-weighted last-strike payoff ESS is2905–3276, with the top
ten contributions about0.52–0.58%. Tail/moment/negative-price diagnostics remain
saved; no population convergence or new Lean theorem is asserted.

- [Canonical call-IV comparison](visualizations/forward-floor021/vix-all-strikes.png)
- [Native IVs using each model forward](visualizations/forward-floor021/vix-native-IV.png)
- [Actual native put prices](visualizations/forward-floor021/native-put-prices.png)
- [Right wings](visualizations/forward-floor021/vix-right-wing.png)
- [Fresh forward errors and the curve/floor](visualizations/forward-floor021/forwards-and-curve.png)

PNG and SVG versions are saved. Pointwise price-sampling bands exclude systematic
integration error and frozen-pilot uncertainty. The current exploratory VIX
candidate is `experimental_model_tight_forwards_free_floor.json` in this run.
Original selected_model.json remains unchanged; no SPX paths were read or priced.
Possible next research would address maturity-specific quiet-state behaviour
and the27-day wing/variance budget. No such extension is implemented here, and
raising the global floor further has not been tested.
