# Option snapshot preprocessing plan

Status: implementation plan, revised 30 September 2026 after consistency review. This document specifies work to do; it does not assert that the planned pipeline or its output already exists. The companion [calibration plan](CALIBRATION_PLAN.md) consumes the contracts defined here. Instructions in the model PDFs are research proposals; the decisions below reflect the present data and code.

## Objective and known inputs

Build a reproducible, auditable snapshot of market targets for the exogenous-reservoir M2 calibration. Preserve raw rows, quote intervals, instrument identities, and rejection reasons. Produce option prices in forward-value units, a conservative forward-variance input, four VIX forward-price targets, and an audit report. Do not fit model parameters in this pipeline.

| Input | Observed contents | Confirmed convention |
| --- | --- | --- |
| `data/raw/SPX_data_19_02_2026.csv` | 4,715 rows, one reference date, 11 expiries | `bid`/`ask` are annualised decimal implied volatilities. |
| `data/raw/VIX_data_19_02_2026.csv` | 520 rows, one reference date, four expiries | `forward` holds VIX forward prices; `bid`/`ask` are annualised decimal implied volatilities. |

Both files have `reference,expiry,callput,strike,exercise,bid,ask,volume_bid,volume_ask,spot,forward,moneyness,expiry_yf`. Reference serial 46072 corresponds to 2026-02-19 under the repository's Excel date convention. `callput` is +1 for a call and -1 for a put. The user confirmed the decimal-IV convention for both files and the VIX forward-price convention. The source and meaning of `volume_bid`/`volume_ask` are unverified; do not interpret them as executed volume or use them for liquidity weights until documented.

Observed two-sided positive, non-crossed rows: 3,550 SPX and 393 VIX. Missing bids: 1,165 SPX and 127 VIX; asks are present. There were no duplicate `(reference,expiry,callput,strike)` keys in the inspected files. This inventory is a baseline to reproduce, not a guarantee about future inputs. The historical index data used by the physical-measure plan end on 2023-11-14, so this snapshot alone cannot supply the 2026-02-19 latent reservoir state.

## Fixed interpretation rules

1. Raw CSVs are immutable. Load paths from configuration or CLI arguments. Do not use the obsolete `/mnt/user-data/uploads/...` constants in `CrossSectionalData/calibrate_smile.py`.
2. `forward` is a maturity-specific SPX forward in the SPX file and a maturity-specific VIX forward price in the VIX file. Neither is the spot index level. The user believes the VIX forwards are traded futures quotes. The CSV does not record the source or venue, so preserve that as a likely but unverified provenance item.
3. Keep original strikes and absolute price units. Compute `moneyness_check = strike / forward`; use the stored moneyness only after verifying agreement. Keep exact provided `expiry_yf` for Black conversion; record that it currently equals integer calendar days divided by 365. Do not claim minute-accurate expiry time.
4. Use *undiscounted* option prices throughout the calibration interface: $P^U=P^{cash}/D(0,T)$. Black-76 with the supplied forward and (T) maps a confirmed decimal IV (s) to an undiscounted call or put. An independent discount curve is unnecessary for comparing these units, but is required to reconstruct cash premiums. Document whether the raw vendor's IV used the same forward and day count.
5. Convert bid IV and ask IV **separately**. Because Black price increases with IV, `price_bid_U = Black76(F,K,T,iv_bid)` and `price_ask_U = Black76(F,K,T,iv_ask)` for the same option type. Keep `price_mid_U = (price_bid_U + price_ask_U)/2` only as a diagnostic; never convert the average IV and call that a price midpoint.
6. A missing bid is missing information, not zero. The primary calibration uses valid two-sided quotes. A separate sensitivity run may use ask-only rows as one-sided upper bounds. Missing, zero, crossed, non-finite, and negative values get distinct reason codes.
7. Do not pool puts and calls by averaging implied volatilities at the same strike. They are two observations of one terminal law, with put-call parity and quote intervals as their consistency check.
8. Apply the same separate bid/ask Black conversion to VIX options, using each row's VIX forward in index points. Retain an explicit quote-unit field so future data cannot silently change this convention.

Black-76 forward-value formulas, for $d_1=[\log(F/K)+s^2T/2]/(s\sqrt T)$ and $d_2=d_1-s\sqrt T$, are $C^U=F\Phi(d_1)-K\Phi(d_2)$ and $P^U=K\Phi(-d_2)-F\Phi(-d_1)$. Implement stable limiting values for tiny $s\sqrt T$, but flag nonpositive (F,K,T) instead of hiding them. For European options, test bid/ask put-call parity in forward units using $C^U-P^U=F-K$; parity is a validation relation, not a licence to average inconsistent quotes.

## Implementation order and deliverables

### 1. Inventory and source contract

Create a pure loader and validator in `matteo/preprocessing/options.py` (or a small `matteo/preprocessing/options/` package if it materially improves clarity). Add a CLI such as `python -m matteo.preprocessing.options --spx ... --vix ... --output data/processed/options/2026-02-19`. The CLI must accept both source paths and an explicit quote-unit enum, defaulting to the confirmed `iv_decimal` for this snapshot. Any later `cash_price` support must be a separate tested conversion, never an inference from magnitudes.

Validate required columns and types; one reference date per run; `reference < expiry`; one option type per row; European exercise; positive finite spot, forward, strike and maturity; consistent forward and spot within each `(instrument,expiry)`; no duplicate contract key; and `expiry_yf` agreement with the declared day count within `1e-10`. Compare `moneyness` with (K/F) within `1e-10`. Halt on schema or key failures. Preserve questionable quotes as excluded rows with reasons; do not overwrite them or silently drop them.

Write SHA-256 hashes, input byte counts, preprocessing version, run time, selected conventions, date conversion, and configuration to `manifest.json`. Because `data/raw` and `data/processed` may be ignored by Git, hashes and a reproducible command are essential for handoff.

### 2. Quote conversion and quote-level audit

Create one output row per raw contract with raw columns and derived fields: `instrument,reference_date,expiry_date,T_years,option_type,strike,spot,forward,K_over_F,source_bid,source_ask,quote_unit,price_bid_U,price_ask_U,price_mid_U,validity,reason_codes`. Keep a stable raw row identifier. Do not overwrite original source values in the derived table.

Validation sequence: missingness; finite values; positive bid/ask in the stated unit; `ask >= bid`; valid Black conversion; `price_bid_U <= price_ask_U`; elementary European price bounds and intrinsic lower bounds when applicable. At very short expiries, retain a distinct `short_expiry` flag. Calculate quote width both in source units and forward price units. An extreme width is a diagnostic, not an automatic rejection unless the rule is predeclared in the configuration.

Use the existing `CrossSectionalData/market_vix.py` Black function as a numerical cross-check, but put the new source-independent conversion in the preprocessing module. Unit tests should cover puts, calls, intrinsic limits, zero-width quotes, missing bids, crossed quotes, and extreme strikes. Do not make importing the new module execute a data run.

### 3. Cross-contract and static-arbitrage checks

For each `(instrument,expiry)`, evaluate every valid matched put/call pair. Convert the put interval to an equivalent call interval by adding (F-K). Record whether that interval intersects the observed call interval and the size of any gap. There are some non-overlapping SPX pairs in the supplied snapshot, so report them by expiry and moneyness; do not replace the original intervals with an averaged IV.

Build an **OTM quote view** for primary calibration: put for (K<F), call for (K>F), and at (K=F) use the intersecting put/call interval if both are present. If (F) lies between listed strikes, retain the two surrounding contracts separately. Keep ITM rows as parity and sensitivity diagnostics. Do not assume the option type by moneyness when reading a row; use `callput`.

For each expiry, transform all OTM put intervals into equivalent call intervals and ask whether a decreasing, convex call curve with discrete slopes in ([-1,0]) can pass through every interval. Use linear feasibility on the actual uneven strike grid, with no invented strike prices. If infeasible, report the minimal inconsistent subset or a reproducible ranked list of interval violations and mark that expiry `interval_surface_infeasible`. Do not widen quotes until they fit. A later optional repair may solve a weighted projection and output both original and repaired prices plus adjustment sizes. Call this a model-consistent projection, not market data.

Default downstream policy: retain individually valid quotes from an infeasible expiry in a separately flagged fit block, but do not claim that full interval coverage is possible. Compute the minimum squared interval violation over arbitrage-free discrete call curves using the calibration plan's fixed price scales and block weights, and report this lower bound beside the model loss. Run this check separately on training quotes and on the complete diagnostic view. Exclude expiries whose training-only surface is infeasible from direct variance-curve stripping in the primary construction; use the documented nuisance-curve fallback if too few reliable training terms remain. A failure involving only held-out quotes does not alter training eligibility. A sensitivity may exclude the entire flagged expiry from option fitting, with that choice frozen before optimisation. Neither the lower-bound curve nor an optional projection replaces original market quotes. These finite-grid checks are necessary consistency checks, not proof of a globally valid cross-maturity surface.

Check calendar comparisons only after putting expiries on compatible forward and discount conventions; raw call monotonicity across expiries with changing forwards is not a valid generic test. Keep actual wings visible in diagnostics. If using a declared core moneyness region for pilot fitting (for example SPX $0.70\le K/F\le1.40$, VIX $0.40\le K/F\le3.00$, matching the older script), preserve excluded wings for independent reporting and sensitivity runs. Do not interpret these example bands as liquidity rules.

### 4. Forward variance and VIX benchmarks

From cleaned SPX OTM prices, calculate a **Cboe-style diagnostic** per expiry using observed strikes, strike spacings and supplied (F,T):

\[
\sigma^2(T)=\frac{2}{T}\sum_i\frac{\Delta K_i}{K_i^2}Q_i^U
 -\frac1T\left(\frac{F}{K_0}-1\right)^2.
\]

Here (K_0) is the largest listed strike less than or equal to (F), and (Q_i^U) is the relevant OTM forward price. At (K_0), use the mean of the converted put and call forward prices when both have valid bids and asks, as in the Cboe rule; otherwise flag the diagnostic as incomplete. Use puts below (K_0) and calls above it. Specify the strike-cutoff rule in code and the manifest. Missing IV bids do not prove that the corresponding *cash-price* bid was zero, so the official two-consecutive-zero-bid rule cannot be reconstructed exactly. Use the actual price intervals to produce lower and upper *sensitivity calculations*, noting that independent extremes need not form an arbitrage-consistent surface. Record wing truncation, strike gaps and the number of included quotes.

**Do not label the result an exact official VIX replication.** The official calculation uses minute-accurate maturities and a specific eligible-expiry/strike algorithm. This snapshot has date-level maturity, and its 29-day SPX expiry is followed by a 57-day expiry, so the exact two-term official spot calculation cannot be reconstructed from the available fields. Compare a documented approximate 30-day result with the stored 19.34 VIX spot as a diagnostic only. See the [Cboe VIX methodology](https://cdn.cboe.com/resources/vix/VIX_Methodology.pdf).

For the M2 input, form total implied variance $W(T)=T\sigma^2(T)$ or a documented variance-swap equivalent at reliable **training-only** SPX expiries and strikes. Assign the splits in Section 5 before this calculation; the all-quote diagnostic above is a separate artifact. Exclude both option types at any held-out strike, including a put/call pair needed at (K_0); if this makes the Cboe-style strip incomplete, use a documented training-only log-contract integration or the nuisance-curve fallback. Fit an absolutely continuous term structure with $W(0)=0$ and an evaluatable derivative $\xi_0(t)\ge\underline v+\varepsilon_v$. Pilot defaults are $\underline v=10^{-4}$ and $\varepsilon_v=10^{-8}$ in annualised variance units. Mere positivity or monotonicity of (W) is insufficient for M2. Enforce the bound on entire interpolation intervals and extrapolated tails, not just at knots; specify one-sided derivative values at any knots.

The curve must cover $[0,H]$, where $H=\max(\max T_{SPX},\max(T_{VIX}+30/365))$ over every maturity actually priced, including diagnostics. Record knot maturities, smoothness penalty, extrapolation rule, floor/margin and integration error. Check that integrating $\xi_0$ recovers fitted (W). A single date cannot identify highly flexible interpolation, especially past the last reliable maturity. Use conservative tails and sensitivity variants. If reconstruction cannot be justified from available training strikes, keep $\xi_0$ as a documented low-dimensional nuisance curve and profile it in calibration; record its basis, bounds and penalty, and enforce the same floor and horizon requirements. Do not silently substitute a flat 20% curve or reintroduce held-out quotes.

Extract exactly one VIX forward target for each of the four VIX option expiries: 27, 55, 118 and 153 calendar days; current values are 20.10, 20.73, 21.19 and 21.72. Validate that the forward is constant within expiry and positive. Mark it `likely_traded_future_user_reported` until source provenance is verified. It is a level constraint, never repeated once per option row. If it later proves derived from the same VIX option quotes, do not call it an independent market observation. VIX spot 19.34 remains a snapshot diagnostic.

**Cross-market second-moment gate.** Under the calibration plan's continuous-window VIX definition, for $\Delta=30/365$,

\[
E_0^Q[\mathrm{VIX}_T^2]=\frac{10^4}{\Delta}[W(T+\Delta)-W(T)],
\qquad (F^{VIX}_{model}(T))^2\le E_0^Q[\mathrm{VIX}_T^2].
\]

VIX options further constrain this moment. For any positive pivot (k), with $X=\mathrm{VIX}_T$ and model mean (f),

\[
E[X^2]=2kf-k^2+2\int_0^k P^U(K)\,dK+2\int_k^\infty C^U(K)\,dK.
\]

Use training bid intervals and payoff monotonicity to obtain conservative finite-strike lower bounds on the integrals: on successive put intervals use the left endpoint's bid; on call intervals use the right endpoint's bid. Integrate only covered segments (including the segment to the pivot when justified); omitted tails contribute a nonnegative amount. Do not call a midpoint trapezoidal estimate a rigorous bound, or a finite-strike ask integral an upper bound without tail assumptions. Report the result at the quoted forward and across a declared forward-error sensitivity band; if the model mean differs from the quote, keep the $2kf-k^2$ term rather than substituting $f^2$ with a mismatched pivot.

Compare these bounds with each candidate curve's window variance before freezing it. Separate bid/ask, missing-wing, curve-interpolation, forward-source and settlement-basis uncertainty from numerical error. A conflict requires a documented curve sensitivity or nuisance-curve profile; it is not evidence by itself that the M2 readout is inadequate. If no admissible curve resolves it under declared tolerances, report incompatible inputs/approximation and do not claim full joint feasibility. Nested Monte Carlo and Gaussian quadrature obey the same second-moment identity and cannot remove this conflict. Full-quote checks are diagnostic only and must not tune a training-only curve or select a model using final-test quotes.

### 5. Calibration views and audit report

Write deterministic CSV/JSON outputs under a versioned `data/processed/options/2026-02-19/` directory:

| File | Contract |
| --- | --- |
| `quotes_all.csv` | One row per source contract, all derived values and reason codes. |
| `quotes_calibration.csv` | One row per selected OTM contract, valid price interval, maturity group, surface-status flag and deterministic train/validation/final-test assignment; no duplicate call/put strike observations. |
| `vix_forwards.csv` | One row per VIX expiry, forward in VIX points and source provenance. |
| `spx_variance_terms.csv` | Term variance estimates, bounds/sensitivities, source strike coverage and flags. |
| `forward_variance.csv` | Evaluatable $\xi_0(t)\ge\underline v+\varepsilon_v$ knots or coefficients, or a nuisance-family specification, with horizon, interpolation/extrapolation and training-source metadata. |
| `cross_market_checks.csv` | VIX-window second moments, option-implied lower bounds, forward sensitivities, quote split and compatibility status for each curve variant. |
| `audit.json` | Counts by instrument and expiry, every rejection reason, parity gaps, arbitrage checks and minimum-violation benchmarks, spread/wing distributions, diagnostic VIX comparison. |
| `manifest.json` | Input hashes, code/config versions, conventions and exact run command. |

Before stripping the training curve, assign a deterministic split by zero-based sorted eligible OTM strike rank within each expiry: ranks congruent to 4 modulo 10 are model-selection validation; ranks congruent to 9 modulo 10 are final-test strikes; all others train. This partitions the previous every-fifth-strike holdout into two disjoint roles. Store the exact rule (no randomness is needed), eligibility rules and source row identities; apply each strike's role to both put and call source rows for curve construction. The calibration document fixes the maturity split, with the 119-day SPX expiry wholly reserved for final maturity testing. Keep all strikes in a separate arbitrage audit, but do not use final-test diagnostics to tune the fit. `quotes_calibration.csv` retains original market prices even if a projected curve is optionally produced.

The primary curve excludes validation and final-test source quotes. A separately labelled sensitivity may use a full-data curve, but both strike and maturity scores then constitute **conditional validation**. Freeze model and hyperparameter selection before evaluating final-test quotes; if those results trigger changes, subsequent scores on that set are exploratory. With one date, even an untouched split tests cross-sectional interpolation only, not time stability.

## Acceptance checks and handoff

- Re-running on identical raw files and configuration produces byte-identical data rows and numerically identical audit summaries apart from explicitly labelled run metadata.
- Row-count reconciliation is exact for each source and expiry: retained + excluded = source rows. Every excluded row has a reason. No quote with a missing bid acquires a fabricated midpoint.
- Confirmed SPX IV quotes convert to price intervals in the same units as the stated forward. Round-trip Black inversion recovers the original IV inside a declared tolerance. Parity and convexity diagnostics report, rather than conceal, conflicts.
- Every fitted $\xi_0(t)$ stays above the declared variance floor by at least $\varepsilon_v$ over the full required horizon, has $W(0)=0$, and integrates to the stored term curve within a stated numerical tolerance. Reject positive-but-below-floor curves. Quote and strike uncertainty are kept distinct from numerical integration error.
- Training curve provenance excludes validation/final-test quotes; infeasible-expiry handling and minimum-violation benchmarks are explicit. Cross-market moment gates either pass under declared tolerances or produce an unresolved compatibility status, never an unconditional readiness claim.
- VIX forwards occur once per expiry, in index points. Both VIX IV bounds are converted separately to Black forward-value price bounds, and the confirmed convention is in the manifest.
- The final audit lists all gates and failures; it does not treat an exact spot-VIX replication or a date-specific latent-state estimate as accomplished.

Implement these as targeted unit tests using small synthetic quote tables plus one integration test on a miniature snapshot. The supplied CSVs are data for the final reconciliation, not golden outputs to hard-code. Keep the existing A2 scripts as historical references; use the new interface for the M2 project.
