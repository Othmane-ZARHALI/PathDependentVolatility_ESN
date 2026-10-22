# Maturity-adjusted central SPX band and coverage — 22 October 2026

User requests a maturity-adaptive criterion and meaningful coverage, suggesting
that including less than 50% would be unsatisfactory. This replaces the fixed
5% reporting band as the priority-membership rule. All quotes and maturities
remain in the full scope; no model fit or pricing simulation is run here.

## Definition and source correspondence

Define `z = abs(log(K/F(T))) / (sigma_ATM_proxy(T)*sqrt(T))` and prioritize **z <= 2**.
Time is in years, using the stored dated T_years. The scale comes from the
**nearest absolute log-moneyness training quote** at each expiry, tie-broken by
strike then option type. Invert its undiscounted mid-price using the existing
Black76 routine. Anchors lie only 0.016–0.077% from the forward, so this is a
near-ATM proxy, not an exactly ATM observation. Retain bid/ask-IV sensitivity.
Do not pool put/call quotes or alter preprocessing. No model-simulated price,
fitted volatility, control performance or held-out quote enters the rule.

The factor sqrt(T) makes the strike band widen with maturity. ATM market IV
also varies with maturity (about 16.07%, 16.01%, 17.07%, 17.57% here). Two units
is a round **central-band convention** chosen to meet the user's coverage
objective. It is broader than strict near ATM. It is not a theorem about where
prices converge, a delta definition, or a claim of Gaussian probability coverage.
A downside smile and non-Gaussian model tails are not represented by this scale.

## Training quote counts

| Maturity | Training quotes | Included, z <= 2 | Coverage | K/F limits |
| --- | ---: | ---: | ---: | --- |
| 29 days | 328 | 179 | 54.6% | 0.9134–1.0948 |
| 57 days | 320 | 196 | 61.3% | 0.8811–1.1349 |
| 211 days | 164 | 88 | 53.7% | 0.7714–1.2963 |
| 302 days | 168 | 103 | 61.3% | 0.7264–1.3767 |
| **All training maturities** | **980** | **566** | **57.8%** | maturity-specific |
| Short maturities, 29/57 days | 648 | 375 | 57.9% | maturity-specific |

This counts option contracts/quotes, not trading volume, open interest, price
weights or all vendor rows. It covers more than half at **every** training
maturity, avoiding an aggregate pass that hides a poorly covered expiry.
Validation/final-test quotes were not used for anchor choice or coverage tuning.

| Cutoff c | All-training included / 980 | Coverage |
| --- | ---: | ---: |
| 0.5 | 186 | 19.0% |
| 1.0 | 350 | 35.7% |
| 1.5 | 468 | 47.8% |
| 1.8 | 530 | 54.1% |
| **2.0** | **566** | **57.8%** |
| 2.5 | 642 | 65.5% |
| 3.0 | 707 | 72.1% |

At c=1.8 the 29-day expiry only just reaches 50%, and drops below it with the
bid-scale anchor. At c=2, using anchor bid/ask IV gives **565–568 / 980** globally,
with over 50% at each expiry for both scales. Bid/ask sensitivity is not a
statistical confidence interval for the true ATM volatility.

## Correction to earlier denominator shorthand

The full SPX **training** surface has **980** quotes at 29/57/211/302 days. The
**648**-quote count used in the recent short-horizon investigation refers to
29/57 days only; it must not be called the full SPX training surface. The
12-contract simulation set and four central representatives were much smaller
samples of those short maturities. None of those experiments priced all 566
quotes in the new priority band, or qualified the two longer SPX maturities.
This distinction is recorded as a scope/reporting correction, without changing
historical results or claiming new convergence.

## Operational interpretation and verification

Use this market-only membership rule to define the next central SPX work.
The four previously prioritized contracts remain representatives inside the
band. Next sampling diagnostics should broaden strike coverage within the band;
29/57-day cost/convergence evidence cannot simply be transferred to 211/302 days.
Keep all wings as secondary diagnostics and all four VIX maturities active.
Negative controlled estimates remain unusable as prices and retained as signed
research output; this selection rule does not guarantee positivity or accuracy.

Run `python -m matteo.model_comparison.option_pricing.spx_scaled_moneyness --output
<fresh-directory>` with the project Python environment. Saved run:
`runs/spx_scaled_moneyness_001`. It contains contract-level training membership,
market-IV anchors, all cutoff counts and input/source hashes. Tests verify
maturity scaling and unit invariance, and check that an adversarial held-out
row leaves anchors and training membership unchanged. Outputs
replay byte-for-byte. This is numerical/source evidence, not a formal theorem.

## Requested wider cutoffs: 2.5, 3 and 3.5

Using identical training quotes and ATM-IV anchors, the coverage is:

| Maturity (training count) | c=2.5 | c=3 | c=3.5 |
| --- | ---: | ---: | ---: |
| 29 days (328) | 202 (61.6%) | 222 (67.7%) | 236 (72.0%) |
| 57 days (320) | 218 (68.1%) | 238 (74.4%) | 247 (77.2%) |
| 211 days (164) | 102 (62.2%) | 113 (68.9%) | 123 (75.0%) |
| 302 days (168) | 120 (71.4%) | 134 (79.8%) | 144 (85.7%) |
| **All (980)** | **642 (65.5%)** | **707 (72.1%)** | **750 (76.5%)** |
| Short 29/57 days (648) | 420 (64.8%) | 460 (71.0%) | 483 (74.5%) |

Run `runs/spx_scaled_moneyness_002` extends the grid through c=3.5. Anchor and
membership files are byte-identical to run001 (the priority flag still refers
to c=2); all previous cutoff rows reproduce exactly. Coverage is monotone in c.
Bid/ask-IV scaling gives 641–643, 705–708 and 750 quotes at these cutoffs.
No priority cutoff was changed by this comparison. Higher coverage includes
more wing quotes and does not establish their numerical pricing accuracy.
