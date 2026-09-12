# ESN exponential–quadratic model: experiment comparison and current conclusion

## Executive conclusion

The modified M2 architecture can jointly reproduce the tested daily-frequency physical-measure stylised facts on synthetic paths. The successful region is not one lucky parameter vector: in the untouched Stage-3 design, M2 passed the declared joint proxy in 62.0% of all rows, 75.8% of configurations passed at least three of five independent cache seeds, 69.5% passed at least four, and five configurations passed all five.

The decisive architecture change was not a quintic readout or a nonlinear echo layer. It was a slower, slow-mode-weighted quadratic spike lift combined with a low-curvature, high-shift feedback quadratic. This remains an exponential–quadratic model descended from ESN-A2 and preserves path reuse.

M2 is the recommended baseline. The standardized orthogonal-energy extension (M2+O) is optional: it improves roughness and raises joint acceptance modestly, but slightly weakens leverage and Zumbach. M3 is not recommended; consequently M4 is not recommended either.

## Selection and validation chronology

| Stage | Role | Independent design | Main result |
|---|---|---|---|
| Pilot 1 | Broad exploration | 3 cache seeds, 16 readout draws, 3 scenarios | Broad proxy pass around 29%, but strict joint pass 0%; tails and tight roughness were weak. |
| Stage 2 | Untouched test of pilot-selected high-curvature box | 5 new seeds, 32 new draws, 4 scenarios, 16 paths | 0% strict joint pass. Tails became strong, but mean H fell near 0.04 and leverage/Zumbach remained seed-sensitive. |
| Mechanism probes | Training only | Reused Stage-2 seeds and common shocks | Identified slow spike rates, slow-mode weighting, and a high-shift/low-curvature leverage ridge. |
| Stage 3 | Untouched confirmation | 5 new seeds, 32 new draws, 4 new scenarios, 24 paths | 62.0% M2 and 64.7% M2+O joint acceptance; stable configurations found across all five seeds. |

No Stage-3 seed, Sobol point, or structural scenario was used in the mechanism probes. The primary Stage-2 proxy gates were retained unchanged in Stage 3.

## What failed in Stage 2

The aggressive Stage-2 box increased feedback and spike curvatures. That solved finite-sample tail thickness—M2 mean excess kurtosis rose from 2.37 in Pilot 1 to 8.50—but it placed too much variance energy at fast scales:

| Model | Joint pass | Rough H | Rank leverage | Rank Zumbach | Excess kurtosis |
|---|---:|---:|---:|---:|---:|
| Stage-2 M2 | 0.0% | 0.0408 | -0.00582 | 0.0152 | 8.50 |
| Stage-2 M2+O | 0.0% | 0.0482 | -0.00531 | 0.0149 | 8.75 |

The standardized orthogonal energy was no longer inert and raised H by about 0.0074, but it could not repair the box. The sampled echo effect was effectively zero. These results reject “increase every curvature” as a general calibration strategy.

## Why the Stage-3 architecture works

The feedback term can be expanded as

\[
q_Z(Y_t-\beta)^2=q_ZY_t^2-2q_Z\beta Y_t+q_Z\beta^2.
\]

This exposes two useful roles. A small \(q_Z\) preserves slow roughness; a larger \(\beta\) creates a strong negative linear feedback coefficient \(-2q_Z\beta\) for leverage. The constant is largely removed by Q normalisation when the cap is inactive.

The original spike grid, 20–2000 inverse years, was too fast for the 5/10/20-day Zumbach windows. Stage 3 uses four geometric rates from 10 to 500 and projection power -0.5, which weights slower modes more strongly. The quadratic spike channel then supports a large positive time-reversal asymmetry and heavy tails while the shifted feedback channel supplies leverage.

The frozen Stage-3 readout box was:

| Quantity | Range |
|---|---:|
| Cluster loading | 0.20–0.30 |
| Feedback curvature | 0.060–0.085 |
| Feedback shift | 2.70–3.30 |
| Spike curvature | 0.18–0.26 |
| Orthogonal curvature, optional | 0.02–0.25, log sampled |
| Feedback/equity correlation | 0.970–0.999 |
| Spike/equity correlation | 0.940–0.995 |

Risk-premium parameters remained fixed; the experiment varied only physical-measure correlations and cheap readout parameters.

## Stage-3 results

### Fixed-design means across five independent cache seeds

| Model | Joint pass | Rough H | Volatility ACF | Rank leverage | Rank Zumbach | Excess kurtosis | Taylor gap | Max return ACF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| M2 | 62.0% | 0.0754 | 0.209 | -0.0175 | 0.0416 | 9.84 | 0.0290 | 0.0754 |
| M2+O | 64.7% | 0.0777 | 0.218 | -0.0162 | 0.0395 | 9.90 | 0.0308 | 0.0757 |

All M2 rows passed the tail, Taylor, and return-ACF gates. The volatility-ACF gate passed 98.8%, leverage 95.8%, Zumbach 97.0%, and roughness 66.7%. Thus roughness, especially on one adverse seed, is the remaining bottleneck.

### Same-configuration seed stability

| Model | Configurations | Pass ≥3/5 seeds | Pass ≥4/5 | Pass 5/5 |
|---|---:|---:|---:|---:|
| M2 | 128 | 75.8% | 69.5% | 3.9% (5 configurations) |
| M2+O | 128 | 78.9% | 70.3% | 5.5% (7 configurations) |

One untouched all-seed M2 configuration had

\[
c_C=0.2333,\quad q_Z=0.06940,\quad \beta=2.9359,\quad q_J=0.2496,
\]

with feedback/equity correlation 0.9963 and spike/equity correlation 0.9791. Its five-seed means were H 0.0845, rank leverage -0.0151, rank Zumbach 0.0481, excess kurtosis 10.28, volatility ACF 0.207, and max return ACF 0.0776. Every individual seed passed the declared joint proxy; the worst H was 0.0701 and the weakest Zumbach 0.0332.

### Strong-Zumbach sensitivity

For M2, 75.2% of rows had rank Zumbach above 0.03 and 28.9% exceeded 0.05. Requiring all other joint gates plus Zumbach above 0.03 retained 54.1% of rows; requiring above 0.05 retained 26.2%. This is materially stronger than merely obtaining the correct sign, although the between-seed variation remains non-negligible.

### Optional orthogonal block

The paired M2+O effect on acceptance was +2.66 percentage points, with a five-seed Student interval of +0.13 to +5.19 points. It raised mean H by 0.00229 and the volatility ACF by 0.0086, but changed rank leverage by +0.00128 (less negative) and rank Zumbach by -0.00206. The improvement is real but small; M2+O is best viewed as an optional roughness stabilizer, not the default architecture.

## What is and is not established

Established within this synthetic daily experiment:

- simultaneous finite-band roughness, volatility clustering, negligible return ACF, negative leverage, positive and sometimes strong Zumbach asymmetry, fat returns, and a positive Taylor gap;
- a non-trivial local parameter volume rather than only one optimized point;
- exact common-random-number path reuse for cheap readout changes and deterministic P/Q shifts;
- no need for a quintic readout or the nonlinear echo layer for these P-measure results.

Not established:

- empirical fit to SPX returns, realised variance, or intraday horizons;
- SPX/VIX option-surface calibration under Q, since no option dataset or pricing loss was available;
- accurate conditional VIX pricing—the current Q normalizer is cross-sectional Monte Carlo;
- an out-of-sample estimate of Q-normalisation error. The recorded forward-mean error is algebraically near zero because the same paths estimate and verify the normalizer;
- asymptotic power-law tails or genuine zero-scale roughness;
- formal identification of the latent lift, whose original calibration still hits several search bounds.

## Recommendation

Use the slower-spike/high-shift M2 as the next baseline. Keep standardized M2+O behind a feature flag for roughness-sensitive experiments. Leave M3 off and do not call the combined model M4 unless a future out-of-sample study finds a material echo benefit.

The next scientifically useful test is frequency robustness on intraday simulations and empirical daily targets. After option data become available, retain the same exogenous path caches and calibrate deterministic P-to-Q risk prices and Q readout parameters on a training surface, then evaluate SPX and VIX smiles jointly on untouched dates.
