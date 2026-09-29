# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository is a research/validation codebase for **ESN A2**, an Echo State Network
(reservoir-computing) stochastic-volatility model — the downstream/validation half of a broader
research program on ESNs and randomized signatures for path-dependent volatility (see "Research
context" below). It has no application code, no tests, and no package structure — it is five
largely independent, self-contained Python scripts ("protocols"), each of which:

1. simulates a known data-generating process (DGP) that exhibits one stylized fact of financial-market
   volatility,
2. calibrates the ESN (and one or two baseline models) against that DGP through a shared 11-term
   scoring function,
3. produces diagnostic matplotlib figures (PNG) and a LaTeX write-up (`_protocol.tex`, optionally
   compiled to `.pdf` if `pdflatex` is on `PATH`).

Each protocol folder is self-contained: `<script>.py`, a `_protocol.tex` write-up, a `_doc.pdf`
(pre-compiled write-up), and the PNG figures / JSON results the script produces — **all of these
outputs are checked into git**. Running a script with a prefix that collides with an existing
`--out` value will overwrite tracked files.

| Folder | Script | Stylized fact under test | DGP (ground truth) |
|---|---|---|---|
| `GridSearch_HyperParameters_Tuning/` | `esn_hyperparam_search.py` | tunes the ESN's own 7 shared hyperparameters + `(Nr, Nz)` across 4 DGPs at once | log-normal SV fBM (rough + persistent variants), MRW, rough Bergomi |
| `Roughness/` | `roughness_evaluation.py` | Hurst roughness of log-volatility | log-normal SV with fractional Brownian motion |
| `Taylor_Effect/` | `mrw_taylor_comparison.py` | Taylor effect: ACF(\|r\|) > ACF(r²) | Multifractal Random Walk (MRW) |
| `Leverage_Effect/` | `heston_leverage_comparison.py` | leverage: Corr(r_t, V_{t+L}) ≪ 0 | Heston model, strong negative correlation |
| `Zumbach_Effect/` | `zumbach_pdvgl.py` | Zumbach / time-reversal asymmetry | Guyon–Lekeufack path-dependent vol (PDV-GL), 4-factor Markovian form |

Every script except the grid-search one compares three models against the DGP:
- **ESN A2-981003** — the reservoir model under study (this repo's subject).
- **QRH** — quadratic rough Heston, Euler–Volterra scheme with a 252-step ring buffer.
- **PDV-GL** — Guyon–Lekeufack two-factor power-law-kernel path-dependent-vol model. Used as a
  candidate model in most protocols; used as the ground-truth DGP itself in `Zumbach_Effect/`, where
  it is therefore *excluded* from the model comparison (`MODEL_NAMES = ["ESN A2-981003", "QRH"]` there).

## Research context

This codebase is the downstream/validation half of a larger research program: whether Echo State
Networks (ESNs) / randomized signatures can serve as a tractable model for **path-dependent
volatility (PDV)**, targeting joint SPX/VIX calibration and the standard stylized facts (roughness,
exploding ATM skew, spot/vol leverage, the Zumbach effect). The theory/audit trail behind it —
manuscript drafts, peer-review passes, and literature-mapping conversations — lives in
`matteo/web_claude/` (a claude.ai Project export, not code) and is **not duplicated here**; see
`matteo/web_claude/CLAUDE.md` for full navigation, document lineage, and per-conversation summaries.

A few facts from that research are worth carrying over explicitly, since they're easy to
misremember or conflate:

- **The Zumbach effect is not what fits the VIX call smile.** VIX calls need a right-skewed,
  fat-tailed forward-variance law from vol-of-vol/variance-dependent diffusion loading; roughness
  governs the short-maturity SPX skew; Zumbach is a ℙ-measure time-reversal statistic relevant to
  *joint* SPX/VIX coupling and skew *dynamics*, not the static ℚ cross-section of VIX strikes at one
  date. `Zumbach_Effect/zumbach_pdvgl.py` validates the stylized fact in isolation — its existence
  isn't a claim that this alone would fit VIX options.
- A separate, two-layer "Deep ESN" architecture is under active theoretical audit in
  `matteo/web_claude/Analyzing second moment in ESN documents/` (current reference document:
  `deep_esn_master_document.md`). It is **related to but not verified identical to** the ESN
  A2-981003 model in this repo. As of the last recorded audit conversation (2026-07-27), Layer 2
  there was fixed drift-only ("Case A"), which forces an even-component amendment to its activation
  (a purely odd activation provably can't sign-fold — "Prop Z"); three decisions (the target
  even-activation form, a target tail index, and `σ_min`) were still open at that point. Treat this
  as a snapshot that may now be stale — confirm current status before relying on it — and don't
  assume this repo's z-layer (tanh over pairwise products of r-layer taps) has actually been checked
  against those propositions.
- The founding motivation for reservoir/signature methods over directly simulating a DGP (e.g.
  Guyon-style PDV) is a **"simulate once, cheaply re-fit"** property: avoid re-simulating an entire
  log-return path from scratch whenever one downstream parameter changes.

## Working style

Standing instructions from Matteo for this research program (originally given as claude.ai Project
custom instructions; apply them here too):

- **State assumptions explicitly.** Don't silently pick a default for an ambiguous analysis/coding
  decision — say what you're assuming.
- **Ask rather than guess** when a request is underspecified.
- **Literature first**, with correct, checkable references — including non-finance fields (ML,
  neuroscience, rough-path theory) where relevant, not just quantitative finance.
- **Brainstorm before planning, and plan before coding.**
- **OOP + reuse existing, tested code** rather than rewriting from scratch.

Recurring patterns from how this research program is actually run, which carry over to this repo:

- When asked to brainstorm alternatives, don't anchor on the existing approach or document —
  genuinely consider other constructions before settling.
- In any write-up (e.g. a protocol's `_protocol.tex`), distinguish claims established in the
  literature from this project's own derivations explicitly, at the level of individual results —
  not just in a references section at the end.
- Back analytical/closed-form claims with a numerical check where practical — this repo's own DGP
  sanity-check assertions inside each `run()` (see "No tests, no linter" below) are exactly this
  pattern; keep it up in any new script.
- Before introducing a divergence between sibling protocols (a new tolerance, weight, or starting
  point), check whether it's a deliberate choice already covered under "Per-protocol score tweaks"
  below rather than assuming it's an oversight.
- If a `.tex`/`.md` file has a Markdown table with LaTeX math in a cell, escape literal `|`
  characters (`\lvert...\rvert`, `\Vert...\Vert`) — bare pipes break table parsing. This has bitten
  this research program's manuscripts more than once.

## Running the protocols

There is no build step, package manager, or lockfile in this repo. Each script is invoked directly
with `python3`. Dependencies are not pinned anywhere — install manually:

```bash
pip install numpy scipy matplotlib
```

LaTeX (`pdflatex`) is optional; if it isn't on `PATH`, the script still writes the `.tex` file and
just skips the PDF-compile step.

Full run (default budget — can take minutes to hours depending on `n_dgp`/`n_sim`/`T`):
```bash
python3 Roughness/roughness_evaluation.py
python3 Taylor_Effect/mrw_taylor_comparison.py
python3 Leverage_Effect/heston_leverage_comparison.py
python3 Zumbach_Effect/zumbach_pdvgl.py
python3 GridSearch_HyperParameters_Tuning/esn_hyperparam_search.py
```

Every script accepts `--fast` for a quick, reduced-budget smoke test (small `n_dgp`/`n_sim`/`T`,
fewer grid points) — use this to sanity-check a change before committing to a full run:
```bash
python3 Taylor_Effect/mrw_taylor_comparison.py --fast
```

Common CLI flags (all scripts): `--out <prefix>` (output file prefix), `--n_dgp`/`--n_sim`/`--n_cal`
(paths for the DGP reference / final model comparison / inner calibration loop), `--T`, `--burn`,
`--dt`. Each script also exposes a flag for its own swept parameter: `--H`/`--lam` (Roughness),
`--lam2` (Taylor), `--rho` (Leverage), `--H` meaning PDV-GL's `alpha_r` (Zumbach). Note that the
figures/tex/json checked into each folder were produced with a **non-default `--out` prefix**
matching that specific historical run (e.g. `roughness_h05_h10_out_*`, `taylor_lam01_03_*`,
`heston_leverage_extreme_*`, `zumbach_final_out_*`) — running a script with its default prefix will
not clobber them.

`esn_hyperparam_search.py` has extra modes beyond a plain CMA-ES run:
```bash
python3 GridSearch_HyperParameters_Tuning/esn_hyperparam_search.py --grid_matrix_only    # re-evaluate a saved/default hyperparameter set over the full (Nr,Nz) grid only
python3 GridSearch_HyperParameters_Tuning/esn_hyperparam_search.py --full_grid_report    # like above, plus regenerates the overview/deviation figures
```

## No tests, no linter

There is no test suite, CI config, linter, or formatter in this repo. The only correctness checks
available are: running a script (ideally with `--fast`) and inspecting the console summary table it
prints per DGP cell (`H_hat`, vol%, the stylized-fact-specific columns, and score `S` per model —
lower is better, `S=0` is a perfect match); the sanity-check assertions each `run()` embeds on the
DGP itself before calibrating anything against it (e.g. `Leverage_Effect` asserts the Heston DGP's
leverage is negative; `Taylor_Effect` asserts the MRW DGP's Taylor gap/fraction sit in the expected
direction); and the generated PNG figures.

## Architecture

### The pattern repeated in every protocol script

Each of the 5 scripts is built from the same numbered sections (see the `# === N. ... ===` headers
inside each file) and, critically, **each script embeds its own copy of the shared machinery below
instead of importing it from a common module**. Docstrings often claim a function is "identical to"
or "unchanged from" another script's — true for the core of `compute_statistics`, but score
weights/tolerances and calibration starting points have been independently tuned per protocol (see
"Per-protocol score tweaks" below). Never assume the copies are byte-identical; diff before porting a
fix from one file to another.

1. **Constants** — `TRADING_DAYS`, `TV_ANN`/`TV_DAY` (target annualised/daily vol, 20%), `SIG_MIN`,
   and an `_ARCH`/`OPTIMAL_ARCH` dict holding the ESN's 7 shared hyperparameters + reservoir sizes
   (`n_r`/`Nr`, `n_z`/`Nz`). Reservoir sizes intentionally differ by file: the grid-search script's
   `OPTIMAL_ARCH` uses Nr=64, Nz=8 (the reported best cell of its own sweep); the four downstream
   comparison protocols use Nr=96, Nz=12. This is not an inconsistency to "fix".
2. **`compute_statistics(daily_x, daily_var)`** — the shared 11-term stylized-fact vector computed
   from any simulated (or DGP) path: `H_hat` (Hurst exponent, via dyadic log-vol variogram OLS
   slope), `mean_vol_ann`/`q995_vol_ann`/`max_vol_ann`, `mean_vol_acf`, `taylor_gap`/`taylor_frac`
   (ACF(|r|) vs ACF(r²)), `zumbach` (windowed past-magnitude→future-vol vs past-vol→future-magnitude
   asymmetry), `leverage` (Corr(r_t, V_{t+L})), `kurtosis`, `max_ret_acf`.
3. **`make_score_ref(dgp_stats)` / `score_fn(st, ref)`** — turns the DGP's own cross-path mean/std
   into a data-adaptive calibration target, then scores any candidate path against it. `S >= 0`, is
   **minimised**, and `S=0` ⇔ a perfect match on all 11 terms with no stress event. Uses a **smooth
   Gaussian-kernel proximity** `f(x,c,s) = exp(-0.5*((x-c)/s)^2)` (see `esn_hyperparam_search.py`'s
   `score_fn_tent` vs `score_fn_smooth` docstrings: the original piecewise-linear "tent" kernel hits
   exactly zero gradient beyond one tolerance width, which is why Nelder-Mead calibration used to
   stall well short of `S=0`). A hard `+5.0` "stress" penalty fires on blow-up/collapse
   (`max_vol_ann` too high, or `mean_vol_ann` outside `[5%, 150%]`).
4. **DGP simulator(s)** — one per protocol (`dgp_lnsv`, `dgp_mrw`, `dgp_heston`, `dgp_pdvgl`, or the
   4-scenario `sim_dgp` dispatcher in the grid-search script). fBM-driven DGPs build a Cholesky factor
   of the fractional covariance matrix (cached in a module-level dict keyed by `(T, H, dt)`), and
   **the time argument must be annualised** (`t * TIME_DT`, `TIME_DT = 1/252`). Using the raw day
   index instead is a real historical bug, documented at length in `esn_hyperparam_search.py`'s
   module docstring: it made `log sigma_t`'s Itô correction term blow up for `H` not close to 0,
   silently producing NaN/garbage downstream statistics and a meaningless calibration target. If you
   add a new fBM-based DGP, follow the existing `_get_fbm_chol`/`TIME_DT` pattern, not raw day counts.
5. **Model simulators** — `_sim_esn`/`_sim_esn_with_params` (the ESN, see below), `_sim_qrh`
   (Euler–Volterra rough Heston with a ring buffer), `_sim_gl`/`_sim_gl_4factor` (PDV-GL, power-law or
   4-factor-Markovian OU-mixture form).
6. **`quick_calibrate(dgp_sts, ...)`** — builds the score reference from the DGP, then fits each model
   to it via **multi-start Nelder-Mead** (`scipy.optimize.minimize(method="Nelder-Mead")`), keeping
   the best of several hand-picked starting points (`ESN_INNER_STARTS`/`INNER_NM_STARTS`,
   `GL_STARTS`, `QRH_STARTS`). A single fixed start systematically misses persistent/slow-reservoir
   optima — see the rationale comment above `INNER_NM_STARTS` in `esn_hyperparam_search.py`.
7. **Figures** — matplotlib, `Agg` backend, saved as `{save_prefix}_<fig>.png`, one function per
   figure (e.g. `plot_L1_leverage_curves`, `plot_T1_taylor_curves`, `plot_hurst_estimation`,
   `plot_Z1_curves`).
8. **`write_latex(path)` + `run(...)`** — assembles the `_protocol.tex` write-up and, if `pdflatex`
   is on `PATH`, compiles it twice (for cross-references) to PDF.
9. **CLI (`if __name__ == "__main__":`)** — `argparse`, with a `--fast` flag that shrinks every
   budget knob for a smoke test.

### The ESN model ("ESN A2-981003")

Two coupled internal state banks, built by `_build_esn`/`build_esn_matrices` and stepped by
`_sim_esn`/`_sim_esn_with_params` (see `GridSearch_HyperParameters_Tuning/esn_persist_protocol.tex`
§1 for the authoritative write-up):

- **r-layer**: `Nr` discrete-time OU processes sharing one scalar innovation per step, mean-reversion
  rates geometrically spaced over calibrated bounds `[lam_lo, lam_hi]`. Projected onto a "rough
  factor" direction `q = lambdas^(0.5-H) / ‖·‖` — `H` here is the ESN's *own* internal memory-kernel
  exponent, calibrated so the simulated `H_hat` matches the DGP's measured Hurst exponent; it is a
  different object from the DGP's true Hurst parameter, despite the shared name.
- **z-layer**: `Nz` slower nonlinear feedback modes (`tanh` activation) reading pairwise products of
  r-layer taps, driving multiscale/Zumbach-type effects.
- **Read-out**: `eta_t = b0 + rough_orientation*rough_scale*(q·r_t) + sum_j w_{j,z}*z_{j,t} + r_tᵀ Q r_t`,
  with `Q = (1/Nr)*(m1*I + m2*q qᵀ)` a symmetric (not necessarily positive-definite) quadratic term.
  `sigma_t = sqrt(SIG_MIN² + softplus(eta_t)²)`.
- **Hyperparameters vs. calibrated parameters**: 7 hyperparameters (`z_strength`, `even_strength`,
  `linear_strength`, `gamma_norm`, `local_z_strength`, `zz_scale`, `sign_prob_neg`) plus `Nr`/`Nz` are
  shared across every DGP scenario and found by the outer CMA-ES search in `esn_hyperparam_search.py`;
  12 parameters (`H`, `lam_lo`, `lam_hi`, `az_lo`, `az_hi`, `b0_delta`, `scale`, `rough_scale`,
  `zr_lo`, `zr_hi`, `m1`, `m2`) are calibrated **per DGP scenario** by the inner Nelder-Mead loop
  (`inner_calibrate`/`quick_calibrate`).
- The outer optimiser (`CMAESOptimiser` in `esn_hyperparam_search.py`) is a small hand-rolled
  implementation, not the `cma` PyPI package — no extra dependency needed for it.
- `rough_orientation` is fixed to `-1.0` by an admissibility constraint (`kappa0 < 0`, required for a
  negative leverage effect); scripts assert this immediately after building `_ARCH`.

### Per-protocol score tweaks

Although `score_fn`'s formula is shared, each downstream protocol re-tunes the weight/tolerance of
the term it's investigating, so the target a model is chasing is not identical across files:
- `Leverage_Effect`: leverage weight raised 1.1 → 3.5, its tolerance tightened to 0.15× (from 0.30×)
  the DGP-centre.
- `Taylor_Effect`: `taylor_gap`/`taylor_frac` weights raised to 3.5/2.5 (from 1.0/0.8), tolerances
  tightened to 0.15×.
- `Zumbach_Effect`: the single `zumbach` difference term is split into two independently-scored
  channels (`zumbach_term1`, `zumbach_term2`, weight 1.75 each) so a model can't match the
  *difference* while getting both underlying channels wrong.

Check whether a protocol already has its own deliberate deviation before porting a score change from
one file into another.


# IMPORTANT
- Always share with me if you are making assumptions and what you are assuming. Do not silent any choice you are taking for analysis/coding, and so on.
- Don't be shy: always ask me questions if you have doubts or I haven't specified something well enough.
- Give priority to scientific literature and always report the correct references
- Try to brainstorm before writing a plan (which then could be used to write code).
- FOR CODING: let's try to stay compliant with OOP and to reuse code that we have already written/tested.
- Try to use ASD-STE-100 when you speak to the operator, but you're welcome to use specific or technical lingo when needed, just introduce/define it the first time you use it. Same holds true for acronyms: remember to define them at least the first time you use them.