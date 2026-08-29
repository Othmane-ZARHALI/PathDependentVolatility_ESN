"""
final_grid_trf_search.py
======================================
THE METHODOLOGY (final, consolidated):

  1. HYPERPARAMETERS (architecture: Nr, Nz, and 7 shared z-bank
     hyperparameters) are determined by a MINIMISATION ALGORITHM --
     bounded Levenberg-Marquardt (scipy.optimize.least_squares,
     method="trf"), which supports native box bounds and any
     residual/parameter count.

  2. PARAMETERS (p): H and rough_scale are explored via a fixed,
     ADMISSIBLE GRID, restricted to the range realistic for the 5
     market assets being fit (H in [0.03,0.09], rough_scale in
     [0.3,0.6]) rather than the model's full mathematical bounds --
     averaging over the full range mixes calm and extremely turbulent
     simulated regimes and distorts the tail statistics; the realistic
     grid avoids this. The other 10 entries of p are held at neutral
     defaults.

  3. SURROGATE (SMOOTH) SCORE: the 11 stylised facts are computed with
     smooth surrogates in place of every hard order statistic or
     indicator (smooth_stats.py): a log-spaced Hurst regression, a
     Harrell-Davis quantile, two log-sum-exp soft-maxes, a
     sigmoid-relaxed fraction, and a tanh-winsorised kurtosis. This
     eases the optimisation landscape TRF has to navigate.

  For a given hyperparameter candidate: at every point of the
  (H, rough_scale) grid, the ESN is simulated and its smooth stylised
  facts computed; these are AVERAGED across the grid (never
  minimised) into one representative vector; that vector is compared
  against each of 5 real assets to build a 55-entry weighted residual
  vector (5 assets x 11 stylised facts, r_k = sqrt(w_k)*(x_k-c_k)/s_k);
  TRF minimises the sum of squares of this vector over the 9
  hyperparameters.
"""
import sys, time, json, pickle, os, math
sys.path.insert(0, ".")
import numpy as np
from scipy.optimize import least_squares
import esn_base as E
import market_data as M
import smooth_stats as SM

cache = M.build_market_regime_cache(SM.compute_statistics_smooth, verbose=False)
refs = {name: E.make_score_ref(sts) for name, sts in cache.items()}
ASSETS = list(refs.keys())
STAT_KEYS = [k for k, _ in E.SCORE_TERM_WEIGHTS]
WEIGHTS = dict(E.SCORE_TERM_WEIGHTS)
ONE_SIDED = {"q995_vol_ann", "max_vol_ann"}

N_CAL_PATHS = 10    # per explicit request
T_CAL = 1342         # = the market data's own window length: each asset's
                      # full history (8,054 trading days, identical across
                      # all 5 assets) split into 6 non-overlapping windows
                      # of 8054//6=1342 days each (see market_data.py) --
                      # this makes a simulated path directly comparable in
                      # length to the real windows the target references
                      # (c_k, s_k) are themselves built from.
BURN_CAL = 150        # scaled up in proportion to T_CAL

# ALL 12 entries of p are gridded (2 levels each), via a randomised
# 2-level fractional design (the full 2^12=4096-point factorial is
# intractable). WIDER GRID per explicit request: H, lam_lo, lam_hi,
# rough_scale, zr_lo, zr_hi now span a substantially wider admissible
# range than the narrow-bound fix established earlier in this note
# family (which used H=[0.03,0.09], rough_scale=[0.35,0.55], etc.).
# az_lo, az_hi, m1, m2, b0_delta, scale (not specified in this request)
# are kept at the previously-established narrow bounds.
P_BOUNDS_LOHI = dict(
    H=(0.01, 0.3), lam_lo=(1e-5, 1e-1), lam_hi=(1.0, 5.0),
    az_lo=(1/450, 1/350), az_hi=(1/9, 1/6), rough_scale=(0.2, 0.7),
    zr_lo=(0.01, 0.1), zr_hi=(0.1, 0.4), m1=(-0.1, 0.1), m2=(-0.1, 0.1),
    b0_delta=(-0.02, 0.02), scale=(0.98, 1.02),
)
P_KEYS = ["H", "lam_lo", "lam_hi", "az_lo", "az_hi", "rough_scale",
          "zr_lo", "zr_hi", "m1", "m2", "b0_delta", "scale"]
N_GRID_POINTS = 24
_grid_rng = np.random.default_rng(2026)
P_GRID_POINTS = [
    {k: P_BOUNDS_LOHI[k][int(bit)] for k, bit in zip(P_KEYS, row)}
    for row in _grid_rng.integers(0, 2, size=(N_GRID_POINTS, len(P_KEYS)))
]

# Hyperparameter bounds: Nr, Nz continuous over a wide integer range
# (log-encoded), the 7 shared hyperparameters over their own documented
# admissible ranges (esn_base.HP_SPEC).
NR_BOUNDS = (8.0, 200.0)
NZ_BOUNDS = (2.0, 60.0)
HP_KEYS = [name for name, *_ in E.HP_SPEC]         # the 7 shared hyperparameters
HP_BOUNDS = {name: (lo, hi) for name, default, lo, hi, transform in E.HP_SPEC}
HP_ISLOG = {name: (transform == "log") for name, default, lo, hi, transform in E.HP_SPEC}

# Full 9-dim theta = [log(Nr), log(Nz), 7 hyperparameters in THEIR OWN
# units directly, since TRF handles bounds natively -- no need for the
# log/linear HP_SPEC transform used by the CMA-ES version elsewhere;
# TRF's bounds argument enforces the box directly].
THETA_KEYS = ["Nr", "Nz"] + HP_KEYS
THETA_LO = np.array([math.log(NR_BOUNDS[0]), math.log(NZ_BOUNDS[0])] +
                     [HP_BOUNDS[k][0] for k in HP_KEYS])
THETA_HI = np.array([math.log(NR_BOUNDS[1]), math.log(NZ_BOUNDS[1])] +
                     [HP_BOUNDS[k][1] for k in HP_KEYS])


def decode_theta(theta):
    """theta[0:2] = log(Nr), log(Nz); theta[2:] = the 7 shared
    hyperparameters directly in their own units (TRF's bounds keep them
    admissible, so no further transform is needed)."""
    nr = int(round(math.exp(theta[0])))
    nz = int(round(math.exp(theta[1])))
    nr = int(np.clip(nr, *NR_BOUNDS))
    nz = int(np.clip(nz, *NZ_BOUNDS))
    arch = dict(rough_orientation=-1.0, Nr=nr, Nz=nz, matrix_seed=E.MATRIX_SEED)
    for i, k in enumerate(HP_KEYS):
        lo, hi = HP_BOUNDS[k]
        arch[k] = float(np.clip(theta[2 + i], lo, hi))
    return arch


def stats_at_p_point(arch, p_point, seed_base):
    ip = dict(p_point)
    b0_delta = ip.pop("b0_delta")
    scale = ip.pop("scale")
    mat = E.build_esn_matrices(arch, ip)
    if mat["kappa0"] >= 0:
        return None
    stats_accum = {k: [] for k in STAT_KEYS}
    for q in range(N_CAL_PATHS):
        x, v = E._sim_esn_with_params(seed_base + q, T_CAL, arch, mat,
                                       b0_delta=b0_delta, scale=scale)
        st = SM.compute_statistics_smooth(x[BURN_CAL:], v[BURN_CAL:])
        for k in STAT_KEYS:
            stats_accum[k].append(st[k])
    mean_st = {k: float(np.mean(v)) for k, v in stats_accum.items()}
    if any(not np.isfinite(mean_st[k]) for k in STAT_KEYS):
        return None
    return mean_st


_eval_counter = [0]


def grid_averaged_smooth_stats(arch):
    """Grid over ALL 12 entries of p (a randomised 2-level fractional
    design, N_GRID_POINTS combinations); average the SMOOTH stylised
    facts across the grid points -- an average of an average, never a
    minimum."""
    _eval_counter[0] += 1
    seed_base = 800000 + _eval_counter[0] * 100
    per_point_stats = []
    for p_point in P_GRID_POINTS:
        st = stats_at_p_point(arch, p_point, seed_base)
        seed_base += N_CAL_PATHS
        if st is not None:
            per_point_stats.append(st)
    if not per_point_stats:
        return None
    return {k: float(np.mean([s[k] for s in per_point_stats])) for k in STAT_KEYS}


def residual_vector(theta, return_avg_st=False):
    """55-entry weighted residual vector: 5 assets x 11 stylised facts,
    r_k = sqrt(w_k)*(x_k-c_k)/s_k (one-sided terms clipped at 0)."""
    arch = decode_theta(theta)
    avg_st = grid_averaged_smooth_stats(arch)
    if avg_st is None:
        r = np.full(len(STAT_KEYS) * len(ASSETS), math.sqrt(E.SCORE_MAX))
        return (r, None) if return_avg_st else r
    r_full = []
    for asset in ASSETS:
        ref = refs[asset]
        for k in STAT_KEYS:
            w = WEIGHTS[k]
            s = max(ref[k + "_s"], 1e-8)
            z = (avg_st[k] - ref[k + "_c"]) / s
            if k in ONE_SIDED:
                z = max(0.0, z)
            r_full.append(math.sqrt(w) * z)
    r = np.array(r_full)
    return (r, avg_st) if return_avg_st else r


def s_avg_from_stats(avg_st):
    """The Gaussian-kernel S_avg -- the SAME metric plotted in the
    (Nr,Nz) heat map and the 5-seed stability table -- computed from an
    already-averaged stylised-fact vector (no extra simulation)."""
    if avg_st is None:
        return E.SCORE_MAX
    per_scenario = {r: E.score_fn_smooth(avg_st, refs[r]) for r in ASSETS}
    return float(np.mean(list(per_scenario.values())))


def run_trf(theta0, max_nfev=300, diff_step=0.1, ckpt_path=None, max_seconds=None,
            incremental_ckpt=None, base_history_q=None, base_history_s=None):
    """A thin wrapper around least_squares(method='trf') that also logs
    every evaluation -- BOTH Q (the quadratic sum-of-squares LM
    actually minimises) and S_avg (the Gaussian-kernel score used
    everywhere else in this note, e.g.\ the heat map) -- for a
    convergence plot on the heat map's own scale, and supports simple
    checkpointing across multiple invocations. If incremental_ckpt is
    given, the checkpoint is rewritten after EVERY evaluation (not
    just at the end of a chunk), so no progress is lost even if the
    process is killed mid-chunk."""
    history_q = list(base_history_q) if base_history_q else []
    history_s = list(base_history_s) if base_history_s else []
    t_start = time.time()

    def logged_residual(theta):
        r, avg_st = residual_vector(theta, return_avg_st=True)
        q = float(np.sum(r ** 2))
        s_avg = s_avg_from_stats(avg_st)
        history_q.append(q)
        history_s.append(s_avg)
        print(f"  eval {len(history_q):3d}  Q={q:.3f}  S_avg={s_avg:.3f}  "
              f"({time.time()-t_start:.0f}s)", flush=True)
        if incremental_ckpt:
            if q <= min(history_q):   # only overwrite the saved theta if this is the best seen
                with open(incremental_ckpt, "w") as f:
                    json.dump({"theta": theta.tolist(), "history_q": history_q,
                               "history_s": history_s}, f)
            else:
                with open(incremental_ckpt) as f:
                    prev = json.load(f)
                prev["history_q"] = history_q
                prev["history_s"] = history_s
                with open(incremental_ckpt, "w") as f:
                    json.dump(prev, f)
        return r

    res = least_squares(logged_residual, theta0, bounds=(THETA_LO, THETA_HI),
                         method="trf", max_nfev=max_nfev,
                         xtol=1e-12, ftol=1e-12, gtol=1e-12, diff_step=diff_step)
    elapsed = time.time() - t_start
    return res, history_q, history_s, elapsed


def run_trf_chunked(theta0, chunk_nfev, ckpt_path, diff_step=0.1):
    """Warm-restart wrapper: each call runs at most chunk_nfev
    evaluations; the NEXT call re-initialises bounded LM at the last
    theta evaluated so far (a practical approximation to one
    continuous run -- the trust-region radius itself is not
    preserved across calls, but the search direction and history
    are). The checkpoint is rewritten after EVERY evaluation (not
    just at the end of a chunk), since one evaluation now costs ~35s
    (N_CAL_PATHS=10, T_CAL=1342) and can exceed a chunk's own time
    budget -- losing progress mid-chunk would otherwise be silent."""
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            ck = json.load(f)
        theta_cur = np.array(ck["theta"])
        history_q = ck["history_q"]
        history_s = ck["history_s"]
        # BUGFIX: _eval_counter resets to 0 every process start, so
        # without this offset the seeds (and hence every simulated
        # result) exactly repeat the previous chunk's Jacobian probe
        # around the same theta_cur -- the search would silently
        # stagnate, re-evaluating identical points forever. Offsetting
        # by the number of evaluations already done keeps every new
        # evaluation's random seeds genuinely fresh.
        _eval_counter[0] = len(history_q)
        print(f"  resumed from checkpoint: {len(history_q)} evaluations so far, "
              f"best S_avg={min(history_s):.4f}", flush=True)
    else:
        theta_cur = np.array(theta0)
        history_q, history_s = [], []

    res, history_q, history_s, elapsed = run_trf(
        theta_cur, max_nfev=chunk_nfev, diff_step=diff_step,
        incremental_ckpt=ckpt_path, base_history_q=history_q, base_history_s=history_s)

    print(f"  chunk done in {elapsed:.0f}s; {len(history_q)} total evaluations; "
          f"checkpoint at {ckpt_path}", flush=True)
    return res, history_q, history_s


STARTS = [
    # (Nr, Nz, z_strength, even_strength, linear_strength, gamma_norm,
    #  local_z_strength, zz_scale, sign_prob_neg)
    (48, 20, 0.34, 1.50, 0.25, 1.00, 0.030, 0.080, 0.22),   # A2-981003 defaults
    (20, 10, 0.34, 1.50, 0.25, 1.00, 0.030, 0.080, 0.22),   # small reservoir
    (120, 40, 0.34, 1.50, 0.25, 1.00, 0.030, 0.080, 0.22),  # large reservoir
    (48, 20, 0.10, 3.50, 0.60, 0.30, 0.150, 0.020, 0.45),   # strong z-bank, weak gamma
    (48, 20, 1.20, 0.30, 0.05, 2.50, 0.005, 0.300, 0.08),   # weak z-bank, strong gamma
]


def make_theta0(start_row):
    nr, nz = start_row[0], start_row[1]
    hp_vals = dict(zip(HP_KEYS, start_row[2:]))
    return np.concatenate([[math.log(nr), math.log(nz)],
                            [hp_vals[k] for k in HP_KEYS]])


MAX_TOTAL_EVALS = 200   # per explicit request: Figure 1 stops at 200 iterations


if __name__ == "__main__":
    chunk_nfev = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    ckpt_path = sys.argv[2] if len(sys.argv) > 2 else "ckpt_final_trf_widegrid.json"

    # Start at the previously reported A2-981003 hyperparameters and
    # Nr=48, Nz=20 (same starting point as every earlier version of
    # this search, for direct comparability).
    theta0 = np.concatenate([[math.log(48.0), math.log(20.0)],
                              [E.OPTIMAL_ARCH[k] for k in HP_KEYS]])

    # Cap this chunk so the run never overshoots MAX_TOTAL_EVALS.
    n_so_far = 0
    if os.path.exists(ckpt_path):
        with open(ckpt_path) as f:
            n_so_far = len(json.load(f)["history_q"])
    chunk_nfev = min(chunk_nfev, max(0, MAX_TOTAL_EVALS - n_so_far))
    if chunk_nfev == 0:
        print(f"Already at {n_so_far} >= {MAX_TOTAL_EVALS} evaluations -- nothing to do.")
        sys.exit(0)

    res, history_q, history_s = run_trf_chunked(theta0, chunk_nfev, ckpt_path)

    best_arch = decode_theta(res.x)
    final_avg_st = grid_averaged_smooth_stats(best_arch)
    final_r = residual_vector(res.x)
    final_score = float(np.sum(final_r ** 2))
    final_s_avg = s_avg_from_stats(final_avg_st)

    print(f"\nSo far: {len(history_q)}/{MAX_TOTAL_EVALS} evaluations, "
          f"latest status={res.status} ({res.message})")
    print(f"Q(theta) [PRIMARY METRIC, quadratic least-squares]: "
          f"start={history_q[0]:.4f}  best-so-far={min(history_q):.4f}")
    print(f"S_avg (Gaussian, reference only): start={history_s[0]:.4f}  "
          f"best-so-far={min(history_s):.4f}  latest-chunk-end={final_s_avg:.4f}")
    print(f"Latest architecture: Nr={best_arch['Nr']}, Nz={best_arch['Nz']}")
    for k in HP_KEYS:
        print(f"  {k:20s} = {best_arch[k]}")

    json.dump({"history_q": history_q, "history_s": history_s,
               "best_arch": {k: v for k, v in best_arch.items() if k != "matrix_seed"},
               "final_score": final_score, "final_s_avg": final_s_avg,
               "final_avg_stats": final_avg_st,
               "status": res.status, "message": res.message,
               "p_grid_points": P_GRID_POINTS, "n_grid_points": N_GRID_POINTS,
               "n_cal_paths": N_CAL_PATHS, "t_cal": T_CAL,
               "p_bounds": P_BOUNDS_LOHI},
              open("market_final_trf_widegrid.json", "w"), indent=2)
    print("Saved market_final_trf_widegrid.json (checkpoint at", ckpt_path,
          "-- rerun same command to continue)")
