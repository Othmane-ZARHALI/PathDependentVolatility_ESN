"""Read-only reproduction of the checked-in ESN option engine.

Run from the repository root with ``python -m
matteo.model_comparison.option_pricing.benchmark --help``. No fitting or
source repairs are performed. Legacy conventions are deliberately preserved.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import importlib
import json
from pathlib import Path
import platform
import sys
import zipfile

import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[3]
LEGACY = ROOT / "CrossSectionalData"
# These are rounded printed values, NOT an archived final optimizer vector.
REPORT_THETA = [0.0369, 0.200, np.log10(1.1e-3), np.log10(1.7e-2),
                0.260, -1.477, 0.803, -0.437]
REFERENCE_THETA = [0.011746584544205775, 0.4127761798792011,
                   -3.7219875698226055, -0.07996888792183689,
                   0.28164092850989053, 0.0, 0.999999990007693,
                   -0.5817249325806491]


@lru_cache(maxsize=1)
def legacy_modules():
    """Load the exact legacy files; reject ambiguous absolute-import collisions."""
    names = ("esn_base", "joint_calibration", "market_vix", "vectorized_sim",
             "model_pricing", "calibrate_smile")
    for name in names:
        if name in sys.modules:
            origin = Path(sys.modules[name].__file__).resolve()
            if origin != LEGACY / f"{name}.py":
                raise RuntimeError(f"Conflicting legacy module {name}: {origin}")
    sys.path.insert(0, str(LEGACY))
    try:
        modules = {name: importlib.import_module(name) for name in names}
    finally:
        sys.path.pop(0)
    return modules


def source_manifest():
    """Compare uploaded archive bytes to current engine bytes without extracting."""
    hashes, archive_match = {}, {}
    with zipfile.ZipFile(LEGACY / "esn_market_test_code.zip") as archive:
        for path in sorted(LEGACY.glob("*.py")):
            raw = path.read_bytes()
            hashes[str(path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
            member = f"code/{path.name}"
            archive_match[path.name] = (archive.read(member) == raw
                                       if member in archive.namelist() else None)
    return {"sha256": hashes, "archive_byte_match": archive_match}


def scalar_batch_check():
    """Same innovations, single path: a stronger oracle than independent moments."""
    m = legacy_modules()
    c = m["calibrate_smile"]
    mat = c.build_matrices_raw(0.386, 0.672, 1e-5, 4.8809, 0.1745, -0.2)
    maximum = {"return_max_abs_error": 0.0, "variance_max_abs_error": 0.0}
    for alpha in (0.22, 1.0):
        for bias in (0.0, -0.437):
            x, v, _ = m["esn_base"]._sim_esn_with_params(
                31415, 128, c.OPT_ARCH, mat, alpha=alpha, b0_delta=bias)
            xb, vb, _ = m["vectorized_sim"].simulate_batch(
                c.OPT_ARCH, mat, alpha, 128, 1, seed=31415, b0_delta=bias)
            maximum["return_max_abs_error"] = max(
                maximum["return_max_abs_error"], float(np.max(np.abs(x - xb[:, 0]))))
            maximum["variance_max_abs_error"] = max(
                maximum["variance_max_abs_error"], float(np.max(np.abs(v - vb[:, 0]))))
    maximum["passed_at_1e-12"] = all(v < 1e-12 for v in maximum.values())
    return maximum


def fit_metrics(model_iv, market_iv, weights):
    """Expose failures and their weight instead of reporting only surviving IVs."""
    valid = np.isfinite(model_iv)
    difference = model_iv[valid] - market_iv[valid]
    mass = float(weights[valid].sum())
    return {
        "n_targets": len(model_iv), "n_finite_model_iv": int(valid.sum()),
        "n_failed_model_iv": int((~valid).sum()),
        "failed_target_weight": float(weights[~valid].sum()),
        "finite_only_weighted_rmse": (float(np.sqrt(
            np.sum(weights[valid] * difference**2) / mass)) if mass else None),
        "finite_only_rmse": float(np.sqrt(np.mean(difference**2))) if valid.any() else None,
        "complete_iv_comparison": bool(valid.all()),
    }


def run(args):
    modules = legacy_modules()
    c = modules["calibrate_smile"]
    # Only replace obsolete input locations, never quote processing or equations.
    c.SPX_CSV, c.VIX_CSV = str(args.spx.resolve()), str(args.vix.resolve())
    calib = c.CombinedCalibrator(K_spx=args.spx_paths, K_vix_outer=args.outer,
                                K_vix_inner=args.inner, burn=args.burn,
                                seed_base=args.seed)
    theta = np.array(REPORT_THETA if args.parameters == "report-rounded" else REFERENCE_THETA)
    p = calib.unpack(theta)
    mat = c.build_matrices_raw(p["H"], p["rough_scale"], p["lam_lo"],
                               p["lam_hi"], p["m1"], p["m2"])
    inputs = {kind: pd.read_csv(path) for kind, path in (("spx", args.spx), ("vix", args.vix))}
    rows, summaries = [], {}
    for spec in calib.specs:
        if spec["kind"] == "spx":
            dx, _, _ = modules["vectorized_sim"].simulate_batch(
                c.OPT_ARCH, mat, p["alpha"], args.burn + spec["n_days"],
                args.spx_paths, seed=args.seed, b0_delta=p["b0_delta"])
            terminal = np.exp(dx[args.burn:].sum(axis=0))
            forward = 1.0
            draws = terminal
            iv_upper_bound = 6.0
        else:
            draws = c.vix_ensemble(mat, p["alpha"], p["b0_delta"], spec["n_days"],
                                   args.burn, args.outer, args.inner, calib.vix_grid, args.seed)
            forward = float(draws.mean())
            iv_upper_bound = 8.0
        model_iv, prices = [], []
        for money in spec["m"]:
            strike = money * forward
            payoff = np.maximum(draws - strike, 0) if money >= 1 else np.maximum(strike - draws, 0)
            price = float(payoff.mean())
            prices.append(price)
            model_iv.append(c.implied_vol_from_price(
                forward, strike, spec["T_years"], price, money >= 1, hi=iv_upper_bound))
        model_iv = np.array(model_iv)
        summary = fit_metrics(model_iv, spec["v"], spec["w"])
        raw = inputs[spec["kind"]]
        expiry = raw.reference.iloc[0] + spec["cal_days"]
        raw = raw[raw.expiry == expiry].copy()
        eligible = raw[(raw.bid.notna()) & (raw.ask.notna())]
        low, high = getattr(c, f'{spec["kind"].upper()}_MONEYNESS_RANGE')
        money_raw = eligible.strike / eligible.forward
        summary.update(
            calendar_days=spec["cal_days"], trading_steps=spec["n_days"],
            diffusion_years=spec["n_days"] / 252,
            black_years=spec["T_years"],
            legacy_moneyness_range=[low, high],
            raw_rows=len(raw), rows_with_bid_and_ask=len(eligible),
            rows_outside_legacy_range=int(((money_raw < low) | (money_raw > high)).sum()),
            model_forward=forward, market_forward=float(raw.forward.iloc[0]),
            forward_comparison_units=("normalized SPX ratio; not absolute index points"
                                      if spec["kind"] == "spx" else "VIX index points"),
            draw_mean=float(draws.mean()),
            draw_mean_standard_error=float(draws.std(ddof=1) / np.sqrt(len(draws))),
        )
        summaries[spec["label"]] = summary
        for money, target, weight, price, iv in zip(spec["m"], spec["v"], spec["w"], prices, model_iv):
            rows.append(dict(label=spec["label"], moneyness=float(money),
                             target_iv=float(target), model_iv=float(iv) if np.isfinite(iv) else None,
                             weight=float(weight), model_price_U=price,
                             model_strike=float(money * forward),
                             iv_failed=not bool(np.isfinite(iv))))
    result = dict(
        status="diagnostic only; no optimizer run or fit reproduction claimed",
        parameter_source=args.parameters,
        parameter_caveat=("rounded report values, not exact archived optimum"
                          if args.parameters == "report-rounded" else "exact theta0 from legacy entry point"),
        theta=theta.tolist(), parameters=p, source=source_manifest(),
        input_sha256={kind: hashlib.sha256(path.read_bytes()).hexdigest()
                      for kind, path in (("spx", args.spx), ("vix", args.vix))},
        versions=dict(python=platform.python_version(), numpy=np.__version__,
                      scipy=scipy.__version__, pandas=pd.__version__),
        budget=dict(spx_paths=args.spx_paths, outer=args.outer, inner=args.inner,
                    burn=args.burn, seed=args.seed),
        scalar_batch_check=scalar_batch_check(), expiries=summaries,
        limitations=["Legacy burn-in state ensemble, not a dated filtered state.",
                     "VIX draws include finite inner-MC and truncated-strike-grid error.",
                     "VIX forward SE measures outer variation at fixed inner count, not total error.",
                     "No legacy scaled/robust optimizer objective reconstructed by this diagnostic.",
                     "Finite-only RMSE is incomplete when any IV inversion fails."],
    )
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "run.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    pd.DataFrame(rows).to_csv(args.output / "points.csv", index=False)
    print(json.dumps({label: {k: v for k, v in s.items() if k in
          ("n_targets", "n_failed_model_iv", "finite_only_weighted_rmse", "model_forward")}
          for label, s in summaries.items()}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spx", type=Path, default=ROOT / "data/raw/SPX_data_19_02_2026.csv")
    parser.add_argument("--vix", type=Path, default=ROOT / "data/raw/VIX_data_19_02_2026.csv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parameters", choices=["report-rounded", "reference"], default="report-rounded")
    for flag, default in (("spx-paths", 512), ("outer", 32), ("inner", 64), ("burn", 500), ("seed", 90000000)):
        parser.add_argument(f"--{flag}", type=int, default=default)
    args = parser.parse_args()
    if min(args.spx_paths, args.outer, args.inner) < 2 or args.burn < 0:
        parser.error("path counts must be >=2 and burn must be nonnegative")
    if args.output.exists():
        parser.error("output already exists; choose a new run directory")
    run(args)


if __name__ == "__main__":
    main()
