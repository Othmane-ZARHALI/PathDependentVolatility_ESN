"""Unfitted dated-M2 SPX pricing diagnostic on common training price targets."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from .provenance import capture_sources, source_hashes, check_sources, save_sources

ROOT = Path(__file__).resolve().parents[3]
M2 = ROOT / "matteo/exogenous_reservoir_eq_m2"
sys.path.insert(0, str(M2 / "src"))
from esn_eq.config import ArchitectureConfig, ReadoutParameters, RiskPremium
from esn_eq.market_targets import MarketTargets, interval_residuals
from esn_eq.pricing import DatedM2Pricer, coarsen_paths, dated_grid, payoff_statistics
from esn_eq.variance_reduction import (deterministic_stock_controls, fit_payoff_controls,
                                      controlled_payoff_statistics)


def frozen_stage3():
    """Load the existing candidate-18/scenario-0 anchor, without fitting or default drift."""
    manifest_path = M2 / "experiments/stage3_001/manifest.json"
    candidate_path = M2 / "experiments/stage3_001/top_configurations.csv"
    manifest = json.loads(manifest_path.read_text())
    table = pd.read_csv(candidate_path)
    row = table[(table.variant == "M2") & (table.candidate_id == 18)
                & (table.structural_scenario_id == 0)].iloc[0]
    base = ArchitectureConfig.from_mapping(manifest["architecture"])
    architecture = replace(base,
        feedback=replace(base.feedback, asset_correlation=float(row.feedback_asset_correlation)),
        spike=replace(base.spike, asset_correlation=float(row.spike_asset_correlation)),
        echo=replace(base.echo, enabled=False))
    parameters = ReadoutParameters(**{key: float(row[key]) for key in (
        "cluster_loading", "feedback_curvature", "feedback_shift", "spike_curvature",
        "cap_level", "cap_sharpness")}, orthogonal_curvature=0., echo_loading=0.)
    premium = RiskPremium(equity=float(row.equity_risk_price),
        feedback_idiosyncratic=float(row.feedback_idiosyncratic_risk_price),
        spike_idiosyncratic=float(row.spike_idiosyncratic_risk_price))
    hashes = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in (manifest_path, candidate_path)}
    return architecture, parameters, premium, hashes


def run(args):
    captured = capture_sources([Path(__file__), *sorted((M2/"src/esn_eq").glob("*.py"))])
    targets = MarketTargets.load(args.snapshot)
    quotes = targets.select("SPX", "train", args.expiry_days)
    if quotes.empty or set(quotes.expiry_days) != set(args.expiry_days):
        raise ValueError("Every requested maturity needs training quotes")
    curve = targets.curve(args.curve_id)
    architecture, parameters, premium, frozen_hashes = frozen_stage3()
    pricer = DatedM2Pricer(architecture, premium)
    initial, prior_penalty = pricer.initial_state(args.initial_projections)
    maturities = np.sort(quotes.T_years.unique())
    events = np.r_[maturities, [k for k in curve.ends if k < maturities[-1]]]
    coarse_grid = dated_grid(events, 1 / (252 * args.steps_per_day))
    grid = np.r_[0., np.concatenate([np.linspace(a, b, args.refinement_factor + 1)[1:]
                                     for a, b in zip(coarse_grid[:-1], coarse_grid[1:])])]
    cache = pricer.prepare(grid, initial, args.paths, args.seed, checkpoints=maturities)
    lower = pricer.spx(cache, parameters, curve, maturities, order=args.quadrature_order)
    result = pricer.spx(cache, parameters, curve, maturities, order=args.check_order)
    pilot_seed = args.seed+100000 if args.pilot_seed is None else args.pilot_seed
    if args.pilot_paths:
        if pilot_seed == args.seed:
            raise ValueError("Control pilot seed must differ from evaluation seed")
        pilot = pricer.prepare(grid, initial, args.pilot_paths, pilot_seed)
        pilot_result = pricer.spx(pilot, parameters, curve, maturities, order=args.check_order)
        reference_pilot = deterministic_stock_controls(pilot, curve, maturities)
        reference_evaluation = deterministic_stock_controls(cache, curve, maturities)
    timestep = (pricer.spx(coarsen_paths(cache, coarse_grid), parameters, curve,
                          maturities, order=args.check_order) if args.refinement_factor > 1 else None)
    frames, summaries = [], []
    for maturity, group in quotes.groupby("T_years"):
        forward = float(group.forward.iloc[0])
        strikes = group.strike.to_numpy(float) / forward
        kind = group.option_type.to_numpy()
        draws = result.terminal_ratios[float(maturity)]
        prices, se = payoff_statistics(draws, strikes, kind, forward_scale=forward)
        coarse, _ = payoff_statistics(lower.terminal_ratios[float(maturity)], strikes, kind, forward_scale=forward)
        frame = group.copy()
        frame["model_price_U"], frame["mc_se_U"] = prices, se
        frame["quadrature_price_change_U"] = prices - coarse
        frame["interval_residual"] = interval_residuals(prices, frame.price_bid_U, frame.price_ask_U, frame.price_scale_U)
        frame["inside_bid_ask"] = frame.interval_residual == 0
        frame["error_to_price_mid_U"] = prices - frame.price_mid_U
        frame["outside_interval_U"] = np.abs(frame.interval_residual) * frame.price_scale_U
        control_summary = {}
        if args.pilot_paths:
            total_variance = curve.integral(0., float(maturity))
            coefficients = fit_payoff_controls(pilot_result.terminal_ratios[float(maturity)],
                reference_pilot[float(maturity)], strikes, kind, total_variance)
            controlled, controlled_se = controlled_payoff_statistics(draws,
                reference_evaluation[float(maturity)], strikes, kind, total_variance,
                coefficients, forward_scale=forward)
            frame["control_price_estimate_U"], frame["control_mc_se_U"] = controlled, controlled_se
            frame["control_beta_stock"], frame["control_beta_reference_option"] = coefficients.T
            frame["control_estimate_negative"] = controlled < 0
            # Never clip a negative unbiased statistical estimate to a model price.
            control_summary = dict(
                control_negative_estimates=int(np.sum(controlled < 0)),
                control_maximum_mc_se_U=float(controlled_se.max()),
                control_median_variance_reduction_factor=float(np.median(se**2/np.maximum(controlled_se**2, 1e-30))),
                control_fraction_quotes_passing_mc_gate=float(np.mean(
                    (2*controlled_se <= frame.price_scale_U/4) & (controlled >= 0))),
                control_inside_bid_ask_fraction=float(np.mean(
                    (controlled >= frame.price_bid_U) & (controlled <= frame.price_ask_U))))
        timestep_summary = {}
        if timestep is not None:
            coarse_draws = timestep.terminal_ratios[float(maturity)]
            paired_mean, paired_se = [], []
            for strike, option_type in zip(strikes, kind):
                sign = 1 if option_type == "call" else -1
                difference = forward * (np.maximum(sign * (draws - strike), 0)
                                         - np.maximum(sign * (coarse_draws - strike), 0))
                paired_mean.append(difference.mean())
                paired_se.append(difference.std(ddof=1) / np.sqrt(len(difference)))
            frame["paired_timestep_change_U"] = paired_mean
            frame["paired_timestep_change_se_U"] = paired_se
            timestep_summary = dict(
                maximum_abs_paired_timestep_change_U=float(np.max(np.abs(paired_mean))),
                maximum_paired_timestep_change_se_U=float(np.max(paired_se)),
                maximum_timestep_change_in_price_scales=float(np.max(np.abs(paired_mean) / frame.price_scale_U)))
            if args.pilot_paths:
                # The deterministic-volatility reference is identical on coupled grids;
                # its payoff cancels. Use the same independently fitted stock coefficient.
                adjusted_means, adjusted_errors = [], []
                for strike, option_type, beta in zip(strikes, kind, coefficients):
                    sign = 1 if option_type == "call" else -1
                    difference = forward*(np.maximum(sign*(draws-strike),0)
                        -np.maximum(sign*(coarse_draws-strike),0)-beta[0]*(draws-coarse_draws))
                    adjusted_means.append(difference.mean())
                    adjusted_errors.append(difference.std(ddof=1)/np.sqrt(len(difference)))
                frame["control_paired_timestep_change_U"] = adjusted_means
                frame["control_paired_timestep_change_se_U"] = adjusted_errors
                control_summary.update(
                    control_maximum_abs_timestep_change_U=float(np.max(np.abs(adjusted_means))),
                    control_maximum_timestep_change_se_U=float(np.max(adjusted_errors)))
        frames.append(frame)
        summaries.append(dict(
            expiry_days=int(group.expiry_days.iloc[0]), T_years=float(maturity),
            targets=len(group), inside_bid_ask_fraction=float(frame.inside_bid_ask.mean()),
            mean_squared_scaled_interval_error=float(np.mean(frame.interval_residual**2)),
            midpoint_price_rmse_U=float(np.sqrt(np.mean(frame.error_to_price_mid_U**2))),
            maximum_price_mc_se_U=float(se.max()),
            maximum_quadrature_price_change_U=float(np.max(np.abs(prices - coarse))),
            maximum_quadrature_change_in_price_scales=float(np.max(np.abs(prices - coarse) / frame.price_scale_U)),
            fraction_quotes_passing_mc_gate=float(np.mean(2 * se <= frame.price_scale_U / 4)),
            ratio_forward_mean=float(draws.mean()),
            ratio_forward_se=float(draws.std(ddof=1) / np.sqrt(len(draws))), **timestep_summary, **control_summary))
    result_record = dict(
        status="unfitted M2 diagnostic; numerical gates and input compatibility unresolved",
        command=sys.argv, reference_date=targets.manifest["reference_date"],
        calibration_performed=False, quote_split="train_only",
        budget=dict(paths=args.paths, seed=args.seed,
                    max_step_years=float(np.max(np.diff(grid))),
                    coarse_max_step_years=1/(252*args.steps_per_day),
                    refinement_factor=args.refinement_factor,
                    grid_steps=len(grid)-1, quadrature_order=args.check_order,
                    comparison_order=args.quadrature_order),
        control_pilot=dict(paths=args.pilot_paths, seed=pilot_seed if args.pilot_paths else None,
                           coefficients="independent pilot; frozen per strike"),
        initial_state_convention="fixed full state, minimum prior norm for supplied C/Y/F projections; not filtered",
        initial_projections_C_Y_F=args.initial_projections, initial_state=initial.tolist(),
        initial_state_mahalanobis_squared=prior_penalty,
        architecture=asdict(architecture), readout=asdict(parameters), premium=asdict(premium),
        curve_id=args.curve_id, variance_floor=1e-4,
        input_compatibility=targets.audit["variance_audit"],
        target_artifact_sha256=targets.artifact_sha256, frozen_source_sha256=frozen_hashes,
        implementation_sha256=source_hashes(captured, ROOT),
        source_snapshot="source/; captured before calculation, checked unchanged at completion",
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__),
        normalizer_max_relative_order_change=float(np.max(np.abs(result.normalizer / lower.normalizer - 1))),
        cap_fraction=result.cap_fraction, expiries=summaries,
        limitations=["No optimization or relative-model performance claim.",
                     "Pricing diagnostics use a curve flagged by preprocessing as cross-market inconsistent.",
                     "Monte Carlo SE is conditional on quadrature and discretization choices.",
                     "Two quadrature orders alone do not certify normalization accuracy.",
                     "Time-step convergence and more independent seeds remain required.",
                     "Forward scaling uses quoted deterministic carry; missing settlement details remain."])
    check_sources(captured)
    args.output.mkdir(parents=True, exist_ok=False)
    save_sources(args.output, captured, ROOT)
    (args.output / "run.json").write_text(json.dumps(result_record, indent=2, allow_nan=False) + "\n")
    pd.concat(frames).to_csv(args.output / "prices.csv", index=False)
    pd.DataFrame(dict(time_years=grid[:-1], xi=curve(grid[:-1]),
        normalizer=result.normalizer, comparison_normalizer=lower.normalizer,
        variance_mean=result.variance_mean, variance_mean_se=result.variance_mean_se)).to_csv(
            args.output / "normalizer.csv", index=False)
    print(json.dumps(dict(normalizer_max_relative_order_change=result_record["normalizer_max_relative_order_change"],
                          expiries=summaries), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=ROOT / "data/processed/options/2026-02-19")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--paths", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--pilot-paths", type=int, default=0)
    parser.add_argument("--pilot-seed", type=int)
    parser.add_argument("--steps-per-day", type=int, default=2)
    parser.add_argument("--refinement-factor", type=int, default=1)
    parser.add_argument("--quadrature-order", type=int, default=20)
    parser.add_argument("--check-order", type=int, default=32)
    parser.add_argument("--expiry-days", nargs="+", type=int, default=[29, 57])
    parser.add_argument("--initial-projections", nargs=3, type=float, default=[0., 0., 0.])
    parser.add_argument("--curve-id", default="training_mid_piecewise_v1")
    args = parser.parse_args()
    if args.paths < 2 or args.steps_per_day < 1 or args.refinement_factor < 1 or args.quadrature_order < 2 or args.check_order <= args.quadrature_order:
        parser.error("Require paths>=2, steps-per-day>=1, and check-order>quadrature-order>=2")
    if args.pilot_paths not in (0,) and args.pilot_paths < 2:
        parser.error("pilot-paths must be zero (disabled) or >=2")
    if args.output.exists():
        parser.error("Output exists; choose a new run directory")
    run(args)


if __name__ == "__main__":
    main()
