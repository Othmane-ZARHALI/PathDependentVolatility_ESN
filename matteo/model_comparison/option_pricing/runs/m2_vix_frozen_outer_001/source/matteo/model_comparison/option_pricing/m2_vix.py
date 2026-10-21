"""Unfitted conditional M2 VIX diagnostics on actual training strikes/forwards."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import pandas as pd
import scipy

from .m2_spx import ROOT, M2, frozen_stage3
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.market_targets import MarketTargets, interval_residuals
from esn_eq.pricing import DatedM2Pricer, payoff_statistics
from esn_eq.vix_pricing import ConditionalVIXPricer


def statistics(values):
    return float(np.mean(values)), float(np.std(values, ddof=1)/np.sqrt(len(values)))


def run(args):
    captured = capture_sources([Path(__file__), Path(__file__).with_name("m2_spx.py"),
                               *sorted((M2/"src/esn_eq").glob("*.py"))])
    targets = MarketTargets.load(args.snapshot)
    quotes = targets.select("VIX", "train", args.expiry_days)
    if quotes.empty or set(quotes.expiry_days) != set(args.expiry_days):
        raise ValueError("Every requested VIX expiry needs training quotes")
    architecture, parameters, premium, frozen_hashes = frozen_stage3()
    dated = DatedM2Pricer(architecture, premium)
    initial, _ = dated.initial_state([0., 0., 0.])
    curve = targets.curve(args.curve_id)
    pricer = ConditionalVIXPricer(dated, parameters, curve, initial, normalizer_order=args.normalizer_order)
    summaries, all_prices, all_states = [], [], []
    for maturity, group in quotes.groupby("T_years"):
        maturity = float(maturity)
        days = int(group.expiry_days.iloc[0])
        outer_seed = args.seed+days
        states = dated.prepare([0., maturity], initial, args.outer, outer_seed,
                               checkpoints=[maturity]).checkpoint_states[maturity]
        window = pricer.window(maturity, max_step=1/(252*args.window_steps_per_day))
        evaluations, timings = {}, {}
        for name, options in (
            ("quadrature_low", dict(order=args.order)),
            ("quadrature", dict(order=args.check_order)),
            ("nested_low", dict(engine="nested_mc", inner_paths=args.inner, seed=outer_seed+100000)),
            ("nested", dict(engine="nested_mc", inner_paths=args.refined_inner, seed=outer_seed+100000))):
            start = time.perf_counter()
            evaluations[name] = pricer.evaluate(states, window, **options)
            timings[name] = time.perf_counter()-start
        values = evaluations["quadrature"].vix
        frame = group.copy()
        prices, errors = payoff_statistics(values, group.strike.to_numpy(float), group.option_type.to_numpy())
        frame["model_price_U"], frame["outer_mc_se_U"] = prices, errors
        hits = np.array([np.sum(values > strike) if kind == "call" else np.sum(values < strike)
                         for strike, kind in zip(group.strike, group.option_type)])
        frame["positive_payoff_outer_count"] = hits
        frame["interval_residual"] = interval_residuals(prices, frame.price_bid_U, frame.price_ask_U, frame.price_scale_U)
        for name, result in evaluations.items():
            candidate, _ = payoff_statistics(result.vix, group.strike.to_numpy(float), group.option_type.to_numpy())
            frame[f"{name}_price_U"] = candidate
        all_prices.append(frame)
        state_table = pd.DataFrame({name: value.vix for name, value in evaluations.items()})
        state_table["nested_integral_se"] = evaluations["nested"].integrated_variance_se
        state_table["expiry_days"], state_table["outer_index"] = days, np.arange(args.outer)
        all_states.append(state_table)
        mean, se = statistics(values)
        moment, moment_se = statistics(values**2)
        quad_change, quad_change_se = statistics(values-evaluations["quadrature_low"].vix)
        nested_change, nested_change_se = statistics(evaluations["nested"].vix-evaluations["nested_low"].vix)
        discrepancy, discrepancy_se = statistics(evaluations["nested"].vix-values)
        forward = float(targets.vix_forwards.loc[targets.vix_forwards.expiry == group.expiry.iloc[0], "forward_vix_points"].iloc[0])
        summaries.append(dict(expiry_days=days, T_years=maturity, targets=len(group),
            model_forward=mean, outer_forward_se=se, market_forward=forward,
            model_second_moment=moment, second_moment_outer_se=moment_se,
            curve_second_moment=pricer.moment_target(window),
            forward_quadrature_order_change=quad_change,
            forward_quadrature_order_change_outer_se=quad_change_se,
            forward_inner_refinement_change=nested_change,
            forward_inner_refinement_change_outer_se=nested_change_se,
            nested_minus_quadrature_forward=discrepancy,
            nested_minus_quadrature_forward_outer_se=discrepancy_se,
            max_option_quadrature_order_change=float(np.max(np.abs(frame.quadrature_price_U-frame.quadrature_low_price_U))),
            max_option_inner_refinement_change=float(np.max(np.abs(frame.nested_price_U-frame.nested_low_price_U))),
            max_option_outer_se=float(errors.max()),
            fraction_quotes_passing_nominal_se_screen=float(np.mean(2*errors <= frame.price_scale_U/4)),
            zero_payoff_sample_quotes=int(np.sum(hits == 0)),
            fraction_quotes_passing_mc_gate=float(np.mean((2*errors <= frame.price_scale_U/4) & (hits >= 2))),
            inside_bid_ask_fraction=float(np.mean(frame.interval_residual == 0)),
            window_nodes=len(window.times), wall_seconds=timings))
    record = dict(status="unfitted VIX diagnostic; input and numerical gates unresolved",
        command=sys.argv, quote_split="train_only", calibration_performed=False,
        budget=dict(outer=args.outer, inner=args.inner, refined_inner=args.refined_inner,
                    conditional_orders=[args.order, args.check_order], normalizer_order=args.normalizer_order,
                    window_max_step=1/(252*args.window_steps_per_day), seed=args.seed),
        initial_state_convention="fixed projected C/Y/F=0 reference state, not filtered",
        vix_definition="continuous 30/365-year conditional integrated variance, not settlement replication",
        normalizer_convention="date-zero conditional law fixed across every outer/inner state",
        architecture=asdict(architecture), parameters=asdict(parameters), premium=asdict(premium),
        curve_id=args.curve_id, input_compatibility=targets.audit["variance_audit"],
        target_artifact_sha256=targets.artifact_sha256, frozen_source_sha256=frozen_hashes,
        implementation_sha256=source_hashes(captured, ROOT),
        source_snapshot="source/; captured before calculation, checked unchanged at completion",
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__),
        expiries=summaries,
        limitations=["Shared analytic OU law and readout; engines independently integrate conditional expectation, not independent model specifications.",
                     "Outer SE excludes normalizer, conditional quadrature and time-integration error.",
                     "Nested estimates include finite-inner nonlinear bias; larger outer counts alone cannot remove it.",
                     "Conditional quadrature, inner count, outer count and integration step require separate refinement.",
                     "A zero-hit payoff has zero sample SE but unresolved tail expectation; it cannot pass the MC screen.",
                     "The screen requiring at least two positive payoffs is a minimal diagnostic, not a tail-error guarantee.",
                     "Underlying diagnostic curve remains cross-market inconsistent; no fit or ranking claim."])
    check_sources(captured)
    args.output.mkdir(parents=True, exist_ok=False)
    save_sources(args.output, captured, ROOT)
    (args.output/"run.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")
    pd.concat(all_prices).to_csv(args.output/"prices.csv", index=False)
    pd.concat(all_states).to_csv(args.output/"outer_values.csv", index=False)
    print(json.dumps(summaries, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=ROOT/"data/processed/options/2026-02-19")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--outer", type=int, default=64)
    parser.add_argument("--inner", type=int, default=128)
    parser.add_argument("--refined-inner", type=int, default=512)
    parser.add_argument("--order", type=int, default=12)
    parser.add_argument("--check-order", type=int, default=20)
    parser.add_argument("--normalizer-order", type=int, default=48)
    parser.add_argument("--window-steps-per-day", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--expiry-days", nargs="+", type=int, default=[27,55,118,153])
    parser.add_argument("--curve-id", default="training_mid_piecewise_v1")
    args = parser.parse_args()
    if (min(args.outer,args.inner,args.order,args.normalizer_order) < 2
            or args.refined_inner <= args.inner or args.check_order <= args.order or args.window_steps_per_day < 1):
        parser.error("Require valid counts and strictly increasing inner/quadrature checks")
    if args.output.exists():
        parser.error("Output exists; choose a new directory")
    run(args)


if __name__ == "__main__":
    main()
