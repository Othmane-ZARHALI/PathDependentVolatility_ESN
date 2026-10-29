"""Score the unfitted M2 anchor with archived SPX and fresh conditional VIX.

This is a restricted-maturity plumbing and integration diagnostic. Heterogeneous
numerical budgets remain visible; it is neither a fitted model nor a ranking.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd

from .m2_spx import ROOT, M2, frozen_stage3
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.calibration_loss import JointPriceObjective, LossConfig
from esn_eq.market_targets import MarketTargets
from esn_eq.pricing import DatedM2Pricer, payoff_statistics
from esn_eq.vix_pricing import ConditionalVIXPricer


def run(args):
    start = time.perf_counter()
    sources = capture_sources([Path(__file__), Path(__file__).with_name("m2_spx.py"),
                               *sorted((M2/"src/esn_eq").glob("*.py"))])
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    architecture, parameters, premium, hashes = frozen_stage3()
    dated = DatedM2Pricer(architecture, premium)
    initial, _ = dated.initial_state([0., 0., 0.])
    # Reuse expensive high-resolution SPX data only after checking its context.
    archive = Path(__file__).with_name("runs")/"m2_spx_128_256_001"
    previous = json.loads((archive/"run.json").read_text())
    expected = json.loads(json.dumps(dict(architecture=asdict(architecture), parameters=asdict(parameters),
        premium=asdict(premium), initial_state=initial.tolist(), target_sha256=targets.artifact_sha256)))
    if any(previous[key] != value for key, value in expected.items()):
        raise ValueError("Archived SPX model, state or target snapshot changed")
    # Only VIX time integration and scoring changed; retain the archived SPX implementation.
    for name in ("pricing.py", "variance_reduction.py", "states.py", "features.py", "config.py"):
        path = M2/"src/esn_eq"/name
        if previous["implementation_sha256"][str(path.relative_to(ROOT))] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError(f"Archived SPX pricing dependency changed: {name}")
    spx = targets.select("SPX", "train", [29, 57])
    keys = ["instrument", "expiry", "strike", "option_type"]
    old = pd.read_csv(archive/"prices.csv").set_index(keys)
    contracts = pd.MultiIndex.from_frame(spx[keys])
    if len(old) != len(spx) or not old.index.is_unique:
        raise ValueError("Archived SPX contracts differ")
    old = old.loc[contracts].reset_index()
    for name in ("price_bid_U", "price_ask_U", "price_scale_U", "forward", "T_years"):
        np.testing.assert_allclose(old[name], spx[name], rtol=0, atol=1e-10)
    prices = [old.fine_control_price_U.to_numpy()]
    spx["mc_se_U"] = old.fine_control_se_U.to_numpy()
    spx["positive_payoff_count"] = old.fine_positive_payoff_count.to_numpy()
    spx["paired_grid_change_U"] = old.control_paired_change_U.to_numpy()
    spx["paired_grid_change_se_U"] = old.control_paired_change_se_U.to_numpy()
    spx["numerical_source"] = "archived SPX 256/day; paired 128/day; 4096 eval + 2048 independent pilot"
    tables, reports, state_tables = [spx], [], []
    forwards = targets.vix_forwards.sort_values("expiry").copy()
    predicted_forwards = {}
    curve = targets.curve("training_mid_piecewise_v1")
    engine = ConditionalVIXPricer(dated, parameters, curve, initial, normalizer_order=48)
    for maturity, group in targets.select("VIX").groupby("T_years"):
        tick = time.perf_counter()
        maturity, expiry = float(maturity), int(group.expiry.iloc[0])
        seed = args.seed+int(group.expiry_days.iloc[0])
        states = dated.prepare([0., maturity], initial, args.outer, seed,
                              checkpoints=[maturity]).checkpoint_states[maturity]
        low_time = engine.gaussian_window(maturity, time_order=8)
        high_time = engine.gaussian_window(maturity, time_order=16)
        low_order = engine.evaluate(states, low_time, order=20).vix
        base = engine.evaluate(states, low_time, order=32).vix
        high = engine.evaluate(states, high_time, order=32).vix
        strikes, kinds = group.strike.to_numpy(), group.option_type.to_numpy()
        price, se = payoff_statistics(base, strikes, kinds)
        time_price, _ = payoff_statistics(high, strikes, kinds)
        order_price, _ = payoff_statistics(low_order, strikes, kinds)
        group = group.copy()
        group["mc_se_U"] = se
        group["positive_payoff_count"] = [np.sum(base > k) if kind == "call" else np.sum(base < k)
                                           for k, kind in zip(strikes, kinds)]
        group["time_order_change_U"] = time_price-price
        group["conditional_order_change_U"] = price-order_price
        group["numerical_source"] = f"fresh VIX Gaussian time 8/16; conditional 20/32; iid outer {args.outer}"
        tables.append(group); prices.append(price)
        predicted_forwards[expiry] = float(base.mean())
        report = dict(expiry=expiry, expiry_days=int(group.expiry_days.iloc[0]), seed=seed,
            outer=args.outer,time_nodes=[len(low_time.times),len(high_time.times)],
            model_forward=float(base.mean()),forward_se=float(base.std(ddof=1)/np.sqrt(args.outer)),
            forward_time_change=float(np.mean(high-base)),
            forward_conditional_order_change=float(np.mean(base-low_order)),
            max_state_time_change=float(np.max(np.abs(high-base))),
            max_option_time_change=float(np.max(np.abs(time_price-price))),
            max_option_conditional_order_change=float(np.max(np.abs(price-order_price))),
            max_option_mc_se=float(se.max()),zero_hit_quotes=int(np.sum(group.positive_payoff_count == 0)),
            wall_seconds=time.perf_counter()-tick)
        reports.append(report)
        state_tables.append(pd.DataFrame(dict(expiry=expiry,outer_index=np.arange(args.outer),
            base_vix=base,time_refined_vix=high,conditional_low_order_vix=low_order)))
        print(json.dumps(report), flush=True)
    q = pd.concat(tables, ignore_index=True)
    config = LossConfig()
    result = JointPriceObjective(q, forwards, config).evaluate(np.concatenate(prices),
        [predicted_forwards[int(e)] for e in forwards.expiry])
    record = dict(status="unfitted restricted-maturity diagnostic; numerical and input gates unresolved",
        calibration_performed=False,quote_split="train_only",spx_expiry_days=[29,57],vix_expiry_days=[27,55,118,153],
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        initial_state=initial.tolist(),state_convention="fixed zero projections, not filtered",
        curve_id="training_mid_piecewise_v1",input_compatibility=targets.audit["variance_audit"],
        loss_config=asdict(config),loss_components=result["components"],total_loss=result["total"],
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(sources, ROOT),
        archived_spx=str(archive.relative_to(ROOT)),
        archived_spx_sha256={name:hashlib.sha256((archive/name).read_bytes()).hexdigest()
                             for name in ("run.json","prices.csv")},
        vix=reports,normalizer_order=48,conditional_orders=[20,32],time_orders=[8,16],
        wall_seconds=time.perf_counter()-start,
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),
        limitations=["Only 29/57-day SPX; 211/302-day training expiries are not scored here.",
            "Loss values are noisy descriptive scores, not fit, feasibility or comparative ranking evidence.",
            "VIX MC SE excludes shared normalizer, quadrature and time errors; tail expectations remain unresolved.",
            "State/curve penalties are zero because their parameters are fixed, not profiled.",
            "SPX and VIX use different path budgets; numerical errors are reported separately and never reweight quotes."])
    check_sources(sources)
    args.output.mkdir(parents=True, exist_ok=False)
    save_sources(args.output, sources, ROOT)
    (args.output/"run.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")
    result["quotes"].to_csv(args.output/"price_residuals.csv", index=False)
    result["forwards"].to_csv(args.output/"forward_residuals.csv", index=False)
    result["blocks"].to_csv(args.output/"blocks.csv", index=False)
    pd.concat(state_tables, ignore_index=True).to_csv(args.output/"vix_outer_values.csv", index=False)
    print(json.dumps(dict(total=result["total"],components=result["components"])), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--outer", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20261021)
    args = parser.parse_args()
    if args.output.exists() or args.outer < 2:
        parser.error("Require a new output directory and at least two outer states")
    run(args)
