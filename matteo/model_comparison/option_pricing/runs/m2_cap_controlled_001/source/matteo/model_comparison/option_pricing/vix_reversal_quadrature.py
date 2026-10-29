"""Matched-state quadrature attribution at the reversal's 153-day VIX block."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .m2_spx import ROOT,M2,frozen_stage3
from .provenance import capture_sources,source_hashes,check_sources,save_sources
from esn_eq.market_targets import MarketTargets,interval_residuals
from esn_eq.market_calibration import NuisanceCurveFamily
from esn_eq.outer_sampling import GaussianOuterSampler
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def run(args):
    sources = capture_sources([Path(__file__),Path(__file__).with_name("m2_spx.py"),
                               *sorted((M2/"src/esn_eq").glob("*.py"))])
    fit_path = Path(__file__).with_name("runs")/"m2_vix_curve_fit_001"
    fit = json.loads((fit_path/"run.json").read_text())
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    architecture,parameters,premium,hashes = frozen_stage3()
    expected = json.loads(json.dumps(dict(architecture=asdict(architecture),parameters=asdict(parameters),
        premium=asdict(premium),target_sha256=targets.artifact_sha256)))
    if any(fit[k] != value for k,value in expected.items()):
        raise ValueError("Model or target snapshot changed")
    dated = DatedM2Pricer(architecture,premium)
    initial,_ = dated.initial_state([0.,0.,0.])
    family = NuisanceCurveFamily.from_targets(targets)
    template = family.build([.04,.05,fit["selected_levels"][2]])
    low = ConditionalVIXPricer(dated,parameters,template,initial,normalizer_order=48)
    high = ConditionalVIXPricer(dated,parameters,template,initial,normalizer_order=64)
    quotes = targets.select("VIX",expiry_days=[153]).reset_index(drop=True)
    expiry,maturity = int(quotes.expiry.iloc[0]),float(quotes.T_years.iloc[0])
    windows = [p.gaussian_window(maturity,time_order=8) for p in (low,high)]
    np.testing.assert_array_equal(windows[0].times,windows[1].times)
    sampler = GaussianOuterSampler(dated,maturity,initial)
    state_groups = []
    for role in (0,1):
        seeds = [row["seed"] for row in fit["seed_records"] if row["role"] == role and row["expiry"] == expiry]
        state_groups.append(np.concatenate([sampler.draw(fit["arguments"]["power"],seed) for seed in seeds]))
    states = np.concatenate(state_groups)
    arrays = [np.zeros((len(states),len(family.ends))) for _ in windows]
    segments = np.minimum(np.searchsorted(family.ends,windows[0].times,side="right"),len(family.ends)-1)
    for node,(time,weight,segment) in enumerate(zip(windows[0].times,windows[0].weights,segments)):
        mean,cov = dated.gaussian_projection_law(float(time-maturity),states)
        expected_g = dated.gaussian_expectation(mean,cov,parameters,48)
        for coeff,window in zip(arrays,windows):
            coeff[:,segment] += weight*expected_g/window.normalizer[node]
        if (node+1) % 20 == 0:
            print(f"quadrature nodes {node+1}/{len(segments)}",flush=True)
    archive = np.load(fit_path/"exposures.npz")
    rows = []
    count = len(state_groups[0])
    for role_id,role in enumerate(("train","check")):
        coeffs = [archive[f"role{role_id}_expiry{expiry}"],
                  *[c[role_id*count:(role_id+1)*count] for c in arrays]]
        for rule,coeff,window in zip(("GH32_M48","GH48_M48","GH48_M64"),coeffs,
                                    (windows[0],windows[0],windows[1])):
            for candidate in ("starting","selected"):
                levels = np.array(fit[f"{candidate}_levels"])
                values = 100*np.sqrt((low.variance_floor*window.weights.sum()+coeff @ (levels-low.variance_floor))/window.delta)
                sign = np.where(quotes.option_type.to_numpy() == "call",1.,-1.)
                prices = np.maximum(sign*(values[:,None]-quotes.strike.to_numpy()),0.).mean(axis=0)
                residual = interval_residuals(prices,quotes.price_bid_U,quotes.price_ask_U,quotes.price_scale_U)
                forward = float(values.mean())
                forward_residual = (forward-float(quotes.forward.iloc[0]))/.5
                rows.append(dict(role=role,rule=rule,candidate=candidate,
                    weighted_option_loss=float(np.mean(residual**2)/4),
                    weighted_forward_loss=forward_residual**2/4,forward=forward))
    table = pd.DataFrame(rows)
    changes = []
    for (role,rule),group in table.groupby(["role","rule"]):
        selected = group[group.candidate == "selected"].iloc[0]
        starting = group[group.candidate == "starting"].iloc[0]
        changes.append(dict(role=role,rule=rule,selected_minus_starting_block_loss=float(
            selected.weighted_option_loss+selected.weighted_forward_loss-
            starting.weighted_option_loss-starting.weighted_forward_loss)))
    check_sources(sources)
    args.output.mkdir(parents=True,exist_ok=False)
    save_sources(args.output,sources,ROOT)
    record = dict(status="one-block fixed-state quadrature attribution; no refitting or full-model convergence claim",
        expiry_days=153,role_states=count,fit_sha256=hashlib.sha256((fit_path/"run.json").read_bytes()).hexdigest(),
        cache_sha256=hashlib.sha256((fit_path/"exposures.npz").read_bytes()).hexdigest(),
        target_sha256=targets.artifact_sha256,implementation_sha256=source_hashes(sources,ROOT),changes=changes)
    (args.output/"run.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    table.to_csv(args.output/"block_losses.csv",index=False)
    print(json.dumps(changes,indent=2),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists")
    run(args)
