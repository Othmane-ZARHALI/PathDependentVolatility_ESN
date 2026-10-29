"""Refine Gaussian time integration against archived 128/day fixed-state VIX."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform
import time

import numpy as np

from .m2_spx import ROOT, M2, frozen_stage3
from .provenance import capture_sources, source_hashes, check_sources, save_sources
from esn_eq.market_targets import MarketTargets
from esn_eq.pricing import DatedM2Pricer
from esn_eq.vix_pricing import ConditionalVIXPricer


def run(args):
    files = [Path(__file__), Path(__file__).with_name("m2_spx.py"),
             *sorted((M2/"src/esn_eq").glob("*.py"))]
    sources = capture_sources(files)
    targets = MarketTargets.load(ROOT/"data/processed/options/2026-02-19")
    architecture, parameters, premium, hashes = frozen_stage3()
    dated = DatedM2Pricer(architecture, premium)
    initial, _ = dated.initial_state([0., 0., 0.])
    curve = targets.curve("training_mid_piecewise_v1")
    engine = ConditionalVIXPricer(dated, parameters, curve, initial, normalizer_order=48)
    projections = [[0., 0., 0.], [0., -1., -1.], [1., -2., -2.]]
    states = np.array([dated.initial_state(u)[0] for u in projections])
    reference_path = Path(__file__).with_name("runs")/"m2_vix_64_128_states_001/run.json"
    reference = json.loads(reference_path.read_text())
    if (reference["architecture"] != asdict(architecture)
            or reference["parameters"] != asdict(parameters)
            or reference["premium"] != asdict(premium)
            or reference["target_sha256"] != targets.artifact_sha256):
        # JSON converts tuples to lists, so compare canonical JSON representations.
        current = json.loads(json.dumps(dict(architecture=asdict(architecture),
            parameters=asdict(parameters),premium=asdict(premium),target_sha256=targets.artifact_sha256)))
        if any(reference[k] != value for k, value in current.items()):
            raise ValueError("Archived fixed-model reference does not match")
    fine = reference["vix"]["levels"][-1]
    if fine["steps_per_day"] != 128 or fine["state_projections_C_Y_F"] != projections:
        raise ValueError("Unexpected archived reference states or grid")
    baseline = np.array(fine["quadrature_vix"])
    levels = []
    for order in (4, 8, 16):
        start = time.perf_counter()
        window = engine.gaussian_window(55/365, time_order=order)
        values = engine.evaluate(states, window, order=32).vix
        levels.append(dict(time_order=order, nodes=len(window.times), vix=values.tolist(),
            minus_midpoint_128=(values-baseline).tolist(), wall_seconds=time.perf_counter()-start))
        print(json.dumps(levels[-1]), flush=True)
    record = dict(status="fixed-state time-rule diagnostic; no outer accuracy or calibration claim",
        architecture=asdict(architecture),parameters=asdict(parameters),premium=asdict(premium),
        target_sha256=targets.artifact_sha256,frozen_model_sha256=hashes,
        implementation_sha256=source_hashes(sources, ROOT),
        reference=str(reference_path.relative_to(ROOT)),
        reference_sha256=__import__("hashlib").sha256(reference_path.read_bytes()).hexdigest(),
        reference_nodes=fine["nodes"],reference_vix=baseline.tolist(),
        state_projections_C_Y_F=projections,normalizer_order=48,conditional_order=32,
        expiry_days=55,delta=30/365,curve_id="training_mid_piecewise_v1",levels=levels,
        input_compatibility=targets.audit["variance_audit"],
        versions=dict(python=platform.python_version(),numpy=np.__version__))
    check_sources(sources)
    args.output.mkdir(parents=True, exist_ok=False)
    save_sources(args.output, sources, ROOT)
    (args.output/"run.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists")
    run(args)
