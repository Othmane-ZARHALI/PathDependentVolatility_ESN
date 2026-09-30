"""Compare historical and simulated roughness on identical lag bands."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from compare import ROOT, load_legacy

sys.path.insert(0, str(ROOT))
from matteo.preprocessing import hurst, market_data  # noqa: E402


HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simulation", type=Path,
                        default=HERE / "nine_assets_fresh_001.json")
    parser.add_argument("--output", type=Path,
                        default=HERE / "roughness_lag_sensitivity_001.json")
    args = parser.parse_args()
    simulation = json.loads(args.simulation.read_text())
    legacy, _ = load_legacy()
    frames = market_data.load_bench9()
    result = {}
    for asset, frame in frames.items():
        log_vol = 0.5 * np.log(market_data.garman_klass_variance(frame))
        real_1_40 = hurst.estimate_hurst(log_vol, lags=range(1, 41))
        real_2_40 = hurst.estimate_hurst(log_vol, lags=range(2, 41))
        if not np.isclose(real_2_40, simulation["assets"][asset]["target"]["Hhat"],
                          rtol=0, atol=1e-12):
            raise AssertionError(f"{asset}: rebuilt real H differs from saved target")
        models = {}
        for name, record in simulation["assets"][asset]["models"].items():
            model_1_40 = record["statistics"]["Hhat"]
            model_2_40 = record["Hhat_common_lags_2_40"]
            non_h_score = (record["weighted_squared_error_total"] -
                           record["weighted_squared_error_by_block"]["Hhat"])
            tolerance_h = legacy.TOLERANCES["Hhat"]
            models[name] = {
                "H_lags_1_40": model_1_40,
                "H_lags_2_40": model_2_40,
                "absolute_error_common_1_40": abs(model_1_40 - real_1_40),
                "absolute_error_common_2_40": abs(model_2_40 - real_2_40),
                "full_sse_common_1_40": non_h_score +
                    (3 * (model_1_40 - real_1_40) / tolerance_h) ** 2,
                "full_sse_common_2_40": non_h_score +
                    (3 * (model_2_40 - real_2_40) / tolerance_h) ** 2,
            }
        result[asset] = {
            "historical_H_lags_1_40": real_1_40,
            "historical_H_lags_2_40": real_2_40,
            "models": models,
        }
    output = {"source": "same historical GK variance and fresh simulated paths",
              "assets": result}
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
