"""Measure the effect of Stage-2's omitted fitted-alpha argument.

Uses archived Stage-1-like coordinates from the saved final fit as a fixed
comparison point. The saved final m1 has uncertain provenance; this only
isolates the alpha argument under common simulation seeds.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import joint_calibration_reviewed as J
from compare import LEGACY, ROOT


HERE = Path(__file__).resolve().parent


def kurtosis_at(fit, m1, alpha, n_paths, days, seed_offset):
    params = dict(J.P_FIXED, H=fit["H_fit"], rough_scale=fit["rough_scale_fit"],
                  lam_lo=fit["lam_lo_fit"], lam_hi=fit["lam_hi_fit"], m1=m1)
    matrices = J.E.build_esn_matrices(J.OPT_ARCH, params)
    values = []
    for q in range(n_paths):
        returns, _, _ = J.E._sim_esn_with_params(
            16_600_000 + seed_offset + q, days, J.OPT_ARCH, matrices,
            alpha=alpha,
        )
        values.append(J.excess_kurtosis(returns))
    return float(np.mean(values))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--asset", default="SP500", choices=J.NAMES9)
    parser.add_argument("--paths", type=int, default=12)
    parser.add_argument("--days", type=int, default=2500)
    parser.add_argument("--output", type=Path, default=HERE / "stage2_alpha_sp500_001.json")
    args = parser.parse_args()
    fit = json.loads((LEGACY / "joint_calibration_final.json").read_text())[args.asset]
    target = json.loads((ROOT / "data" / "processed" /
                         "bench9_real_kurtosis.json").read_text())[args.asset]
    seed_offset = J.NAMES9.index(args.asset) * 1000
    inputs = (fit["H_fit"], fit["rough_scale_fit"], fit["lam_lo_fit"],
              fit["lam_hi_fit"], target)
    m1_fitted, k_fitted = J.bisect_m1_for_kurtosis(
        *inputs, n_paths=args.paths, T_eval=args.days,
        seed_offset=seed_offset, alpha=fit["alpha_fit"])
    m1_default, k_default = J.bisect_m1_for_kurtosis(
        *inputs, n_paths=args.paths, T_eval=args.days,
        seed_offset=seed_offset, alpha=1.0)
    k_default_at_fitted = kurtosis_at(
        fit, m1_default, fit["alpha_fit"], args.paths, args.days, seed_offset)
    result = {
        "status": "controlled_stage2_alpha_sensitivity_not_full_refit",
        "asset": args.asset, "paths": args.paths, "days": args.days,
        "seed_first": 16_600_000 + seed_offset,
        "target_excess_kurtosis": target,
        "archived_fitted_alpha": fit["alpha_fit"],
        "archived_final_m1": fit["m1_fit"],
        "m1_with_fitted_alpha": m1_fitted,
        "kurtosis_with_fitted_alpha": k_fitted,
        "m1_with_omitted_alpha_default_1": m1_default,
        "kurtosis_at_default_alpha_1": k_default,
        "kurtosis_of_default_m1_at_fitted_alpha": k_default_at_fitted,
        "m1_difference_default_minus_fitted": m1_default - m1_fitted,
        "kurtosis_error_default_m1_at_fitted_alpha": k_default_at_fitted - target,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
