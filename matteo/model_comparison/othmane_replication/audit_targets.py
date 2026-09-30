"""Rebuild Othmane's nine-asset target JSONs in memory and check exact equality."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from matteo.preprocessing import hurst, kurtosis, leverage, market_data, taylor_acf, zumbach


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    frames = market_data.load_bench9()
    absolute, squared = taylor_acf.compute_targets(frames)
    rebuilt = {
        "bench9_real_hurst": hurst.compute_targets(frames),
        "bench9_real_leverage_wide": leverage.compute_targets(frames),
        "bench9_real_abs_acf_40": absolute,
        "bench9_real_sq_acf_40": squared,
        "bench9_real_zumbach_profile": zumbach.compute_targets(frames),
        "bench9_real_kurtosis": kurtosis.compute_targets(frames),
    }
    checks = {}
    for stem, values in rebuilt.items():
        path = market_data.PROCESSED_DIR / (stem + ".json")
        saved = json.loads(path.read_text())
        checks[stem] = {"exact_match": saved == values, "saved_sha256": _hash(path)}
    result = {
        "raw_files": {
            "stock": {"path": str(market_data.STOCK_FILE.relative_to(ROOT)),
                      "sha256": _hash(market_data.STOCK_FILE)},
            "index": {"path": str(market_data.INDEX_FILE.relative_to(ROOT)),
                      "sha256": _hash(market_data.INDEX_FILE)},
        },
        "sample_rows": {name: len(frame) for name, frame in frames.items()},
        "checks": checks,
        "all_exact": all(check["exact_match"] for check in checks.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(args.output)
    print(f"All six target files match: {result['all_exact']}")
    if not result["all_exact"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
