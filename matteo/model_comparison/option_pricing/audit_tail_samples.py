"""Append a zero-payoff-sample audit without rewriting original VIX run records."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def audit(directory):
    directory = Path(directory)
    output = directory/"tail_diagnostics.json"
    if output.exists():
        raise FileExistsError(output)
    prices = pd.read_csv(directory/"prices.csv")
    values = pd.read_csv(directory/"outer_values.csv")
    records = []
    for days, group in prices.groupby("expiry_days"):
        draws = values.loc[values.expiry_days == days, "quadrature"].to_numpy()
        hits = np.array([np.sum(draws > strike) if kind == "call" else np.sum(draws < strike)
                         for strike, kind in zip(group.strike, group.option_type)])
        nominal = 2*group.outer_mc_se_U.to_numpy() <= group.price_scale_U.to_numpy()/4
        records.append(dict(expiry_days=int(days), targets=len(group),
            zero_payoff_sample_quotes=int(np.sum(hits == 0)),
            nominal_se_screen_count=int(nominal.sum()),
            nominal_passes_without_two_payoff_hits=int(np.sum(nominal & (hits < 2))),
            minimal_nonzero_hit_screen_count=int(np.sum(nominal & (hits >= 2))),
            hit_counts=hits.tolist()))
    record = dict(status="diagnostic screen only; zero hits cannot bound rare-payoff expectation",
        input_sha256={name: hashlib.sha256((directory/name).read_bytes()).hexdigest()
                      for name in ("run.json", "prices.csv", "outer_values.csv")},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note="Original records remain intact. Requiring two hits is not a confidence guarantee.",
        expiries=records)
    output.write_text(json.dumps(record, indent=2)+"\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    audit(parser.parse_args().directory)
