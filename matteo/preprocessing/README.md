# Nine-asset benchmark preprocessing

Place the historical OHLC pickle and index workbook in `data/raw/` at the
repository root. Install the local dependencies with:

```bash
python -m pip install -r matteo/preprocessing/requirements.txt
```

From the repository root, generate the six JSON targets with:

```bash
python -m matteo.preprocessing.leverage
python -m matteo.preprocessing.taylor_acf
python -m matteo.preprocessing.hurst
python -m matteo.preprocessing.zumbach
python -m matteo.preprocessing.kurtosis
```

The targets are written to `data/processed/`, which is ignored by Git. The
joint calibration scripts read them from there. The stock sample follows
Othmane's seed-2026 selection; the index sample starts in 1994.

# Option snapshot preprocessing

Run from the repository root with the dependencies in `requirements.txt`:

```sh
python -m matteo.preprocessing.options \
  --spx data/raw/SPX_data_19_02_2026.csv \
  --vix data/raw/VIX_data_19_02_2026.csv \
  --output data/processed/options/2026-02-19 \
  --quote-unit iv_decimal
```

The command preserves every source row in `quotes_all.csv`, converts the two
IV bounds separately to undiscounted Black-76 prices, and records excluded
quotes with reason codes. `quotes_calibration.csv` contains one OTM interval
per eligible strike and deterministic train, validation, and final-test roles.
The 119-day SPX maturity is a final maturity test. The original wings and ITM
quotes remain in `quotes_all.csv` and the audit.

`spx_variance_terms.csv` separates the all-quote Cboe-style diagnostic from
training-only terms. The calculation uses date-level maturities and all valid
observed strikes; it cannot reconstruct the official VIX strike-cutoff rule.
`forward_variance.csv` stores piecewise-constant forward variance fitted to
training-only SPX terms when all pilot terms are complete and training surfaces
are feasible. The curve is right-continuous at knots, uses the last fitted
slope out to the required pricing horizon, and needs tail sensitivity analysis.
Bid and ask curve variants are sensitivity calculations: simultaneous quote
extremes need not form an arbitrage-free surface. If a fit is unsupported or
fails a cross-market moment gate, the file also describes an unfitted
low-dimensional nuisance family for later profiling.

`cross_market_checks.csv` compares each curve's exact 30-day variance window
moment with finite-strike VIX option lower bounds at the quoted forward and
at ±0.5 VIX points. A failed necessary gate is reported; it does not certify a
model failure. Quote uncertainty, omitted wings, forward provenance, and
settlement basis remain separate from the exact piecewise integration error.
`audit.json` and `manifest.json` contain counts, parity gaps, static-arbitrage
checks, source hashes, conventions, and the reproducible command. The vendor's
IV forward and day-count inputs, VIX forward venue, and a 2026 latent state are
not established by these files.
