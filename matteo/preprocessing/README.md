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
