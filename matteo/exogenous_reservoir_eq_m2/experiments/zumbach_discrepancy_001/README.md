# Pearson/rank Strong-Zumbach discrepancy study

This diagnostic experiment rebuilds the five untouched Stage-4 validation seeds for the frozen
primary and the preceding M2 anchor. No parameter is fitted or selected.

At 5-, 10-, and 20-day windows, it separately records

- the forward correlation between past squared returns and future realised variance;
- the backward correlation between past realised variance and future squared returns;
- their difference under Pearson, Spearman rank, Gaussian-rank, and Kendall dependence;
- Pearson results after upper-tail winsorisation or removal;
- upper-tail covariance concentration and conditional decile profiles;
- both log-return and matched Ito-relative-return definitions.

Run from the package root:

```bash
python experiments/zumbach_discrepancy_001/run.py
```
