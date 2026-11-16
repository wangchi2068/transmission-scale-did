# M0-OI variance: G/(G-1) factor and the Windmeijer correction

Reviewer M6 noted that Section 2.5 of the manuscript reports the `G/(G-1)` factor and cites
Windmeijer (2005) for two-step estimators, whereas `qdgmm.fit_two_step` uses neither. This directory
records the corrected variance and its effect on the reported application numbers.

## What was changed

`qdgmm.fit_two_step` gains an optional `variance` argument:

- `variance="plain"` reproduces the previous implementation, `V = (Gamma' W Gamma)^{-1} / G`.
- `variance="windmeijer"` (default) applies the first-order correction of Windmeijer (2005),

  `V = (D'WD)^{-1} (D'W Omega W D + D' W (dGamma) (I - W Gamma (D'WD)^{-1} Gamma' W) (dGamma)' W D) (D'WD)^{-1} / G`,

  with the second-derivative term `dGamma` computed by central differences in the beta block.
- The `G/(G-1)` small-sample factor enters through `qdgmm.small_sample_factor(G)`.

The J statistic is unchanged: the correction affects the variance of the estimator, not the
over-identification test.

## Run

```
cd experiment/code && python -I ../../verification/application/m0oi_variance.py
```

The script reports, for the application panel, the plain and corrected standard errors of
`tau`, `kappa` and `phi`, the correction ratio, and the J statistic, and writes them to
`m0oi_variance.json`.
