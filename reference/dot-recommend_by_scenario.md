# Best estimator under each collection scenario (internal)

Maps each instrument/NC collection scenario to the estimators it makes
available, then picks the most robust eligible one from the
`iconic_recommend` ranking (robustness = per-estimand max\|bias\| over
the Phase 3 rho x omega surface). An estimator appears only under
scenarios that supply its required data:

- UNADJ: none (always available)

- COCA: W (negative controls only; no instrument)

- DIRECT, IV2SLS, PGC: G + W

- IV2SLS2: G + Gm (optional path-specific W1/W2 augmentation)

- PGC2: G + W1 + W2

- PGC2Gm: G + Gm + W1 + W2

## Usage

``` r
.recommend_by_scenario(rec)
```

## Arguments

- rec:

  An `iconic_recommendation` object.

## Value

A data frame with columns `scenario`, `estimator`, and `robustness_NDE`.
