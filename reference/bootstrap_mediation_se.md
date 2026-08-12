# Bootstrap SE for a mediation estimator

Resamples the data `n_boot` times and returns the SD of the bootstrap
NDE/NIE distribution. This is an opt-in alternative to the delta-method
SE for users who want a non-parametric SE ( bootstrapped reasonably").
Slower than the delta method.

## Usage

``` r
bootstrap_mediation_se(estimator_fn, n, n_boot = 500)
```

## Arguments

- estimator_fn:

  A closure `function(idx)` returning a fit list with `NDE` and `NIE`
  (and optionally `NDE_se`, `NIE_se`).

- n:

  Sample size (length of the resampling index).

- n_boot:

  Number of bootstrap resamples. Default 500.

## Value

A list with NDE_boot_se, NIE_boot_se, NDE_boot_dist, NIE_boot_dist.

## Details

The estimator is supplied as a closure `estimator_fn(idx)` that captures
all model-specific data (outcome, exposure, mediator, instruments,
negative controls, covariates) and subsets each by `idx` internally.
This guarantees every resampled draw uses a synchronised bootstrap
sample across all variables — critical for the instrumented estimators
(IV2SLS, PGC, PGC2Gm) whose instruments and NC panels must be resampled
in lockstep with y, X, and M.
