# Benchmark estimators across confounding scenarios on synthetic data

For each scenario in the grid, generates `n_iter` synthetic datasets
with
[`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
(using the trained generator for texture), runs every estimator, and
summarises bias / RMSE / power. This is the multi-confounder,
negative-control-aware generalisation of the package's parameter sweeps.

## Usage

``` r
gan_sensitivity(
  trained_gan = NULL,
  conf_grid = c(0.2, 0.5, 0.8),
  coverage_grid = c(0.3, 0.7, 1),
  k_grid = 1,
  nc_model = "proxy",
  n_iter = 50,
  n_samples = 500,
  n_features = 20,
  beta_X = 0.1,
  alpha_M = 0.5,
  beta_M = 0.3,
  effect_size = NULL,
  base_seed = 700,
  n_cores = 1,
  outcome_type = c("continuous", "survival"),
  effect_scale = c("loghr", "rmst"),
  surv_h0 = 0.1,
  surv_event_frac = 0.6,
  surv_censor_rate = NULL
)
```

## Arguments

- trained_gan:

  An `iconic_gan` (or `NULL` to use default texture).

- conf_grid:

  Confounding-strength values to sweep. Default `c(0.2, 0.5, 0.8)`.

- coverage_grid:

  Negative-control coverage values in `[0,1]`. Default `c(0.3, 0.7, 1)`.

- k_grid:

  Numbers of latent confounders to sweep. Default `1`.

- nc_model:

  Negative-control model (function or name). Default `"proxy"`.

- n_iter:

  Replicates per scenario. Default 50.

- n_samples:

  Samples per replicate. Default 500.

- n_features:

  Features per replicate. Default 20.

- beta_X, alpha_M, beta_M:

  Causal paths (ground truth). Defaults 0.10 / 0.50 / 0.30.

- effect_size:

  Optional pure-direct total effect override (see
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)).

- base_seed:

  Base RNG seed. Default 700.

- n_cores:

  Parallel workers across replicates. Default 1.

- outcome_type:

  `"continuous"` (default) or `"survival"` When survival, the DGP
  generates time-to-event outcomes and estimation uses the Cox / RMST
  survival drivers via
  [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).

- effect_scale:

  `"loghr"` (default) or `"rmst"`. Only used when
  `outcome_type = "survival"`.

- surv_h0:

  Baseline hazard for the survival DGP. Default 0.1.

- surv_event_frac:

  Target fraction of observed events. Default 0.6.

- surv_censor_rate:

  Explicit censoring rate. Default NULL.

## Value

A list with `summary` (one row per scenario x method, with
`conf_strength`, `coverage`, `k`, `true_total` and the columns from
`summarise_results`) and `grid` (the scenario grid).

## Examples

``` r
sens <- gan_sensitivity(NULL, conf_grid = 0.8, coverage_grid = 0.7,
  n_iter = 2, n_samples = 100, n_features = 5)
head(sens$summary)
#>   conf_strength coverage k true_total method       mean     median         sd
#> 1           0.8      0.7 1       0.25  UNADJ  0.6288594  0.6377314 0.07535704
#> 2           0.8      0.7 1       0.25 DIRECT  0.3671867  0.3533574 0.11724910
#> 3           0.8      0.7 1       0.25   COCA -0.6149253 -0.5515430 0.22603691
#> 4           0.8      0.7 1       0.25 IV2SLS  0.2203744  0.2251649 0.09778570
#> 5           0.8      0.7 1       0.25    PGC  0.2784006  0.2934452 0.09297978
#>          bias   abs_bias       rmse power  n
#> 1  0.37885941 0.37885941 0.38554541   1.0 10
#> 2  0.11718666 0.11718666 0.16157144   0.9 10
#> 3 -0.86492528 0.86492528 0.89111120   1.0 10
#> 4 -0.02962563 0.02962563 0.09738335   0.5 10
#> 5  0.02840064 0.02840064 0.09266774   1.0 10
```
