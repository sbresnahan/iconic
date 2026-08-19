# Benchmark mediation estimators across confounding scenarios

For each scenario in the grid, generates `n_iter` synthetic datasets
with mediator-outcome confounding (via `mo_confounding`), runs every
mediation estimator, and summarises NDE/NIE bias / RMSE / Type I error.
This is the mediation analogue of
[`gan_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_sensitivity.md).

## Usage

``` r
gan_mediation_sensitivity(
  trained_gan = NULL,
  conf_grid = c(0.2, 0.5, 0.8),
  coverage_grid = c(0.3, 0.7, 1),
  k_grid = 1,
  mo_confounding = 0.8,
  phi = 0,
  rho_G1 = 0,
  rho_G2 = 0,
  rho_pop = 0,
  lambda_XM = NULL,
  lambda_MY = NULL,
  omega_1 = NULL,
  omega_2 = NULL,
  nc_model = "proxy",
  n_iter = 50,
  n_samples = 500,
  n_features = 20,
  beta_X = 0.1,
  alpha_M = 0.5,
  beta_M = 0.3,
  base_seed = 750,
  n_cores = 1,
  outcome_type = c("continuous", "survival", "binary"),
  effect_scale = c("loghr", "rmst", "logor", "riskdiff"),
  surv_h0 = 0.1,
  surv_event_frac = 0.6,
  surv_censor_rate = NULL,
  bin_prev = 0.5
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

- mo_confounding:

  Strength of U1 -\> M (mediator-outcome confounding). Default 0.80.

- phi:

  Strength of the mediator instrument Gm -\> M. 0 = no mediator
  instrument (five estimators). \> 0 = generates Gm and includes the
  2-stage MR estimator (IV2SLS2). Default 0.

- rho_G1:

  Correlation of G1 with conf_XM. Default 0.

- rho_G2:

  Correlation of G2 with conf_MY. Default 0.

- rho_pop:

  Shared population structure. Default 0.

- lambda_XM:

  Optional per-path confounder loading vector (X-\>M path).

- lambda_MY:

  Optional per-path confounder loading vector (M-\>Y path).

- omega_1:

  Coverage of conf_XM by W1. NULL = use `coverage`.

- omega_2:

  Coverage of conf_MY by W2. NULL = use `coverage`.

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

- base_seed:

  Base RNG seed. Default 750.

- n_cores:

  Parallel workers across replicates. Default 1.

- outcome_type:

  `"continuous"` (default), `"survival"`, or `"binary"`. When survival,
  the DGP generates time-to-event outcomes and estimation uses the Cox /
  RMST survival mediation drivers via
  [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).
  When binary, the DGP generates a 0/1 outcome and estimation uses the
  logistic / linear-probability-model binary mediation drivers via
  [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).

- effect_scale:

  `"loghr"` (default), `"rmst"`, `"logor"`, or `"riskdiff"`.
  `"loghr"`/`"rmst"` apply to survival outcomes; `"logor"`/`"riskdiff"`
  apply to binary outcomes (an inert default is remapped with a
  message).

- surv_h0:

  Baseline hazard for survival DGP. See
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md).

- surv_event_frac:

  Target event fraction for survival DGP.

- surv_censor_rate:

  Censoring rate for survival DGP.

- bin_prev:

  Target prevalence for the binary DGP. Default 0.5.

## Value

A list with `summary` (one row per scenario x method, with
`conf_strength`, `coverage`, `k`, `mo_confounding`, `phi`, `true_NDE`,
`true_NIE` and NDE/NIE bias/RMSE/Type I columns) and `grid`.

## Details

When `phi > 0`, a mediator-specific genetic instrument (Gm) is generated
and the 2-stage MR estimator (IV2SLS2) is included in the results,
enabling point identification of NDE/NIE under M-O confounding.

## Examples

``` r
sens <- gan_mediation_sensitivity(NULL, conf_grid = 0.8,
  coverage_grid = 0.7, mo_confounding = 0.8,
  n_iter = 2, n_samples = 100, n_features = 5)
head(sens$summary)
#>   conf_strength coverage k mo_confounding phi rho_G1 rho_G2 rho_pop true_NDE
#> 1           0.8      0.7 1            0.8   0      0      0       0      0.1
#> 2           0.8      0.7 1            0.8   0      0      0       0      0.1
#> 3           0.8      0.7 1            0.8   0      0      0       0      0.1
#> 4           0.8      0.7 1            0.8   0      0      0       0      0.1
#> 5           0.8      0.7 1            0.8   0      0      0       0      0.1
#>   true_NIE method     NDE_mean    NDE_bias NDE_pct_bias     NDE_sd  NDE_rmse
#> 1     0.15  UNADJ -0.042487429 -0.14248743   -1.4248743 0.05561129 0.1519409
#> 2     0.15 DIRECT  0.008833256 -0.09116674   -0.9116674 0.11552854 0.1425606
#> 3     0.15   COCA  0.114220263  0.01422026    0.1422026 0.13341666 0.1273665
#> 4     0.15 IV2SLS -0.038700536 -0.13870054   -1.3870054 0.05206206 0.1472320
#> 5     0.15    PGC -0.015543555 -0.11554356   -1.1554356 0.06698824 0.1318673
#>   NDE_mean_se NDE_coverage    NIE_mean   NIE_bias NIE_pct_bias     NIE_sd
#> 1   0.1217711            1  0.71981233  0.5698123     3.798749 0.08766402
#> 2   0.1476969            1  0.32691101  0.1769110     1.179407 0.05025958
#> 3   0.2826143            1 -0.09316795 -0.2431679    -1.621120 0.06778980
#> 4   0.1117838            1  0.41817424  0.2681742     1.787828 0.08425477
#> 5   0.1116669            1  0.30099204  0.1509920     1.006614 0.05776075
#>    NIE_rmse NIE_mean_se NIE_coverage NIE_type1 NDE_type1 n_NDE n_NIE lambda_XM
#> 1 0.5758494  0.11546967          0.0         1         0    10    10          
#> 2 0.1832237  0.11062555          0.9         1         0    10    10          
#> 3 0.2515285  0.17767340          0.8         0         0    10    10          
#> 4 0.2798328  0.09712824          0.3         1         0    10    10          
#> 5 0.1606277  0.09010352          0.7         1         0    10    10          
#>   lambda_MY
#> 1          
#> 2          
#> 3          
#> 4          
#> 5          
```
