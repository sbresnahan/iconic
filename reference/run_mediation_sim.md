# Run repeated mediation simulations for a single parameter configuration

Generates `n_iter` synthetic datasets with mediator-outcome confounding
(via `mo_confounding`), runs all five (or six, when `phi > 0`) mediation
estimators, and summarises bias / RMSE / Type I error for both NDE and
NIE.

## Usage

``` r
run_mediation_sim(
  n_iter = 100,
  n_samples = 500,
  n_features = 20,
  n_mediators = 1,
  beta_X = 0.1,
  alpha_M = 0.5,
  beta_M = 0.3,
  conf_str = 0.8,
  w_signal = 0.7,
  mo_confounding = 0.8,
  phi = 0,
  rho_G1 = 0,
  rho_G2 = 0,
  rho_pop = 0,
  lambda_XM = NULL,
  lambda_MY = NULL,
  omega_1 = NULL,
  omega_2 = NULL,
  feat_cor = 0,
  base_seed = 100,
  n_cores = 1
)
```

## Arguments

- n_iter:

  Number of simulation replicates. Default 100.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Number of outcome and negative-control features. Default 20.

- n_mediators:

  Number of independent mediators. When \> 1, each mediator has its own
  genetic instrument Gm and contributes additively to Y. Default 1
  (single mediator, backward compatible).

- beta_X:

  Direct effect of X on Y (true NDE). Default 0.10.

- alpha_M:

  Effect of X on mediator. Default 0.50.

- beta_M:

  Effect of mediator on Y (per-mediator true NIE = alpha_M \* beta_M).
  Default 0.30.

- conf_str:

  Confounding strength delta. Default 0.80.

- w_signal:

  Proxy quality omega. Default 0.70.

- mo_confounding:

  Strength of U1 -\> M (mediator-outcome confounding). Default 0.80.

- phi:

  Strength of the mediator instrument Gm -\> M. 0 = no mediator
  instrument (five estimators, backward compatible). \> 0 = generates Gm
  and includes the 2-stage MR estimator (IV2SLS2). Default 0.

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

  Coverage of conf_XM by W1. NULL = use w_signal.

- omega_2:

  Coverage of conf_MY by W2. NULL = use w_signal.

- feat_cor:

  Within-module correlation for block-diagonal co-expression modules in
  Y and W. 0 = independent features. Default 0.

- base_seed:

  Starting seed; replicate i uses base_seed + i. Default 100.

- n_cores:

  Number of parallel workers. Default 1.

## Value

A list with `summary` (one row per method with NDE/NIE bias, RMSE, Type
I), `raw` (full per-feature results), `true_NDE`, `true_NIE`, and
`params`.

## Details

When are non-default (rho_G1, rho_G2, rho_pop, lambda_XM, lambda_MY,
omega_1, omega_2), the two-stage proximal estimators (PGC2, PGC2Gm) are
also included.

## Examples

``` r
res <- run_mediation_sim(n_iter = 3, n_samples = 100,
                         mo_confounding = 0.8)
res$summary
#>   method   NDE_mean   NDE_bias NDE_pct_bias     NDE_sd  NDE_rmse NDE_mean_se
#> 1  UNADJ -0.4692764 -0.5692764    -5.692764 0.13046612 0.5837922  0.06182445
#> 2 DIRECT -0.4026100 -0.5026100    -5.026100 0.14829814 0.5236818  0.11644642
#> 3   COCA  2.3579413  2.2579413    22.579413 4.42753239 4.9269858  5.55311186
#> 4 IV2SLS -0.1893373 -0.2893373    -2.893373 0.07783502 0.2994552  0.05814444
#> 5    PGC -0.4413775 -0.5413775    -5.413775 0.13261134 0.5571197  0.07335239
#>   NDE_coverage   NIE_mean   NIE_bias NIE_pct_bias    NIE_sd  NIE_rmse
#> 1   0.00000000  1.1309836  0.9809836     6.539891 0.1913857 0.9991730
#> 2   0.01666667  0.8873532  0.7373532     4.915688 0.1684348 0.7560339
#> 3   0.97826087 -2.1607016 -2.3107016   -15.404678 4.4217077 4.9462924
#> 4   0.01666667  0.6710078  0.5210078     3.473386 0.1071212 0.5317263
#> 5   0.00000000  0.7834623  0.6334623     4.223082 0.1400766 0.6485128
#>   NIE_mean_se NIE_coverage  NIE_type1  NDE_type1 n_NDE n_NIE
#> 1  0.07189956    0.0000000 1.00000000 1.00000000    60    60
#> 2  0.11620532    0.0000000 1.00000000 0.86666667    60    60
#> 3  5.53584617    0.9782609 0.02173913 0.06521739    46    46
#> 4  0.06855147    0.0000000 1.00000000 0.81666667    60    60
#> 5  0.07627316    0.0000000 1.00000000 0.98333333    60    60
```
