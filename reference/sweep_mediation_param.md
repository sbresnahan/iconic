# Sweep a single mediation simulation parameter across a grid

Sweep a single mediation simulation parameter across a grid

## Usage

``` r
sweep_mediation_param(
  param,
  param_grid,
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
  u_strength = NULL,
  w_coverage_profile = NULL,
  base_seed = 0,
  n_cores = 1
)
```

## Arguments

- param:

  Parameter to vary: one of "beta_X", "conf_str", "w_signal", "alpha_M",
  "beta_M", "n_samples", "mo_confounding", "phi", "rho_G1", "rho_G2",
  "rho_pop", "omega_1", "omega_2", "feat_cor".

- param_grid:

  Numeric vector of values to sweep.

- n_iter:

  Replicates per grid point. Default 100.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Features per replicate. Default 20.

- n_mediators:

  Number of independent mediators. Default 1.

- beta_X:

  Baseline direct effect. Default 0.10.

- alpha_M:

  Baseline mediator path. Default 0.50.

- beta_M:

  Baseline mediator effect. Default 0.30.

- conf_str:

  Baseline confounding strength. Default 0.80.

- w_signal:

  Baseline proxy quality. Default 0.70.

- mo_confounding:

  Baseline M-O confounding. Default 0.80.

- phi:

  Baseline mediator-instrument strength. 0 = no mediator instrument.
  Default 0.

- rho_G1:

  Baseline G1-conf_XM correlation. Default 0.

- rho_G2:

  Baseline G2-conf_MY correlation. Default 0.

- rho_pop:

  Baseline population structure. Default 0.

- lambda_XM:

  Optional per-path confounder loading vector (X-\>M path).

- lambda_MY:

  Optional per-path confounder loading vector (M-\>Y path).

- omega_1:

  Baseline W1 coverage. NULL = use w_signal.

- omega_2:

  Baseline W2 coverage. NULL = use w_signal.

- feat_cor:

  Baseline within-module feature correlation. Default 0.

- u_strength:

  Numeric vector: per-confounder strength scaling. `NULL` = uniform
  (legacy behavior). See
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md).

- w_coverage_profile:

  A list with `w1`/`w2` per-control coverage vectors. `NULL` = uniform
  coverage. See
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md).

- base_seed:

  Seed offset. Default 0.

- n_cores:

  Parallel workers. Default 1.

## Value

A list with `summary` (data frame) and `iter_bias` (per-iteration
NDE/NIE bias).

## Examples

``` r
res <- sweep_mediation_param("conf_str", c(0.2, 0.8), n_iter = 3,
                             n_samples = 100)
res$summary
#>       param param_value true_NDE true_NIE method     NDE_mean   NDE_bias
#> 1  conf_str         0.2      0.1     0.15  UNADJ -0.052175951 -0.1521760
#> 2  conf_str         0.2      0.1     0.15 DIRECT -0.050195520 -0.1501955
#> 3  conf_str         0.2      0.1     0.15   COCA -0.161803415 -0.2618034
#> 4  conf_str         0.2      0.1     0.15 IV2SLS -0.006311205 -0.1063112
#> 5  conf_str         0.2      0.1     0.15    PGC -0.051572633 -0.1515726
#> 6  conf_str         0.8      0.1     0.15  UNADJ -0.478802481 -0.5788025
#> 7  conf_str         0.8      0.1     0.15 DIRECT -0.425246117 -0.5252461
#> 8  conf_str         0.8      0.1     0.15   COCA  2.077246902  1.9772469
#> 9  conf_str         0.8      0.1     0.15 IV2SLS -0.163416385 -0.2634164
#> 10 conf_str         0.8      0.1     0.15    PGC -0.464073478 -0.5640735
#>    NDE_pct_bias     NDE_sd  NDE_rmse NDE_mean_se NDE_coverage   NIE_mean
#> 1     -1.521760 0.04375792 0.1582415  0.03656377   0.05000000  0.3548806
#> 2     -1.501955 0.08898512 0.1741984  0.08530117   0.58333333  0.3127152
#> 3     -2.618034 5.34036888 5.2742092  5.18411700   0.89189189  0.4025030
#> 4     -1.063112 0.04214475 0.1142307  0.04339490   0.28333333  0.3179926
#> 5     -1.515726 0.04720124 0.1586351  0.03922785   0.05000000  0.3020370
#> 6     -5.788025 0.13297456 0.5936327  0.06968797   0.00000000  1.1751432
#> 7     -5.252461 0.16258928 0.5494344  0.13405549   0.03333333  0.9492700
#> 8     19.772469 4.71939996 5.0682632  5.21090364   0.93333333 -1.8130655
#> 9     -2.634164 0.09403182 0.2794330  0.06466726   0.08333333  0.6960097
#> 10    -5.640735 0.14306818 0.5816410  0.08164359   0.00000000  0.8572318
#>      NIE_bias NIE_pct_bias     NIE_sd  NIE_rmse NIE_mean_se NIE_coverage
#> 1   0.2048806     1.365871 0.04068571 0.2088153  0.03913716   0.00000000
#> 2   0.1627152     1.084768 0.08882045 0.1850238  0.08004840   0.43333333
#> 3   0.2525030     1.683354 5.34428361 5.2776127  5.18165205   0.89189189
#> 4   0.1679926     1.119951 0.03776188 0.1721154  0.04933904   0.00000000
#> 5   0.1520370     1.013580 0.04543845 0.1585733  0.03897774   0.03333333
#> 6   1.0251432     6.834288 0.20035075 1.0442174  0.07718812   0.00000000
#> 7   0.7992700     5.328467 0.19274012 0.8218042  0.13245578   0.00000000
#> 8  -1.9630655   -13.087103 4.71682131 5.0603973  5.19160868   0.93333333
#> 9   0.5460097     3.640064 0.10536109 0.5559159  0.06928902   0.00000000
#> 10  0.7072318     4.714879 0.15933812 0.7246671  0.08396233   0.00000000
#>     NIE_type1  NDE_type1 n_NDE n_NIE
#> 1  1.00000000 0.33333333    60    60
#> 2  0.95000000 0.13333333    60    60
#> 3  0.10810811 0.16216216    37    37
#> 4  1.00000000 0.05000000    60    60
#> 5  1.00000000 0.33333333    60    60
#> 6  1.00000000 1.00000000    60    60
#> 7  1.00000000 0.85000000    60    60
#> 8  0.06666667 0.08888889    45    45
#> 9  1.00000000 0.61666667    60    60
#> 10 1.00000000 0.98333333    60    60
```
