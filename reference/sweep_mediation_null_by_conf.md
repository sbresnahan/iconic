# Sweep mediation Type I error across confounding strength levels

Sweep mediation Type I error across confounding strength levels

## Usage

``` r
sweep_mediation_null_by_conf(
  conf_grid = c(0.2, 0.4, 0.6, 0.8, 1),
  n_iter = 100,
  n_samples = 500,
  n_features = 20,
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
  base_seed = 900,
  n_cores = 1,
  alpha = 0.05
)
```

## Arguments

- conf_grid:

  Numeric vector of confounding strength values. Default c(0.2, 0.4,
  0.6, 0.8, 1.0).

- n_iter:

  Replicates per conf_str value. Default 100.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Features per replicate. Default 20.

- w_signal:

  Proxy quality omega. Default 0.70.

- mo_confounding:

  Strength of U1 -\> M. Default 0.80.

- phi:

  Strength of the mediator instrument Gm -\> M. 0 = no mediator
  instrument. Default 0.

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

  Within-module feature correlation. Default 0.

- base_seed:

  Seed offset. Default 900.

- n_cores:

  Parallel workers. Default 1.

- alpha:

  Significance threshold. Default 0.05.

## Value

A data frame with columns: conf_str, method, NIE_type1, NDE_type1.

## Examples

``` r
t1e <- sweep_mediation_null_by_conf(c(0.2, 0.8), n_iter = 3,
                                    n_samples = 100)
head(t1e)
#>   conf_str method NIE_type1  NDE_type1
#> 1      0.2  UNADJ 0.7666667 0.01666667
#> 2      0.2 DIRECT 0.0000000 0.08333333
#> 3      0.2   COCA 0.0000000 0.00000000
#> 4      0.2 IV2SLS 0.1833333 0.05000000
#> 5      0.2    PGC 0.0000000 0.06666667
#> 6      0.8  UNADJ 1.0000000 0.05000000
```
