# Run null mediation simulations to estimate Type I error rates

Sets `alpha_M = 0` and `beta_M = 0` (no true NIE) while keeping `beta_X`
active (NDE still present). Reports Type I error for both NIE (should be
0 under the null) and NDE (should reflect power for the direct effect).

## Usage

``` r
run_null_mediation_sim(
  n_iter = 200,
  n_samples = 500,
  n_features = 20,
  beta_X = 0.1,
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
  base_seed = 300,
  n_cores = 1,
  alpha = 0.05
)
```

## Arguments

- n_iter:

  Number of replicates. Default 200.

- n_samples:

  Observations per replicate. Default 500.

- n_features:

  Features per replicate. Default 20.

- beta_X:

  Direct effect (NDE, still active under null NIE). Default 0.10.

- conf_str:

  Confounding strength delta. Default 0.80.

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

  Within-module correlation for block-diagonal co-expression modules in
  Y and W. 0 = independent features. Default 0.

- base_seed:

  Seed offset. Default 300.

- n_cores:

  Parallel workers. Default 1.

- alpha:

  Significance threshold. Default 0.05.

## Value

A list with `rates` (data frame: method, NIE_type1, NDE_type1) and `raw`
(full results).

## Examples

``` r
null <- run_null_mediation_sim(n_iter = 2, n_samples = 100,
  n_features = 5, mo_confounding = 0.8, phi = 0.8)
null$rates
#>    method NIE_type1 NDE_type1
#> 1   UNADJ       0.1       1.0
#> 2  DIRECT       0.0       1.0
#> 3    COCA       0.0       0.0
#> 4  IV2SLS       0.0       0.5
#> 5     PGC       0.0       1.0
#> 6 IV2SLS2       0.0       0.5
```
