# Scenario manifest: truth and parameter ranges for a simulation

Returns a reader-orienting summary of a simulation scenario: the
ground-truth estimands (NDE, NIE, total effect), the modifiable
parameters with their swept ranges, and the fixed parameters with their
values. Intended to be rendered as a table preceding simulation results
in the manuscript (, \#582: "state the truth and the parameter ranges up
front to orient the reader").

## Usage

``` r
scenario_manifest(
  dat_or_params,
  conf_grid = NULL,
  coverage_grid = NULL,
  mo_confounding_grid = NULL,
  phi_grid = NULL,
  rho_G1_grid = NULL,
  rho_G2_grid = NULL,
  rho_pop_grid = NULL,
  omega_1_grid = NULL,
  omega_2_grid = NULL
)
```

## Arguments

- dat_or_params:

  A `dat` list from
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md)
  or a named list of DGP parameters.

- conf_grid:

  Swept values of `conf_str` (confounding strength). Default NULL.

- coverage_grid:

  Swept values of NC coverage (`w_signal` / `omega`). Default NULL.

- mo_confounding_grid:

  Swept values of mediator-outcome confounding strength. Default NULL.

- phi_grid:

  Swept values of the mediator-instrument strength `phi`. Default NULL.

- rho_G1_grid:

  Swept values of the G1-conf_XM correlation. Default NULL.

- rho_G2_grid:

  Swept values of the G2-conf_MY correlation. Default NULL.

- rho_pop_grid:

  Swept values of the population-stratification correlation. Default
  NULL.

- omega_1_grid:

  Swept values of `omega_1` (coverage of conf_XM by W1). Default NULL.

- omega_2_grid:

  Swept values of `omega_2` (coverage of conf_MY by W2). Default NULL.

## Value

A named list with elements:

- `estimands`:

  Named numeric vector: `NDE`, `NIE`, `total`.

- `modifiable_parameters`:

  Data frame with columns `parameter`, `value`, `swept_range`. `value`
  is the scalar used (or NA when only a grid was supplied);
  `swept_range` is a comma-separated string of grid values, or NA when
  the parameter was held fixed.

- `fixed_parameters`:

  Data frame with columns `parameter`, `value`.

## Details

`dat_or_params` may be either:

- a list returned by
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md),
  in which case the estimands and fixed parameters are read from the
  object (`true_NDE`, `true_NIE`, `true_total`, `n`, `n_features`,
  `n_mediators`, `lambda_XM`, `lambda_MY`, `feat_cor`); or

- a bare named list of DGP parameters (`beta_X`, `alpha_M`, `beta_M`,
  `n_mediators`, `n`, `n_features`, `lambda_XM`, `lambda_MY`,
  `feat_cor`, ...), in which case the estimands are recomputed as
  `NDE = beta_X`, `NIE = n_mediators * alpha_M * beta_M`,
  `total = NDE + NIE`.

The `*_grid` arguments record the swept ranges for the modifiable
parameters. When a grid argument is omitted or NULL, the corresponding
modifiable parameter is reported with its scalar value only (no range).

## Examples

``` r
dat <- generate_toy_data(n = 200, phi = 0.8, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
omega_1 = 0.7, omega_2 = 0.7, seed = 42)
scenario_manifest(dat, conf_grid = c(0.3, 0.8),
coverage_grid = c(0.3, 0.7, 1))
#> $estimands
#>   NDE   NIE total 
#>  0.10  0.15  0.25 
#> 
#> $modifiable_parameters
#>          parameter value   swept_range
#> 1         conf_str   0.8      0.3, 0.8
#> 2 w_signal / omega   0.7 0.3, 0.7, 1.0
#> 3   mo_confounding   0.0          <NA>
#> 4              phi   0.8          <NA>
#> 5           rho_G1   0.0          <NA>
#> 6           rho_G2   0.0          <NA>
#> 7          rho_pop   0.0          <NA>
#> 8          omega_1   0.7          <NA>
#> 9          omega_2   0.7          <NA>
#> 
#> $fixed_parameters
#>              parameter                  value
#> 1                    n                    200
#> 2           n_features                     20
#> 3          n_mediators                      1
#> 4 confounder_structure path-specific loadings
#> 5             feat_cor                      0
#> 
```
