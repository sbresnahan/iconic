# Check negative-control validity across confounding scenarios

Focuses on whether the negative-control estimators (COCA, PGC) hold as
the controls' coverage of the confounder subspace drops and as the
number of latent confounders grows. Reports their bias/RMSE alongside
IV2SLS, and flags scenarios where identification is not credible — in
particular when the number of latent confounders exceeds the effective
number of valid controls (the proximal-inference completeness
condition).

## Usage

``` r
nc_validity_check(
  trained_gan = NULL,
  coverage_grid = c(0.2, 0.5, 0.8, 1),
  k_grid = c(1, 2, 3),
  conf_strength = 0.8,
  n_valid_controls = 1,
  n_iter = 50,
  n_samples = 500,
  n_features = 20,
  nc_model = "proxy",
  base_seed = 800,
  n_cores = 1
)
```

## Arguments

- trained_gan:

  An `iconic_gan` (or `NULL`).

- coverage_grid:

  Coverage values to sweep. Default `c(0.2, 0.5, 0.8, 1)`.

- k_grid:

  Numbers of confounders. Default `c(1, 2, 3)`.

- conf_strength:

  Fixed confounding strength. Default 0.8.

- n_valid_controls:

  Number of *distinct valid* controls the design provides (for the
  identifiability flag). Default 1.

- n_iter, n_samples, n_features, nc_model, base_seed, n_cores:

  As in
  [`gan_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_sensitivity.md).

## Value

A list with `summary` (COCA/PGC/IV2SLS bias & RMSE per scenario, with an
`identified` flag) and `verdict` (short per-scenario diagnosis).

## Details

The matrix-bridge PGC is the estimator for which the completeness
condition is binding: when `k > n_valid_controls`, PGC is
under-identified while IV2SLS remains unbiased (it does not depend on NC
completeness).

## Examples

``` r
chk <- nc_validity_check(NULL, coverage_grid = 0.7, k_grid = 1,
  n_iter = 2, n_samples = 100, n_features = 5)
chk$verdict
#>       coverage k identified nc_abs_bias                 diagnosis
#> 0.7.1      0.7 1       TRUE       0.307 negative controls holding
```
