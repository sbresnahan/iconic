# CpG-predicted-expression negative-control model (SCENIC case)

Simulates spatially correlated CpG methylation whose signal is partly
driven by the captured confounders, then forms each negative control as
a linear prediction from the methylation sites ("CpG-predicted
expression"). The controls therefore carry confounder information only
to the extent the methylation does, mediated through a realistic spatial
methylation layer.

## Usage

``` r
nc_cpg(U, covariates, params)
```

## Arguments

- U:

  `n x k` confounder matrix.

- covariates:

  Covariate data frame (unused; kept for the NC contract).

- params:

  List with `n_features`, and optionally `coverage`
  (confounder-\>methylation strength, default 0.7), `captured`
  (confounder indices, default all), `n_cpg` (methylation sites, default
  60), `rho` (AR(1) spatial correlation across sites, default 0.6),
  `MMCpG` (methylation-confounding multiplier, default 1), `MMCon`
  (default 1), and `noise_cor` (a `p x p` correlation matrix for
  correlated idiosyncratic noise, or `NULL` for independent noise).

## Value

`n x n_features` matrix of CpG-predicted negative controls.

## Details

When `noise_cor` is supplied (a `p x p` correlation matrix), the
idiosyncratic noise added to each control is drawn from a multivariate
normal with that correlation structure, so the controls retain realistic
cross-feature correlations conditional on the confounder.

## Examples

``` r
U <- matrix(rnorm(100), 100, 1)
W <- nc_cpg(U, covariates = NULL,
  params = list(n_features = 10, coverage = 0.7, n_cpg = 20))
dim(W)
#> [1] 100  10
```
