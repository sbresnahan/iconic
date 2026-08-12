# COCA estimator: Negative-Control Outcome Correction via ratio

Implements the Correlated Outcome Control Approach (COCA). Fits
`w ~ y + X + covars` and recovers the causal effect as \\\hat\beta =
-\hat\beta_X / \hat\beta_Y\\. Standard errors are obtained via the delta
method.

## Usage

``` r
fit_coca(y, X, w, covars = NULL, ratio_cap = 10, se_cap = 5)
```

## Arguments

- y:

  Numeric primary outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- w:

  Numeric negative-control outcome vector (length n). Recommended: pass
  `rowMeans(W_matrix)` for stability.

- covars:

  Optional data frame of additional covariates (n rows).

- ratio_cap:

  Maximum absolute value of the ratio estimate before flagging as
  unstable and returning `NA`. Default 10.

- se_cap:

  Maximum SE before flagging as unstable. Default 5.

## Value

Named list: `beta`, `se`, `pvalue`. Returns
`list(beta=NA, se=NA, pvalue=NA)` if estimation is unstable (near-zero
\\\hat\beta_Y\\ or extreme ratio).

## Details

The negative-control W should be an outcome that shares the same
unmeasured confounders as Y but has no direct causal path from X.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, seed = 1)
fit_coca(dat$Y[, 1], dat$X, rowMeans(dat$W))
#> $beta
#>          X 
#> 0.07120217 
#> 
#> $se
#> [1] 0.05005562
#> 
#> $pvalue
#>         X 
#> 0.1548925 
#> 
```
