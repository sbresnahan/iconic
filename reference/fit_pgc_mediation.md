# PGC mediation estimator: bridge-function-adjusted natural effects (matrix bridge)

Extends
[`fit_pgc`](https://seantbresnahan.com/iconic/reference/fit_pgc.md)
(matrix bridge) to the mediation setting by constructing a confounding
proxy \\\hat W\\ from the full W matrix and including it in both the
mediator and outcome regressions.

## Usage

``` r
fit_pgc_mediation(y, X, M, g, W, covars = NULL)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- M:

  Numeric mediator vector (length n).

- g:

  Numeric instrument vector (length n).

- W:

  Numeric negative-control matrix (n x q) or vector (length n). If a
  matrix, the bridge uses all q columns.

- covars:

  Optional data frame of additional covariates (n rows).

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`.

## Details

Steps:

1.  Residualise X on G -\> X_resid (U-driven residual).

2.  Bridge X_resid on the FULL W matrix -\> W_hat (proxy for U).

3.  `M ~ X + W_hat` -\> alpha_M (adjusted for confounding proxy).

4.  `Y ~ X + M + W_hat` -\> NDE = beta_X, beta_M (adjusted).

The matrix bridge requires `ncol(W) >= k` (proximal completeness) for
the bridge to span the confounder subspace.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 200, mo_confounding = 0.8, seed = 1)
fit_pgc_mediation(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$W)
#> $NDE
#> [1] -0.3886088
#> 
#> $NDE_se
#> [1] 0.04901454
#> 
#> $NDE_p
#> [1] 1.647096e-13
#> 
#> $NIE
#> [1] 0.7579615
#> 
#> $NIE_se
#> [1] 0.05033667
#> 
#> $NIE_p
#> [1] 3.066427e-51
#> 
#> $alpha_M
#> [1] 0.5826071
#> 
#> $alpha_se
#> [1] 0.01642378
#> 
#> $beta_M
#> [1] 1.300982
#> 
#> $beta_M_se
#> [1] 0.07822876
#> 
```
