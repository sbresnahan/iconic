# PGC-2 mediation estimator: two-stage proximal mediation with path-specific bridges

Extends the proximal-inference approach to the two-linked-DAG mediation
setting. Uses **path-specific** negative controls — W1 capturing the
X-\>M confounder conf_XM and W2 capturing the M-\>Y confounder conf_MY —
to purge confounding at both stages via bridge functions.

## Usage

``` r
fit_pgc_mediation2(y, X, M, g, W1, W2, gm = NULL, covars = NULL, min_f = 10)
```

## Arguments

- y:

  Numeric outcome vector (length n).

- X:

  Numeric exposure vector (length n).

- M:

  Numeric mediator vector (length n).

- g:

  Numeric instrument for X (length n).

- W1:

  Numeric negative-control matrix (n x q) or vector for the X-\>M path
  (captures conf_XM). If a matrix, the bridge uses all q columns.

- W2:

  Numeric negative-control matrix (n x q) or vector for the M-\>Y path
  (captures conf_MY). If a matrix, the bridge uses all q columns.

- gm:

  Optional numeric mediator instrument vector (length n). When `NULL`
  (default), Stage 2 uses pure NC identification. When supplied, Gm
  helps isolate conf_MY before bridging — robust to Gm-U correlation.

- covars:

  Optional data frame of additional covariates (n rows).

- min_f:

  Minimum acceptable partial F-statistic for the excluded instrument G1.
  Default 10.

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`. Returns
all-`NA` if the first-stage partial F for G1 is below `min_f`.

## Details

When `gm` is `NULL` (no mediator instrument), Stage 2 residualises M on
the cleaned exposure X_hat to isolate the conf_MY-driven component, then
bridges on W2. This is pure negative-control identification at both
stages — no mediator instrument is required.

When `gm` is supplied (mediator instrument present but its exogeneity
may be in doubt), Stage 2 uses Gm to help isolate conf_MY's effect on M
before bridging on W2. The bridge W_hat_M does the confounding removal,
so the estimator is robust to Gm being correlated with conf_MY — the
residual bias from Gm-U correlation is smaller than IV2SLS2's.

Strategy (three stages):

1.  **Bridge for X** (purge conf_XM from X):
    `X_resid = residuals(X ~ G1 + C)`; `W_hat_X = bridge(X_resid ~ W1)`
    (fitted values, proxy for conf_XM);
    `X_hat = fitted(X ~ G1 + W_hat_X + C)`. Weak-IV check: partial F for
    G1 \>= `min_f`.

2.  **Bridge for M** (purge conf_MY from M):
    `M_resid = residuals(M ~ X_hat + C)` (`gm = NULL`), or
    `M_resid = residuals(M ~ Gm + C)` (`gm` supplied);
    `W_hat_M = bridge(M_resid ~ W2)` (fitted values, proxy for conf_MY);
    `M_hat = fitted(M ~ X_hat + W_hat_M + C)` (`gm = NULL`), or
    `M_hat = fitted(M ~ X_hat + Gm + W_hat_M + C)` (`gm` supplied);
    `alpha_M = coefficient on X_hat`.

3.  **Outcome**: `Y ~ X_hat + M_hat + W_hat_X + W_hat_M + C`;
    `NDE = coefficient on X_hat`, `beta_M = coefficient on M_hat`.

NIE = alpha_M \* beta_M (delta-method SE).

The key advantage over
[`fit_iv2sls_mediation2`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2.md):
PGC-2 does not require instrument exogeneity. When the mediator
instrument Gm is correlated with the confounder conf_MY (rho_G2 \> 0),
IV2SLS2 is biased but PGC-2's bridge absorbs conf_MY regardless of the
instrument violation. The tipping-point simulation maps where PGC-2 bias
crosses below IV2SLS2 bias.

## References

Miao, W., Geng, Z., & Tchetgen Tchetgen, E. (2018). Identifying causal
effects with proxy variables of an unmeasured confounder. *Biometrika*,
105(4), 987-993.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 500, mo_confounding = 0.8,
rho_G2 = 0.3, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
omega_1 = 0.7, omega_2 = 0.7, seed = 1)
# Without mediator instrument (pure NC identification)
fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1, dat$W1, dat$W2)
#> $NDE
#> [1] -0.5101239
#> 
#> $NDE_se
#> [1] 0.02827378
#> 
#> $NDE_p
#> [1] 2.55428e-56
#> 
#> $NIE
#> [1] 0.739548
#> 
#> $NIE_se
#> [1] 0.02236263
#> 
#> $NIE_p
#> [1] 7.839686e-240
#> 
#> $alpha_M
#> [1] 0.5002583
#> 
#> $alpha_se
#> [1] 0.0116492
#> 
#> $beta_M
#> [1] 1.478332
#> 
#> $beta_M_se
#> [1] 0.02851672
#> 
# With (possibly imperfect) mediator instrument
fit_pgc_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G1, dat$W1, dat$W2,
gm = dat$Gm)
#> $NDE
#> [1] -0.5101239
#> 
#> $NDE_se
#> [1] 0.02827378
#> 
#> $NDE_p
#> [1] 2.55428e-56
#> 
#> $NIE
#> [1] 0.739548
#> 
#> $NIE_se
#> [1] 0.02236263
#> 
#> $NIE_p
#> [1] 7.839686e-240
#> 
#> $alpha_M
#> [1] 0.5002583
#> 
#> $alpha_se
#> [1] 0.0116492
#> 
#> $beta_M
#> [1] 1.478332
#> 
#> $beta_M_se
#> [1] 0.02851672
#> 
```
