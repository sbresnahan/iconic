# IV2SLS2 mediation estimator: 2-stage MR with instruments for both X and M

Uses TWO genetic instruments: G for the exposure X and Gm for the
mediator M. This is the key extension: by instrumenting both endogenous
variables, NDE and NIE become **point-identified** even under
mediator-outcome (M-O) confounding, resolving the identification failure
that limits the single-instrument estimators.

## Usage

``` r
fit_iv2sls_mediation2(
  y,
  X,
  M,
  g,
  gm,
  covars = NULL,
  min_f = 10,
  W1 = NULL,
  W2 = NULL,
  w = NULL
)
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

- gm:

  Numeric instrument for M (length n).

- covars:

  Optional data frame of additional covariates (n rows).

- min_f:

  Minimum acceptable partial F-statistic for each excluded instrument.
  Default 10.

- W1:

  Optional negative-control panel (vector length n or matrix n x q)
  proxying the exposure-mediator confounder (X-\>M path); added to stage
  1 only. Default `NULL` (stage 1 unaugmented).

- W2:

  Optional negative-control panel (vector length n or matrix n x q)
  proxying the mediator-outcome confounder (M-\>Y path); added to stages
  2 and 3. Default `NULL` (stages 2–3 unaugmented).

- w:

  Defunct. The pooled single-panel argument was removed because
  conditioning on a pooled panel in all three stages opens a collider
  path under multi-confounder designs. Use `W1` and/or `W2` instead.

## Value

Named list: `NDE`, `NDE_se`, `NDE_p`, `NIE`, `NIE_se`, `NIE_p`. Returns
all-`NA` if either first-stage partial F is below `min_f`.

## Details

The motivating design uses two distinct genetic instruments: a mediator
instrument (e.g. a cis-eQTL) instruments the mediator (M), while a
separate exposure instrument (e.g. a polygenic score) instruments the
exposure (X). The mediator set must be restricted to features for which
a mediator instrument is available.

Strategy (sequential 2SLS, three OLS stages), with optional
**path-specific** negative-control (NC) augmentation:

1.  `X ~ G (+ W1) + covars` -\> X_hat (purge U1 from X). Weak-IV check:
    partial F for G \>= `min_f`.

2.  `M ~ X_hat + Gm (+ W2) + covars` -\> M_hat, alpha_M = coef on X_hat.
    Weak-IV check: partial F for Gm \>= `min_f`.

3.  `Y ~ X_hat + M_hat (+ W2) + covars` -\> NDE = beta_X_hat, beta_M =
    coef on M_hat.

NIE = alpha_M \* beta_M (delta-method SE).

**Path-specific NC augmentation (optional).** `W1` proxies the
exposure-mediator confounder (U1, the X-\>M path) and is added to stage
1 only; `W2` proxies the mediator-outcome confounder (U2, the M-\>Y
path) and is added to stages 2 and 3. Either panel may be omitted: with
`W1 = NULL` stage 1 is unaugmented, with `W2 = NULL` stages 2–3 are
unaugmented, and with both `NULL` the estimator reduces to plain
two-instrument 2-stage MR. Conditioning on a *pooled* panel
`(W1 + W2) / 2` in all three stages is a collider under multi-confounder
designs (the pooled panel is a common child of the independent
confounders U1 and U2) and is therefore not supported: if `W1` and `W2`
are identical they are treated as absent (pure MR). Pass the two
path-specific panels separately to obtain the coverage-improves-accuracy
behaviour.

Path-specific augmentation is only defined when `W1` and `W2` proxy
*distinct* confounders. When the two panels are distinct-noise proxies
of the *same* latent composite (the single-confounder / shared-loading
design), their column spaces are near-collinear and stage-1 `W1` would
inject the shared M-\>Y confounder into `X_hat`. The estimator detects
this (leading-PC correlation between the panels) and falls back to plain
two-instrument 2-stage MR, exactly as for an identical pooled panel.
Genuinely distinct panels are unaffected.

Including X_hat in the M first-stage (stage 2) is essential: M
structurally depends on X, so the X -\> M path must be captured for
alpha_M and the NIE to be correctly estimated. Gm provides the exogenous
variation that identifies the M -\> Y effect net of U1 confounding.

The sequential 2SLS implementation matches the pattern of
[`fit_iv2sls_mediation`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation.md)
(three separate OLS stages). Standard errors in the outcome stage do not
account for the two-stage estimation of M_hat, consistent with the
existing estimators; the simulation benchmarks bias, not SE accuracy. A
unit test cross-validates this estimator's NDE and beta_M against
[`AER::ivreg`](https://rdrr.io/pkg/AER/man/ivreg.html) on the
just-identified system (`Y ~ X + M | G + Gm`, pure-MR form).

## References

Rudolph, K. E., et al. (2024). Natural direct and indirect effects with
an instrumental variable. *Biometrics*.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 500, mo_confounding = 0.8,
phi = 0.8, lambda_XM = c(1, 0), lambda_MY = c(0, 1),
omega_1 = 0.7, omega_2 = 0.7, seed = 1)
fit_iv2sls_mediation2(dat$Y[, 1], dat$X, dat$M, dat$G[, 1], dat$Gm,
W1 = dat$W1, W2 = dat$W2)
#> $NDE
#> [1] 0.1074552
#> 
#> $NDE_se
#> [1] 0.01497136
#> 
#> $NDE_p
#> [1] 2.738294e-12
#> 
#> $NIE
#> [1] 0.1524762
#> 
#> $NIE_se
#> [1] 0.008039873
#> 
#> $NIE_p
#> [1] 3.320173e-80
#> 
#> $alpha_M
#> [1] 0.5039503
#> 
#> $alpha_se
#> [1] 0.01205152
#> 
#> $beta_M
#> [1] 0.302562
#> 
#> $beta_M_se
#> [1] 0.01421859
#> 
```
