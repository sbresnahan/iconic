# IV2SLS2 survival mediation estimator: 2-stage MR with Cox / RMST outcome

Two-stage Mendelian-randomization mediation with a survival outcome
stage (2SPS), with optional **path-specific** negative-control (NC)
augmentation:

1.  OLS: `X ~ g (+ W1) + covars` -\> X_hat (purge U1 from X).

2.  OLS: `M ~ X_hat + gm (+ W2) + covars` -\> M_hat, alpha_M.

3.  Cox / RMST: `Surv(t,e) ~ X_hat + M_hat (+ W2) + covars` -\> NDE
    (coef on X_hat), beta_M (coef on M_hat).

`NIE = alpha_M * beta_M`. Weak-instrument gates (partial F for G and Gm)
apply to the OLS first stages.

## Usage

``` r
fit_iv2sls_mediation2_surv(
  time,
  event,
  X,
  M,
  g,
  gm,
  covars = NULL,
  min_f = 10,
  W1 = NULL,
  W2 = NULL,
  w = NULL,
  effect_scale = c("loghr", "rmst"),
  tau = NULL
)
```

## Arguments

- time:

  Numeric follow-up time vector (length n).

- event:

  Numeric 0/1 event indicator (length n).

- X:

  Numeric exposure vector (length n).

- M:

  Numeric mediator vector (length n).

- g:

  Numeric instrument for X (length n).

- gm:

  Numeric instrument for M (length n).

- covars:

  Optional data frame of covariates (n rows).

- min_f:

  Minimum partial F for each excluded instrument. Default 10.

- W1:

  Optional NC panel (vector length n or matrix n x q) proxying the
  exposure-mediator confounder (X-\>M path); added to stage 1 only.
  Default `NULL`.

- W2:

  Optional NC panel (vector length n or matrix n x q) proxying the
  mediator-outcome confounder (M-\>Y path); added to stages 2 and 3.
  Default `NULL`.

- w:

  Defunct. The pooled single-panel argument was removed (collider under
  multi-confounder designs). Use `W1` and/or `W2` instead.

- effect_scale:

  Character: `"loghr"` or `"rmst"`.

- tau:

  RMST horizon (rmst only).

## Value

Named list (same fields as
[`fit_unadj_mediation_surv`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_surv.md)).
Returns all-NA if either first-stage partial F is below `min_f`.

## Details

**Path-specific NC augmentation (optional).** `W1` proxies the
exposure-mediator confounder (U1, X-\>M path) and is added to stage 1
only; `W2` proxies the mediator-outcome confounder (U2, M-\>Y path) and
is added to stages 2 and 3. Either panel may be omitted; with both
`NULL` the estimator reduces to plain two-instrument 2-stage MR.
Conditioning on a *pooled* panel in all three stages is a collider under
multi-confounder designs and is not supported: identical `W1`/`W2` are
treated as absent (pure MR). When `W1` and `W2` are distinct-noise
proxies of the *same* latent composite (single-confounder design), their
column spaces are near-collinear and the estimator likewise falls back
to plain two-instrument 2-stage MR.

## References

Rudolph, K. E., et al. (2024). Natural direct and indirect effects with
an instrumental variable. *Biometrics*.

## Examples

``` r
set.seed(1)
dat <- generate_toy_data(n = 500, outcome_type = "survival",
mo_confounding = 0.8, phi = 0.8, lambda_XM = c(1, 0),
lambda_MY = c(0, 1), omega_1 = 0.7, omega_2 = 0.7, seed = 1)
fit_iv2sls_mediation2_surv(dat$surv_time, dat$surv_event, dat$X, dat$M,
dat$G[, 1], dat$Gm, W1 = dat$W1, W2 = dat$W2)
#> $NDE
#> [1] 0.0954067
#> 
#> $NDE_se
#> [1] 1.100106
#> 
#> $NDE_p
#> [1] 0.2662602
#> 
#> $NIE
#> [1] 0.137888
#> 
#> $NIE_se
#> [1] 0.6625555
#> 
#> $NIE_p
#> [1] 0.8351389
#> 
#> $alpha_M
#> [1] 0.5039503
#> 
#> $alpha_se
#> [1] 0.01205152
#> 
#> $beta_M
#> [1] 0.2736143
#> 
#> $beta_M_se
#> [1] 1.314708
#> 
```
