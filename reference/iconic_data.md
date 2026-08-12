# Construct a standardized data object for ICONIC model selection

Creates an `iconic_data` S3 object from the user's real data,
standardizing vectors and matrices into a consistent format that all
downstream model selection functions consume.

## Usage

``` r
iconic_data(
  X,
  Y = NULL,
  M = NULL,
  G = NULL,
  Gm = NULL,
  W = NULL,
  W1 = NULL,
  W2 = NULL,
  covariates = NULL,
  feature_names = NULL,
  mediator_names = NULL,
  trained_gan = NULL,
  outcome_type = c("continuous", "survival"),
  surv_time = NULL,
  surv_event = NULL,
  scale = TRUE,
  recycle_lone_panel = FALSE,
  Z = NULL
)
```

## Arguments

- X:

  Exposure: numeric vector (length n) or features x samples matrix. If a
  matrix, column means are taken and scaled (one exposure per sample).

- Y:

  Outcome: numeric vector (length n) or features x samples matrix. When
  a matrix, estimation runs per-feature.

- M:

  Optional mediator: numeric vector (length n) or features x samples
  matrix. When a matrix, estimation runs per-mediator x per-outcome.

- G:

  Optional exposure instrument: numeric vector (length n) or n x
  n_features matrix (as returned by
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md)).
  If a matrix, the first column is extracted. E.g., a polygenic risk
  score.

- Gm:

  Optional mediator instrument: numeric vector (length n) or matrix (n x
  n_mediators). When a matrix, each column is the instrument for the
  corresponding mediator (e.g., per-isoform cis-eQTLs).

- W:

  Optional negative-control panel: features x samples matrix.
  Single-panel NCs used for COCA, PGC.

- W1:

  Optional path-specific NCs for the X-\>M path: features x samples
  matrix. Captures conf_XM. When W1/W2 are absent but W is present, W1 =
  W2 = W.

- W2:

  Optional path-specific NCs for the M-\>Y path: features x samples
  matrix. Captures conf_MY.

- covariates:

  Optional data frame of sample-level covariates (n rows). Recognized
  columns `sex`, `GA`, `mother_ethnicity` are encoded; other numeric
  columns are z-scored. Names colliding with estimator-reserved tokens
  are renamed.

- feature_names:

  Optional character vector of outcome feature names.

- mediator_names:

  Optional character vector of mediator names.

- trained_gan:

  Optional `iconic_gan` from
  [`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md).
  When supplied, the GAN is attached to the data object and reused by
  [`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
  and
  [`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md)
  instead of auto-training a new one. This avoids retraining when the
  same data is used across multiple workflow steps.

- outcome_type:

  Character: `"continuous"` (default, backward-compatible) or
  `"survival"`. When `"survival"`, `Y` is not required; instead supply
  `surv_time` and `surv_event`. Estimation uses Cox proportional-hazards
  ([`coxph`](https://rdrr.io/pkg/survival/man/coxph.html)) or RMST
  pseudo-observation OLS (see `effect_scale` in
  [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md)).

- surv_time:

  Numeric follow-up time vector (length n). Required when
  `outcome_type = "survival"`; ignored otherwise.

- surv_event:

  Numeric 0/1 event indicator (length n; 1 = event observed, 0 =
  censored). Required when `outcome_type = "survival"`; ignored
  otherwise.

- scale:

  Logical: center and scale all continuous inputs (X, Y, M, G, Gm, W,
  W1, W2, and numeric covariates) to mean 0 / sd 1. Default `TRUE`.
  Scaling parameters are recorded in `$scaling` for back-transformation.
  Set `FALSE` to preserve the original scale.

- recycle_lone_panel:

  Logical: when exactly one of `W1` / `W2` is supplied (no pooled `W`),
  use that lone panel as BOTH path-specific bridges (`W1 = W2`), making
  the two-bridge estimators PGC2 / PGC2Gm eligible. Default `FALSE`: a
  lone panel is retained for IV2SLS2's path-specific augmentation and
  used to derive the pooled `W` (so DIRECT / COCA / PGC and the NC
  validity screens run), but `has_path_nc` stays `FALSE` and PGC2 /
  PGC2Gm remain ineligible. Set to `TRUE` only when the single panel is
  assumed complete for BOTH path confounder composites (the shared-panel
  special case); IV2SLS2 is the more defensible primary estimator when
  coverage of the other path's composite is in doubt. A warning is
  emitted when the recycle is activated.

- Z:

  Defunct. Renamed to `X`; passing a value errors with a message
  pointing to `X`. Retained in the signature only to catch and redirect
  old calls.

## Value

An `iconic_data` S3 object: a named list with `$X`, `$Y`, `$M`, `$G`,
`$Gm`, `$W`, `$W1`, `$W2`, `$covariates`, `$n`, `$n_features`,
`$n_mediators`, `$has_instrument`, `$has_mediator_instrument`,
`$has_nc`, `$has_path_nc`, `$is_mediation`, `$feature_names`,
`$mediator_names`, `$trained_gan`, `$outcome_type`, and (when survival)
`$surv_time`, `$surv_event`.

## Examples

``` r
# Total-effect only (no mediation)
data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100))

# Full mediation with instruments and NCs
data <- iconic_data(
X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
W = matrix(rnorm(100 * 10), 10, 100)
)
print(data)
#> <iconic_data> 100 samples, 10 outcome features, 1 mediator(s)
#>  Available: G (exposure instrument), Gm (mediator instrument), W (negative controls), W1/W2 (path-specific NCs) 
#>  Mode: mediation 

# Survival outcome
data <- iconic_data(
X = rnorm(100), outcome_type = "survival",
surv_time = rexp(100), surv_event = rbinom(100, 1, 0.6),
M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
W = matrix(rnorm(100 * 10), 10, 100)
)
print(data)
#> <iconic_data> 100 samples, survival outcome (57 events / 100), 1 mediator(s)
#>  Available: G (exposure instrument), Gm (mediator instrument), W (negative controls), W1/W2 (path-specific NCs) 
#>  Mode: mediation 
```
