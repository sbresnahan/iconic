# Estimate mediation effects for a single feature with specified methods (internal)

Format-agnostic per-feature mediation estimation: takes explicit vectors
and covariates (no dependency on the `dat` list format). Called by both
the simulation driver (.analyze_mediation_feature) and the real-data
driver (iconic_estimate).

## Usage

``` r
.estimate_mediation_feature(
  X,
  y,
  M_vec,
  g = NULL,
  gm = NULL,
  w = NULL,
  W_mat = NULL,
  W1_mat = NULL,
  W2_mat = NULL,
  W_avg = NULL,
  covars = NULL,
  methods = NULL,
  feature_idx = 1L,
  min_f = 10,
  se_method = c("delta", "bootstrap", "composite"),
  n_boot = 500
)
```

## Arguments

- X:

  Numeric exposure vector (length n).

- y:

  Numeric outcome vector (length n).

- M_vec:

  Numeric mediator vector (length n).

- g:

  Numeric instrument vector (length n), or NULL.

- gm:

  Numeric mediator instrument vector (length n), or NULL.

- w:

  Numeric NC vector (length n) or matrix (n x q), or NULL (for DIRECT,
  IV2SLS, IV2SLS2 — full panel as covariates).

- W_mat:

  Numeric NC matrix (n x q), or NULL (for PGC matrix bridge).

- W1_mat:

  Numeric path-specific NC matrix (n x q) for X-\>M, or NULL.

- W2_mat:

  Numeric path-specific NC matrix (n x q) for M-\>Y, or NULL.

- W_avg:

  Numeric vector (length n): row means of NC panel for COCA. If NULL but
  W_mat present, computed inline.

- covars:

  Optional data frame of covariates (n rows).

- methods:

  Character vector of methods to run. Default: all applicable.

- feature_idx:

  Integer or character: feature identifier. Default 1L.

- se_method:

  "delta" (default), "bootstrap", or "composite". When "bootstrap",
  NDE_se/NIE_se are replaced by the SD of `n_boot` nonparametric
  bootstrap resamples. Point estimates (NDE/NIE) and p-values are
  unchanged; only the SE column is swapped. When "composite", NIE_p is
  replaced by the Huang (2019) JT-comp composite null p-value
  (post-processed across features by the driver).

- n_boot:

  Number of bootstrap resamples when `se_method="bootstrap"`.

## Value

Data frame: `feature`, `method`, `NDE`, `NDE_se`, `NDE_p`, `NIE`,
`NIE_se`, `NIE_p`. Returns NULL if \< 20 complete cases. When
`se_method="composite"`, also includes `alpha_M`, `alpha_se`, `beta_M`,
`beta_M_se` for downstream composite p-value computation.
