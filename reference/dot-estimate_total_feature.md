# Estimate total effect for a single feature with specified methods (internal)

Format-agnostic per-feature estimation: takes explicit vectors and
covariates (no dependency on the `dat` list format). Called by both the
simulation driver (.analyze_feature) and the real-data driver
(iconic_estimate).

## Usage

``` r
.estimate_total_feature(
  X,
  y,
  g = NULL,
  w = NULL,
  W_mat = NULL,
  W_avg = NULL,
  covars = NULL,
  methods = .methods_all,
  feature_idx = 1L,
  min_f = 10
)
```

## Arguments

- X:

  Numeric exposure vector (length n).

- y:

  Numeric outcome vector (length n).

- g:

  Numeric instrument vector (length n), or NULL.

- w:

  Numeric NC vector (length n) or matrix (n x q), or NULL. Used for
  DIRECT, IV2SLS (full panel as covariates), and COCA (via W_avg
  scalar).

- W_mat:

  Numeric NC matrix (n x q), or NULL. Used for PGC (matrix bridge). If
  NULL, PGC is skipped.

- W_avg:

  Numeric vector (length n): row means of the NC panel for COCA. If NULL
  but W_mat is present, computed inline.

- covars:

  Optional data frame of covariates (n rows).

- methods:

  Character vector of methods to run. Default: all five. Methods whose
  required inputs are missing are silently skipped.

- feature_idx:

  Integer or character: feature identifier for the output `feature`
  column. Default 1L.

## Value

Data frame: `feature`, `method`, `beta`, `se`, `pvalue`. Returns NULL if
fewer than 20 complete cases.
