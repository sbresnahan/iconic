# Apply all mediation estimators to a single feature (internal)

Thin wrapper around .estimate_mediation_feature() that extracts vectors
from the `dat` list format used by the simulation pipeline.

## Usage

``` r
.analyze_mediation_feature(
  dat,
  f,
  W_avg,
  W_valid = NULL,
  se_method = "delta",
  n_boot = 500
)
```

## Arguments

- dat:

  Dataset list (from generate_toy_data / run_single_iteration).

- f:

  Feature (column) index.

- W_avg:

  Row means of the full negative-control panel (for COCA).

- W_valid:

  Optional: validity-screened W matrix for matrix-bridge PGC.

- se_method:

  "delta" (default) or "bootstrap".

- n_boot:

  Number of bootstrap resamples when `se_method="bootstrap"`.

## Value

Data frame of five to eight rows (one per method) or NULL. When `dat$M`
is a matrix (n_mediators \> 1), results are returned for each mediator
separately with a `mediator` column.

## Details

When `dat$Gm` is present (i.e. a mediator instrument was supplied), the
2-stage MR estimator
[`fit_iv2sls_mediation2`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2.md)
is also run, producing a sixth row (method = "IV2SLS2"). When `dat$W1`
and `dat$W2` are present, the two-stage proximal estimator
[`fit_pgc_mediation2`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2.md)
is run as "PGC2" (without Gm) and "PGC2Gm" (with Gm, when Gm is also
present).
