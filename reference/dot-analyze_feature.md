# Analyse a single outcome feature with all estimators (internal)

Thin wrapper around .estimate_total_feature() that extracts vectors from
the `dat` list format used by the simulation pipeline.

## Usage

``` r
.analyze_feature(dat, f, W_avg, W_valid = NULL)
```

## Arguments

- dat:

  Dataset list (from generate_toy_data / run_single_iteration).

- f:

  Feature (column) index into `dat$Y`, `dat$W`, `dat$G`.

- W_avg:

  Row means of the full negative-control panel (for COCA).

- W_valid:

  Optional: validity-screened W matrix (n x q) for the matrix-bridge
  PGC. If NULL, uses the full W panel.

## Value

Data frame of five rows (one per method) or `NULL` if too few cases.
