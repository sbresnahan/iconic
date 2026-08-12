# Apply all estimators across features (internal)

Apply all estimators across features (internal)

## Usage

``` r
run_methods(dat, n_features = ncol(dat$Y), W_valid = NULL, n_cores = 1)
```

## Arguments

- dat:

  List returned by
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md)
  /
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md).

- n_features:

  Number of outcome columns to process.

- W_valid:

  Optional: validity-screened W matrix for matrix-bridge PGC.

- n_cores:

  Number of parallel workers. Default 1 (sequential).

## Value

Data frame with columns: feature, method, beta, se, pvalue.
