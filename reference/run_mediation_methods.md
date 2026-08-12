# Apply all mediation estimators across features (internal)

Apply all mediation estimators across features (internal)

## Usage

``` r
run_mediation_methods(
  dat,
  n_features = ncol(dat$Y),
  W_valid = NULL,
  n_cores = 1,
  se_method = "delta",
  n_boot = 500
)
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

- se_method:

  "delta" (default), "bootstrap", or "composite".

- n_boot:

  Number of bootstrap resamples when `se_method="bootstrap"`.

## Value

Data frame with columns: feature, method, NDE, NDE_se, NDE_p, NIE,
NIE_se, NIE_p. When `se_method="composite"`, also includes alpha_M,
alpha_se, beta_M, beta_M_se, var_a, var_b.
