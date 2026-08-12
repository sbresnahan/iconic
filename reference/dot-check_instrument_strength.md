# Check instrument strength via first-stage partial F (internal)

Computes partial F-statistics for G (in X ~ G + W + covars) and Gm (in M
~ X_hat + Gm + W + covars), matching the first-stage regressions used by
fit_iv2sls() and fit_iv2sls_mediation2().

## Usage

``` r
.check_instrument_strength(
  data,
  min_f = 10,
  g_threshold = NULL,
  gm_threshold = NULL,
  n_cores = 1
)
```

## Arguments

- data:

  iconic_data object.

- min_f:

  Minimum acceptable partial F. Default 10.

## Value

List: F_G, F_Gm, weak_G (logical), weak_Gm (logical).
