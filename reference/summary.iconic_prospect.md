# Summary method for iconic_prospect objects

Prints the full prospective summary (same as
[`print()`](https://rdrr.io/r/base/print.html)).

## Usage

``` r
# S3 method for class 'iconic_prospect'
summary(object, ...)
```

## Arguments

- object:

  An `iconic_prospect` object.

- ...:

  Unused.

## Value

Invisibly returns `object`.

## Examples

``` r
if (check_torch_setup()) {
  data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    M = rnorm(100))
  result <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
    gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
  summary(result)
}
```
