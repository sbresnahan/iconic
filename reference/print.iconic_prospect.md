# Print method for iconic_prospect objects

Print method for iconic_prospect objects

## Usage

``` r
# S3 method for class 'iconic_prospect'
print(x, ...)
```

## Arguments

- x:

  An `iconic_prospect` object.

- ...:

  Unused.

## Value

Invisibly returns `x` (the `iconic_prospect` object); called for its
side effect of printing a human-readable summary.

## Examples

``` r
if (check_torch_setup()) {
  data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    M = rnorm(100))
  result <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
    gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
  print(result)
}
```
