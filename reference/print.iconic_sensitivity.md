# Print method for iconic_sensitivity objects

Print method for iconic_sensitivity objects

## Usage

``` r
# S3 method for class 'iconic_sensitivity'
print(x, ...)
```

## Arguments

- x:

  An `iconic_sensitivity` object.

- ...:

  Unused.

## Value

Invisibly returns `x` (the `iconic_sensitivity` object); called for its
side effect of printing a human-readable summary.

## Examples

``` r
if (check_torch_setup()) {
  data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
    W = matrix(rnorm(100 * 10), 10, 100))
  sens <- iconic_sensitivity(data, n_iter = 2, gan_epochs = 5,
    rho_G1_grid = c(0, 0.2), rho_G2_grid = c(0, 0.2))
  print(sens)
}
```
