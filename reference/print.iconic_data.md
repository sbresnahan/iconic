# Print method for iconic_data objects

Print method for iconic_data objects

## Usage

``` r
# S3 method for class 'iconic_data'
print(x, ...)
```

## Arguments

- x:

  An `iconic_data` object.

- ...:

  Unused.

## Value

Invisibly returns `x` (the `iconic_data` object); called for its side
effect of printing a human-readable summary.

## Examples

``` r
data <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 5), 5, 50),
  G = rnorm(50), W = matrix(rnorm(50 * 5), 5, 50))
print(data)
#> <iconic_data> 50 samples, 5 outcome features
#>  Available: G (exposure instrument), W (negative controls), W1/W2 (path-specific NCs) 
#>  Mode: total effect 
```
