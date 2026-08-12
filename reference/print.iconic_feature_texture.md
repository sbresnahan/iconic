# Print method for iconic_feature_texture objects

Print method for iconic_feature_texture objects

## Usage

``` r
# S3 method for class 'iconic_feature_texture'
print(x, ...)
```

## Arguments

- x:

  An `iconic_feature_texture` object.

- ...:

  Unused.

## Value

Invisibly returns `x` (the `iconic_feature_texture` object); called for
its side effect of printing a human-readable summary.

## Examples

``` r
ft <- train_feature_texture(matrix(rnorm(20 * 100), 20, 100))
print(ft)
#> <iconic_feature_texture>
#>  features : 20 
#>  samples : 100 
#>  method : auto 
#>  marginals : normal (20) 
#>  copula cor: mean=-0.000, range=[-0.235, 0.238]
```
