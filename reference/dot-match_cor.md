# Match a correlation matrix to the target feature count (internal)

If `cor_mat` is a valid `p_target x p_target` matrix, return it as-is.
If it's a different size, return NULL (the simulation falls back to
independent noise). This guards against dimension mismatches when the
simulation uses a different `n_features` than the training data.

## Usage

``` r
.match_cor(cor_mat, p_target)
```

## Arguments

- cor_mat:

  A correlation matrix or NULL.

- p_target:

  Target number of features.

## Value

A `p_target x p_target` correlation matrix, or NULL.
