# Build a block-diagonal correlation matrix (internal)

Creates a `p x p` block-diagonal correlation matrix modelling
co-expression modules: features within the same block have pairwise
correlation `rho`, features in different blocks are uncorrelated. The
number of blocks is `ceiling(sqrt(p))`, with block sizes as equal as
possible.

## Usage

``` r
.block_diag_cor(p, rho)
```

## Arguments

- p:

  Number of features.

- rho:

  Within-block pairwise correlation (off-diagonal).

## Value

A `p x p` numeric matrix with unit diagonal.
