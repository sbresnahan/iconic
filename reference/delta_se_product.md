# Delta-method SE for a product alpha \* beta (internal)

Delta-method SE for a product alpha \* beta (internal)

## Usage

``` r
delta_se_product(alpha, alpha_se, beta, beta_se, cov_ab = 0)
```

## Arguments

- alpha:

  Point estimate of alpha.

- alpha_se:

  SE of alpha.

- beta:

  Point estimate of beta.

- beta_se:

  SE of beta.

- cov_ab:

  Estimated covariance of alpha and beta. Default 0.

## Value

Scalar SE of the product alpha \* beta.
