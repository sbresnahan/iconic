# Composite null hypothesis test p-value

Computes the Huang (2019) JT-comp composite p-value for testing H0:
alpha \* beta = 0 against H1: alpha \* beta != 0.

## Usage

``` r
composite_p_value(a, b, var_a = 1, var_b = 1)
```

## Arguments

- a:

  Standardized z-statistic for alpha (alpha_hat / SE).

- b:

  Standardized z-statistic for beta (beta_hat / SE).

- var_a:

  Variance of the z-statistic for alpha across the collection of tests.
  Under the point null this is ~1; under H0(2)/H0(3) it is \>1. Default
  1 (conservative).

- var_b:

  Variance of the z-statistic for beta across the collection of tests.
  Default 1.

## Value

Numeric scalar in \\\[0, 1\]\\.

## Details

The composite null decomposes into three cases (alpha=beta=0,
alpha!=0/beta=0, alpha=0/beta!=0). The Wald/Sobel test is conservative
under the first case because the product of two independent normals
follows a normal product distribution, not a normal. This function
computes a closed-form p-value that accounts for all three cases without
estimating their proportions.

## References

Huang, Y.-T. (2019). Genome-wide analyses of sparse mediation effects
under composite null hypotheses. *Annals of Applied Statistics*, 13(1),
60-84.

## Examples

``` r
composite_p_value(0, 0) # null: p = 1
#> [1] 1
composite_p_value(2, 2) # strong signal: small p
#> [1] 0.006459626
composite_p_value(2, 0) # H0(2): p = 1 (no mediation)
#> [1] 1
```
