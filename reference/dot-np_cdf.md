# Normal product distribution CDF

Computes F(z) = 2 \* \\integral\_{\|z\|}^{Inf}\\ K0(x)/pi dx, the
two-sided tail probability of the standard normal product distribution.
Under H0(1), if Z1, Z2 ~ N(0,1) independent, then P(\|Z1\*Z2\| \>=
\|z\|) = F(z). K0 is the modified Bessel function (besselK(x, nu=0)).

## Usage

``` r
.np_cdf(z)
```

## Arguments

- z:

  Numeric scalar or vector.

## Value

Numeric scalar or vector of CDF values in \\\[0, 1\]\\.
