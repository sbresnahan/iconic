# Check first-stage instrument strength (partial F)

Regresses the exposure `X` on the instrument `G` (plus optional
covariates) and reports the partial F statistic and partial R2 for the
instrument — the weak-instrument diagnostic. The partial F (not the
overall model F) is the relevant quantity when covariates such as
ancestry PCs explain much of the exposure.

## Usage

``` r
check_instrument_strength(G, X, covariates = NULL, min_f = 10)
```

## Arguments

- G:

  Numeric vector (length n): the exposure instrument (e.g. a polygenic
  score from
  [`build_prs_ldpred2()`](https://seantbresnahan.com/iconic/reference/build_prs_ldpred2.md)
  or
  [`score_pgs_panel()`](https://seantbresnahan.com/iconic/reference/score_pgs_panel.md)).

- X:

  Numeric vector (length n): the exposure.

- covariates:

  Optional data.frame or matrix (n rows) of covariates to partial out
  (e.g. ancestry PCs).

- min_f:

  Numeric: weak-instrument threshold on the partial F. Default 10
  (Staiger-Stock rule of thumb).

## Value

A list with `F` (partial F statistic), `df1`, `df2`, `pvalue`,
`partial_r2`, `n`, and `weak` (`TRUE` when F \< `min_f`).

## Examples

``` r
n <- 300
G <- rnorm(n)
X <- 0.3 * G + rnorm(n)
pcs <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("PC", 1:3)))
check_instrument_strength(G, X, covariates = pcs)
#> $F
#> [1] 38.16507
#> 
#> $df1
#> [1] 1
#> 
#> $df2
#> [1] 295
#> 
#> $pvalue
#> [1] 2.152696e-09
#> 
#> $partial_r2
#> [1] 0.114553
#> 
#> $n
#> [1] 300
#> 
#> $weak
#> [1] FALSE
#> 
```
