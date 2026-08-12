# Screen negative controls for exposure dependence (W *\|* X \| C)

Regresses each negative-control feature on the exposure X (plus observed
covariates C) and flags controls whose association with X survives
Benjamini-Hochberg FDR control at the requested level.

## Usage

``` r
nc_validity_screen(
  dat,
  fdr_level = 0.1,
  alpha = 0.05,
  n_cores = 1,
  criterion = c("both", "fdr", "magnitude"),
  magnitude_threshold = 0.1,
  u_proxy = NULL
)
```

## Arguments

- dat:

  Dataset list from
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  or
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md),
  containing `W`, `X`, and `synthetic_data` (covariates).

- fdr_level:

  Target false-discovery rate for BH correction. Default 0.10.

- alpha:

  Per-test significance level used before FDR adjustment (informational
  only). Default 0.05.

- n_cores:

  Number of parallel workers. Default 1 (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

- criterion:

  Which criterion to use: `"fdr"` (legacy), `"magnitude"`, or `"both"`
  (default).

- magnitude_threshold:

  \|partial_r(W,X\|C)\| cutoff for the magnitude branch. Default 0.10.

- u_proxy:

  Optional outcome-derived proxy for U (e.g. the leading
  residualized-outcome principal component) used to compute
  `relative_effect` = \|partial_r(W,X\|C)\| /
  \|partial_r(W,U_proxy\|C)\|. When `NULL` (default), `relative_effect`
  is `NA` and only the FDR/magnitude branches are used.

## Value

A data frame with one row per control feature: `feature`, `p_value`,
`fdr`, `partial_r`, `relative_effect`, `significant` (composite, per
`criterion`), `verdict` (composite), `verdict_fdr`, `verdict_magnitude`.

## Details

A control that is significantly associated with X after adjusting for C
violates the negative-control assumption (A1): it is either Such
controls should be dropped before running COCA / PGC.

## Assumption

A1 states that negative-control outcomes are independent of the exposure
**conditional on observed covariates C and the unmeasured confounder U**
(U is unobserved). NC outcome proxies may cause or be caused by the
outcome or instruments; the screen tests the observable marginal
association with X given C.

## Criterion

`criterion = "fdr"` (the legacy behavior) flags controls whose
association with X survives BH-FDR. `criterion = "magnitude"` flags
controls by partial-correlation *size* (\|partial_r(W,X\|C)\| \>
`magnitude_threshold`), which is less vulnerable to the "always
significant via U" problem: because W shares U with X by construction,
the FDR test can flag valid controls as violations purely from the
intended confounder-sharing signal. `criterion = "both"` (default)
requires a control to pass both branches to be deemed valid.

Distinguishing "W downstream of X" (a true A1 violation) from "W shares
a cause with X" (the intended NC behavior) remains an open problem
without a clean test statistic. The `relative_effect` column
(\|partial_r(W,X\|C)\| / \|partial_r(W,U_proxy\|C)\|, when `u_proxy` is
supplied) is exposed as a *diagnostic*: a control downstream of X tends
to have a larger W-X association relative to its W-U association, but
the separation is not always clean. Analysts should inspect `partial_r`
and `relative_effect` together rather than relying on a single gate.

## Examples

``` r
dat <- run_single_iteration(n_features = 10, seed = 1)
nc_validity_screen(dat)
#> NC validity screen: 10 tasks (sequential)
#>  NC validity screen: 10% (1/10) [0s]
#>  NC validity screen: 20% (2/10) [0s]
#>  NC validity screen: 30% (3/10) [0s]
#>  NC validity screen: 40% (4/10) [0s]
#>  NC validity screen: 50% (5/10) [0s]
#>  NC validity screen: 60% (6/10) [0s]
#>  NC validity screen: 70% (7/10) [0s]
#>  NC validity screen: 80% (8/10) [0s]
#>  NC validity screen: 90% (9/10) [0s]
#>  NC validity screen: 100% (10/10) [0s]
#>    feature      p_value          fdr partial_r relative_effect significant
#> 1        1 5.569730e-73 1.113946e-72 0.6946748              NA        TRUE
#> 2        2 3.163889e-71 3.954861e-71 0.6885451              NA        TRUE
#> 3        3 1.420593e-71 2.029418e-71 0.6897724              NA        TRUE
#> 4        4 7.347684e-75 1.836921e-74 0.7010751              NA        TRUE
#> 5        5 3.236826e-75 1.370757e-74 0.7022684              NA        TRUE
#> 6        6 7.438239e-71 7.438239e-71 0.6872282              NA        TRUE
#> 7        7 4.112272e-75 1.370757e-74 0.7019206              NA        TRUE
#> 8        8 7.555540e-77 7.555540e-76 0.7076628              NA        TRUE
#> 9        9 4.484837e-71 4.983152e-71 0.6880085              NA        TRUE
#> 10      10 1.929682e-72 3.216137e-72 0.6928056              NA        TRUE
#>                    verdict                   verdict_fdr
#> 1  drop: associated with X drop: associated with X (FDR)
#> 2  drop: associated with X drop: associated with X (FDR)
#> 3  drop: associated with X drop: associated with X (FDR)
#> 4  drop: associated with X drop: associated with X (FDR)
#> 5  drop: associated with X drop: associated with X (FDR)
#> 6  drop: associated with X drop: associated with X (FDR)
#> 7  drop: associated with X drop: associated with X (FDR)
#> 8  drop: associated with X drop: associated with X (FDR)
#> 9  drop: associated with X drop: associated with X (FDR)
#> 10 drop: associated with X drop: associated with X (FDR)
#>                      verdict_magnitude
#> 1  drop: associated with X (magnitude)
#> 2  drop: associated with X (magnitude)
#> 3  drop: associated with X (magnitude)
#> 4  drop: associated with X (magnitude)
#> 5  drop: associated with X (magnitude)
#> 6  drop: associated with X (magnitude)
#> 7  drop: associated with X (magnitude)
#> 8  drop: associated with X (magnitude)
#> 9  drop: associated with X (magnitude)
#> 10 drop: associated with X (magnitude)
```
