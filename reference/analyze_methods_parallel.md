# Parallel version of [`analyze_methods_robust()`](https://seantbresnahan.com/iconic/reference/analyze_methods_robust.md)

Distributes the per-feature analysis across workers via
[`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html)
(falling back to sequential on Windows / one core).

## Usage

``` r
analyze_methods_parallel(
  iteration_data,
  test_features = NULL,
  alpha = 0.05,
  debug = FALSE,
  n_cores = 1
)
```

## Arguments

- iteration_data:

  Dataset list from
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  (or
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md)).

- test_features:

  Optional integer indices of outcome features to test. Default `NULL`
  (all features).

- alpha:

  Significance threshold for the `significant` flag. Default 0.05.

- debug:

  If `TRUE`, message per-feature progress. Default `FALSE`.

- n_cores:

  Number of workers. Default 1.

## Value

Data frame: `feature`, `method`, `beta`, `se`, `pvalue`, `significant`.

## Examples

``` r
dat <- run_single_iteration(NULL, n_synthetic_samples = 100,
  n_features = 5, n_confounders = 1, seed = 1)
analyze_methods_parallel(dat, n_cores = 1)
#> analyze_methods_parallel: 5 tasks (sequential)
#>  analyze_methods_parallel: 20% (1/5) [0.7s]
#>  analyze_methods_parallel: 40% (2/5) [0.7s]
#>  analyze_methods_parallel: 60% (3/5) [0.8s]
#>  analyze_methods_parallel: 80% (4/5) [0.8s]
#>  analyze_methods_parallel: 100% (5/5) [0.8s]
#>         feature method       beta         se       pvalue significant
#> UNADJ         1  UNADJ  0.5261014 0.05357861 2.987041e-16        TRUE
#> DIRECT        1 DIRECT  0.3610414 0.10312635 7.234765e-04        TRUE
#> COCA          1   COCA -0.5410075 0.22242994 1.500493e-02        TRUE
#> IV2SLS        1 IV2SLS  0.2176178 0.09301922 2.149797e-02        TRUE
#> PGC           1    PGC  0.2755532 0.06758716 9.474808e-05        TRUE
#> UNADJ1        2  UNADJ  0.5302118 0.05187722 4.022556e-17        TRUE
#> DIRECT1       2 DIRECT  0.3705756 0.09715792 2.501072e-04        TRUE
#> COCA1         2   COCA -0.4275786 0.18451456 2.048649e-02        TRUE
#> IV2SLS1       2 IV2SLS  0.2108103 0.08800189 1.864500e-02        TRUE
#> PGC1          2    PGC  0.2764705 0.06419232 4.030979e-05        TRUE
#> UNADJ2        3  UNADJ  0.4984669 0.05121875 4.617412e-16        TRUE
#> DIRECT2       3 DIRECT  0.4048192 0.09989213 1.073455e-04        TRUE
#> COCA2         3   COCA -0.5739017 0.23948376 1.655648e-02        TRUE
#> IV2SLS2       3 IV2SLS  0.1931068 0.09134674 3.725080e-02        TRUE
#> PGC2          3    PGC  0.2821838 0.06580514 4.326252e-05        TRUE
#> UNADJ3        4  UNADJ  0.4327221 0.04984945 8.773302e-14        TRUE
#> DIRECT3       4 DIRECT  0.3741117 0.10202212 4.149712e-04        TRUE
#> COCA3         4   COCA -0.9309171 0.39474764 1.836088e-02        TRUE
#> IV2SLS3       4 IV2SLS  0.1793920 0.09287439 5.652621e-02       FALSE
#> PGC3          4    PGC  0.2642272 0.06745048 1.686514e-04        TRUE
#> UNADJ4        5  UNADJ  0.5350279 0.05522612 5.758312e-16        TRUE
#> DIRECT4       5 DIRECT  0.4092637 0.10640144 2.234730e-04        TRUE
#> COCA4         5   COCA -0.5686515 0.23588898 1.592305e-02        TRUE
#> IV2SLS4       5 IV2SLS  0.2056640 0.09687073 3.646193e-02        TRUE
#> PGC4          5    PGC  0.2938509 0.06987364 5.897800e-05        TRUE
```
