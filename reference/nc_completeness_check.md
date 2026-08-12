# Check negative-control completeness (dimensional + covariance-capture)

Reports two components of the proximal-inference completeness condition:

## Usage

``` r
nc_completeness_check(
  dat,
  n_valid_controls = NULL,
  fdr_level = 0.1,
  n_cores = 1,
  outcome = "Y",
  n_perm = 1000,
  capture_thresholds = list(strong = 0.3, weak = 0.1)
)
```

## Arguments

- dat:

  Dataset list from
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  or
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md).

- n_valid_controls:

  Optional override: the number of valid controls known from the study
  design. If `NULL` (default), the function runs both empirical screens
  and counts controls that pass both.

- fdr_level:

  FDR level for the empirical screens (used only when `n_valid_controls`
  is `NULL`).

- n_cores:

  Number of parallel workers for the empirical screens. Default 1
  (sequential).

- outcome:

  Outcome block for the capture test: `"Y"` (default) or `"M"`
  (mediator). When `NULL`, the capture component is skipped and only the
  dimensional verdict is returned.

- n_perm:

  Number of permutations for the capture-test null. Default 1000.

- capture_thresholds:

  Named list with `strong` and `weak` R^2 cutoffs for the capture
  verdict. Default `list(strong = 0.3, weak = 0.1)`.

## Value

A list with: `n_valid_controls` (count), `k` (number of confounders),
`dim_W`, `dimensional` ("satisfied", "borderline", "under-identified"),
`capture` (output of
[`nc_completeness_capture()`](https://seantbresnahan.com/iconic/reference/nc_completeness_capture.md),
or NULL), `completeness` (composite: "satisfied", "borderline",
"under-identified", or "weak-capture"), `screen_X` (A1 screen results,
if run), `screen_G` (A2 screen results, if run).

## Details

\(1\) **Dimensional** (the legacy count-based check): the number of
valid negative-control features vs the number of latent confounders k.
Bridge-function estimators (COCA, PGC) require at least as many valid
controls as confounders (Miao, Geng & Tchetgen Tchetgen, 2018). When
`dim(W_valid) < k`, no estimator built on those controls can recover the
causal effect, regardless of sample size.

\(2\) **Covariance-capture**: whether the controls actually *capture the
confounder covariance* — i.e., whether adding W reduces U's contribution
to the outcome toward zero, operationalized as the incremental R^2 of W
for the outcome above covariates alone, with a permutation null. This
addresses the concern that completeness is about covariance captured,
not proxy number.

The composite `completeness` verdict requires the dimensional component
to pass (satisfied/borderline) AND the capture component to be non-
negligible (strong/weak). When the dimensional component passes but
capture is negligible, the verdict is "weak-capture".

## Examples

``` r
dat <- run_single_iteration(n_features = 10, n_confounders = 1, seed = 1)
nc_completeness_check(dat, n_perm = 50)
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
#> NC independence (G): 10 tasks (sequential)
#>  NC independence (G): 10% (1/10) [0s]
#>  NC independence (G): 20% (2/10) [0s]
#>  NC independence (G): 30% (3/10) [0s]
#>  NC independence (G): 40% (4/10) [0s]
#>  NC independence (G): 50% (5/10) [0s]
#>  NC independence (G): 60% (6/10) [0s]
#>  NC independence (G): 70% (7/10) [0s]
#>  NC independence (G): 80% (8/10) [0s]
#>  NC independence (G): 90% (9/10) [0s]
#>  NC independence (G): 100% (10/10) [0s]
#> NC capture null: 50 tasks (sequential)
#>  NC capture null: 10% (5/50) [0.1s]
#>  NC capture null: 20% (10/50) [0.2s]
#>  NC capture null: 30% (15/50) [0.3s]
#>  NC capture null: 40% (20/50) [0.4s]
#>  NC capture null: 50% (25/50) [0.5s]
#>  NC capture null: 60% (30/50) [0.7s]
#>  NC capture null: 70% (35/50) [0.8s]
#>  NC capture null: 80% (40/50) [0.9s]
#>  NC capture null: 90% (45/50) [1s]
#>  NC capture null: 100% (50/50) [1.1s]
#> $n_valid_controls
#> [1] 0
#> 
#> $k
#> [1] 1
#> 
#> $dim_W
#> [1] 10
#> 
#> $dimensional
#> [1] "under-identified"
#> 
#> $capture
#> $capture$capture_R2
#> [1] 0.616692
#> 
#> $capture$capture_pvalue
#> [1] 0
#> 
#> $capture$capture_verdict
#> [1] "strong"
#> 
#> $capture$null_distribution
#>  [1] 0.012081537 0.022761326 0.029848473 0.012926332 0.030538428 0.026354722
#>  [7] 0.014845387 0.015766020 0.021626695 0.014062661 0.016607813 0.020639402
#> [13] 0.019293768 0.018038274 0.028841584 0.021431851 0.008042623 0.008051202
#> [19] 0.021791833 0.025587591 0.023736943 0.016190043 0.024885332 0.024825202
#> [25] 0.015914368 0.009267005 0.021967634 0.033023687 0.013065615 0.013956668
#> [31] 0.025094305 0.006760867 0.015833784 0.011581334 0.007824060 0.036081752
#> [37] 0.017785892 0.009041063 0.021326898 0.023749467 0.012056004 0.013645462
#> [43] 0.021877932 0.030514910 0.014532153 0.014014346 0.024894856 0.027739813
#> [49] 0.023882658 0.028064650
#> 
#> $capture$n_features
#> [1] 10
#> 
#> 
#> $completeness
#> [1] "under-identified"
#> 
#> $screen_X
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
#> 
#> $screen_G
#>    feature    partial_r   p_value       fdr significant verdict
#> 1        1 -0.012333113 0.7834542 0.9349553       FALSE   valid
#> 2        2 -0.024678277 0.5823385 0.9349553       FALSE   valid
#> 3        3 -0.027609123 0.5383489 0.9349553       FALSE   valid
#> 4        4 -0.017619234 0.6945948 0.9349553       FALSE   valid
#> 5        5  0.003662625 0.9349553 0.9349553       FALSE   valid
#> 6        6 -0.013731240 0.7596215 0.9349553       FALSE   valid
#> 7        7 -0.009543787 0.8315904 0.9349553       FALSE   valid
#> 8        8 -0.007105246 0.8742022 0.9349553       FALSE   valid
#> 9        9 -0.021990376 0.6240949 0.9349553       FALSE   valid
#> 10      10 -0.014680248 0.7435717 0.9349553       FALSE   valid
#> 
```
