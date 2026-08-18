# Summary method for iconic_recommendation objects

Prints the full recommendation summary (same as
[`print()`](https://rdrr.io/r/base/print.html)).

## Usage

``` r
# S3 method for class 'iconic_recommendation'
summary(object, ...)
```

## Arguments

- object:

  An `iconic_recommendation` object.

- ...:

  Unused.

## Value

Invisibly returns `object`.

## Examples

``` r
data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
  G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
diag <- iconic_diagnose(data)
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
#> Estimating features: 10 tasks (sequential)
#>  Estimating features: 10% (1/10) [0s]
#>  Estimating features: 20% (2/10) [0s]
#>  Estimating features: 30% (3/10) [0s]
#>  Estimating features: 40% (4/10) [0s]
#>  Estimating features: 50% (5/10) [0.1s]
#>  Estimating features: 60% (6/10) [0.1s]
#>  Estimating features: 70% (7/10) [0.1s]
#>  Estimating features: 80% (8/10) [0.1s]
#>  Estimating features: 90% (9/10) [0.1s]
#>  Estimating features: 100% (10/10) [0.1s]
#> omega_1 NC coverage: 10 tasks (sequential)
#>  omega_1 NC coverage: 10% (1/10) [0s]
#>  omega_1 NC coverage: 20% (2/10) [0s]
#>  omega_1 NC coverage: 30% (3/10) [0s]
#>  omega_1 NC coverage: 40% (4/10) [0s]
#>  omega_1 NC coverage: 50% (5/10) [0s]
#>  omega_1 NC coverage: 60% (6/10) [0.1s]
#>  omega_1 NC coverage: 70% (7/10) [0.1s]
#>  omega_1 NC coverage: 80% (8/10) [0.1s]
#>  omega_1 NC coverage: 90% (9/10) [0.1s]
#>  omega_1 NC coverage: 100% (10/10) [0.1s]
#> omega_2 NC coverage: 10 tasks (sequential)
#>  omega_2 NC coverage: 10% (1/10) [0s]
#>  omega_2 NC coverage: 20% (2/10) [0s]
#>  omega_2 NC coverage: 30% (3/10) [0s]
#>  omega_2 NC coverage: 40% (4/10) [0s]
#>  omega_2 NC coverage: 50% (5/10) [0s]
#>  omega_2 NC coverage: 60% (6/10) [0.1s]
#>  omega_2 NC coverage: 70% (7/10) [0.1s]
#>  omega_2 NC coverage: 80% (8/10) [0.1s]
#>  omega_2 NC coverage: 90% (9/10) [0.1s]
#>  omega_2 NC coverage: 100% (10/10) [0.1s]
#> k permutation analysis: 100 tasks (sequential)
#>  k permutation analysis: 10% (10/100) [0s]
#>  k permutation analysis: 20% (20/100) [0s]
#>  k permutation analysis: 30% (30/100) [0s]
#>  k permutation analysis: 40% (40/100) [0s]
#>  k permutation analysis: 50% (50/100) [0s]
#>  k permutation analysis: 60% (60/100) [0s]
#>  k permutation analysis: 70% (70/100) [0s]
#>  k permutation analysis: 80% (80/100) [0s]
#>  k permutation analysis: 90% (90/100) [0s]
#>  k permutation analysis: 100% (100/100) [0s]
#> NC capture null: 200 tasks (sequential)
#>  NC capture null: 10% (20/200) [0.4s]
#>  NC capture null: 20% (40/200) [0.8s]
#>  NC capture null: 30% (60/200) [1.2s]
#>  NC capture null: 40% (80/200) [1.6s]
#>  NC capture null: 50% (100/200) [2s]
#>  NC capture null: 60% (120/200) [2.4s]
#>  NC capture null: 70% (140/200) [2.8s]
#>  NC capture null: 80% (160/200) [3.2s]
#>  NC capture null: 90% (180/200) [3.6s]
#>  NC capture null: 100% (200/200) [4s]
#> iconic_diagnose complete. Call summary() or print() on the result for the full diagnosis.
est <- iconic_estimate(data, diagnosis = diag)
#> Estimating features: 10 tasks (sequential)
#>  Estimating features: 10% (1/10) [0s]
#>  Estimating features: 20% (2/10) [0s]
#>  Estimating features: 30% (3/10) [0s]
#>  Estimating features: 40% (4/10) [0s]
#>  Estimating features: 50% (5/10) [0s]
#>  Estimating features: 60% (6/10) [0s]
#>  Estimating features: 70% (7/10) [0s]
#>  Estimating features: 80% (8/10) [0s]
#>  Estimating features: 90% (9/10) [0s]
#>  Estimating features: 100% (10/10) [0s]
rec <- iconic_recommend(data, diagnosis = diag, estimate = est,
  auto_sensitivity = FALSE)
summary(rec)
#> <iconic_recommendation>
#>  Recommended (composite): COCA
#>  Estimated total effect (mean): 0.0901
#> 
#>  Eligible alternatives: UNADJ, DIRECT 
#> 
#>  Full ranking:
#>  1. COCA -- requires: valid NCs (A1), completeness (A2 not required)
#>  2. UNADJ -- requires: no assumptions (naive OLS)
#>  3. DIRECT -- requires: G + W as covariates (no causal identification)
#> 
#>  Ineligible:
#>  IV2SLS -- ineligible -- requires G + F_G>=10 (F_G=0.7)
#>  PGC -- ineligible -- requires G + W + F_G>=10 + completeness (completeness: satisfied)
#>  IV2SLS2 -- ineligible -- requires mediation data (supply M)
#>  PGC2 -- ineligible -- requires mediation data (supply M)
#>  PGC2Gm -- ineligible -- requires mediation data (supply M)
```
