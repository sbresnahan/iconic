# Print method for iconic_diagnosis objects

Print method for iconic_diagnosis objects

## Usage

``` r
# S3 method for class 'iconic_diagnosis'
print(x, ...)
```

## Arguments

- x:

  An `iconic_diagnosis` object.

- ...:

  Unused.

## Value

Invisibly returns `x` (the `iconic_diagnosis` object); called for its
side effect of printing a human-readable summary.

## Examples

``` r
set.seed(1)
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
#>  Estimating features: 50% (5/10) [0s]
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
#>  omega_1 NC coverage: 60% (6/10) [0s]
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
#>  omega_2 NC coverage: 60% (6/10) [0s]
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
#>  NC capture null: 20% (40/200) [0.7s]
#>  NC capture null: 30% (60/200) [1.1s]
#>  NC capture null: 40% (80/200) [1.4s]
#>  NC capture null: 50% (100/200) [1.8s]
#>  NC capture null: 60% (120/200) [2.1s]
#>  NC capture null: 70% (140/200) [2.5s]
#>  NC capture null: 80% (160/200) [2.8s]
#>  NC capture null: 90% (180/200) [3.2s]
#>  NC capture null: 100% (200/200) [3.5s]
#> iconic_diagnose complete. Call summary() or print() on the result for the full diagnosis.
print(diag)
#> <iconic_diagnosis>
#> Diagnostic summary:
#>  G (exposure instrument): partial F = 0.7 (WEAK)
#>  NC validity (A1): 10/10 controls valid (0 flagged)
#>  NC independence (A2): 10/10 controls valid (0 flagged)
#>  Completeness: 10 valid NCs vs k=3 (inferred) -> satisfied
#>   Capture: incremental R^2 = 0.103 (p = 0.455) -> weak
#>   Support: R^2(U~|W) = 0.097 -> narrow coverage (0/10 controls add unique coverage)
#> 
#>  Eligible estimators: 3/8
#>   UNADJ, DIRECT, COCA 
```
