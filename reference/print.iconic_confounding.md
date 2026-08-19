# Print method for iconic_confounding objects

Print method for iconic_confounding objects

## Usage

``` r
# S3 method for class 'iconic_confounding'
print(x, ...)
```

## Arguments

- x:

  An `iconic_confounding` object.

- ...:

  Unused.

## Value

Invisibly returns `x` (the `iconic_confounding` object); called for its
side effect of printing a human-readable summary.

## Examples

``` r
set.seed(1)
data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
  G = rnorm(100), W = matrix(rnorm(100 * 10), 10, 100))
cf <- infer_confounding(data, max_infer_tasks = 5)
#> Estimating features: 5 tasks (sequential)
#>  Estimating features: 20% (1/5) [0s]
#>  Estimating features: 40% (2/5) [0s]
#>  Estimating features: 60% (3/5) [0s]
#>  Estimating features: 80% (4/5) [0s]
#>  Estimating features: 100% (5/5) [0s]
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
print(cf)
#> <iconic_confounding>
#>  conf_strength    default (0.8) -- weak instrument
#>  mo_confounding   default (0.8) -- no mediation estimates
#>  omega_1          0.199 (sqrt(R^2) of W on Y residualized on X+C)
#>  warning: composite: coverage x confounder strength, not pure coverage
#>  omega_2          0.199 (sqrt(R^2) of W on Y residualized on X+C)
#>  warning: composite: coverage x confounder strength, not pure coverage
#>  k                3 [CI: 2, 4] (parallel analysis (Horn, 1965))
#> 
#>  Unavailable: rho_G1, rho_G2, conf_strength, mo_confounding 
#> 
#>  Warnings:
#>   conf_strength: weak instrument (F_G=NA < 10), inference unreliable, using default 0.8. 
#>   mo_confounding: no mediation estimates available, using default 0.8. 
```
