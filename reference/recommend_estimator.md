# Recommend the preferred estimator from a sensitivity sweep

Ranks methods by RMSE within each scenario, and identifies an overall
pick that is robust across scenarios (smallest worst-case RMSE). UNADJ
is excluded from recommendations (it is a bias reference).

## Usage

``` r
recommend_estimator(sens, exclude = c("UNADJ"))
```

## Arguments

- sens:

  Object returned by
  [`gan_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_sensitivity.md).

- exclude:

  Methods to exclude from the recommendation. Default `c("UNADJ")`.

## Value

A list with `per_scenario` (best method + RMSE per scenario),
`worst_case` (max RMSE per method across scenarios), and `overall` (the
method minimising worst-case RMSE).

## Examples

``` r
sens <- gan_sensitivity(NULL, conf_grid = 0.8, coverage_grid = 0.7,
  n_iter = 2, n_samples = 100, n_features = 5)
recommend_estimator(sens)
#> $per_scenario
#>           conf_strength coverage k best_method       rmse
#> 0.8.0.7.1           0.8      0.7 1         PGC 0.09266774
#> 
#> $worst_case
#>   method worst_rmse
#> 1   COCA 0.89111120
#> 2 DIRECT 0.16157144
#> 3 IV2SLS 0.09738335
#> 4    PGC 0.09266774
#> 
#> $overall
#> [1] "PGC"
#> 
```
