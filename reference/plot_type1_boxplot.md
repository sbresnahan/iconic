# Type I error boxplot per method

Type I error boxplot per method

## Usage

``` r
plot_type1_boxplot(
  null_result,
  methods = iconic_method_order,
  conf_str = 0.8,
  alpha = 0.05
)
```

## Arguments

- null_result:

  Object returned by run_null_sim().

- methods:

  Methods to include. Default: all five.

- conf_str:

  Confounding strength used (for the title). Default 0.80.

- alpha:

  Nominal significance level. Default 0.05.

## Value

Invisibly returns `NULL`; called for the side effect of drawing a
boxplot of Type I error rates.

## Examples

``` r
null <- run_null_sim(n_iter = 2, n_samples = 100, n_features = 5)
plot_type1_boxplot(null)
```
