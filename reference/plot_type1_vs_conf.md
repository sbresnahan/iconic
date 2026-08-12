# Type I error rate vs confounding strength

Type I error rate vs confounding strength

## Usage

``` r
plot_type1_vs_conf(
  t1e_df,
  methods = iconic_method_order,
  alpha = 0.05,
  title = "Type I Error Rate vs Confounding Strength"
)
```

## Arguments

- t1e_df:

  Data frame returned by sweep_null_by_conf().

- methods:

  Methods to plot. Default: all five.

- alpha:

  Nominal significance level. Default 0.05.

- title:

  Plot title.

## Value

Invisibly returns `NULL`; called for the side effect of drawing Type I
error vs confounding strength.

## Examples

``` r
t1e <- sweep_null_by_conf(c(0.5, 0.8), n_iter = 2, n_samples = 100,
  n_features = 5)
plot_type1_vs_conf(t1e)
```
