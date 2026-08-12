# Bar chart of Type I error rates

Bar chart of Type I error rates

## Usage

``` r
plot_type1_error(null_df, title = "Type I Error Rate by Method")
```

## Arguments

- null_df:

  Data frame: either run_null_sim()\$rates or a data frame with columns
  method and type1_error.

- title:

  Plot title.

## Value

Invisibly returns `NULL`; called for the side effect of plotting Type I
error rates by method.

## Examples

``` r
t1e <- sweep_null_by_conf(0.8, n_iter = 2, n_samples = 100,
  n_features = 5)
plot_type1_error(t1e)
```
