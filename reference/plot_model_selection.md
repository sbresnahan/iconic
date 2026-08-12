# Model selection workflow figure

Produces a 3-panel figure showing the ICONIC model selection workflow on
example data. Panels: (A) eligibility table from
[`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md),
(B) NDE/NIE forest plot from
[`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md),
(C) degradation surface heatmap from
[`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md).

## Usage

``` r
plot_model_selection(
  diagnosis,
  estimate,
  sensitivity,
  recommendation,
  file = NULL,
  width = 12,
  height = 14
)
```

## Arguments

- diagnosis:

  Result of
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md).

- estimate:

  Result of
  [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).

- sensitivity:

  Result of
  [`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md).

- recommendation:

  Result of
  [`iconic_recommend()`](https://seantbresnahan.com/iconic/reference/iconic_recommend.md).

- file:

  Optional file path to save the figure.

- width:

  Figure width in inches.

- height:

  Figure height in inches.

## Value

A `patchwork` ggplot object.

## Examples

``` r
if (check_torch_setup()) {
  data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
    W = matrix(rnorm(100 * 10), 10, 100))
  diag <- iconic_diagnose(data)
  est <- iconic_estimate(data, diagnosis = diag)
  sens <- iconic_sensitivity(data, n_iter = 2, gan_epochs = 5,
    rho_G1_grid = c(0, 0.2), rho_G2_grid = c(0, 0.2))
  rec <- iconic_recommend(data, diagnosis = diag, estimate = est,
    sensitivity = sens, auto_sensitivity = FALSE)
  plot_model_selection(diag, est, sens, rec)
}
```
