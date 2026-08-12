# Negative-control support/range check

Diagnostic for whether the negative-control panel captures the full
support of the confounder, or only part of it. Uses a confounder proxy
`U_tilde = resid(X ~ G + C)` (the instrument-purged exposure residual,
which carries the confounding signal) and asks how much of `U_tilde` the
NC panel explains, and whether each individual NC covers a distinct
share of that signal.

## Usage

``` r
nc_support_check(dat, fdr_level = 0.1)
```

## Arguments

- dat:

  Dataset list from
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md),
  [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md),
  or the `.to_nc_dat()` bridge in
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md),
  containing `W`, `X`, and `synthetic_data` (covariates). `G` is used
  when present.

- fdr_level:

  FDR level for flagging controls whose unique contribution is
  indistinguishable from zero. Default 0.10.

## Value

A list with: `R2_utilde_given_W` (multivariate R^2 of the confounder
proxy on the full NC panel), `support` (data frame with per-NC
`support_ratio` = partial correlation of the control with `U_tilde`
given the other controls, `p_value`, and `adds_coverage` flag),
`n_controls`, and `verdict` ("broad" when `R2_utilde_given_W >= 0.5`,
"partial" when \>= 0.2, else "narrow").

## Details

A panel can pass the count-based completeness check
(`dim(W_valid) >= k`) and the covariance-capture test while still
covering only part of the confounder support — for example when every
control loads on the same single confounder direction. This diagnostic
reports the multivariate `R^2(U_tilde | W)` (how much of the confounder
proxy the panel explains) and a per-NC `support_ratio` (the partial
correlation of each control with `U_tilde` given the other controls),
which flags controls that add no unique coverage.

This is a diagnostic, not a gate: `U_tilde` is itself an imperfect proxy
(it mixes the confounder with exposure noise), so the values are
interpretable only comparatively across controls and panels.

## Examples

``` r
dat <- run_single_iteration(n_features = 10, n_confounders = 1, seed = 1)
nc_support_check(dat)
#> $R2_utilde_given_W
#> [1] 0.7298259
#> 
#> $support
#>    control support_ratio    unique_R2    p_value     p_adj adds_coverage
#> 1        1   0.019987145 3.994860e-04 0.39653132 0.6072618         FALSE
#> 2        2   0.015904803 2.529628e-04 0.49983121 0.6247890         FALSE
#> 3        3   0.031000982 9.610609e-04 0.18872998 0.3774600         FALSE
#> 4        4   0.053481583 2.860280e-03 0.02360514 0.1733597         FALSE
#> 5        5   0.018802935 3.535504e-04 0.42508328 0.6072618         FALSE
#> 6        6   0.002752147 7.574314e-06 0.90703018 0.9070302         FALSE
#> 7        7   0.034561783 1.194517e-03 0.14292175 0.3669097         FALSE
#> 8        8   0.049889052 2.488918e-03 0.03467195 0.1733597         FALSE
#> 9        9   0.006411452 4.110671e-05 0.78557838 0.8728649         FALSE
#> 10      10   0.034232210 1.171844e-03 0.14676388 0.3669097         FALSE
#> 
#> $n_controls
#> [1] 10
#> 
#> $verdict
#> [1] "broad"
#> 
```
