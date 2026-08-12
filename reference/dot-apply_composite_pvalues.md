# Apply composite null p-values to a mediation results data frame (internal)

Two-pass post-processing for `se_method = "composite"`.

## Usage

``` r
.apply_composite_pvalues(res)
```

## Arguments

- res:

  Data frame from
  [`.estimate_mediation_feature()`](https://seantbresnahan.com/iconic/reference/dot-estimate_mediation_feature.md)
  / `.estimate_mediation_driver()` /
  [`run_mediation_methods()`](https://seantbresnahan.com/iconic/reference/run_mediation_methods.md).
  Must contain columns: `method`, `NIE`, `NIE_p`, `alpha_M`, `alpha_se`,
  `beta_M`, `beta_M_se`.

## Value

The same data frame with `NIE_p` replaced by composite p-values. Also
adds columns `var_a`, `var_b` for transparency.

## Details

The JT-comp test (Huang 2019) requires Var(a) and Var(b), the variances
of the standardized z-statistics `a = alpha_M / alpha_se` and
`b = beta_M / beta_M_se` across the collection of tests sharing the same
estimator. This function:

1.  Groups rows by `method` (and `mediator` when present).

2.  For each group, computes the z-statistics from the delta-method
    estimates already stored in the data frame.

3.  Estimates Var(a) / Var(b) as the sample variance of the
    z-statistics. When fewer than 5 tests are available in a group,
    falls back to Var = 1 (the point-null value, conservative).

4.  Replaces `NIE_p` with the JT-comp composite p-value.

NDE_p is left unchanged (the NDE is a single coefficient, not a product,
so the composite null does not apply). NIE_se is left as the
delta-method SE (used for CI construction); only the p-value is
replaced.
