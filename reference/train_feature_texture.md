# Train a feature-level texture model for the mediator panel

Learns the joint distribution of a features x samples mediator matrix
using a Gaussian copula with flexible marginal distributions. For each
feature, the marginal is fitted as either an empirical CDF (default) or
the best parametric family (normal, log-normal, gamma, beta) selected by
AIC with a Kolmogorov-Smirnov goodness-of-fit check. Cross-feature
dependence is captured by the Gaussian copula correlation matrix.

## Usage

``` r
train_feature_texture(M_matrix, marginal_method = "auto")
```

## Arguments

- M_matrix:

  Mediator panel, `features x samples` (e.g. gene or transcript
  expression, one row per feature, one column per sample).

- marginal_method:

  Marginal fitting method: `"auto"` (default, uses parametric if KS test
  passes at p \> 0.05, otherwise empirical), `"empirical"` (always use
  empirical CDF), or `"parametric"` (always use best parametric fit by
  AIC).

## Value

An `iconic_feature_texture` S3 object: a named list with `marginals`
(list of per-feature marginal specs), `copula_cor` (p x p correlation
matrix), `n_features`, `n_samples`, `marginal_method`, and
`marginal_types` (summary of how many features used each method).

## Details

The resulting model can be sampled via
[`sample_feature_texture()`](https://seantbresnahan.com/iconic/reference/sample_feature_texture.md)
to produce realistic synthetic feature vectors that preserve the
marginal shapes and correlation structure of the user's mediator panel.

## Examples

``` r
M <- matrix(rnorm(30 * 200), 30, 200) # 30 transcripts, 200 samples
ft <- train_feature_texture(M)
draws <- sample_feature_texture(ft, 100) # 30 x 100 synthetic draws
```
