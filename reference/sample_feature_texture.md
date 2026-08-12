# Draw synthetic feature vectors from a trained feature texture model

Samples from the Gaussian copula: draws n_samples from MVN(0,
copula_cor), then transforms each feature through its inverse marginal
CDF. The result is centered and scaled to zero mean and unit variance
per feature, so the texture acts as noise that the structural signal in
[`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
provides the mean for.

## Usage

``` r
sample_feature_texture(feature_texture, n_samples, n_features = NULL)
```

## Arguments

- feature_texture:

  An `iconic_feature_texture` object from
  [`train_feature_texture()`](https://seantbresnahan.com/iconic/reference/train_feature_texture.md).

- n_samples:

  Number of synthetic samples to draw.

- n_features:

  Target number of features. If NULL, uses the number of features in the
  training data. If larger, additional features are drawn by sampling
  existing columns with replacement and adding independent noise. If
  smaller, the first n_features are used.

## Value

A `n_features x n_samples` matrix of synthetic texture values, centered
and scaled per feature.

## Examples

``` r
M <- matrix(rnorm(30 * 200), 30, 200) # 30 transcripts, 200 samples
ft <- train_feature_texture(M)
draws <- sample_feature_texture(ft, 500, n_features = 20)
```
