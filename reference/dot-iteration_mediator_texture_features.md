# Exogenous, mediator-only texture component for one iteration (internal)

Returns a p x n matrix of realistic mediator feature vectors drawn from
the feature-level copula texture model. Each row is a feature
(transcript), each column a sample. The draws preserve the marginal
distributions and cross-feature correlation structure of the user's
mediator panel. The texture is drawn independently of U and never fed
into X, so it cannot open a confounding backdoor.

## Usage

``` r
.iteration_mediator_texture_features(trained_gan, n, p)
```

## Arguments

- trained_gan:

  An iconic_gan object (may be NULL).

- n:

  Number of synthetic samples.

- p:

  Number of features (target n_features for this iteration).

## Value

A p x n matrix, or NULL.

## Details

When no feature_texture is available (trained_gan is NULL or lacks
\$feature_texture), returns NULL and the caller falls back to scalar
noise.
