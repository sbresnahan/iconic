# Train a generative texture model on real data

Fits an adversarial network (torch GAN) to the tidy training frame.
Requires the torch package; if torch is not available, an error is
raised with installation instructions.

## Usage

``` r
train_gan_on_real_data(
  real_data,
  feature_correlations = NULL,
  feature_texture = NULL,
  epochs = 300,
  batch_size = 32,
  lr = 2e-04,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- real_data:

  Data frame of numeric variables (e.g. the `gan_training_data` element
  from
  [`load_real_input_data()`](https://seantbresnahan.com/iconic/reference/load_real_input_data.md)).
  Complete cases are used.

- feature_correlations:

  Optional list with elements `Y`, `W`, each a p x p residual
  correlation matrix from
  [`load_real_input_data()`](https://seantbresnahan.com/iconic/reference/load_real_input_data.md).
  When supplied, the matrices are stored in the returned `iconic_gan`
  object and used by
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  to generate correlated noise. Default `NULL`.

- feature_texture:

  Optional `iconic_feature_texture` object from
  [`train_feature_texture()`](https://seantbresnahan.com/iconic/reference/train_feature_texture.md).
  When supplied, it is stored in the returned `iconic_gan` object and
  used by
  [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  to inject realistic mediator texture. Default `NULL`.

- epochs:

  GAN training epochs. Default 300.

- batch_size:

  Mini-batch size. Default 32.

- lr:

  Adam learning rate. Default 2e-4.

- seed:

  Optional RNG seed (sets both R and torch seeds).

- verbose:

  Print progress. Default TRUE.

## Value

An `iconic_gan` object: a list with `model_type` (`"gan"`), `columns`,
`norm` (per-column centre/scale), `binary_cols` (names of 0/1 columns),
`onehot_groups` (list of one-hot dummy groups), `feature_correlations`
(list or NULL), `feature_texture` (object or NULL), training statistics,
and the trained networks + loss history.

## Details

Binary columns (values in \\{0, 1}\\, e.g. encoded sex or one-hot
ethnicity dummies) are detected automatically and stored in the returned
object.
[`sample_texture()`](https://seantbresnahan.com/iconic/reference/sample_texture.md)
rounds them back to 0/1 and enforces one-hot mutual exclusivity, so
synthetic draws respect the discrete structure.

Feature-level residual correlation matrices (for the Y and W panels) can
be attached via the `feature_correlations` argument. When present,
[`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
uses them to inject correlated noise into the simulated outcome and
negative-control panels.

A feature-level copula texture model for the mediator (M) panel can be
attached via the `feature_texture` argument (an `iconic_feature_texture`
object from
[`train_feature_texture()`](https://seantbresnahan.com/iconic/reference/train_feature_texture.md)).
When present,
[`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
uses it to draw realistic mediator feature vectors that preserve the
marginal distributions and cross-feature correlation structure of the
user's mediator panel.

## Examples

``` r
if (check_torch_setup()) {
  dat <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(dat$gan_training_data,
    feature_correlations = dat$feature_correlations,
    feature_texture = dat$feature_texture,
    epochs = 5, verbose = FALSE)
  head(sample_texture(gan, 5))
}
```
