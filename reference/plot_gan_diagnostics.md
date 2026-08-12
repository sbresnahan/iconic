# Compare real vs synthetic marginals from a trained generator

For continuous variables, overlays kernel densities of the real training
variables against draws from the trained generator. For binary (0/1)
variables – detected at training time and stored in the `iconic_gan`
object – uses side-by-side bar charts of the proportion of 1s, which is
the appropriate visualisation for discrete data. One panel per variable.

## Usage

``` r
plot_gan_diagnostics(trained_gan, real_data, n_draw = NULL, M_matrix = NULL)
```

## Arguments

- trained_gan:

  An `iconic_gan` from
  [`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md).

- real_data:

  The data frame the generator was trained on.

- n_draw:

  Synthetic rows to draw for the overlay. Default: number of real rows.

- M_matrix:

  Optional mediator panel (features x samples) for feature-level
  marginal comparison. When supplied and the GAN carries a
  feature_texture, an additional panel of density overlays compares real
  vs copula-sampled mediator marginals.

## Value

Called for its side effects; invisibly returns `NULL`.

## Examples

``` r
if (check_torch_setup()) {
  dat <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(dat$gan_training_data, epochs = 5,
    verbose = FALSE)
  plot_gan_diagnostics(gan, dat$gan_training_data)
}
```
