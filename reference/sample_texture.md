# Draw synthetic base rows from a trained texture model

Binary columns (detected at training time) are rounded to \\{0, 1}\\
after de-normalisation, and one-hot dummy groups are made mutually
exclusive (the column with the highest pre-rounding value wins per row).

## Usage

``` r
sample_texture(trained_gan, n)
```

## Arguments

- trained_gan:

  An `iconic_gan` object from
  [`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md).

- n:

  Number of rows to draw.

## Value

A data frame of `n` rows with the trained columns, on the original
(de-normalised) scale. Binary columns contain only 0/1 values.

## Examples

``` r
if (check_torch_setup()) {
  dat <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(dat$gan_training_data, epochs = 5,
    verbose = FALSE)
  head(sample_texture(gan, 5))
}
```
