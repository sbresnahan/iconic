# Print method for iconic_gan objects

Print method for iconic_gan objects

## Usage

``` r
# S3 method for class 'iconic_gan'
print(x, ...)
```

## Arguments

- x:

  An `iconic_gan` object.

- ...:

  Unused.

## Value

Invisibly returns `x` (the `iconic_gan` object); called for its side
effect of printing a human-readable summary.

## Examples

``` r
if (check_torch_setup()) {
  dat <- load_real_input_data(example = TRUE)
  gan <- train_gan_on_real_data(dat$gan_training_data, epochs = 5,
    verbose = FALSE)
  print(gan)
}
```
