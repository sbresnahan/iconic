# Resolve the trained GAN for a simulation function (internal)

Priority: explicit `trained_gan` argument \> `data$trained_gan`
(attached at iconic_data() construction) \> auto-train from `data`.

## Usage

``` r
.resolve_gan(trained_gan, data, epochs = 100)
```

## Arguments

- trained_gan:

  Explicit argument (may be NULL).

- data:

  An `iconic_data` object.

- epochs:

  Epochs for auto-training. Default 100.

## Value

A list with `gan` (the resolved iconic_gan) and `source` (a string:
"user-supplied GAN", "data-attached GAN", or "auto-trained from data").
