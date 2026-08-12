# Auto-train a texture model from an iconic_data object (internal)

Converts the user's iconic_data to the load_real_input_data() format and
trains a GAN (or multivariate-normal fallback) on the resulting tidy
frame. The trained model supplies realistic covariate and outcome
texture to run_single_iteration(), and carries feature-level residual
correlation matrices for the Y, M, and W panels so the simulation can
inject correlated noise.

## Usage

``` r
.auto_train_gan(data, epochs = 100)
```

## Arguments

- data:

  An `iconic_data` object.

- epochs:

  GAN training epochs. Default 100.

## Value

An `iconic_gan` object.
