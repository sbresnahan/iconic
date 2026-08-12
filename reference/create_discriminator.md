# Build the GAN discriminator network (internal to the torch path)

MLP `input_dim -> 1024 -> 512 -> 256 -> 1` with leaky-ReLU and dropout,
returning a real-vs-fake logit (no sigmoid; paired with a
BCE-with-logits loss).

## Usage

``` r
create_discriminator(input_dim)
```

## Arguments

- input_dim:

  Number of input variables.

## Value

A torch `nn_module` instance. Requires torch.
