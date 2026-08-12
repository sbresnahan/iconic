# Build the GAN generator network (internal to the torch path)

MLP mapping a noise vector to a synthetic row:
`noise_dim -> 256 -> 512 -> 1024 -> output_dim`, ReLU + batch-norm +
dropout, linear output.

## Usage

``` r
create_generator(output_dim, noise_dim = 100)
```

## Arguments

- output_dim:

  Number of output variables (columns of the training frame).

- noise_dim:

  Latent noise dimension. Default 100.

## Value

A torch `nn_module` instance. Requires torch.
