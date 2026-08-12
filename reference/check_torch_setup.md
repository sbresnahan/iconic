# Check whether a working torch installation is available

Verifies the package is installed, libtorch is present, and a basic
tensor op succeeds. Used to gate the GAN path; a `FALSE` result triggers
the multivariate-normal fallback.

## Usage

``` r
check_torch_setup()
```

## Value

`TRUE` if torch is usable, otherwise `FALSE`.

## Examples

``` r
check_torch_setup()
#> [1] FALSE
```
