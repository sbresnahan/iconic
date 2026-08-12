# Simulate a single genetic instrument

Draws one standardised genetic instrument for the exposure, as an
additive score over `n_snps` biallelic SNPs (dosages
`Binomial(2, maf)`). Mirrors the single-instrument (PRS) design used
throughout the SCENIC framework.

## Usage

``` r
simulate_single_genetic_instrument(n, n_snps = 20, maf = 0.3, seed = NULL)
```

## Arguments

- n:

  Number of samples.

- n_snps:

  SNPs contributing to the score. Default 20.

- maf:

  Minor-allele frequency. Default 0.3.

- seed:

  Optional RNG seed.

## Value

List with `G` (standardised instrument, length `n`), `dosages`
(`n x n_snps`), and `maf`.

## Examples

``` r
g <- simulate_single_genetic_instrument(100, n_snps = 10, seed = 1)
length(g$G)
#> [1] 100
```
