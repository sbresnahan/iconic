# Build an exposure polygenic score with LDpred2-auto

Runs the LDpred2-auto Bayesian shrinkage pipeline (Privé et al.) via the
bigsnpr package: match summary statistics to an LD reference, estimate
SNP heritability by LD score regression, run a grid of auto chains,
filter chains by convergence, and average the per-chain effects.
Optionally scores a target genotype set.

## Usage

``` r
build_prs_ldpred2(
  sumstats,
  ld_ref,
  genotypes = NULL,
  hapmap3 = NULL,
  max_variants = NULL,
  n_chains = 30,
  p_range = c(1e-04, 0.2),
  burn_in = 500,
  num_iter = 500,
  shrink_corr = 0.95,
  allow_jump_sign = FALSE,
  use_mle = FALSE,
  chain_filter_quantile = 0.95,
  impute_mean = TRUE,
  ncores = 1,
  verbose = TRUE
)
```

## Arguments

- sumstats:

  A data.frame of GWAS summary statistics, ideally pre-processed by
  [`qc_gwas_sumstats()`](https://seantbresnahan.com/iconic/reference/qc_gwas_sumstats.md).
  Must contain (after standardization) `chr`, `pos`, `a0`, `a1`, `beta`,
  `beta_se`, `n_eff`; `freq` is used for the SD-ratio QC when present.

- ld_ref:

  Either a path to a bigsnpr `.rds` backing file (from
  [`bigsnpr::snp_readBed()`](https://privefl.github.io/bigsnpr/reference/snp_readBed.html))
  or an attached `bigSNP` object. Should match the ancestry of the GWAS.

- genotypes:

  Optional: path to a bigsnpr `.rds` or an attached `bigSNP` object for
  the target cohort to be scored. When `NULL` (default), only the
  LDpred2 effect sizes are returned.

- hapmap3:

  Optional character vector of variant IDs (rsIDs or `chr:pos:a0:a1`
  keys) to restrict the analysis to (e.g. HapMap3). Restriction improves
  chain stability; `NULL` (default) uses all matched variants.

- max_variants:

  Optional integer: if the matched variant count exceeds this,
  systematically thin to approximately `max_variants` variants (evenly
  spaced along the genome) before computing LD. Default `NULL` (no
  thinning). Use e.g. `250000` to bound memory.

- n_chains:

  Integer: number of LDpred2-auto chains. Default 30.

- p_range:

  Numeric length-2: range of the initial per-chain
  proportion-of-causal-variants grid, log-spaced. Default
  `c(1e-4, 0.2)`.

- burn_in:

  Integer: burn-in iterations per chain. Default 500.

- num_iter:

  Integer: sampling iterations per chain (after burn-in). Default 500.

- shrink_corr:

  Numeric: LD shrinkage parameter. Default 0.95; reduce to ~0.4 when the
  SD-ratio QC flags \>10% of variants.

- allow_jump_sign:

  Logical: allow sign jumps across chains. Default `FALSE` (recommended
  for auto mode).

- use_mle:

  Logical: use MLE for the hyperparameter updates. Default `FALSE`
  (appropriate for moderate GWAS sample sizes).

- chain_filter_quantile:

  Numeric: chains are kept when the range of their posterior `corr_est`
  exceeds `chain_filter_quantile`-quantile of all chain ranges times
  0.95. Default 0.95.

- impute_mean:

  Logical: mean-impute missing target genotypes before scoring. Default
  `TRUE`.

- ncores:

  Integer: threads for bigsnpr. Default 1.

- verbose:

  Logical: print progress. Default `TRUE`.

## Value

A list with:

- beta:

  data.frame of LDpred2 effect sizes (`chr`, `pos`, `a0`, `a1`, `beta`)
  in LD-reference orientation.

- score:

  numeric vector of per-sample polygenic scores for `genotypes` (scaled
  to unit variance), or `NULL` when `genotypes` was not supplied.

- h2:

  LD-score-regression SNP heritability estimate used to initialize the
  chains.

- chains_kept:

  integer indices of chains passing the convergence filter.

- qc:

  list with variant counts at each stage.

## Details

This function requires the bigsnpr package (listed under `Suggests`) and
a PLINK/bed-format LD reference panel converted with
[`bigsnpr::snp_readBed()`](https://privefl.github.io/bigsnpr/reference/snp_readBed.html).

## Examples

``` r
if (requireNamespace("bigsnpr", quietly = TRUE)) {
  # Tiny fake genotype panel standing in for an LD reference
  # (snp_fake genotypes are all-missing mocks; fill with random dosages)
  fake <- bigsnpr::snp_fake(100, 500)
  fake$genotypes[] <- rbinom(100 * 500, 2, 0.3)
  fake$map$chromosome <- 1L
  fake$map$physical.pos <- sort(sample(1:1e6, 500))
  ss <- data.frame(
    chr = fake$map$chromosome, pos = fake$map$physical.pos,
    a0 = fake$map$allele2, a1 = fake$map$allele1,
    beta = rnorm(500, 0, 0.05), beta_se = 0.02, n_eff = 20000
  )
  prs <- build_prs_ldpred2(ss, ld_ref = fake, n_chains = 3,
                           burn_in = 50, num_iter = 50, verbose = FALSE)
  str(prs$beta)
}
#> 500 variants to be matched.
#> 0 ambiguous SNPs have been removed.
#> 500 variants have been matched; 0 were flipped and 0 were reversed.
#> Loading required package: foreach
#> Loading required package: rngtools
#> Warning: No LDpred2-auto chains passed the convergence filter; returning the average across all chains. Inspect chain diagnostics before trusting these effects (consider HapMap3 restriction, larger burn_in/num_iter, or an ancestry-matched LD reference).
#> 'data.frame':    500 obs. of  5 variables:
#>  $ chr : chr  "chr1" "chr1" "chr1" "chr1" ...
#>  $ pos : int  311 420 5032 6090 6833 8070 8580 9729 9991 10298 ...
#>  $ a0  : chr  "C" "C" "C" "T" ...
#>  $ a1  : chr  "T" "T" "T" "C" ...
#>  $ beta: num  NA NA NA NA NA NA NA NA NA NA ...
```
