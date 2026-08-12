# Standardize and QC GWAS summary statistics

Renames common GWAS summary-statistic column conventions to the bigsnpr
convention (`chr`, `pos`, `a0`, `a1`, `beta`, `beta_se`, `n_eff`,
`freq`) and applies the standard LDpred2 QC filters: missing/invalid
values, extreme effect sizes, ambiguous-strand (A/T and C/G) SNPs, and —
when effect allele frequencies are available — the standard-deviation
ratio check from the LDpred2 tutorial.

## Usage

``` r
qc_gwas_sumstats(
  sumstats,
  n_eff = NULL,
  drop_ambiguous = TRUE,
  beta_iqr_mult = 10,
  sd_ratio_range = c(0.5, 2),
  sd_y = NULL
)
```

## Arguments

- sumstats:

  A data.frame of GWAS summary statistics. Recognized column aliases:
  `chr`/`chromosome`/`CHR`; `pos`/`base_pair_location`/ `BP`;
  `a1`/`effect_allele`/`ALT`; `a0`/`other_allele`/`REF`;
  `beta`/`BETA`/`effect`/`logOR`; `beta_se`/`standard_error`/`se`/`SE`;
  `freq`/`effect_allele_frequency`/`EAF`/`FRQ`; `n_eff`/`n`/`N`;
  `p`/`p_value`/`P`; `rsid`/`rs_id`/`SNP`/`ID`.

- n_eff:

  Optional scalar: effective sample size, used when the input has no
  sample-size column.

- drop_ambiguous:

  Logical: drop A/T and C/G SNPs whose strand cannot be resolved.
  Default `TRUE`.

- beta_iqr_mult:

  Numeric: drop SNPs with \\\|\beta - \mathrm{median}(\beta)\| \>\\
  `beta_iqr_mult` \\\times \mathrm{IQR}(\beta)\\. Default 10. Set to
  `Inf` to skip.

- sd_ratio_range:

  Numeric length-2: acceptable range for the ratio of the
  summary-statistic-implied genotype SD to the frequency-implied SD,
  \\\sqrt{2 f (1-f)}\\. SNPs outside the range are dropped. Default
  `c(0.5, 2)` (LDpred2 tutorial). Requires `freq`; skipped with a
  message when frequencies are absent.

- sd_y:

  Optional scalar: phenotypic SD used in the SD-ratio check. When `NULL`
  (default), estimated as the 1st percentile of \\\sqrt{0.5 (n\_{eff}
  \cdot se^2 + \beta^2)}\\ across SNPs, the continuous-trait
  approximation used in the LDpred2 tutorial.

## Value

A list with:

- sumstats:

  data.frame of QC-passing variants with standardized columns `chr`,
  `pos`, `a0`, `a1`, `beta`, `beta_se`, and (when available) `n_eff`,
  `freq`, `p`, `rsid`.

- qc:

  data.frame with counts of input variants and variants dropped by each
  filter.

- sd_ratio:

  numeric vector of SD ratios for retained variants, or `NULL` when the
  check was skipped.

## Examples

``` r
ss <- data.frame(
  chromosome = 1, base_pair_location = 1000 + 0:9 * 1000L,
  effect_allele = rep(c("A", "C"), 5), other_allele = rep(c("G", "T"), 5),
  beta = rnorm(10, 0, 0.05), standard_error = 0.01,
  effect_allele_frequency = runif(10, 0.1, 0.5)
)
out <- qc_gwas_sumstats(ss, n_eff = 50000)
str(out$sumstats)
#> 'data.frame':    10 obs. of  8 variables:
#>  $ chr    : chr  "1" "1" "1" "1" ...
#>  $ pos    : int  1000 2000 3000 4000 5000 6000 7000 8000 9000 10000
#>  $ a1     : chr  "A" "C" "A" "C" ...
#>  $ a0     : chr  "G" "T" "G" "T" ...
#>  $ beta   : num  0.0253 0.0591 -0.0423 -0.0497 -0.0457 ...
#>  $ beta_se: num  0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01
#>  $ freq   : num  0.211 0.429 0.246 0.388 0.185 ...
#>  $ n_eff  : num  50000 50000 50000 50000 50000 50000 50000 50000 50000 50000
out$qc
#>   n_input dropped_missing dropped_extreme_beta dropped_ambiguous
#> 1      10               0                    0                 0
#>   dropped_sd_ratio n_output
#> 1                0       10
```
