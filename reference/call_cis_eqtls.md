# Scan for cis-eQTLs of each gene

For each gene, tests every cis SNP (within `cis_dist` of the gene start,
optionally restricted to a supplied candidate set) for association with
expression, after residualizing both expression and dosages on
`covariates` (Frisch-Waugh-Lovell). A lightweight, pure-R alternative to
MatrixEQTL for preparing mediator-instrument candidate sets.

## Usage

``` r
call_cis_eqtls(
  expression,
  genotypes,
  gene_pos,
  snp_pos,
  covariates = NULL,
  cis_snps = NULL,
  cis_dist = 1e+06,
  min_maf = 0.01,
  fdr = 0.05,
  n_cores = 1
)
```

## Arguments

- expression:

  Numeric matrix, genes x samples (e.g. inverse-normal-transformed TPM).

- genotypes:

  Numeric dosage matrix, variants x samples (0-2).

- gene_pos:

  data.frame with columns `gene`, `chr`, and `tss` (or `start`; the gene
  start is used as the cis-window anchor).

- snp_pos:

  data.frame with columns `snp`, `chr`, `pos`.

- covariates:

  Optional data.frame or matrix (n rows) of covariates to partial out
  (e.g. ancestry PCs, technical factors).

- cis_snps:

  Optional named list: gene name -\> character vector of candidate SNP
  IDs (e.g. significant GTEx pairs). When supplied, only these SNPs are
  tested for each gene (still intersected with the cis window).

- cis_dist:

  Numeric: cis window half-width around the gene start, in bp. Default
  1e6.

- min_maf:

  Numeric: minimum minor allele frequency (computed from mean dosage).
  Default 0.01.

- fdr:

  Numeric: BH FDR level for the per-gene `pass` flag in the `best`
  table. Default 0.05.

- n_cores:

  Integer: parallel workers via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html).
  Default 1.

## Value

A list with:

- all:

  data.frame of all tested gene-SNP pairs: `gene`, `snp`, `beta`, `se`,
  `t`, `p`.

- best:

  data.frame with the top SNP per gene (`gene`, `snp`, `beta`, `p`, `q`
  (BH within gene), `pass`).

## Examples

``` r
n <- 80
dos <- matrix(rbinom(200 * n, 2, 0.3), nrow = 200,
              dimnames = list(paste0("1:", 1:200, ":A:G"), NULL))
expr <- matrix(rnorm(20 * n), nrow = 20,
               dimnames = list(paste0("Gene", 1:20), NULL))
expr[1, ] <- expr[1, ] + 0.5 * dos[5, ]   # one true cis-eQTL
gp <- data.frame(gene = paste0("Gene", 1:20), chr = "1",
                 tss = seq(1, 191, by = 10))
sp <- data.frame(snp = rownames(dos), chr = "1", pos = 1:200)
hits <- call_cis_eqtls(expr, dos, gp, sp)
head(hits$best)
#>     gene       snp       beta        se         t            p         q  pass
#> 1  Gene1 1:134:A:G  0.4403331 0.1796693  2.450798 0.0164921069 0.6190431 FALSE
#> 2 Gene10  1:58:A:G  0.4700744 0.1442743  3.258199 0.0016617558 0.2424813 FALSE
#> 3 Gene11   1:1:A:G  0.4160499 0.1590465  2.615900 0.0106810254 0.8235129 FALSE
#> 4 Gene12 1:135:A:G -0.5341922 0.1539292 -3.470377 0.0008499441 0.1699888 FALSE
#> 5 Gene13  1:84:A:G  0.6162287 0.2008581  3.067981 0.0029627583 0.5443979 FALSE
#> 6 Gene14 1:138:A:G -0.3848913 0.1576141 -2.441985 0.0168701271 0.9901600 FALSE
```
