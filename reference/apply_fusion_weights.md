# Apply FUSION/TWAS eQTL weights to a dosage matrix

Predicts genetically regulated expression by applying published FUSION
per-gene weight files (`.wgt.RDat`) to a genotype dosage matrix. For
each gene, the model with the best cross-validation performance
(`cv.performance["rsq", ]`) is used. The resulting predicted-expression
matrix is germline-determined and — after PCA via
[`build_w_pcs()`](https://seantbresnahan.com/iconic/reference/build_w_pcs.md)
— serves as a negative-control panel that captures confounding structure
without direct exposure mediation.

## Usage

``` r
apply_fusion_weights(
  genotypes,
  pos,
  weights_dir = NULL,
  weights = NULL,
  min_snps = 2,
  n_cores = 1
)
```

## Arguments

- genotypes:

  Numeric dosage matrix, variants x samples, with rownames matching the
  SNP IDs used in the weight files (usually rsIDs).

- pos:

  A data.frame with (at least) columns `ID` (gene ID) and `WGT`
  (weight-file name), as read from a FUSION `.pos` file; or a path to
  such a file.

- weights_dir:

  Directory containing the `.wgt.RDat` files named in `pos$WGT`. Ignored
  when `weights` is supplied.

- weights:

  Optional named list of pre-loaded weight objects, each a list with
  elements `wgt.matrix`, `snps`, and `cv.performance` (as stored in
  FUSION `.wgt.RDat` files). Names must match `pos$ID`. Supplying this
  avoids repeated file I/O and makes the function usable without the
  original files.

- min_snps:

  Integer: minimum number of matched SNPs required to predict a gene.
  Default 2.

- n_cores:

  Integer: parallel workers via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html).
  Default 1.

## Value

A list with:

- predicted:

  numeric matrix, genes x samples: predicted expression.

- info:

  data.frame per predicted gene: `gene`, `n_snps`, `best_method`, `rsq`
  (cross-validation R2 of the selected model).

## Examples

``` r
# Simulate two FUSION-style weight objects in memory
dos <- matrix(rbinom(10 * 50, 2, 0.3), nrow = 10,
              dimnames = list(paste0("rs", 1:10), paste0("S", 1:50)))
mk_wgt <- function(snps, w) {
  list(wgt.matrix = matrix(w, ncol = 1, dimnames = list(NULL, "enet")),
       snps = data.frame(V2 = snps),
       cv.performance = matrix(0.2, nrow = 1,
                               dimnames = list("rsq", "enet")))
}
wlist <- list(GeneA = mk_wgt(paste0("rs", 1:5), rep(0.1, 5)),
              GeneB = mk_wgt(paste0("rs", 6:10), rep(-0.2, 5)))
pos <- data.frame(ID = c("GeneA", "GeneB"), WGT = c("a.wgt.RDat", "b.wgt.RDat"))
out <- apply_fusion_weights(dos, pos = pos, weights = wlist)
dim(out$predicted)
#> [1]  2 50
```
