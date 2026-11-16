# ============================================================
# genomic_ranges.R: GenomicRanges-based coordinate handling for
# the mediator-instrument builders (call_cis_eqtls(),
# build_mediator_instruments()).
#
# Gene and variant coordinates may be supplied as data.frames
# (legacy interface), as GRanges, or implicitly through the
# rowRanges() of a RangedSummarizedExperiment. All inputs are
# normalised to GRanges and cis windows are resolved with
# GenomicRanges::findOverlaps().
# ============================================================


#' Normalise chromosome names to a common style (internal)
#'
#' Strips a leading "chr" (case-insensitive) and maps the
#' mitochondrial aliases "M"/"MT" to "MT", so "chr1" and "1" (or
#' "chrM" and "MT") are treated as the same sequence.
#' @param x Character or factor vector of chromosome names.
#' @return Character vector.
#' @keywords internal
#' @noRd
.norm_chr <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  x[toupper(x) %in% c("M", "MT")] <- "MT"
  x
}


#' Rebuild a GRanges with normalised seqnames (internal)
#' @noRd
.restyle_granges <- function(gr) {
  GenomicRanges::GRanges(
    seqnames = .norm_chr(GenomicRanges::seqnames(gr)),
    ranges   = IRanges::ranges(gr),
    strand   = GenomicRanges::strand(gr),
    S4Vectors::mcols(gr)
  )
}


#' Coerce gene coordinates to a 1-bp TSS GRanges (internal)
#'
#' @param gene_pos A data.frame with columns `gene`, `chr`, and `tss`
#'   (or `start`), or a GRanges whose gene IDs are `names()` or a
#'   `gene` metadata column. For GRanges the TSS is the strand-aware
#'   5' end; unstranded ranges use `start()`.
#' @param dup How to treat genes listed more than once: `"first"`
#'   keeps the first record, `"drop"` drops the gene.
#' @return A GRanges of width 1 at each gene's TSS, with metadata
#'   column `gene`, normalised seqnames, and no duplicated genes.
#' @keywords internal
#' @noRd
.as_gene_granges <- function(gene_pos, dup = c("first", "drop")) {
  dup <- match.arg(dup)
  if (methods::is(gene_pos, "GRanges")) {
    ids <- S4Vectors::mcols(gene_pos)$gene
    if (is.null(ids)) ids <- names(gene_pos)
    if (is.null(ids))
      stop("gene_pos GRanges needs gene IDs as names() or a `gene` ",
           "metadata column.")
    tss_gr <- GenomicRanges::resize(gene_pos, width = 1L, fix = "start")
    gr <- GenomicRanges::GRanges(
      seqnames = .norm_chr(GenomicRanges::seqnames(tss_gr)),
      ranges   = IRanges::IRanges(GenomicRanges::start(tss_gr), width = 1L),
      gene     = as.character(ids))
  } else {
    gene_pos <- as.data.frame(gene_pos)
    if (!all(c("gene", "chr") %in% colnames(gene_pos)))
      stop("gene_pos needs columns `gene` and `chr` (and `tss` or `start`), ",
           "or must be a GRanges.")
    if (!"tss" %in% colnames(gene_pos)) {
      if (!"start" %in% colnames(gene_pos))
        stop("gene_pos needs a `tss` or `start` column.")
      gene_pos$tss <- gene_pos$start
    }
    gr <- GenomicRanges::GRanges(
      seqnames = .norm_chr(gene_pos$chr),
      ranges   = IRanges::IRanges(as.integer(gene_pos$tss), width = 1L),
      gene     = as.character(gene_pos$gene))
  }
  g <- S4Vectors::mcols(gr)$gene
  if (dup == "first") {
    gr <- gr[!duplicated(g)]
  } else {
    gr <- gr[!(g %in% g[duplicated(g)])]
  }
  gr
}


#' Coerce variant coordinates to a 1-bp GRanges (internal)
#'
#' @param snp_pos A data.frame with columns `snp`, `chr`, `pos`, or a
#'   GRanges whose variant IDs are `names()` or an `snp` metadata
#'   column (position = `start()`).
#' @return A GRanges of width 1 with metadata column `snp` and
#'   normalised seqnames, in input order.
#' @keywords internal
#' @noRd
.as_snp_granges <- function(snp_pos) {
  if (methods::is(snp_pos, "GRanges")) {
    ids <- S4Vectors::mcols(snp_pos)$snp
    if (is.null(ids)) ids <- names(snp_pos)
    if (is.null(ids))
      stop("snp_pos GRanges needs variant IDs as names() or an `snp` ",
           "metadata column.")
    GenomicRanges::GRanges(
      seqnames = .norm_chr(GenomicRanges::seqnames(snp_pos)),
      ranges   = IRanges::IRanges(GenomicRanges::start(snp_pos), width = 1L),
      snp      = as.character(ids))
  } else {
    snp_pos <- as.data.frame(snp_pos)
    if (!all(c("snp", "chr", "pos") %in% colnames(snp_pos)))
      stop("snp_pos needs columns `snp`, `chr`, `pos`, or must be a GRanges.")
    GenomicRanges::GRanges(
      seqnames = .norm_chr(snp_pos$chr),
      ranges   = IRanges::IRanges(as.integer(snp_pos$pos), width = 1L),
      snp      = as.character(snp_pos$snp))
  }
}


#' Resolve cis variants for each gene with findOverlaps (internal)
#'
#' A variant is in cis when it lies on the same chromosome as the
#' gene's TSS and `|pos - tss| <= cis_dist` (a symmetric window).
#' @param gene_gr Output of `.as_gene_granges()`.
#' @param snp_gr Output of `.as_snp_granges()`.
#' @param cis_dist Window half-width in bp.
#' @return Named list: gene ID -> character vector of variant IDs, in
#'   the order the variants appear in `snp_gr`. Genes with no cis
#'   variant are absent.
#' @keywords internal
#' @noRd
.cis_snps_by_gene <- function(gene_gr, snp_gr, cis_dist) {
  half <- as.integer(floor(cis_dist))
  tss <- GenomicRanges::start(gene_gr)
  win <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gene_gr),
    ranges   = IRanges::IRanges(start = pmax(1L, tss - half),
                                end   = tss + half))
  hits <- GenomicRanges::findOverlaps(snp_gr, win, ignore.strand = TRUE)
  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)
  o <- order(s, q)
  genes <- S4Vectors::mcols(gene_gr)$gene
  snps  <- S4Vectors::mcols(snp_gr)$snp
  split(snps[q[o]], factor(genes[s[o]], levels = unique(genes[s[o]])))
}


#' Unpack expression / genotype inputs and coordinates (internal)
#'
#' Accepts plain matrices with explicit `gene_pos` / `snp_pos`, or
#' (Ranged)SummarizedExperiment objects whose `rowRanges()` supply the
#' coordinates when `gene_pos` / `snp_pos` are `NULL`.
#' @return list(expression, genotypes, gene_gr, snp_gr)
#' @keywords internal
#' @noRd
.prepare_cis_inputs <- function(expression, genotypes, gene_pos, snp_pos,
                                expr_assay = 1L, geno_assay = 1L,
                                dup = "first") {
  if (methods::is(expression, "SummarizedExperiment")) {
    if (is.null(gene_pos)) {
      if (!methods::is(expression, "RangedSummarizedExperiment"))
        stop("gene_pos is NULL but `expression` has no rowRanges(); ",
             "supply gene_pos or a RangedSummarizedExperiment.")
      gene_pos <- SummarizedExperiment::rowRanges(expression)
      if (is.null(names(gene_pos)) && is.null(S4Vectors::mcols(gene_pos)$gene))
        names(gene_pos) <- rownames(expression)
    }
    expression <- SummarizedExperiment::assay(expression, expr_assay)
  }
  if (methods::is(genotypes, "SummarizedExperiment")) {
    if (is.null(snp_pos)) {
      if (!methods::is(genotypes, "RangedSummarizedExperiment"))
        stop("snp_pos is NULL but `genotypes` has no rowRanges(); ",
             "supply snp_pos or a RangedSummarizedExperiment.")
      snp_pos <- SummarizedExperiment::rowRanges(genotypes)
      if (is.null(names(snp_pos)) && is.null(S4Vectors::mcols(snp_pos)$snp))
        names(snp_pos) <- rownames(genotypes)
    }
    genotypes <- SummarizedExperiment::assay(genotypes, geno_assay)
  }
  if (is.null(gene_pos)) stop("gene_pos is required.")
  if (is.null(snp_pos)) stop("snp_pos is required.")
  list(expression = as.matrix(expression),
       genotypes  = as.matrix(genotypes),
       gene_gr    = .as_gene_granges(gene_pos, dup = dup),
       snp_gr     = .as_snp_granges(snp_pos))
}
