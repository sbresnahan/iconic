# GenomicRanges / RangedSummarizedExperiment inputs to the mediator
# instrument builders must give results identical to the data.frame
# interface, and cis windows must match the documented definition
# (same chromosome, |pos - tss| <= cis_dist).

.cis_toy <- function(n = 80, seed = 21) {
    set.seed(seed)
    p <- 200
    dos <- matrix(rbinom(p * n, 2, 0.3), nrow = p,
        dimnames = list(paste0("rs", seq_len(p)),
            paste0("S", seq_len(n))))
    expr <- matrix(rnorm(20 * n), nrow = 20,
        dimnames = list(paste0("Gene", 1:20), paste0("S", seq_len(n))))
    expr[1, ] <- expr[1, ] + 0.6 * dos[5, ]
    expr[12, ] <- expr[12, ] + 0.6 * dos[160, ]
    chr_s <- rep(c("1", "2"), each = 100)
    pos_s <- c(seq(1000, 100000, length.out = 100),
        seq(1000, 100000, length.out = 100))
    gp <- data.frame(gene = paste0("Gene", 1:20),
        chr = rep(c("1", "2"), each = 10),
        tss = rep(seq(5000, 95000, length.out = 10), 2))
    sp <- data.frame(snp = rownames(dos), chr = chr_s, pos = round(pos_s))
    list(dos = dos, expr = expr, gp = gp, sp = sp)
}

test_that("cis windows follow |pos - tss| <= cis_dist on the same chr", {
    t <- .cis_toy()
    gr_g <- iconic:::.as_gene_granges(t$gp)
    gr_s <- iconic:::.as_snp_granges(t$sp)
    cm <- iconic:::.cis_snps_by_gene(gr_g, gr_s, cis_dist = 12000)
    for (g in t$gp$gene) {
        gi <- t$gp[t$gp$gene == g, ]
        ref <- t$sp$snp[t$sp$chr == gi$chr & abs(t$sp$pos - gi$tss) <= 12000]
        got <- cm[[g]]
        if (length(ref) == 0) expect_null(got) else expect_identical(got, ref)
    }
})

test_that("call_cis_eqtls gives identical results for GRanges and data.frame input", {
    t <- .cis_toy()
    ref <- call_cis_eqtls(t$expr, t$dos, t$gp, t$sp, cis_dist = 15000)
    gp_gr <- GenomicRanges::GRanges(paste0("chr", t$gp$chr),
        IRanges::IRanges(t$gp$tss, width = 1),
        gene = t$gp$gene)
    sp_gr <- GenomicRanges::GRanges(t$sp$chr,
        IRanges::IRanges(t$sp$pos, width = 1))
    names(sp_gr) <- t$sp$snp
    got <- call_cis_eqtls(t$expr, t$dos, gp_gr, sp_gr, cis_dist = 15000)
    expect_identical(got$all, ref$all)
    expect_identical(got$best, ref$best)
    expect_true("Gene1" %in% ref$best$gene)
})

test_that("strand-aware TSS is taken from GRanges gene models", {
    gr <- GenomicRanges::GRanges(c("1", "1"),
        IRanges::IRanges(c(100, 100), c(500, 500)),
        strand = c("+", "-"), gene = c("A", "B"))
    tss <- iconic:::.as_gene_granges(gr)
    expect_equal(GenomicRanges::start(tss), c(100, 500))
})

test_that("call_cis_eqtls reads coordinates from RangedSummarizedExperiment", {
    t <- .cis_toy()
    ref <- call_cis_eqtls(t$expr, t$dos, t$gp, t$sp, cis_dist = 15000)
    rse_e <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expr = t$expr),
        rowRanges = GenomicRanges::GRanges(t$gp$chr,
            IRanges::IRanges(t$gp$tss, width = 1)))
    rse_g <- SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = t$dos),
        rowRanges = GenomicRanges::GRanges(t$sp$chr,
            IRanges::IRanges(t$sp$pos, width = 1)))
    got <- call_cis_eqtls(rse_e, rse_g, cis_dist = 15000)
    expect_identical(got$best, ref$best)
    expect_error(call_cis_eqtls(t$expr, t$dos, cis_dist = 15000),
        "gene_pos is required")
})

test_that("build_mediator_instruments gives identical results for GRanges input", {
    skip_if_not_installed("glmnet")
    t <- .cis_toy(n = 120)
    ref <- build_mediator_instruments(t$expr, t$dos, t$gp, t$sp,
        cis_dist = 15000, seed = 3,
        verbose = FALSE)
    gp_gr <- GenomicRanges::GRanges(t$gp$chr,
        IRanges::IRanges(t$gp$tss, width = 1))
    names(gp_gr) <- t$gp$gene
    sp_gr <- GenomicRanges::GRanges(paste0("chr", t$sp$chr),
        IRanges::IRanges(t$sp$pos, width = 1),
        snp = t$sp$snp)
    got <- build_mediator_instruments(t$expr, t$dos, gp_gr, sp_gr,
        cis_dist = 15000, seed = 3,
        verbose = FALSE)
    expect_identical(got$qc, ref$qc)
    expect_identical(got$Gm, ref$Gm)
})

test_that("mixed chr naming styles are harmonised", {
    expect_identical(iconic:::.norm_chr(c("chr1", "1", "CHRX", "chrM", "MT")),
        c("1", "1", "X", "MT", "MT"))
})
