#' @references
#' Bresnahan ST, Xiong C, Head ST, Chang Y-H, Bhattacharya A, Huang JY.
#' ICONIC: An R Package for Integrating Instrumental Variable- and
#' Negative-Control-Informed Causal Discovery and Diagnostics in Multiomic
#' Studies. medRxiv (2026). \doi{10.64898/2026.08.26.26361466}
#' @keywords internal
"_PACKAGE"

#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom graphics abline axis barplot box boxplot grid image legend lines mtext par points polygon rect text
#' @importFrom stats aggregate as.formula ave coef complete.cases cov dbeta density fitted integrate lm median p.adjust pnorm prcomp pt quantile rbinom residuals rnorm runif sd setNames vcov
#' @importFrom utils modifyList tail
#' @importFrom methods is
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList experiments sampleMap
#' @importFrom SummarizedExperiment SummarizedExperiment assay colData rowRanges
#' @importFrom S4Vectors DataFrame metadata mcols queryHits subjectHits
#' @importFrom GenomicRanges GRanges seqnames strand start resize findOverlaps
#' @importFrom IRanges IRanges ranges
NULL

# `self` is bound by torch::nn_module() inside the generator/discriminator
# forward methods; declare it to avoid a spurious "no visible binding" note.
utils::globalVariables("self")

# ggplot2 uses non-standard evaluation (NSE) for aes() mappings, so column
# names referenced inside aes() appear as "no visible binding" global
# variables to R CMD check. Declare them here to suppress the NOTE.
# These are all data-frame column names used in aes() across the plotting
# functions in figures.R and plots.R.
utils::globalVariables(c(
    ".data", "bias", "bias_sd", "capture_R2_mean", "capture_R2_sd",
    "completeness", "conf_str", "coverage", "eligible", "estimate",
    "estimand", "estimator", "feat_cor", "frac_sig", "gamma_G", "group",
    "k", "label", "mean_F", "method", "n_valid", "NDE_bias", "NIE_type1",
    "null_R2_mean", "omega", "partial_F", "pgc_bias", "pleio", "rate",
    "recommended", "rejected", "rho_G1", "rho_G2", "se", "support_R2_mean",
    "support_R2_sd", "t1e", "value", "winner", "x", "xg", "xint", "xvar",
    "yint"
))
