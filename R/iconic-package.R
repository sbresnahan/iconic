#' @keywords internal
"_PACKAGE"

#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom graphics abline axis barplot box boxplot grid image legend lines mtext par points polygon rect text
#' @importFrom stats aggregate as.formula ave coef complete.cases cov dbeta density fitted integrate lm median p.adjust pnorm prcomp pt quantile rbinom residuals rnorm runif sd setNames vcov
#' @importFrom utils modifyList tail
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
  ".data", "bias", "bias_sd", "completeness", "conf_str", "coverage",
  "eligible", "estimate", "estimand", "estimator", "feat_cor", "gamma_G", "group",
  "k", "label", "mean_F", "method", "n_valid", "NDE_bias", "NIE_type1",
  "partial_F", "pgc_bias", "pleio", "rate", "recommended", "rejected",
  "rho_G1", "rho_G2", "se", "t1e", "value", "winner", "x", "xg", "xint",
  "xvar", "yint"
))
