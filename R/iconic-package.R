#' @keywords internal
"_PACKAGE"

#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom graphics abline axis barplot box boxplot grid image legend lines
#'   mtext par points polygon rect text
#' @importFrom stats as.formula coef complete.cases cov density fitted lm median
#'   pnorm rbinom residuals rnorm runif sd setNames vcov
#' @importFrom utils modifyList tail
NULL

# `self` is bound by torch::nn_module() inside the generator/discriminator
# forward methods; declare it to avoid a spurious "no visible binding" note.
utils::globalVariables("self")
