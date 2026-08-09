## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5, comment = "#>"
)

set.seed(1)
library(iconic)
has_torch <- check_torch_setup()

## ----workflow-data------------------------------------------------------------
set.seed(1)
dat <- generate_toy_data(n = 200, n_features = 10, seed = 1)

# Assemble an iconic_data object. Y and M are features x samples matrices;
# X is a length-n vector; G is an n x n_snps matrix; W is an n x n_ncs matrix.
idat <- iconic_data(
  X = dat$X, Y = dat$Y, M = dat$M,
  G = dat$G, W = dat$W,
  covariates = dat$synthetic_data
)

## ----workflow-diagnose--------------------------------------------------------
diag <- iconic_diagnose(idat)
diag$eligibility

## ----workflow-estimate--------------------------------------------------------
est <- iconic_estimate(idat)
head(est)

## ----workflow-gan, eval=has_torch---------------------------------------------
# # Train a texture model once (few epochs for the vignette).
# input <- load_real_input_data(
#   X_matrix = matrix(dat$X, nrow = 1),
#   Y_matrix = t(dat$Y), M_matrix = t(dat$M),
#   W_matrix = t(dat$W), covariates_df = dat$synthetic_data
# )
# 
# gan <- train_gan_on_real_data(
#   input$gan_training_data,
#   feature_correlations = input$feature_correlations,
#   feature_texture = input$feature_texture,
#   epochs = 5, seed = 1, verbose = FALSE
# )
# 
# # Attach it: iconic_sensitivity() and iconic_prospect() will now reuse it.
# idat$trained_gan <- gan

## ----workflow-no-torch, eval=!has_torch, echo=FALSE, results="asis"-----------
cat(
"The sensitivity, recommendation, and prospective chunks in this vignette ",
"are skipped because torch is not available. Install torch with ",
"`install.packages(\"torch\")` and run `torch::install_torch()` to run them."
)

## ----workflow-sensitivity, eval=has_torch-------------------------------------
# sens <- iconic_sensitivity(idat, n_iter = 3)
# sens

## ----workflow-confounding-----------------------------------------------------
conf <- infer_confounding(idat, diag, est)
conf

## ----workflow-recommend, eval=has_torch---------------------------------------
# rec <- iconic_recommend(idat, diag, est, sens)
# rec

## ----workflow-prospect, eval=has_torch----------------------------------------
# pros <- iconic_prospect(idat, n_iter = 3)
# pros

## ----workflow-survival--------------------------------------------------------
set.seed(1)
sdat <- generate_toy_data(
  n = 100, n_features = 3,
  outcome_type = "survival", seed = 1
)

# For survival, pass surv_time + surv_event instead of Y.
sdata <- iconic_data(
  X = sdat$X,
  surv_time = sdat$surv_time, surv_event = sdat$surv_event,
  G = sdat$G, W = sdat$W,
  outcome_type = "survival",
  covariates = sdat$synthetic_data
)

sest <- iconic_estimate(sdata)
head(sest)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

