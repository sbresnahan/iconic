## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5, comment = "#>"
)
set.seed(1)
library(iconic)
has_torch <- check_torch_setup()

## ----install-bioc, eval=FALSE-------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("iconic")

## ----install-github, eval=FALSE-----------------------------------------------
# BiocManager::install("sbresnahan/iconic")

## ----sens-data----------------------------------------------------------------
dat <- generate_toy_data(
  n = 200, n_features = 10, seed = 1,
  phi = 0.8, omega_1 = 0.7, omega_2 = 0.7
)

idat <- iconic_data(
  X = dat$X, Y = dat$Y, M = dat$M,
  G = dat$G, Gm = dat$Gm,
  W1 = dat$W1, W2 = dat$W2,
  covariates = dat$synthetic_data
)

## ----sens-gan, eval=has_torch-------------------------------------------------
# # Train a texture model briefly and attach it to the iconic_data object.
# input <- load_real_input_data(
#   X_matrix = matrix(dat$X, nrow = 1),
#   Y_matrix = t(dat$Y), M_matrix = t(dat$M),
#   W_matrix = t(dat$W), covariates_df = dat$synthetic_data
# )
# 
# idat$trained_gan <- train_gan_on_real_data(
#   input$gan_training_data,
#   feature_correlations = input$feature_correlations,
#   feature_texture = input$feature_texture,
#   epochs = 5, seed = 1, verbose = FALSE
# )

## ----sens-no-torch, eval=!has_torch, echo=FALSE, results="asis"---------------
cat(
"The sensitivity-sweep chunks in this vignette are skipped because torch ",
"is not available. Install torch with `install.packages(\"torch\")` and ",
"run `torch::install_torch()` to run them."
)

## ----sens-run, eval=has_torch-------------------------------------------------
# sens <- iconic_sensitivity(idat, n_iter = 3)
# sens

## ----sens-infer---------------------------------------------------------------
diag <- iconic_diagnose(idat)
est <- iconic_estimate(idat, diagnosis = diag)
conf <- infer_confounding(idat, diag, est)
conf

## ----sens-inferred, eval=has_torch--------------------------------------------
# sens_inf <- iconic_sensitivity(idat, n_iter = 3, confounding = "inferred")
# sens_inf

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

