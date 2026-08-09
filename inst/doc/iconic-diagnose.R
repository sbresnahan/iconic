## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5, comment = "#>"
)
set.seed(1)
library(iconic)

## ----diagnose-data------------------------------------------------------------
dat <- generate_toy_data(
  n = 100, n_features = 5, seed = 1,
  phi = 0.8, omega_1 = 0.7, omega_2 = 0.7
)

idat <- iconic_data(
  X = dat$X, Y = dat$Y, M = dat$M,
  G = dat$G, Gm = dat$Gm,
  W1 = dat$W1, W2 = dat$W2,
  covariates = dat$synthetic_data
)

## ----diagnose-run-------------------------------------------------------------
diag <- iconic_diagnose(idat)
diag$eligibility

## ----diagnose-instrument------------------------------------------------------
diag$instrument_strength

## ----diagnose-nc--------------------------------------------------------------
diag$nc_validity

## ----diagnose-no-proxy, eval=FALSE--------------------------------------------
# dat0 <- generate_toy_data(n = 100, n_features = 5, seed = 1)
# idat0 <- iconic_data(X = dat0$X, Y = dat0$Y, covariates = dat0$synthetic_data)
# diag0 <- iconic_diagnose(idat0)

