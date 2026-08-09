## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5, comment = "#>"
)
set.seed(1)
library(iconic)

## ----surv-data----------------------------------------------------------------
sdat <- generate_toy_data(
  n = 100, n_features = 3, seed = 1,
  outcome_type = "survival",
  phi = 0.8, omega_1 = 0.7, omega_2 = 0.7
)

# For survival, pass surv_time + surv_event instead of Y.
sdata <- iconic_data(
  X = sdat$X,
  surv_time = sdat$surv_time, surv_event = sdat$surv_event,
  M = sdat$M,
  G = sdat$G, Gm = sdat$Gm,
  W1 = sdat$W1, W2 = sdat$W2,
  outcome_type = "survival",
  covariates = sdat$synthetic_data
)

## ----surv-loghr---------------------------------------------------------------
est_loghr <- iconic_estimate(sdata, effect_scale = "loghr")
head(est_loghr)

## ----surv-rmst----------------------------------------------------------------
est_rmst <- iconic_estimate(sdata, effect_scale = "rmst")
head(est_rmst)

## ----surv-diagnose------------------------------------------------------------
diag <- iconic_diagnose(sdata)
diag$eligibility

