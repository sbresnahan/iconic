## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, message = FALSE, warning = FALSE,
  fig.width = 7, fig.height = 4.5, comment = "#>"
)

set.seed(1)
library(iconic)

## ----walkthrough-step1--------------------------------------------------------
library(iconic)
set.seed(42)

dat <- generate_toy_data(n = 200, n_features = 10, mo_confounding = 0.8,
                         phi = 0.8, separate_U = TRUE,
                         omega_1 = 0.7, omega_2 = 0.7, seed = 42)
idata <- iconic_data(Z = dat$Z, Y = dat$Y, M = dat$M,
                     G = dat$G[, 1], Gm = dat$Gm,
                     W1 = dat$W1, W2 = dat$W2)
print(idata)

## ----walkthrough-step2--------------------------------------------------------
diag <- iconic_diagnose(idata)
print(diag)

## ----walkthrough-step2b-------------------------------------------------------
# Build a multi-mediator panel where eQTL strength varies
set.seed(7)
n <- 200; nm <- 50
G <- rnorm(n)
Z <- 0.3 * G + rnorm(n, 0, 0.9)
# Strong eQTLs for the first 25 mediators, weak for the rest
Gm_mat <- matrix(NA, nm, n)
M_mat  <- matrix(NA, nm, n)
for (m in 1:nm) {
  gm_m <- if (m <= 25) rnorm(n) else rnorm(n, 0, 0.3)
  Gm_mat[m, ] <- gm_m
  M_mat[m, ]  <- 0.5 * gm_m + 0.2 * Z + rnorm(n, 0, 0.7)
}
Y <- 0.1 * Z + 0.3 * M_mat[1, ] + rnorm(n, 0, 0.8)
W  <- matrix(rnorm(10 * n), 10, n)

idata2 <- iconic_data(Z = Z, Y = matrix(Y, 1), M = M_mat,
                      G = G, Gm = Gm_mat, W = W)

# Legacy scalar rule: median F_Gm vs min_f (default 10)
diag_scalar <- iconic_diagnose(idata2)
cat("Scalar rule:\n")
print(diag_scalar)

# Panel rule: eligible if >= 50% of transcripts have F_Gm >= 10
diag_panel <- iconic_diagnose(idata2,
  gm_threshold = list(E = 0.5, R = 10))
cat("\nPanel rule (E=50%, R=10):\n")
print(diag_panel)

# Inspect the full distribution
fgm <- diag_panel$instrument_strength$F_Gm
cat("\nF_Gm distribution (", length(fgm), " mediators):\n", sep = "")
print(summary(fgm))
for (R in c(3, 5, 10))
  cat(sprintf("  F_Gm >= %2d: %5.1f%%\n", R,
              100 * mean(fgm >= R, na.rm = TRUE)))

## ----walkthrough-step3--------------------------------------------------------
est <- iconic_estimate(idata, diagnosis = diag)
head(est, 8)
cat("True NDE =", dat$true_NDE, "  True NIE =", dat$true_NIE, "\n")

## ----walkthrough-step4--------------------------------------------------------
conf <- infer_confounding(idata, diagnosis = diag, estimate = est)
print(conf)

## ----walkthrough-step5--------------------------------------------------------
sens <- iconic_sensitivity(idata, diagnosis = diag,
                           rho_G1_grid = c(0, 0.3),
                           rho_G2_grid = c(0, 0.3),
                           n_iter = 5, n_features = 3,
                           gan_epochs = 20)
print(sens)

## ----walkthrough-step5b-------------------------------------------------------
sens_inf <- iconic_sensitivity(idata, diagnosis = diag,
                               rho_G1_grid = c(0, 0.3),
                               rho_G2_grid = c(0, 0.3),
                               n_iter = 5, n_features = 3,
                               confounding = "inferred", gan_epochs = 20)
print(sens_inf)
print(sens_inf$inferred_confounding)

## ----walkthrough-step6--------------------------------------------------------
rec <- iconic_recommend(idata, diagnosis = diag, estimate = est,
                        sensitivity = sens)
print(rec)

## ----walkthrough-step7--------------------------------------------------------
bare <- iconic_data(Z = rnorm(200), Y = matrix(rnorm(200 * 5), 5, 200),
                    M = rnorm(200))
prospect <- iconic_prospect(bare, gamma_G_grid = c(0.3, 0.6),
                            n_iter = 5, n_features = 3, gan_epochs = 20)
print(prospect)

## ----mediation-sim------------------------------------------------------------
# Generate one dataset with M-O confounding
dat <- generate_toy_data(n = 300, n_features = 5,
                                  mo_confounding = 0.8, seed = 42)
cat("True NDE =", dat$true_NDE, "  True NIE =", dat$true_NIE, "\n")

# Run all eligible mediation estimators
res <- analyze_mediation_robust(dat)
aggregate(NDE ~ method, res, function(v) round(mean(v, na.rm = TRUE), 3))
aggregate(NIE ~ method, res, function(v) round(mean(v, na.rm = TRUE), 3))

## ----mediation-bench, eval=FALSE----------------------------------------------
# # 100 replicates, 20 features, M-O confounding = 0.8
# bench <- run_mediation_sim(n_iter = 100, n_features = 20,
#                            mo_confounding = 0.8, n_cores = 4)
# bench$summary[, c("method", "NDE_bias", "NIE_bias",
#                      "NDE_rmse", "NIE_rmse")]

## ----mediation-sweep, eval=FALSE----------------------------------------------
# # Sweep confounding strength
# sweep_conf <- sweep_mediation_param("conf_str", c(0.2, 0.4, 0.6, 0.8, 1.0),
#                                     n_iter = 50, mo_confounding = 0.8)
# 
# # Sweep proxy quality
# sweep_omega <- sweep_mediation_param("w_signal",
#                                      c(0.2, 0.4, 0.6, 0.7, 0.8, 0.9),
#                                      n_iter = 50, mo_confounding = 0.8)
# 
# # Sweep sample size (bias should NOT shrink — it's structural)
# sweep_n <- sweep_mediation_param("n_samples", c(100, 200, 500, 1000),
#                                  n_iter = 50, mo_confounding = 0.8)
# 
# # Null Type I error across confounding strength
# null_t1e <- sweep_mediation_null_by_conf(
#     c(0.2, 0.4, 0.6, 0.8, 1.0), n_iter = 50)

## ----iv2sls2-demo-------------------------------------------------------------
# Generate one dataset with M-O confounding AND a mediator instrument
dat <- generate_toy_data(n = 500, n_features = 5,
                                  mo_confounding = 0.8,
                                  phi = 0.8, seed = 42)
cat("True NDE =", dat$true_NDE, "  True NIE =", dat$true_NIE, "\n")
cat("Gm present:", !is.null(dat$Gm), "\n")

# Run all eligible mediation estimators (IV2SLS2 included because Gm is present)
res <- analyze_mediation_robust(dat)
aggregate(NDE ~ method, res, function(v) round(mean(v, na.rm = TRUE), 4))
aggregate(NIE ~ method, res, function(v) round(mean(v, na.rm = TRUE), 4))

## ----iv2sls2-direct, eval=FALSE-----------------------------------------------
# # Direct usage: fit_iv2sls_mediation2() with both instruments
# # w can be a vector (single NC) or a matrix (full NC panel)
# fit <- fit_iv2sls_mediation2(
#   y = dat$Y[, 1], Z = dat$Z, M = dat$M,
#   g = dat$G[, 1], gm = dat$Gm, w = dat$W
# )
# fit$NDE   # natural direct effect
# fit$NIE   # natural indirect effect

## ----iv2sls2-bench, eval=FALSE------------------------------------------------
# # 50 replicates with mediator instrument
# bench <- run_mediation_sim(n_iter = 50, n_features = 10, n_samples = 500,
#                            mo_confounding = 0.8, phi = 0.8, n_cores = 4)
# bench$summary[, c("method", "NDE_bias", "NIE_bias",
#                      "NDE_rmse", "NIE_rmse")]

## ----phi-sweep, eval=FALSE----------------------------------------------------
# # Sweep mediator-instrument strength
# sweep_phi <- sweep_mediation_param("phi",
#                                    c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
#                                    n_iter = 30, n_features = 10,
#                                    n_samples = 500,
#                                    mo_confounding = 0.8, n_cores = 4)
# # IV2SLS2 bias vs phi
# iv2sls2 <- sweep_phi$summary[sweep_phi$summary$method == "IV2SLS2", ]
# iv2sls2[, c("param_value", "NDE_bias", "NIE_bias", "NDE_rmse", "NIE_rmse")]

## ----gm-independence----------------------------------------------------------
# Generate a dataset with a mediator instrument
dat <- generate_toy_data(n = 500, n_features = 10,
                                  mo_confounding = 0.8,
                                  phi = 0.8, seed = 1)

# A2': Are any controls associated with the mediator instrument?
gm_check <- nc_independence_check_gm(dat)
head(gm_check)
table(gm_check$verdict)

## ----gm-independence-absent, eval=FALSE---------------------------------------
# dat_no_gm <- generate_toy_data(n = 500, n_features = 10, seed = 1)
# nc_independence_check_gm(dat_no_gm)  # NULL: no mediator instrument

## ----pgc2-demo----------------------------------------------------------------
# Generate a dataset with imperfect independence, separate confounders
dat <- generate_toy_data(n = 500, n_features = 5,
                                  mo_confounding = 0.8, phi = 0.8,
                                  rho_G1 = 0.3, rho_G2 = 0.3,
                                  separate_U = TRUE,
                                  omega_1 = 0.7, omega_2 = 0.7, seed = 42)
cat("True NDE =", dat$true_NDE, "  True NIE =", dat$true_NIE, "\n")
cat("W1 present:", !is.null(dat$W1),
    "  W2 present:", !is.null(dat$W2), "\n")

# Run all eight estimators
# (PGC2, PGC2Gm included because W1/W2 are present)
res <- analyze_mediation_robust(dat)
aggregate(NDE ~ method, res, function(v) round(mean(v, na.rm = TRUE), 4))
aggregate(NIE ~ method, res, function(v) round(mean(v, na.rm = TRUE), 4))

## ----pgc2-direct--------------------------------------------------------------
# Without mediator instrument (pure NC identification)
fit_nc <- fit_pgc_mediation2(dat$Y[, 1], dat$Z, dat$M, dat$G1,
                             dat$W1, dat$W2)
cat("PGC2  NDE =", round(fit_nc$NDE, 4),
    " NIE =", round(fit_nc$NIE, 4), "\n")

# With (possibly imperfect) mediator instrument
fit_gm <- fit_pgc_mediation2(dat$Y[, 1], dat$Z, dat$M, dat$G1,
                             dat$W1, dat$W2, gm = dat$Gm)
cat("PGC2Gm NDE =", round(fit_gm$NDE, 4),
    " NIE =", round(fit_gm$NIE, 4), "\n")

## ----degradation-surface, eval=FALSE------------------------------------------
# # 2D sweep: both instrument exogeneity violations
# results <- data.frame()
# for (rg1 in c(0, 0.1, 0.2, 0.3, 0.5)) {
#   for (rg2 in c(0, 0.1, 0.2, 0.3, 0.5)) {
#     sweep_res <- sweep_mediation_param(
#       "rho_G2", rg2,
#       n_iter = 50, n_samples = 500, n_features = 10,
#       mo_confounding = 0.8, phi = 0.8, rho_G1 = rg1,
#       separate_U = TRUE, omega_1 = 0.7, omega_2 = 0.7)
#     smry <- sweep_res$summary
#     for (m in c("IV2SLS2", "PGC2Gm")) {
#       row <- smry[smry$method == m, ]
#       results <- rbind(results, data.frame(
#         rho_G1 = rg1, rho_G2 = rg2, method = m,
#         NDE_bias = row$NDE_bias, NIE_bias = row$NIE_bias))
#     }
#   }
# }

## ----train--------------------------------------------------------------------
# 1. Load real data as features x samples matrices + a covariate frame.
#    Here we use the built-in example; swap in your own:
#    load_real_input_data(Z_matrix = my_gdm, Y_matrix = my_expr,
#                         covariates_df = my_meta)
input <- load_real_input_data(example = TRUE)
str(input, max.level = 1)

# 2. Train the generator (requires torch).
gan <- train_gan_on_real_data(input$gan_training_data,
                                 epochs = 100, seed = 1)
gan

## ----diag, fig.height=5-------------------------------------------------------
plot_gan_diagnostics(gan, input$gan_training_data)

## ----iterate------------------------------------------------------------------
dat <- run_single_iteration(
  gan,
  n_synthetic_samples = 500,
  n_features          = 12,
  n_confounders       = 2,       # two latent confounders
  beta_Z = 0.10, alpha_M = 0.50, beta_M = 0.30,
  conf_strength       = 0.8,
  seed                = 42
)
cat("True total effect tau =", dat$true_total, "\n")
# 0.10 + 0.50*0.30 = 0.25

# Output has the shape run_methods()/
# analyze_methods_robust() expect.
res <- analyze_methods_robust(dat)
aggregate(beta ~ method, res, function(v) round(mean(v, na.rm = TRUE), 3))

## ----iterate-mo---------------------------------------------------------------
dat_mo <- run_single_iteration(
  gan,
  n_synthetic_samples = 500,
  n_features          = 12,
  n_confounders       = 1,
  mo_confounding      = 0.8,     # U1 -> M (mediator-outcome confounding)
  seed                = 42
)
cat("True NDE =", dat_mo$true_NDE, "  True NIE =", dat_mo$true_NIE, "\n")

# Run mediation estimators
res_mo <- analyze_mediation_robust(dat_mo)
aggregate(NDE ~ method, res_mo,
           function(v) round(mean(v, na.rm = TRUE), 3))
aggregate(NIE ~ method, res_mo,
           function(v) round(mean(v, na.rm = TRUE), 3))

## ----iterate-pleio------------------------------------------------------------
dat_pl <- run_single_iteration(
  gan,
  n_synthetic_samples = 500,
  n_features          = 12,
  n_confounders       = 1,
  pleio = 0.10,  # direct G -> Y path
                # (violates exclusion restriction)
  seed                = 42
)
# IV2SLS is now biased because G affects Y directly, not only through Z
res_pl <- analyze_methods_robust(dat_pl)
aggregate(beta ~ method, res_pl,
           function(v) round(mean(v, na.rm = TRUE), 3))

## ----nc-----------------------------------------------------------------------
list_nc_models()

# Use the CpG-predicted-expression controls instead of the direct proxy:
dat_cpg <- run_single_iteration(gan, n_features = 12, n_confounders = 1,
                                nc_model = "cpg", coverage = 0.7, seed = 7)
dim(dat_cpg$W)

# Or supply your own mechanism — anything honouring the contract:
my_nc <- function(U, covariates, params) {
  # controls = first confounder + independent noise
  scale(matrix(U[, 1], nrow(U), params$n_features) +
        matrix(rnorm(nrow(U) * params$n_features), nrow(U)))
}
dat_custom <- run_single_iteration(gan, n_features = 8,
                                      nc_model = my_nc, seed = 8)
dim(dat_custom$W)

## ----sens---------------------------------------------------------------------
sens <- gan_sensitivity(
  gan,
  conf_grid     = c(0.2, 0.5, 0.8),
  coverage_grid = c(0.3, 0.6, 1.0),
  k_grid        = c(1, 2),
  nc_model      = "proxy",
  n_iter        = 25,
  n_samples     = 400,
  n_features    = 12
)
head(sens$summary, 8)

## ----heat---------------------------------------------------------------------
plot_sensitivity_heatmap(sens, metric = "rmse", method = "IV2SLS", k = 1)

## ----rec----------------------------------------------------------------------
rec <- recommend_estimator(sens)
rec$overall
rec$worst_case

## ----validity-----------------------------------------------------------------
validity <- nc_validity_check(
  gan,
  coverage_grid    = c(0.2, 0.5, 1.0),
  k_grid           = c(1, 2, 3),
  conf_strength    = 0.8,
  n_valid_controls = 1,  # your design provides
                          # 1 distinct valid control
  n_iter           = 20,
  n_features       = 12
)
validity$verdict

## ----nc-diagnostics-----------------------------------------------------------
# Generate a dataset to illustrate (on real data, pass your fitted dataset)
dat <- run_single_iteration(gan, n_features = 10,
                               n_confounders = 1, seed = 1)

# A1: Are any controls associated with the exposure after adjusting for C?
screen_Z <- nc_validity_screen(dat)
head(screen_Z)
table(screen_Z$verdict)

# A2: Are any controls associated with the
# instrument after adjusting for C?
screen_G <- nc_independence_check(dat)
head(screen_G)
table(screen_G$verdict)

# Completeness: enough valid controls for the confounder dimensionality?
comp <- nc_completeness_check(dat)
comp$completeness
comp$n_valid_controls
comp$k

## ----pleio-sens---------------------------------------------------------------
pleio_sens <- gan_pleiotropy_sensitivity(
  gan,
  pleio_grid  = c(0, 0.05, 0.10),   # direct G -> Y path strengths
  conf_grid   = c(0.2, 0.5, 0.8),
  tau = 0.25,  # true total effect
                # for the alternative arm
  n_iter      = 25,
  n_features  = 10,
  n_samples   = 500
)
head(pleio_sens$summary, 12)

## ----iv-strength, eval=FALSE--------------------------------------------------
# # Demonstrate the weak-IV guard with an irrelevant instrument
# dat_weak <- run_single_iteration(NULL, n_synthetic_samples = 200,
#     n_features = 5, effect_size = 0.25, conf_strength = 0.8, seed = 1)
# dat_weak$G <- matrix(rnorm(1000), 200, 5)   # replace G with pure noise
# res <- analyze_methods_robust(dat_weak)
# # IV2SLS returns NA for all features (partial F < 10)
# res[res$method == "IV2SLS", ]

## ----med-sens, eval=FALSE-----------------------------------------------------
# med_sens <- gan_mediation_sensitivity(
#   gan,
#   conf_grid      = c(0.2, 0.5, 0.8),
#   coverage_grid  = c(0.3, 0.7, 1.0),
#   k_grid         = 1,
#   mo_confounding = 0.8,
#   n_iter         = 50,
#   n_features     = 12
# )
# head(med_sens$summary[, c("conf_strength", "coverage", "method",
#                           "NDE_bias", "NIE_bias", "NIE_type1")])

## ----feat-cor-sweep, eval=FALSE-----------------------------------------------
# # Sweep feature correlation strength
# # rho_G1 = 0 ensures the NC bridge is live (A2 passes, all 8 estimators eligible)
# fc_sweep <- sweep_mediation_param("feat_cor", c(0, 0.2, 0.4, 0.6, 0.8),
#     n_iter = 50, mo_confounding = 0.8, phi = 0.8,
#     rho_G1 = 0, separate_U = TRUE, omega_1 = 0.7, omega_2 = 0.7)
# fc_sweep$summary
# 
# # Plot the sweep (3 panels: bias, RMSE, Type I error)
# plot_feature_correlation_sweep(fc_sweep)

## ----pipeline, eval=FALSE-----------------------------------------------------
# library(iconic)
# 
# # ── model selection workflow ──
# # 1. Standardize data (optionally attach a pre-trained GAN)
# data <- iconic_data(Z = my_exposure, Y = my_outcomes, M = my_mediator,
#                     G = my_prs, Gm = my_eqtls, W1 = my_w1, W2 = my_w2,
#                     covariates = my_covariates)
# 
# # 2. Diagnose
# diag <- iconic_diagnose(data)
# 
# # 3. Estimate
# est <- iconic_estimate(data, diagnosis = diag)
# 
# # 4. Infer confounding parameters (optional)
# conf <- infer_confounding(data, diagnosis = diag, estimate = est)
# print(conf)
# 
# # 5. Sensitivity (auto-trains hybrid GAN + copula from data;
# #    or use confounding = "inferred")
# sens <- iconic_sensitivity(data, diagnosis = diag,
#                            confounding = "inferred", n_iter = 50)
# print(sens)
# 
# # 6. Recommend
# rec <- iconic_recommend(data, diagnosis = diag, estimate = est,
#                         sensitivity = sens)
# print(rec)
# 
# # ── Generative pipeline (standalone hybrid GAN + copula training) ──
# input <- load_real_input_data(Z_matrix = my_gdm,      # features x samples
#                               Y_matrix = my_expr,
#                               covariates_df = my_meta)
# gan   <- train_gan_on_real_data(input$gan_training_data, epochs = 300)
# 
# ## Total-effect sensitivity
# sens  <- gan_sensitivity(gan, conf_grid = c(0.2, 0.5, 0.8),
#                          coverage_grid = c(0.3, 0.6, 1.0),
#                          k_grid = c(1, 2, 3),
#                          nc_model = "cpg", n_iter = 100, n_cores = 4)
# 
# recommend_estimator(sens)$overall
# nc_validity_check(gan, k_grid = c(1, 2, 3), n_valid_controls = 2)$verdict
# 
# ## Empirical NC validity screens (on your real fitted dataset)
# # screen_Z <- nc_validity_screen(my_real_dat)
# # screen_G <- nc_independence_check(my_real_dat)
# # screen_Gm <- nc_independence_check_gm(my_real_dat)  # when phi > 0
# # nc_completeness_check(my_real_dat, n_valid_controls = 2)
# 
# ## Pleiotropy sensitivity
# ## (how robust to exclusion-restriction violations?)
# pleio_sens <- gan_pleiotropy_sensitivity(gan,
#     pleio_grid = c(0, 0.05, 0.10), conf_grid = c(0.2, 0.5, 0.8),
#     n_iter = 50, n_features = 10)
# 
# ## Mediation sensitivity (NDE/NIE under M-O confounding)
# ## Set phi > 0 to include the 2-stage MR estimator (IV2SLS2)
# med_sens <- gan_mediation_sensitivity(gan,
#     conf_grid = c(0.2, 0.5, 0.8), coverage_grid = c(0.3, 0.7, 1.0),
#     mo_confounding = 0.8, phi = 0.8, n_iter = 100, n_cores = 4)
# head(med_sens$summary[, c("method", "NDE_bias", "NIE_bias", "NIE_type1")])
# 
# plot_gan_diagnostics(gan, input$gan_training_data)
# plot_sensitivity_heatmap(sens, metric = "rmse", method = "IV2SLS")

