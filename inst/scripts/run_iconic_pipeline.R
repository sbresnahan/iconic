# ============================================================
# End-to-end iconic pipeline.
#
# Cross-platform parallelization via n_cores argument on all
# computationally intensive functions (mclapply on Unix,
# PSOCK on Windows). Progress milestones via message().
# Model selection workflow (diagnose -> estimate ->
# stress-test -> recommend) + prospective analysis.
# GAN sensitivity-analysis pipeline (train generator ->
# sweep confounding -> recommend -> NC validity ->
# mediation benchmark -> plot).
#
# Uses the built-in example dataset when no real data is supplied;
# swap in your own matrices via load_real_input_data(
# Z_matrix, Y_matrix, covariates_df) and (optionally) a real W_matrix.
#
# Rscript inst/scripts/run_iconic_pipeline.R
#
# To use multiple cores, set N_CORES below or pass --args N.
# ============================================================

library(iconic)

## Parallelization setup --------------------------------------------------
# Auto-detect cores (leaving one for the OS) or use a user-specified count.
# Override with: Rscript inst/scripts/run_iconic_pipeline.R --args 8
n_cores <- 1
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  n_cores <- suppressWarnings(as.integer(args[1]))
  if (is.na(n_cores) || n_cores < 1) n_cores <- 1
}
message(sprintf("ICONIC pipeline: using %d core(s)", n_cores))

## 0. Model selection workflow -----------------------------------
# Standardize data, diagnose, estimate, stress-test, recommend.
# This is the primary user-facing workflow for choosing the right
# causal inference model for your data.

# Load real (or example) data and convert to iconic_data format
input <- load_real_input_data(example = TRUE)
message(sprintf("Loaded %d samples x %d features",
                input$n_samples, input$n_features))

# Convert to iconic_data (adds instruments/NCs if available)
data <- as_iconic_data(input)
print(data)

# Diagnose: which estimators are valid for this data?
diag <- iconic_diagnose(data, n_cores = n_cores)
print(diag)

# Estimate: fit all eligible estimators
est <- iconic_estimate(data, diagnosis = diag, n_cores = n_cores)
print(head(est))

# Stress-test: degradation surface (skip if no mediation data)
if (data$is_mediation) {
  sens <- iconic_sensitivity(data, diagnosis = diag, n_iter = 30,
                             n_cores = n_cores)
  print(sens)
} else {
  sens <- NULL
  message("Skipping sensitivity (no mediation data). ",
          "Use gan_sensitivity() for total-effect sensitivity.")
}

# Recommend: which model should you trust?
rec <- iconic_recommend(data, diagnosis = diag, estimate = est,
                        sensitivity = sens)
print(rec)

## 0a. Infer confounding parameters from data -----------------------------
# Estimate held-fixed confounding parameters (delta, mo_confounding,
# omega, k) from the data itself, for use with confounding = "inferred"
# in iconic_sensitivity() and iconic_prospect().
if (data$is_mediation) {
  conf <- infer_confounding(data, diagnosis = diag, estimate = est,
                            n_cores = n_cores)
  print(conf)
}

## 0b. Prospective analysis (if no instruments/NCs) ----------------------
# If you have Z, M, Y but no instruments or negative controls:
# prospect <- iconic_prospect(iconic_data(Z = my_Z, Y = my_Y, M = my_M),
# n_cores = n_cores)
# print(prospect)

## 1. GAN pipeline: train the generator ------------------------------------
# The GAN pipeline uses the same input data for generative benchmarking.
gan <- train_gan_on_real_data(input$gan_training_data, epochs = 300, seed = 1)
print(gan)

## 2. Total-effect sensitivity sweep --------------------------------------
# Vary confounding strength, negative-control coverage, and the number of
# latent confounders. Swap nc_model = "cpg" for the CpG-predicted-expression
# negative control, or pass your own function(U, covariates, params) -> W.
# Use mode = "distinct" in nc_proxy params for the completeness cliff.
sens <- gan_sensitivity(
  gan,
  conf_grid = c(0.2, 0.5, 0.8),
  coverage_grid = c(0.3, 0.6, 1.0),
  k_grid = c(1, 2),
  nc_model = "proxy",
  n_iter = 50,
  n_samples = 500,
  n_features = 20,
  n_cores = n_cores
)
print(head(sens$summary))

## 3. Which estimator is preferred? ----------------------------------------
# recommend_estimator excludes UNADJ (bias reference) and ranks the remaining
# four (DIRECT, COCA, IV2SLS, PGC) by worst-case RMSE across scenarios.
rec <- recommend_estimator(sens)
message("Robust across scenarios: ", rec$overall)
print(rec$per_scenario)

## 4. Do the negative controls hold under confounding? ---------------------
# n_valid_controls encodes how many distinct valid controls your design
# provides (proximal-inference identifiability needs it to be >= #confounders).
# The matrix-bridge PGC is under-identified when k > n_valid_controls;
# IV2SLS remains valid (it does not depend on NC completeness).
validity <- nc_validity_check(
  gan,
  coverage_grid = c(0.2, 0.5, 0.8, 1.0),
  k_grid = c(1, 2, 3),
  conf_strength = 0.8,
  n_valid_controls = 1,
  n_iter = 50
)
print(validity$verdict)

## 5. Mediation benchmark (NDE/NIE under M-O confounding) ------------------
# The mediation DGP adds U1 -> M (mediator-outcome confounding), creating
# the setting where natural effects are not point-identified by a single
# instrument (Rudolph et al. 2024). This benchmark quantifies the bias.
med_bench <- run_mediation_sim(
  n_iter = 50,
  n_samples = 500,
  n_features = 10,
  mo_confounding = 0.8,
  n_cores = n_cores
)
message("\nMediation benchmark (mo_confounding = 0.8):")
message("True NDE = ", med_bench$true_NDE,
        " True NIE = ", med_bench$true_NIE)
print(med_bench$summary[, c("method", "NDE_bias", "NIE_bias",
                            "NDE_rmse", "NIE_rmse", "NIE_type1")])

## 6. Mediation sensitivity sweep ------------------------------------------
# Sweep mediation estimators across confounding scenarios with M-O confounding.
med_sens <- gan_mediation_sensitivity(
  gan,
  conf_grid = c(0.2, 0.5, 0.8),
  coverage_grid = c(0.3, 0.7, 1.0),
  k_grid = 1,
  mo_confounding = 0.8,
  n_iter = 50,
  n_features = 12,
  n_cores = n_cores
)
print(head(med_sens$summary[, c("conf_strength", "coverage", "method",
                                "NDE_bias", "NIE_bias", "NIE_type1")]))

## 7. Null mediation Type I error ------------------------------------------
null_med <- run_null_mediation_sim(
  n_iter = 50,
  n_features = 10,
  mo_confounding = 0.8,
  n_cores = n_cores
)
message("\nNull NIE/NDE Type I error (mo_confounding = 0.8):")
print(null_med$rates)

## 8. Plots -----------------------------------------------------------------
plot_gan_diagnostics(gan, input$gan_training_data) # fidelity check
plot_sensitivity_heatmap(sens, metric = "rmse", method = "IV2SLS", k = 1)

## 9. Persist --------------------------------------------------------------
saveRDS(list(sensitivity = sens, recommendation = rec, validity = validity,
             mediation_benchmark = med_bench,
             mediation_sensitivity = med_sens,
             mediation_null = null_med),
        "iconic_results.rds")
message("Saved iconic_results.rds")
