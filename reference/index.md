# Package index

## Data input & containers

Standardise exposure, outcome, mediator, instrument, and
negative-control data into the object used across the package, and
simulate ground-truth panels.

- [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
  : Construct a standardized data object for ICONIC model selection
- [`as_iconic_data()`](https://seantbresnahan.com/iconic/reference/as_iconic_data.md)
  : Convert external data containers to iconic_data
- [`load_real_input_data()`](https://seantbresnahan.com/iconic/reference/load_real_input_data.md)
  : Load and standardise a real multi-omic dataset for GAN training
- [`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md)
  : Generate one synthetic dataset (internal)
- [`simulate_single_genetic_instrument()`](https://seantbresnahan.com/iconic/reference/simulate_single_genetic_instrument.md)
  : Simulate a single genetic instrument

## Exposure instruments

Build and validate a genetic instrument for the exposure from GWAS
summary statistics or a published polygenic-score panel.

- [`qc_gwas_sumstats()`](https://seantbresnahan.com/iconic/reference/qc_gwas_sumstats.md)
  : Standardize and QC GWAS summary statistics
- [`build_prs_ldpred2()`](https://seantbresnahan.com/iconic/reference/build_prs_ldpred2.md)
  : Build an exposure polygenic score with LDpred2-auto
- [`score_pgs_panel()`](https://seantbresnahan.com/iconic/reference/score_pgs_panel.md)
  : Score a published polygenic-score panel on a dosage matrix
- [`check_instrument_strength()`](https://seantbresnahan.com/iconic/reference/check_instrument_strength.md)
  : Check first-stage instrument strength (partial F)

## Mediator instruments

Call cis-eQTLs and train genetically-predicted expression instruments
for the mediator panel.

- [`call_cis_eqtls()`](https://seantbresnahan.com/iconic/reference/call_cis_eqtls.md)
  : Scan for cis-eQTLs of each gene
- [`build_mediator_instruments()`](https://seantbresnahan.com/iconic/reference/build_mediator_instruments.md)
  : Build mediator instruments via per-gene elastic net

## Negative controls

Construct negative-control panels and run the validity, independence,
support, and completeness diagnostics.

- [`beta_to_m()`](https://seantbresnahan.com/iconic/reference/beta_to_m.md)
  : Convert methylation beta values to M-values

- [`residualize_matrix()`](https://seantbresnahan.com/iconic/reference/residualize_matrix.md)
  : Residualize a feature matrix on covariates

- [`build_w_pcs()`](https://seantbresnahan.com/iconic/reference/build_w_pcs.md)
  : Build a negative-control panel from principal components

- [`apply_fusion_weights()`](https://seantbresnahan.com/iconic/reference/apply_fusion_weights.md)
  : Apply FUSION/TWAS eQTL weights to a dosage matrix

- [`nc_proxy()`](https://seantbresnahan.com/iconic/reference/nc_proxy.md)
  : Direct-proxy negative-control model

- [`nc_cpg()`](https://seantbresnahan.com/iconic/reference/nc_cpg.md) :
  CpG-predicted-expression negative-control model (SCENIC case)

- [`list_nc_models()`](https://seantbresnahan.com/iconic/reference/list_nc_models.md)
  : List the built-in negative-control models

- [`nc_validity_check()`](https://seantbresnahan.com/iconic/reference/nc_validity_check.md)
  : Check negative-control validity across confounding scenarios

- [`nc_validity_screen()`](https://seantbresnahan.com/iconic/reference/nc_validity_screen.md)
  :

  Screen negative controls for exposure dependence (W *\|* X \| C)

- [`nc_independence_check()`](https://seantbresnahan.com/iconic/reference/nc_independence_check.md)
  :

  Test instrument-independence of negative controls (W *\|* G \| C)

- [`nc_independence_check_gm()`](https://seantbresnahan.com/iconic/reference/nc_independence_check_gm.md)
  :

  Test mediator-instrument independence of negative controls (W *\|* Gm
  \| C)

- [`nc_support_check()`](https://seantbresnahan.com/iconic/reference/nc_support_check.md)
  : Negative-control support/range check

- [`nc_completeness_check()`](https://seantbresnahan.com/iconic/reference/nc_completeness_check.md)
  : Check negative-control completeness (dimensional +
  covariance-capture)

- [`nc_completeness_capture()`](https://seantbresnahan.com/iconic/reference/nc_completeness_capture.md)
  : Covariance-capture completeness test

## Core workflow

The end-to-end model-selection pipeline: diagnose eligibility, estimate
effects, stress-test robustness, and recommend an estimator.

- [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
  : Diagnose data and determine estimator eligibility
- [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md)
  : Fit all eligible estimators on real data
- [`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
  : Sensitivity (degradation) surface and effect-decomposition bias
  sweep
- [`iconic_recommend()`](https://seantbresnahan.com/iconic/reference/iconic_recommend.md)
  : Recommend the best causal estimator for the user's data
- [`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md)
  : Prospective bias-reduction analysis for data without instruments or
  negative controls
- [`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)
  : Infer confounding parameters from the user's data
- [`recommend_estimator()`](https://seantbresnahan.com/iconic/reference/recommend_estimator.md)
  : Recommend the preferred estimator from a sensitivity sweep

## Estimators (low-level)

Individual NDE/NIE estimators. Most users call these through
iconic_estimate() rather than directly.

- [`fit_coca()`](https://seantbresnahan.com/iconic/reference/fit_coca.md)
  : COCA estimator: Negative-Control Outcome Correction via ratio
- [`fit_coca_bin()`](https://seantbresnahan.com/iconic/reference/fit_coca_bin.md)
  : COCA binary estimator: NOT supported (returns NA)
- [`fit_coca_mediation()`](https://seantbresnahan.com/iconic/reference/fit_coca_mediation.md)
  : COCA mediation estimator: negative-control calibration of both
  stages
- [`fit_coca_mediation_bin()`](https://seantbresnahan.com/iconic/reference/fit_coca_mediation_bin.md)
  : COCA binary mediation estimator: NOT supported (returns NA)
- [`fit_coca_mediation_surv()`](https://seantbresnahan.com/iconic/reference/fit_coca_mediation_surv.md)
  : COCA survival mediation estimator: NOT supported (returns NA)
- [`fit_coca_surv()`](https://seantbresnahan.com/iconic/reference/fit_coca_surv.md)
  : COCA survival estimator: NOT supported (returns NA)
- [`fit_direct()`](https://seantbresnahan.com/iconic/reference/fit_direct.md)
  : DIRECT estimator: OLS with instrument and negative-control as
  covariates
- [`fit_direct_bin()`](https://seantbresnahan.com/iconic/reference/fit_direct_bin.md)
  : DIRECT binary estimator: logistic / LPM with instrument and NC
  covariates
- [`fit_direct_mediation()`](https://seantbresnahan.com/iconic/reference/fit_direct_mediation.md)
  : DIRECT mediation estimator: OLS with instrument and NC as covariates
- [`fit_direct_mediation_bin()`](https://seantbresnahan.com/iconic/reference/fit_direct_mediation_bin.md)
  : DIRECT binary mediation estimator: logistic / LPM with G and W
  covariates
- [`fit_direct_mediation_surv()`](https://seantbresnahan.com/iconic/reference/fit_direct_mediation_surv.md)
  : DIRECT survival mediation estimator: Cox / RMST with G and W
  covariates
- [`fit_direct_surv()`](https://seantbresnahan.com/iconic/reference/fit_direct_surv.md)
  : DIRECT survival estimator: Cox / RMST with instrument and NC
  covariates
- [`fit_iv2sls()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls.md)
  : IV2SLS estimator: Two-Stage Least Squares with genetic instrument
- [`fit_iv2sls_bin()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_bin.md)
  : IV2SLS binary estimator: two-stage predictor substitution with
  logistic / LPM
- [`fit_iv2sls_mediation()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation.md)
  : IV2SLS mediation estimator: instrumented exposure in both stages
- [`fit_iv2sls_mediation2()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2.md)
  : IV2SLS2 mediation estimator: 2-stage MR with instruments for both X
  and M
- [`fit_iv2sls_mediation2_bin()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2_bin.md)
  : IV2SLS2 binary mediation estimator: 2-stage MR with logistic / LPM
  outcome
- [`fit_iv2sls_mediation2_surv()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2_surv.md)
  : IV2SLS2 survival mediation estimator: 2-stage MR with Cox / RMST
  outcome
- [`fit_iv2sls_mediation_bin()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation_bin.md)
  : IV2SLS binary mediation estimator: single-instrument 2SPS with
  logistic / LPM
- [`fit_iv2sls_mediation_surv()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation_surv.md)
  : IV2SLS survival mediation estimator: single-instrument 2SPS with Cox
  / RMST
- [`fit_iv2sls_surv()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_surv.md)
  : IV2SLS survival estimator: two-stage predictor substitution with Cox
  / RMST
- [`fit_pgc()`](https://seantbresnahan.com/iconic/reference/fit_pgc.md)
  : PGC estimator: Proxy G-Component Correction (matrix bridge)
- [`fit_pgc_bin()`](https://seantbresnahan.com/iconic/reference/fit_pgc_bin.md)
  : PGC binary estimator: proxy G-component correction with logistic /
  LPM
- [`fit_pgc_mediation()`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation.md)
  : PGC mediation estimator: bridge-function-adjusted natural effects
  (matrix bridge)
- [`fit_pgc_mediation2()`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2.md)
  : PGC-2 mediation estimator: two-stage proximal mediation with
  path-specific bridges
- [`fit_pgc_mediation2_bin()`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2_bin.md)
  : PGC2 / PGC2Gm binary mediation estimator: path-specific bridges with
  logistic / LPM
- [`fit_pgc_mediation2_surv()`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation2_surv.md)
  : PGC2 / PGC2Gm survival mediation estimator: path-specific bridges
  with Cox / RMST
- [`fit_pgc_mediation_bin()`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation_bin.md)
  : PGC binary mediation estimator: single-panel bridge with logistic /
  LPM
- [`fit_pgc_mediation_surv()`](https://seantbresnahan.com/iconic/reference/fit_pgc_mediation_surv.md)
  : PGC survival mediation estimator: single-panel bridge with Cox /
  RMST
- [`fit_pgc_scalar()`](https://seantbresnahan.com/iconic/reference/fit_pgc_scalar.md)
  : PGC estimator: Proxy G-Component Correction (scalar bridge)
- [`fit_pgc_scalar_mediation()`](https://seantbresnahan.com/iconic/reference/fit_pgc_scalar_mediation.md)
  : PGC mediation estimator: bridge-function-adjusted natural effects
  (scalar bridge)
- [`fit_pgc_surv()`](https://seantbresnahan.com/iconic/reference/fit_pgc_surv.md)
  : PGC survival estimator: proxy G-component correction with Cox / RMST
- [`fit_unadj_bin()`](https://seantbresnahan.com/iconic/reference/fit_unadj_bin.md)
  : UNADJ binary estimator: unadjusted logistic / linear-probability
  regression
- [`fit_unadj_mediation()`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation.md)
  : UNADJ mediation estimator: naive Baron-Kenny style
- [`fit_unadj_mediation_bin()`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_bin.md)
  : UNADJ binary mediation estimator: naive Baron-Kenny with logistic /
  LPM
- [`fit_unadj_mediation_surv()`](https://seantbresnahan.com/iconic/reference/fit_unadj_mediation_surv.md)
  : UNADJ survival mediation estimator: naive Baron-Kenny with Cox /
  RMST
- [`fit_unadj_surv()`](https://seantbresnahan.com/iconic/reference/fit_unadj_surv.md)
  : UNADJ survival estimator: unadjusted Cox / RMST regression

## Sensitivity & robustness analysis

Generative-model stress tests of instrument exogeneity and pleiotropy,
and p-value combination helpers.

- [`gan_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_sensitivity.md)
  : Benchmark estimators across confounding scenarios on synthetic data
- [`gan_mediation_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_mediation_sensitivity.md)
  : Benchmark mediation estimators across confounding scenarios
- [`gan_pleiotropy_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_pleiotropy_sensitivity.md)
  : Benchmark estimators across pleiotropy and confounding scenarios
- [`composite_p_value()`](https://seantbresnahan.com/iconic/reference/composite_p_value.md)
  : Composite null hypothesis test p-value

## Simulation & sweeps

Simulation drivers and parameter-sweep machinery used for estimator
validation and benchmarking.

- [`run_simulation()`](https://seantbresnahan.com/iconic/reference/run_simulation.md)
  : Run repeated simulations for a single parameter configuration

- [`run_single_iteration()`](https://seantbresnahan.com/iconic/reference/run_single_iteration.md)
  : Generate one synthetic dataset under the generalised SCM

- [`run_mediation_sim()`](https://seantbresnahan.com/iconic/reference/run_mediation_sim.md)
  : Run repeated mediation simulations for a single parameter
  configuration

- [`run_null_sim()`](https://seantbresnahan.com/iconic/reference/run_null_sim.md)
  : Run null simulations to estimate Type I error rates

- [`run_null_mediation_sim()`](https://seantbresnahan.com/iconic/reference/run_null_mediation_sim.md)
  : Run null mediation simulations to estimate Type I error rates

- [`sweep_instrument_strength()`](https://seantbresnahan.com/iconic/reference/sweep_instrument_strength.md)
  : Sweep instrument strength

- [`sweep_mediation_null_by_conf()`](https://seantbresnahan.com/iconic/reference/sweep_mediation_null_by_conf.md)
  : Sweep mediation Type I error across confounding strength levels

- [`sweep_mediation_param()`](https://seantbresnahan.com/iconic/reference/sweep_mediation_param.md)
  : Sweep a single mediation simulation parameter across a grid

- [`sweep_nc_validity()`](https://seantbresnahan.com/iconic/reference/sweep_nc_validity.md)
  : Sweep negative-control validity diagnostics

- [`sweep_null_by_conf()`](https://seantbresnahan.com/iconic/reference/sweep_null_by_conf.md)
  : Sweep Type I error rate across confounding strength levels

- [`sweep_param()`](https://seantbresnahan.com/iconic/reference/sweep_param.md)
  : Sweep a single simulation parameter across a grid

- [`analyze_methods_robust()`](https://seantbresnahan.com/iconic/reference/analyze_methods_robust.md)
  : Run all five estimators (plus UNADJ) on one synthetic dataset

- [`analyze_methods_parallel()`](https://seantbresnahan.com/iconic/reference/analyze_methods_parallel.md)
  :

  Parallel version of
  [`analyze_methods_robust()`](https://seantbresnahan.com/iconic/reference/analyze_methods_robust.md)

- [`analyze_mediation_robust()`](https://seantbresnahan.com/iconic/reference/analyze_mediation_robust.md)
  : Run all mediation estimators on one synthetic dataset

- [`scenario_manifest()`](https://seantbresnahan.com/iconic/reference/scenario_manifest.md)
  : Scenario manifest: truth and parameter ranges for a simulation

## Generative texture model

Train and sample from the hybrid GAN + Gaussian-copula texture model
that calibrates synthetic data to a cohort.

- [`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md)
  : Train a generative texture model on real data
- [`train_feature_texture()`](https://seantbresnahan.com/iconic/reference/train_feature_texture.md)
  : Train a feature-level texture model for the mediator panel
- [`sample_texture()`](https://seantbresnahan.com/iconic/reference/sample_texture.md)
  : Draw synthetic base rows from a trained texture model
- [`sample_feature_texture()`](https://seantbresnahan.com/iconic/reference/sample_feature_texture.md)
  : Draw synthetic feature vectors from a trained feature texture model
- [`check_torch_setup()`](https://seantbresnahan.com/iconic/reference/check_torch_setup.md)
  : Check whether a working torch installation is available

## Plotting

Visualization functions for diagnostics, sensitivity, and results.

- [`plot_bias()`](https://seantbresnahan.com/iconic/reference/plot_bias.md)
  : Plot absolute bias vs a swept parameter
- [`plot_bias_boxplot()`](https://seantbresnahan.com/iconic/reference/plot_bias_boxplot.md)
  : Grouped boxplots of per-seed bias across a parameter sweep
- [`plot_bias_distribution()`](https://seantbresnahan.com/iconic/reference/plot_bias_distribution.md)
  : Baseline bias distribution (single-setting hero plot)
- [`plot_degradation_surface()`](https://seantbresnahan.com/iconic/reference/plot_degradation_surface.md)
  : Degradation surface figure
- [`plot_estimate_distribution()`](https://seantbresnahan.com/iconic/reference/plot_estimate_distribution.md)
  : Boxplot of estimate distributions from run_simulation()
- [`plot_estimated_vs_true()`](https://seantbresnahan.com/iconic/reference/plot_estimated_vs_true.md)
  : Plot estimated effect vs true effect
- [`plot_estimator_benchmark()`](https://seantbresnahan.com/iconic/reference/plot_estimator_benchmark.md)
  : Estimator benchmark figure
- [`plot_feature_correlation_sweep()`](https://seantbresnahan.com/iconic/reference/plot_feature_correlation_sweep.md)
  : Feature correlation sweep figure
- [`plot_gan_diagnostics()`](https://seantbresnahan.com/iconic/reference/plot_gan_diagnostics.md)
  : Compare real vs synthetic marginals from a trained generator
- [`plot_instrument_strength_sweep()`](https://seantbresnahan.com/iconic/reference/plot_instrument_strength_sweep.md)
  : Instrument strength sweep figure
- [`plot_model_selection()`](https://seantbresnahan.com/iconic/reference/plot_model_selection.md)
  : Model selection workflow figure
- [`plot_nc_coverage_comparison()`](https://seantbresnahan.com/iconic/reference/plot_nc_coverage_comparison.md)
  : NC coverage comparison figure
- [`plot_nc_validity_diagnostics()`](https://seantbresnahan.com/iconic/reference/plot_nc_validity_diagnostics.md)
  : NC validity diagnostics figure (5 panels)
- [`plot_pleiotropy_sweep()`](https://seantbresnahan.com/iconic/reference/plot_pleiotropy_sweep.md)
  : Pleiotropy sweep figure
- [`plot_power()`](https://seantbresnahan.com/iconic/reference/plot_power.md)
  : Plot detection rate (power) vs a swept parameter
- [`plot_prospective_analysis()`](https://seantbresnahan.com/iconic/reference/plot_prospective_analysis.md)
  : Prospective analysis figure
- [`plot_sensitivity_heatmap()`](https://seantbresnahan.com/iconic/reference/plot_sensitivity_heatmap.md)
  : Heatmap of a sensitivity metric across the scenario grid
- [`plot_type1_boxplot()`](https://seantbresnahan.com/iconic/reference/plot_type1_boxplot.md)
  : Type I error boxplot per method
- [`plot_type1_error()`](https://seantbresnahan.com/iconic/reference/plot_type1_error.md)
  : Bar chart of Type I error rates
- [`plot_type1_vs_conf()`](https://seantbresnahan.com/iconic/reference/plot_type1_vs_conf.md)
  : Type I error rate vs confounding strength

## Data constants

Package data constants for estimator display (colours and ordering).

- [`iconic_method_colors`](https://seantbresnahan.com/iconic/reference/iconic_method_colors.md)
  : Colour palette for iconic methods
- [`iconic_method_order`](https://seantbresnahan.com/iconic/reference/iconic_method_order.md)
  : Default method display order

## S3 print & summary methods

Print and summary methods for the result objects returned above.

- [`print(`*`<iconic_confounding>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_confounding.md)
  : Print method for iconic_confounding objects
- [`print(`*`<iconic_data>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_data.md)
  : Print method for iconic_data objects
- [`print(`*`<iconic_diagnosis>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_diagnosis.md)
  : Print method for iconic_diagnosis objects
- [`print(`*`<iconic_feature_texture>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_feature_texture.md)
  : Print method for iconic_feature_texture objects
- [`print(`*`<iconic_gan>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_gan.md)
  : Print method for iconic_gan objects
- [`print(`*`<iconic_prospect>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_prospect.md)
  : Print method for iconic_prospect objects
- [`print(`*`<iconic_recommendation>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_recommendation.md)
  : Print method for iconic_recommendation objects
- [`print(`*`<iconic_sensitivity>`*`)`](https://seantbresnahan.com/iconic/reference/print.iconic_sensitivity.md)
  : Print method for iconic_sensitivity objects
- [`summary(`*`<iconic_diagnosis>`*`)`](https://seantbresnahan.com/iconic/reference/summary.iconic_diagnosis.md)
  : Summary method for iconic_diagnosis objects
- [`summary(`*`<iconic_prospect>`*`)`](https://seantbresnahan.com/iconic/reference/summary.iconic_prospect.md)
  : Summary method for iconic_prospect objects
- [`summary(`*`<iconic_recommendation>`*`)`](https://seantbresnahan.com/iconic/reference/summary.iconic_recommendation.md)
  : Summary method for iconic_recommendation objects
- [`summary(`*`<iconic_sensitivity>`*`)`](https://seantbresnahan.com/iconic/reference/summary.iconic_sensitivity.md)
  : Summary method for iconic_sensitivity objects
