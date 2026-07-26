# iconic 0.9.6

## Documentation cleanup

- **biocViews**: removed `Bayesian` from the `biocViews` field.

- **Description**: tightened the `DESCRIPTION` Description from ~1,300
  words to ~180 words, removing version-history prose and internal
  implementation details. The new Description describes the package's
  purpose, the eight estimators, the diagnose/estimate/sensitivity/
  recommend/prospect workflow, survival outcome support, and the torch
  dependency.

- **Coauthor and version-control markers**: removed all coauthor comment
  markers (`JYH #NNN`), review-process markers (`critique #N of the QED
  review`), `TODO(v1.0)` markers, and inline version-history annotations
  (`v0.x.y`) from all documentation, man pages, roxygen blocks, the
  vignette, README, inst scripts, and test files. Substantive
  documentation content was preserved; only changelog-style commentary
  and internal review references were removed.

- **Vignette**: corrected an inaccurate statement that the package falls
  back to a multivariate-normal texture model when `torch` is not
  installed. `torch` is a hard dependency; the vignette now states this
  explicitly.

- **Tests**: removed the source-level test asserting the presence of
  `TODO(v1.0)` markers (the markers were removed in this release).

# iconic 0.9.5

## CRAN / Bioconductor readiness

- **Documentation**: fixed non-ASCII characters in R source files
  (em-dashes, superscripts, box-drawing characters) that triggered
  `R CMD check` NOTEs. Fixed malformed Rd `\eqn{}` braces and
  cross-references in `composite_p_value`, `fit_pgc_mediation2`,
  `sample_texture`, and related help pages. Documented the previously
  undocumented `feat_cor` argument in `run_simulation`,
  `run_mediation_sim`, and `run_null_mediation_sim`. Removed an
  internal roxygen block for the `%||%` operator that produced an
  invalid `\name{}` entry.

- **Global variables**: declared all ggplot2 non-standard-evaluation
  variables (`.data`, `bias`, `estimate`, `method`, etc.) via
  `utils::globalVariables()` and extended the `stats` import list, so
  `R CMD check` no longer reports undeclared global variables.

- **Removed unused code**: deleted the unused, exported
  `cuda_available_safe()` helper. Removed the top-level
  `generate_manuscript_figures.R` script (preserved in the manuscript
  repository) that triggered a non-standard top-level file NOTE.

- **Tests**: updated stale tests in `test-refactoring.R` and
  `test-v0.9.2-revisions.R` to match function contracts that changed
  in v0.8.4 / v0.9.3 / v0.9.4 (`.analyze_feature()` now passes the full
  W matrix; `nc_completeness_capture()` returns `capture_R2` /
  `capture_pvalue` / `capture_verdict`; Rd alias tests use
  `system.file()` instead of relative paths). All 17 test files now
  pass with 0 failures.

- **Vignette**: rewrote `iconic-walkthrough.Rmd` (1582 -> 264 lines)
  to build in ~1.5 minutes. The new vignette trains the GAN texture
  model once and reuses it across `iconic_sensitivity()` and
  `iconic_prospect()`, and drops the simulation benchmarking section
  (retained in the manuscript). Keeps the model-selection workflow,
  GAN training, sensitivity analysis, confounding inference,
  recommendation, prospective analysis, and survival outcomes.

- **README**: rewrote `README.md` (540 -> 100 lines) as a minimal
  CRAN/Bioconductor-style README with package description,
  installation, a tight getting-started code block, a key-functions
  table, and pointers to the vignette and manuscript.

# iconic 0.9.4.1

## Bug fixes

- **Missing `fit_pgc_mediation_surv()` estimator**: the survival
  mediation driver dispatched both "PGC" and "PGC2" to
  `fit_pgc_mediation2_surv()` (path-specific bridges), making them
  identical. In the continuous case, PGC uses
  `fit_pgc_mediation()` (single W-matrix bridge) and PGC2 uses
  `fit_pgc_mediation2()` (path-specific W1/W2 bridges). Added the
  missing `fit_pgc_mediation_surv()` estimator — the survival analogue
  of `fit_pgc_mediation()` — and updated the dispatch in `estimate.R`
  so "PGC" routes to it. PGC and PGC2 now produce distinct results.

- **`fit_pgc_mediation2_surv()` numerical instability**: the
  path-specific bridge estimator could produce catastrophically large
  coefficient estimates (RMSE 4-150 in benchmarks) when the generated
  regressors (Z_hat, M_hat, W_hat_Z, W_hat_M) became collinear or
  extreme in the Cox partial likelihood. Added a numerical-stability
  guard that rejects estimates whose absolute value exceeds 10 on
  either scale (log-HR > 10 corresponds to HR > 22,000; RMST shift > 10
  time units is equally implausible for standardised data), returning
  NA instead. The same guard is applied to `fit_pgc_mediation_surv()`
  for consistency.

- **`.fit_surv_outcome_stage()` coefficient-name flexibility**: the
  shared outcome-stage helper hardcoded `"Z_hat"` and `"M_hat"` as the
  coefficient names to extract. Generalised to accept `nde_name` and
  `med_name` parameters (defaulting to `"Z_hat"` / `"M_hat"`), so
  estimators that use non-instrumented regressors (e.g.,
  `fit_pgc_mediation_surv()` uses `Z` and `M`) can reuse the helper.

## Benchmark script fixes (generate_manuscript_figures.R)

- **Total-effect truth correction**: `SURV_TRUE_TE_LOGHR` was set to
  `SURV_TAU + SURV_ALPHA_M * SURV_BETA_M` (0.40), but the total-effect
  sweep DGP uses `alpha_M = 0, beta_M = 0` (no mediation path), so the
  true total effect is `beta_Z = 0.25`. Corrected to
  `SURV_TRUE_TE_LOGHR <- SURV_TAU`.

- **Survival mediation DGP alignment**: the survival mediation
  benchmark did not set `separate_U`, `mo_confounding`, or `rho_G1`,
  using defaults (single U, no M-O confounding) instead of the
  continuous benchmark settings (`MED_SEPARATE = TRUE`,
  `MED_MO_CONF = 0.8`, `MED_RHO_G1 = 0.3`). Added these parameters and
  passed `W1`/`W2` to `iconic_data()` so path-specific estimators
  (PGC2, PGC2Gm) exercise their intended two-bridge design. Updated
  pre-computed RMST truth values accordingly (NDE: -0.7100 -> -0.7284,
  NIE: -0.5573 -> -0.6198).

- **PGC2 added to mediation figure**: `surv_med_methods_show` now
  includes "PGC2" alongside "PGC" and "PGC2Gm", since PGC and PGC2 are
  now distinct estimators.

# iconic 0.9.4

## New features

- **Time-to-event (survival) outcome support**: `iconic_data()` now
  accepts `outcome_type = "survival"` with `surv_time` and `surv_event`
  arguments, enabling causal inference with genetic instruments and
  negative controls for time-to-event outcomes (e.g., overall survival
  in cancer cohorts). All estimators are available on two effect scales:

  - **Cox log-hazard-ratio** (`effect_scale = "loghr"`, default): the
    outcome stage uses `survival::coxph()` via two-stage predictor
    substitution (2SPS). First-stage regressions (Z ~ G, M ~ Z_hat + Gm,
    bridges) remain OLS because Z, M, and W are continuous; only the
    outcome stage switches to Cox regression. The log-HR is
    non-collapsible, so the product-of-coefficients NIE = alpha * beta
    on this scale is an approximation.

  - **Restricted mean survival time** (`effect_scale = "rmst"`): the
    outcome stage uses OLS on RMST pseudo-observations (leave-one-out
    jackknife of Kaplan-Meier RMST, Graw et al. 2009). The RMST scale is
    collapsible, so the product-of-coefficients NDE/NIE decomposition is
    exact. Default truncation time `tau` is the 90th percentile of
    observed follow-up.

- **Survival estimators**: `fit_unadj_surv()`, `fit_direct_surv()`,
  `fit_iv2sls_surv()`, `fit_pgc_surv()`, `fit_coca_surv()` (total-effect);
  `fit_unadj_mediation_surv()`, `fit_direct_mediation_surv()`,
  `fit_iv2sls_mediation_surv()`, `fit_iv2sls_mediation2_surv()`,
  `fit_pgc_mediation2_surv()`, `fit_coca_mediation_surv()` (mediation).
  All are exported and documented. COCA is structurally incompatible
  with survival outcomes (it regresses W on Y, placing the outcome on
  the RHS, impossible with a Surv object) and returns NA with an
  informative `reason` attribute.

- **Survival simulation DGP**: `generate_toy_data()` and
  `run_single_iteration()` gain `outcome_type`, `surv_h0`,
  `surv_event_frac`, and `surv_censor_rate` arguments. When
  `outcome_type = "survival"`, the continuous linear predictor is
  converted to time-to-event via an exponential proportional-hazards
  model with independent exponential censoring, targeting ~50-60%
  event fraction. The `true_total` / `true_NDE` / `true_NIE` are then
  on the Cox log-HR scale.

- **Survival sensitivity and prospective analysis**:
  `gan_sensitivity()`, `gan_mediation_sensitivity()`,
  `iconic_sensitivity()`, and `iconic_prospect()` all gain
  `outcome_type`, `effect_scale`, and survival DGP parameters. When
  survival is active, estimation uses `iconic_estimate()` with the
  survival drivers instead of `run_methods()` /
  `run_mediation_methods()`.

## Bug fixes

- **`switch()` factor dispatch in survival mediation driver**: the
  mediation survival driver used `switch(method, ...)` where `method`
  came from `expand.grid()`, which converts strings to factors.
  `switch()` treats factors as integers (positional matching), silently
  dispatching methods to the wrong estimators (e.g., PGC calling
  IV2SLS2, COCA calling IV2SLS). Fixed by wrapping in `as.character()`
  and adding `stringsAsFactors = FALSE` to `expand.grid()`.

- **`.find_tipping_points()` NA/NaN handling**: the tipping-point
  detection in `iconic_sensitivity()` crashed with
  `missing value where TRUE/FALSE needed` when any method's bias was
  entirely NA/NaN (e.g., COCA for survival outcomes). Fixed by guarding
  with `is.finite()` before threshold comparisons. Also fixed
  `max(..., na.rm = TRUE)` returning `-Inf` with a warning for all-NA
  methods, replaced with `.safe_max_abs()` returning `NA`.

- **GAN texture model for survival outcomes**: `.resolve_gan()` attempted
  to auto-train a GAN from `data$Y`, which is absent for survival
  outcomes. Now returns `NULL` (default texture) when
  `data$outcome_type == "survival"`.

# iconic 0.9.3.1

## Bug fixes

- **Panel instrument-strength scalar coercion**: `infer_confounding()`,
  `iconic_sensitivity()`, and `iconic_prospect()` crashed with
  `'length = N' in coercion to 'logical(1)'` when the mediator instrument
  (Gm) was a panel (one F-statistic per gene). The diagnosis object
  stores `F_Gm` as a vector, but `infer_confounding()` and the `phi`
  calibration in `iconic_sensitivity()` used scalar comparisons
  (`is.na(F_Gm) || F_Gm < 10`, `F_gm >= 100`). The prospective summary
  builder in `iconic_prospect()` passed the vector to `sprintf("%.1f",
  ...)`. All three sites now collapse to `F_Gm_median` (a scalar) when
  `length(F_Gm) > 1`, using the pre-existing `F_Gm_median` field from
  `iconic_diagnose()`. The same fix is applied to `F_G` for symmetry.

  Affected files: `R/infer_confounding.R`, `R/model_sensitivity.R`,
  `R/prospect.R`.

# iconic 0.9.2

This release addresses coauthor (JYH) comments on the proximal-completeness
condition, negative-control validity screens, the recommendation criterion,
and simulation reporting. All changes are backward-compatible (new arguments
carry defaults preserving v0.9.1 behavior).

## New features

- **Covariance-capture completeness test** (`nc_completeness_capture()`):
  operationalizes the proximal completeness condition as whether the
  negative-control panel captures the confounder covariance (incremental R^2
  of W for the outcome, with a permutation null), rather than only a proxy
  count. `nc_completeness_check()` now returns both the dimensional verdict
  (dim(W_valid) vs k) and the capture verdict, plus a composite. Addresses
  JYH comments that completeness is about covariance captured, not proxy
  number.

- **Magnitude-based A1 screen**: `nc_validity_screen()` gains a
  `criterion = c("fdr", "magnitude", "both")` argument and a
  `magnitude_threshold` (default 0.10). The magnitude branch flags controls
  by partial-correlation size, not significance, addressing the concern that
  FDR hypothesis tests are underpowered and "always significant via U."
  Returns per-branch verdicts (`verdict_fdr`, `verdict_magnitude`).

- **COCA A2 exemption**: COCA no longer requires the A2 (instrument-
  independence) screen for eligibility, only A1 + completeness. A2 remains
  required for IV2SLS, PGC, IV2SLS2, PGC2, PGC2Gm. The eligibility table
  gains an `a2_required` column.

- **Per-scenario recommendations**: `iconic_recommend()` now returns a
  `$per_scenario` table (top estimator per sensitivity-surface cell) in
  addition to the global recommendation, and a `criterion` argument
  (`"minimax_bias"`, `"ci_coverage"`, `"combined"`). The sensitivity surface
  gains `NDE_coverage` and `NIE_coverage` columns (Wald CI coverage from the
  delta-method SE). Tier-A label changed to "identified under stated
  assumptions."

- **CI coverage, %bias, and mean SE in benchmarks**: the simulation
  summarization now reports `NDE_coverage`, `NIE_coverage`,
  `NDE_pct_bias`, `NIE_pct_bias`, `NDE_mean_se`, `NIE_mean_se`. The
  mediation estimators and `iconic_estimate()` gain `se_method =
  c("delta", "bootstrap")` (default `"delta"`) and `n_boot` (default 500)
  for optional bootstrap SE.

- **Scenario manifest**: new `scenario_manifest()` helper exposes the
  ground-truth estimands and the swept vs. fixed parameter ranges;
  `iconic_sensitivity()` and `gan_mediation_sensitivity()` add a `$manifest`
  field to their return value.

- **U/W strength heterogeneity**: `generate_toy_data()` and
  `run_single_iteration()` gain `u_strength` (per-confounder scaling) and
  `w_coverage_profile` (per-control coverage) arguments, allowing some
  confounders/proxies to be stronger than others. Exposed in the returned
  dataset.

- **Optional U-deconfounding in residual correlation**:
  `.residual_correlation()` and the texture-training entry points gain
  `residualize_on = c("ZC", "ZCW")` (default `"ZC"`). Under `"ZCW"`, the
  residual correlation partials out the U-signature captured by W, reducing
  the double-counting of U's cross-feature signature. A
  `$double_count_estimate` diagnostic is exposed.

- **Defaults enumeration and `allow_no_proxy`**: `infer_confounding()`,
  `iconic_sensitivity()`, `iconic_prospect()` documentation gains a
  `@section Defaults` block. `iconic_prospect()` and `iconic_diagnose()`
  gain `allow_no_proxy` (default TRUE): when FALSE, they error if no IV and
  no NC are supplied.

## Documentation

- `iconic_sensitivity()` and `iconic_prospect()` gain `@aliases`
  (`effect_decomposition_bias_sweep`, `bias_reduction_prospective`) and
  updated `@title`/`@section` blocks for emphasis, without renaming the
  functions (backward compatible).

## Deferred (placeholders for the future supplement / companion paper)

- `# TODO(v1.0)` markers added in `R/gan.R`, `R/feature_texture.R`, and
  `R/model_sensitivity.R` for the future move of generative-pipeline and
  technical-sensitivity detail to the supplement.
- The biobank-scale / cross-ancestry case study is deferred to a companion
  paper; no package change.

# iconic 0.9.1

## New features

- **Panel-level instrument-strength distributions**: `.check_instrument_strength()`
  now computes the partial F for every mediator in the panel (not just the
  first transcript) for both G and Gm. The full distributions are exposed at
  `diag[["instrument_strength"]][["F_G"]]` and `diag[["instrument_strength"]][["F_Gm"]]`
  (numeric vectors, one entry per instrument/mediator). `iconic_diagnose()`
  gains `g_threshold = list(E, R)` and `gm_threshold = list(E, R)` arguments:
  an instrument-dependent method is eligible if at least E fraction of
  transcripts have F >= R. Default NULL preserves the legacy scalar behavior
  (median F vs `min_f`).

- **run_all flag**: `iconic_estimate()` gains `run_all = FALSE` and `min_f`
  arguments. `run_all = TRUE` overrides diagnosis eligibility, running every
  method whose required data exists. `min_f` (default NULL) is inherited from
  `diagnosis$min_f` when a diagnosis is supplied and not explicitly set,
  otherwise defaults to 10. It is passed through to the per-transcript
  weak-instrument gate (previously hardcoded to 10 inside the estimator
  functions and never configurable from the top-level API).

## Bug fixes

- **Parallel progress reporting**: `.parallel_lapply()` now reports
  percentage-complete progress during parallel execution (Unix `mclapply`
  and Windows PSOCK), not just a start message. Tasks are split into ~10
  chunks; a progress message is printed after each chunk completes.
  Previously, parallel runs printed only a single start line and produced
  no intermediate output until completion.

- **Sequential timing bug**: The final "100%" completion line in
  sequential mode now reports total elapsed time from task start, not
  just elapsed time since the last 30-second update. A separate
  `start_time` is retained alongside the throttle timer.

- **Unified progress cadence**: Both sequential and parallel paths now
  report at ~10% increments, replacing the previous 30-second throttle
  for sequential mode. Progress messages include cumulative percentage,
  task count, and elapsed seconds.

# iconic 0.9.0

## Major changes

- **Hybrid GAN + copula texture model**: The generative texture model is now a hybrid of two independently trained components. A sample-level torch GAN learns the joint distribution of per-sample exposure, outcome, mediator level, and encoded covariates. A new feature-level Gaussian copula model (`train_feature_texture()` / `sample_feature_texture()`) learns the full joint distribution of the mediator (M) panel — marginal distributions and cross-feature dependence — directly from the user's data.

- **Gaussian copula for mediator panel**: Each mediator feature's marginal is fitted as an empirical CDF by default, with a parametric fallback (normal, log-normal, gamma, or beta, selected by AIC) used when the parametric fit passes a Kolmogorov-Smirnov goodness-of-fit test at p > 0.05. The user may override this automatic selection via the `marginal_method` parameter ("auto", "empirical", "parametric"). Cross-feature dependence is captured by the Gaussian copula correlation matrix of the normal-transformed scores. During simulation, copula draws replace the parametric Gaussian noise in the mediator structural equations, preserving realistic marginal shapes and correlations from the user's cohort while maintaining the closed-form ground truth.

- **torch is now a hard dependency**: The MVN fallback has been removed entirely. `torch` has moved from Suggests to Imports. The GAN component requires torch; `train_gan_on_real_data()` errors if torch is unavailable.

- **W panel unchanged**: The negative-control (W) panel continues to use the pluggable `nc_model` with residual correlation matrices. The copula replaces the residual correlation approach only for the mediator panel.

## New functions

- `train_feature_texture(M_matrix, marginal_method)` — trains the Gaussian copula + marginals model
- `sample_feature_texture(feature_texture, n_samples, n_features)` — draws from the copula, centered and scaled
- `print.iconic_feature_texture()` — summary print method

## Modified functions

- `load_real_input_data()` — gains `marginal_method` parameter; trains feature_texture from M_matrix
- `train_gan_on_real_data()` — gains `feature_texture` argument; errors without torch
- `run_single_iteration()` — uses copula-based mediator texture injection (`.iteration_mediator_texture_features()`)
- `iconic_sensitivity()` / `iconic_prospect()` — pass `feature_texture` through the auto-training pipeline
- `plot_gan_diagnostics()` — gains `M_matrix` argument for feature-level marginal comparison panels

## Bug fixes

- Fixed `test-refactoring.R` to use `iconic:::` namespace prefix for internal function calls

# iconic 0.8.4

## Major changes

- **Full negative-control panel for all estimators**: DIRECT, IV2SLS, and IV2SLS2 now use the full W panel (all q negative-control features) in their regression formulas, rather than matching a single NC feature to each outcome feature. This improves confounding adjustment and aligns the code with the model specification described in the manuscript. A new internal `.expand_w()` helper handles both vector and matrix W inputs. COCA continues to use `rowMeans(W)` as a scalar summary (required by its ratio-estimator form). PGC and PGC2 already used the full panel via bridge regressions.

- **Per-mediator genetic instruments**: `generate_toy_data()` and `run_single_iteration()` now accept an `n_mediators` parameter (default 1, backward compatible). When `n_mediators > 1`, each mediator M_m has its own independent genetic instrument Gm_m, and both M and Gm are returned as `n_mediators x n` matrices. The mediation estimation pipeline estimates each mediator one at a time (M_vec = M[m, ]), mirroring the `iconic_estimate()` behavior. `run_mediation_sim()` and `sweep_mediation_param()` also accept `n_mediators`.

- **Vignette migration**: The walkthrough document has been moved from the package root to `vignettes/iconic-walkthrough.Rmd` with a proper `VignetteIndexEntry` and `knitr::rmarkdown` engine. The output format is now `rmarkdown::html_vignette`. `VignetteBuilder: knitr` added to DESCRIPTION; `knitr` and `rmarkdown` added to Suggests.

## Minor changes

- Fixed Gm dimension validation in `iconic_data()` and `validate_iconic_data()` to correctly check `nrow(Gm)` against `n_mediators` (previously checked `ncol(Gm)`, which always equaled n after transposition)
- Fixed Gm extraction in `estimate.R` to use `nrow(data$Gm) == nm` instead of `ncol(data$Gm) == nm`
- Updated roxygen documentation throughout to reflect matrix M/Gm support and full W panel usage

# iconic 0.8.3

## Major changes

- **Cross-platform parallelization**: All computationally intensive functions now accept an `n_cores` argument for multi-core execution. Uses `parallel::mclapply` on Unix (Linux, macOS) and PSOCK clusters (`makeCluster`/`parLapply`/`stopCluster`) on Windows. No new dependencies required.

- **Progress milestones**: Parallelized functions emit `message()`-based progress reports — start lines announcing task count and core count, and incremental progress every ~10% of tasks in sequential mode.

## Functions with new `n_cores` argument

- `iconic_estimate()` — parallelizes per-feature and per-mediator × feature estimation
- `iconic_sensitivity()` — parallelizes replicate iterations within each grid cell
- `iconic_prospect()` — parallelizes Phase 1 and Phase 2 replicate iterations
- `iconic_diagnose()` — parallelizes NC validity screens (A1, A2, A2')
- `infer_confounding()` — parallelizes NC-coverage loop (omega) and k permutation analysis
- `sweep_nc_validity()` — parallelizes replicate iterations in all four diagnostic panels
- `sweep_instrument_strength()` — parallelizes replicate iterations per grid point
- `analyze_methods_robust()` — parallelizes per-feature estimation
- `analyze_mediation_robust()` — parallelizes per-feature mediation estimation
- `run_methods()` / `run_mediation_methods()` — parallelizes per-feature estimation (internal drivers)
- `nc_validity_screen()` — parallelizes per-control-feature regression
- `nc_independence_check()` / `nc_independence_check_gm()` — parallelizes per-control-feature partial correlation
- `nc_completeness_check()` — delegates `n_cores` to underlying screens

## Internal changes

- Rewrote `.parallel_lapply()` dispatcher with cross-platform support (mclapply Unix, PSOCK Windows) and optional `progress` parameter for milestone messages
- Flattened nested `for(m) for(f)` mediation estimation loop into single `expand.grid` task list for one parallel call
- Default `n_cores = 1` preserves exact backward compatibility with v0.8.2
- No nested parallelism: simulation functions parallelized at the replicate level call `run_methods()` with `n_cores = 1` inside each worker

# iconic 0.8.2

- Initial CRAN submission
