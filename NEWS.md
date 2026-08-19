# iconic 0.99.3

## New features

- **Binary outcome support.** `outcome_type = "binary"` is now
  available throughout the workflow: `iconic_data()` accepts a 0/1
  outcome vector (validated as dichotomous, stored unscaled as
  `Y_bin`), and `iconic_estimate()`, `iconic_sensitivity()`,
  `iconic_prospect()`, and `run_single_iteration()` all dispatch to
  binary-specific estimators. `generate_toy_data()` and
  `run_single_iteration()` gain a `bin_prev` argument (default 0.5)
  controlling the outcome prevalence of the simulated binary DGP.

- **Two effect scales for binary outcomes** via the `effect_scale`
  argument: `"logor"` (default; logistic two-stage
  predictor-substitution outcome stage, mirroring the Cox log-HR
  convention for survival) and `"riskdiff"` (linear probability model
  outcome stage, giving an exactly collapsible NDE/NIE decomposition).
  Inert scale requests are remapped with an informative message
  (e.g. `effect_scale = "loghr"` on binary data maps to `"logor"`).

- **New binary estimators** (`R/bin_estimators.R`,
  `R/bin_mediation.R`): `fit_unadj_bin()`, `fit_direct_bin()`,
  `fit_iv2sls_bin()`, `fit_pgc_bin()`, and the mediation variants
  `fit_unadj_mediation_bin()`, `fit_direct_mediation_bin()`,
  `fit_iv2sls_mediation_bin()`, `fit_iv2sls_mediation2_bin()`,
  `fit_pgc_mediation_bin()`, `fit_pgc_mediation2_bin()`. First stages
  remain OLS; only the outcome stage changes (binomial GLM or linear
  probability model). The partial-F weak-instrument gate and the
  collinear-panel fallback behave as for continuous and survival
  outcomes.

- **COCA is excluded for binary outcomes** (`fit_coca_bin()` and
  `fit_coca_mediation_bin()` return `NA` with a `reason` attribute, the
  same contract as for survival): COCA's ratio identification assumes a
  linear structural outcome model, which neither the log-OR nor the
  risk-difference binary outcome stage satisfies.

## Behavioural notes

- `infer_confounding()` now errors informatively for binary outcomes:
  gap-based calibration relies on a continuous outcome scale. Use
  `confounding = "manual"` or the default confounding grid in
  sensitivity/prospective analyses of binary data.
- The generative texture model (torch GAN) is unavailable for binary
  outcomes because there is no continuous outcome margin to texture;
  `iconic_sensitivity()` falls back to the default texture with an
  informative message.

## Bug fixes

- **Locale/device-safe plot labels (completes the 0.99.2 fix).**
  Remaining multi-byte glyphs that are drawn when the package examples
  run are now portable across locales and graphics devices. In
  `plot_nc_validity_diagnostics()` the internal colour-scale keys still
  used Unicode arrow and superscript-two string literals; these are now
  ASCII (`->`, `R2`). The displayed legend label was already a plotmath
  expression, so the figure is visually unchanged. In
  `plot_pleiotropy_sweep()` and `plot_instrument_strength_sweep()` the
  axis titles carried a Unicode arrow / minus sign as string literals;
  these are now plotmath expressions, which render the real glyphs
  through the symbol font on any device -- including the base `pdf()`
  device used by `R CMD check` -- without `mbcsToSbcs` conversion
  errors. No code line in `R/figures.R` now contains a non-ASCII or
  `\uXXXX`-escaped character.

# iconic 0.99.2

## Bug fixes

- **Portable plotmath labels in `plot_nc_validity_diagnostics()`.** The
  new panels D/E introduced in 0.99.1 drew omega / Utilde / R2 glyphs as
  Unicode string literals, which fail with `mbcsToSbcs` conversion
  errors when the plot is printed on devices/locales without multi-byte
  glyph support (observed on the macOS check runner). All drawn text in
  those panels now uses plotmath expressions, which render through the
  symbol font on any device. The saved figure (cairo_pdf) is visually
  unchanged.

# iconic 0.99.1

## Figure revision

- **Supplementary Figure S6 revision** (`R/figures.R`): the A3
  dimensional heatmap (PGC bias over the `n_valid` x `k` grid) is no
  longer plotted by default; `plot_nc_validity_diagnostics()` now
  assembles a 5-panel figure — (A) A1 screen, (B) A2 screen, (C) A2'
  screen, (D) A3 covariance-capture versus true negative-control
  coverage omega (with the permutation-null mean and the fraction of
  replicates with permutation p < 0.05 on a secondary axis), and (E) A3
  support R2(Utilde | W) versus omega.
- **`sweep_nc_validity()`** gains `omega_grid` (default 0 to 1),
  `n_perm` (default 200, matching `iconic_diagnose()`),
  `cs_confounders` (default 2), and `include_a3_grid` (default FALSE)
  arguments, and returns a new per-replicate `panel_capture_support`
  data frame (capture R2, permutation p-value, permutation-null mean,
  support R2, and the fraction of controls adding unique coverage). The
  legacy A3 dimensional grid sweep remains available via
  `include_a3_grid = TRUE`.
- The figure is saved with `cairo_pdf` so the omega / Utilde / arrow
  glyphs render correctly; the default figure size is now 10 x 6.5
  inches.

# iconic 0.99.0

## New features

- **Exposure-instrument helpers** (`R/instruments_exposure.R`): a
  GWAS-to-instrument workflow distilled from the case-study scripts.
  `qc_gwas_sumstats()` standardises and QC-filters GWAS summary statistics
  (column aliases, missing/invalid rows, IQR-based extreme-beta filter,
  ambiguous-strand removal, SD-ratio check); `build_prs_ldpred2()` fits
  LDpred2-auto weights from summary statistics plus an LD reference and
  scores a target panel; `score_pgs_panel()` applies a published
  per-variant weight panel (e.g. PGS Catalog) to a dosage matrix with
  automatic allele-flip handling; `check_instrument_strength()` computes
  the first-stage partial F / partial R2 and flags weak instruments.

- **Mediator-instrument helpers** (`R/instruments_mediator.R`):
  `call_cis_eqtls()` calls cis-eQTLs for a gene panel from genotype and
  expression matrices; `build_mediator_instruments()` trains per-gene
  elastic-net cis instruments and returns the genetically predicted
  mediator panel (`Gm`) with QC.

- **Negative-control helpers** (`R/negative_controls.R`): `beta_to_m()`
  logit-transforms methylation beta values with clipping;
  `residualize_matrix()` residualises a feature matrix on covariates
  (chunked); `build_w_pcs()` builds a negative-control panel of principal
  components; `apply_fusion_weights()` applies FUSION/TWAS weights to a
  dosage matrix.

- **SummarizedExperiment interop**: `as_iconic_data()` is now an S3 generic
  with a `SummarizedExperiment` method that pulls the outcome panel from an
  assay and sample-level fields (exposure, instruments, negative controls,
  covariates, survival endpoints) from `colData`, in addition to the
  existing `iconic_data` / `load_real_input_data` paths.

- **New vignette**: `vignette("iconic-instruments")` walks through the
  instrument-construction and negative-control workflow end to end.

## Bioconductor submission preparation

- Version bumped to 0.99.0; `biocViews` added (StatisticalMethod, Genetics,
  MultipleComparison, Regression, Transcriptomics, RNASeq, Survival);
  `BiocStyle` used for all vignettes, each gaining Introduction,
  Installation, and Session information sections.
- `withr` moved to Imports; `bigsnpr`, `bigstatsr`, `glmnet`, `irlba`,
  `SummarizedExperiment`, `S4Vectors`, and `BiocStyle` listed under
  Suggests and used conditionally.
- Code hygiene for CRAN/BiocCheck: `set.seed()` replaced with
  `withr::local_seed()`, `cat()` with `message()`, `seq_len()`-style
  indexing, and signaler calls no longer use `paste()` or message keywords.
- Every exported function now has a runnable `@examples` block and a
  documented `@return`; `\dontrun{}` removed (torch-dependent examples are
  guarded by `check_torch_setup()`) and `\donttest{}` minimised.
- Added testthat coverage for all new functions.

# iconic 0.9.9.3

## Documentation and housekeeping

- **Documentation cleanup.** Removed version-control/changelog language from
  roxygen help pages, comments, and defunct-argument error messages (the
  defunct-argument traps for `w`, `Z`, `Z_matrix`, `beta_Z`, and `separate_U`
  are retained and still error informatively). Generalized case-study-specific
  motivating examples in the help pages to domain-neutral omics language.
  Clarified the `iconic_recommend()` documentation to describe the composite
  robustness rule directly. No functional or API changes.

# iconic 0.9.9.2

## New features

- **Data-driven composite recommendation in `iconic_recommend()`.** The
  headline recommendation is now a single data-driven composite rather than
  the NDE-only ranking. For each estimator the composite is the
  *worst-estimand* robustness (`min(score_NDE, score_NIE)`), so an estimator
  is only as trustworthy as its weakest estimand, multiplied by a *confidence
  factor* derived from the graded diagnostic verdict for the assumptions that
  estimator depends on. Bridge-dependent estimators (DIRECT, COCA, PGC, PGC2,
  PGC2Gm) are discounted when path completeness is borderline or weak-capture;
  instrument-only estimators (IV2SLS, IV2SLS2) are not. Structurally naive
  estimators (UNADJ, DIRECT) are demoted below eligible instrument/NC
  estimators. The discount is exposed via the new `completeness_penalty`
  argument (default `c(satisfied = 1.0, borderline = 0.7, "weak-capture" =
  0.5, "under-identified" = 0)`), and the ranking gains `composite`,
  `confidence_mult`, and `final_score` columns. Per-estimand recommendations
  are retained as `$recommended_NDE` and `$recommended_NIE`.

- **`min_f` / `g_threshold` / `gm_threshold` pass-through.** When
  `diagnosis` is `NULL`, `iconic_recommend()` now runs `iconic_diagnose()`
  with the caller's thresholds instead of silently using the defaults, so a
  non-default `min_f` is honoured.

## Bug fixes

- **Threshold-aware requirement labels.** The rationale text previously
  hardcoded `F>=10` regardless of the instrument-strength gate actually
  applied. Labels now interpolate the `min_f` stored in the diagnosis (or the
  `min_f` argument), so the printed requirement matches the gate.

# iconic 0.9.9.1

## Bug fixes

- **Lone path-specific NC panel now derives the pooled `W`.** Supplying only
  `W2` (or only `W1`) without a pooled `W` previously left `W` unset and
  `has_nc = FALSE`, so the single-panel estimators (DIRECT, COCA, PGC) were
  incorrectly ineligible and the NC validity / completeness screens were
  skipped entirely. `iconic_data()` now derives `W` from the lone panel
  (mirroring the existing `W1 + W2 -> W` derivation), so DIRECT / COCA / PGC
  become eligible and the NC screens run. This applies to both continuous
  and survival outcomes. `has_path_nc` still correctly remains `FALSE`, so
  the two-bridge estimators (PGC2, PGC2Gm) stay gated unless the user opts
  in (below).

- **`iconic_sensitivity(confounding = "inferred")` no longer runs the full
  mediator panel.** It previously called `iconic_estimate()` on every
  mediator and passed the result to `infer_confounding()`, bypassing the
  documented random-subset behaviour. It now lets `infer_confounding()` use
  its `max_infer_tasks` random subset (default 50 mediators/features) for
  the gap-based calibration — an unbiased Monte Carlo estimate at a fraction
  of the cost.

- **User-supplied `omega_1` / `omega_2` sweeps take precedence over inferred
  values.** Under `confounding = "inferred"`, the inferred scalar omegas
  previously overwrote user-supplied sweep vectors, collapsing the omega
  facet of the degradation surface to a single point. Inferred values now
  only fill in `omega_1` / `omega_2` / `mo_confounding` when those arguments
  are left at their defaults; explicit values always win. Same fix applied
  to `iconic_prospect()`.

## New features

- **`recycle_lone_panel` argument to `iconic_data()`.** When exactly one of
  `W1` / `W2` is supplied (no pooled `W`), setting
  `recycle_lone_panel = TRUE` uses that lone panel as BOTH path-specific
  bridges (`W1 = W2`), making PGC2 / PGC2Gm eligible. This is the
  shared-panel special case and assumes the single panel is complete for
  BOTH path confounder composites; a warning is emitted, the
  `recycled_lone_panel` flag is recorded on the object, and the
  `iconic_diagnose()` eligibility table annotates PGC2 / PGC2Gm with
  "(shared recycled panel)". Default `FALSE` keeps the previous gating
  (IV2SLS2 remains the more defensible primary estimator when coverage of
  the other path's composite is in doubt).

- **`iconic_sensitivity()` accepts a precomputed `iconic_confounding`
  object.** Passing the result of a prior `infer_confounding()` call as the
  `confounding` argument reuses it as-is (equivalent to `"inferred"`) and
  skips recomputation, avoiding redundant work when `infer_confounding()`
  was already run separately.

# iconic 0.9.9

## Breaking changes

- **IV2SLS2 negative-control augmentation is now path-specific (collider
  fix).** `fit_iv2sls_mediation2()` and `fit_iv2sls_mediation2_surv()` no
  longer take a single pooled negative-control panel `w`. The pooled panel
  `scale((W1 + W2) / 2)` is a common child of the two independent
  confounders (conf_XM and conf_MY) under multi-confounder designs, so
  conditioning on it in all three 2SLS stages opened a collider path and
  made the NDE bias *worse* as NC coverage increased. The estimators now
  take optional path-specific panels:
  - `W1` (proxies the exposure–mediator confounder, X->M path) is added to
    stage 1 (`X ~ G + W1`) only;
  - `W2` (proxies the mediator–outcome confounder, M->Y path) is added to
    stages 2 and 3 (`M ~ X_hat + Gm + W2`, `Y ~ X_hat + M_hat + W2`).
  Either panel may be omitted; with both `NULL` the estimator reduces to
  plain two-instrument 2-stage MR. If `W1` and `W2` are identical they are
  treated as absent (pure MR), because an identical panel is a pooled panel
  in disguise. Passing the old `w` argument is a hard error redirecting to
  `W1`/`W2`. With path-specific panels, higher NC coverage now *improves*
  the NDE estimate.

- **Shared-composite (single-confounder) fallback.** Path-specific
  augmentation is only defined when `W1` and `W2` proxy *distinct*
  confounders. When the two panels are distinct-noise proxies of the *same*
  latent composite (the single-confounder / shared-loading design), their
  column spaces are near-collinear; augmenting stage 1 with `W1` would then
  inject the shared M->Y confounder into `X_hat` and bias the NDE. The
  estimators detect this (leading-PC correlation between the panels) and
  fall back to plain two-instrument 2-stage MR. Genuinely distinct panels
  are unaffected.

- **`iconic_data()` retains a lone path-specific panel.** Supplying only
  `W2` (or only `W1`) without a pooled `W` previously dropped the panel
  silently. A lone panel is now stored so IV2SLS2 can use it for
  path-specific augmentation, while `has_path_nc` remains `FALSE` so the
  two-bridge estimators (PGC2, PGC2Gm) — which require both panels — stay
  ineligible.

# iconic 0.9.8

## Breaking changes

- **`iconic_prospect()` is now sequential-only.** The `n_cores` argument
  has been removed from `iconic_prospect()`; simulation replicates always
  run sequentially. (Parallel replicate execution remains available in
  `iconic_sensitivity()`, `iconic_estimate()`, `iconic_diagnose()`, and the
  simulation drivers via their own `n_cores` arguments.)

- **Quiet by default.** `iconic_prospect()`, `iconic_recommend()`, and
  `iconic_sensitivity()` gain a `verbose` argument (default `FALSE`) that
  gates all progress messages (grid-cell updates, replicate progress, and
  completion banners). Set `verbose = TRUE` to restore the previous
  behaviour.

- **Exposure argument renamed `Z` -> `X` throughout.** The exposure is now
  consistently `X` across the API (`iconic_data(X = ...)`), the structural
  causal model, and all documentation, matching the manuscript's notation.
  Passing the old `Z` argument to `iconic_data()` triggers a hard error with
  an actionable message (`argument `Z` was renamed to `X`; please use
  `X = ...`); there is no silent alias. Related internal names were renamed
  consistently (e.g. the `sweep_instrument_strength()` parameter
  `pi_GZ_grid` is now `pi_GX_grid`; the exposure coefficient is `beta_X`).

- **`separate_U` removed.** The boolean `separate_U` toggle is replaced by a
  general `k`-dimensional unmeasured-confounder space with per-path loading
  vectors `lambda_XM` and `lambda_MY`. The two paths' confounder composites
  are `conf_XM = U %*% lambda_XM` and `conf_MY = U %*% lambda_MY`; distinct
  unit-vector loadings (e.g. `c(1,0)` / `c(0,1)` at `k = 2`) recover the old
  `separate_U = TRUE` behaviour, while identical loadings recover
  `separate_U = FALSE`. Supplying `separate_U` now errors with a pointer to
  the loading-vector parameterization.

## New features

- **Path-specific negative controls and confounder loadings.**
  `generate_toy_data()` / `run_single_iteration()` accept `n_confounders`
  (dimension `k` of `U`), `lambda_XM` / `lambda_MY` (per-path loadings), and
  path-specific negative-control panels `W1` / `W2` with independent coverage
  `omega_1` / `omega_2`. The path-specific proximal estimators PGC2 and PGC2Gm
  use `W1` for the `X -> M` bridge and `W2` for the `M -> Y` bridge.

- **Negative-control coverage (`omega`) sweeps.** `iconic_sensitivity()` and
  `iconic_prospect()` now sweep `omega_1` / `omega_2` (singly or on a grid) so
  estimator performance can be mapped across NC-coverage scenarios, not just
  at a single assumed coverage.

- **`iconic_prospect()` Phase 3: joint exogeneity + coverage robustness
  sweep.** The prospective analysis previously swept instrument strength
  (`gamma_G`, Phase 1) and NC coverage (`omega`, Phase 1/2) but held the
  instruments perfectly exogenous (`rho_G1 = rho_G2 = 0`). A new Phase 3
  sweeps the instrument-exogeneity violations `rho_G1 x rho_G2` jointly with
  NC coverage `omega` (on the diagonal, `omega_1 == omega_2`) at the target
  instrument strength (new arguments `rho_G1_grid`, `rho_G2_grid`, both
  defaulting to `c(0, 0.1, 0.2, 0.3, 0.5)`, `omega_grid_rho` defaulting to
  `c(0.3, 0.7, 1.0)`, and `run_rho_sweep = TRUE` to toggle). The resulting
  degradation surface is returned as `$rho_surface` and is passed to
  `iconic_recommend()` as its `sensitivity` argument, so the recommended
  estimator is chosen by robustness to both imperfect instruments and
  weakening controls rather than by a single-point or eligibility-only
  ranking.

- **`iconic_prospect()` conditional recommendation by collection scenario.**
  The prospective printout no longer reports a single "recommended estimator
  if instruments collected" — which was internally contradictory when the top
  estimator (e.g. COCA) does not use an instrument. A new
  `$recommendation_by_scenario` table reports the best eligible estimator
  under each collection scenario (G1 only / Gm only / G1+Gm / W only /
  W1+W2 / G1+W / G1+Gm+W / G1+W1+W2 / full), mapping each scenario to the
  estimators its data makes available and ranking those by robustness.

- **`iconic_sensitivity()` sweeps NC coverage by default.** `omega_1` and
  `omega_2` now default to `c(0.3, 0.7, 1.0)` (swept on the diagonal when
  identical) instead of the fixed scalar `0.7`, so the degradation surface
  spans both instrument-exogeneity and NC-coverage violations out of the box.

- **`iconic_recommend()` auto-runs the sensitivity suite.** When
  `sensitivity = NULL` and `auto_sensitivity = TRUE` (the default),
  `iconic_recommend()` now calls `iconic_sensitivity()` internally so the
  recommendation is robustness-based by default, not opt-in. The auto-run is
  guarded by torch availability and falls back to eligibility-only ranking
  (with a message) when the torch backend is unavailable or the run fails.
  New arguments `auto_sensitivity`, `rho_G1_grid`, `rho_G2_grid`, `omega_1`,
  `omega_2`, `n_iter_sens`, `gan_epochs`, and `n_cores` control the auto-run.

- **Testability reframe and completeness in `iconic_diagnose()`.** The
  negative-control assumptions are framed as empirically testable projections
  (A1: `W _|_ X | C,U`; A2: `W _|_ G1 | C,U`; A2': `W _|_ Gm | C,U`) plus a
  two-component completeness condition (A3: dimensional check `dim(W_valid)
  >= k` combined with a covariance-capture test). The completeness assessment
  is wired into the diagnosis and the per-estimator eligibility report.

- **Total-effect output.** Estimation surfaces the total effect alongside the
  NDE/NIE decomposition.

## Bug fixes

- **`iconic_prospect()` recommendation no longer crowns UNADJ by list
  position.** Previously the prospective recommendation called
  `iconic_recommend()` with no sensitivity surface and no estimates, so the
  ranking fell back to an eligibility-only stable sort in which `UNADJ` —
  first in the method list — was returned as "recommended" whenever all
  estimators were eligible (the best-case simulated setting). This
  contradicted the function's own Phase 2 bias numbers. Two changes fix it:
  (1) Phase 3 now supplies a real degradation surface, so the recommendation
  is robustness-based; (2) `iconic_recommend()`'s no-sensitivity fallback now
  demotes the unconfoundedness-based estimators (`UNADJ`, `DIRECT`) below the
  instrument/NC-based estimators instead of ranking by list position.

- **`.extract_per_scenario()` robust to degenerate cells.** Per-scenario
  extraction now skips grid cells in which every estimator returned a
  non-finite bias/coverage metric (previously `which.min()` on an all-`NA`
  distance vector returned `integer(0)`, and the `cbind()` of the cell
  metadata with the length-0 result errored with "arguments imply differing
  number of rows"). Within-cell normalization is now computed over finite
  values only.

- **`generate_toy_data()` DGP fix.** With path-specific loadings, the
  non-covered portion of each negative control is now fresh noise rather than
  the other path's confounder composite; previously `W1` was contaminated with
  the `M -> Y` confounder, biasing the PGC purge.

- **`plot_degradation_surface()` colour scales.** The shared colour limit was
  inflated by a single divergent cell, washing out the informative bias range.
  Each estimator panel now uses its own robust (95% quantile, floored) colour
  cap with `scales::squish()` for outliers; the crossover panel preserves a
  direct `|bias|` comparison.

## Behaviour changes

- **Recommendation tiers removed.** `iconic_recommend()` no longer assigns
  estimators to recommendation tiers; it reports per-estimand robustness
  scores (max `|bias|` and CI-coverage distance across the sensitivity grid)
  and ranks on those.

- **Instrument-variable estimators no longer require negative controls.**
  `IV2SLS` (instrument `G`) and `IV2SLS2` (instruments `G` + `Gm`) are
  identified by the instrument(s) alone — relevance, independence, and the
  exclusion restriction — and do not need a negative-control panel `W`. They
  are now eligible and run whenever the required instrument(s) are present,
  with or without `W`; when `W` is supplied it is used as an optional
  proximal augmentation (improving efficiency) rather than as a requirement.
  The proximal bridge estimators (`PGC`, `PGC2`, `PGC2Gm`) and the
  negative-control estimators (`COCA`, `DIRECT`) still require `W`. This
  fixes the case where a dataset with instruments but no negative controls
  was previously left with only `UNADJ`.

- **`infer_confounding()` infers confounding strength on a random subset.**
  When `estimate = NULL`, the confounding-strength (`delta`, the UNADJ–IV2SLS
  gap) and mediator-outcome (`delta_mo`, the IV2SLS–IV2SLS2 gap) parameters
  are averages across the mediator × feature grid. They are now estimated on
  a random subset of at most `max_infer_tasks` mediators and
  `max_infer_tasks` features (default 50) instead of the full panel — an
  unbiased Monte Carlo estimate of the same quantity. This makes
  `iconic_diagnose(k = NULL)`, `iconic_prospect()`, and
  `iconic_sensitivity()` far faster on wide mediator panels (previously the
  internal `iconic_estimate()` call fit every mediator × feature cell). The
  subset sizes are recorded in the returned object as `$inference_subset` and
  noted in the `conf_strength` / `mo_confounding` method strings. Panels
  smaller than the cap are used in full (no behaviour change).

# iconic 0.9.7

## Examples

- **Fixed example failures under `R CMD check --as-cran --run-donttest`.**
  Several examples errored when run by the check (which executes
  `\donttest{}` blocks) on environments without a libtorch backend or under
  CRAN's core limit:
  - The `iconic_prospect`, `iconic_sensitivity`, `gan_mediation_sensitivity`,
    `gan_pleiotropy_sensitivity`, `train_gan_on_real_data`, and
    `run_single_iteration` examples reached a GAN-training path inside
    `\donttest{}`; switched to `\dontrun{}` so they are not executed by the
    check (they require the optional torch backend).
  - The `infer_confounding` and `iconic_prospect` examples used
    `n_cores = 4`, exceeding CRAN's two-core limit
    (`.check_ncores`); reduced to `n_cores = 2`.
  - The `sample_feature_texture` example referenced an undefined `M`;
    made it self-contained by defining `M` first.
  - The `nc_completeness_capture` and `nc_completeness_check` examples run
    a 1000-permutation loop and exceeded the 5s example-time flag; wrapped
    in `\dontrun{}`.

## Test suite

- **Guarded all torch-dependent tests so `R CMD check` passes without a
  libtorch backend.** Forty `test_that()` blocks across `test-gan.R`,
  `test-model-sensitivity.R`, `test-v06-workflow.R`,
  `test-infer-confounding.R`, and `test-ground-truth-regression.R` reach a
  GAN-training code path (`train_gan_on_real_data()` directly, or
  `iconic_sensitivity()` / `iconic_prospect()` via the internal
  `.auto_train_gan()` step) but lacked a skip guard. On environments where
  the `torch` R package is installed but its libtorch backend is not (e.g.
  the GitHub Actions runners, which do not run `torch::install_torch()`),
  these tests errored with "torch is required for the generative texture
  model" instead of skipping, failing the check. Each block now begins with
  `skip_if_not_installed("torch")` and `skip_if_not(check_torch_setup())`,
  matching the guard convention already used elsewhere in the suite. The
  tests still run (and pass) wherever the libtorch backend is available.

## Packaging and vignettes

- **Removed the shipped pre-trained GAN** (`inst/pretrained_gan.rds`). The
  serialized object used serialization format version 3, which forced an
  implicit `R (>= 3.5.0)` dependency and emitted an `R CMD build` warning.
  Removing it drops that constraint. Texture models are now always trained
  from data, consistent with the package's design that GANs are fit to the
  user's cohort rather than loaded from a shipped artifact.

- **Vignettes train the texture model instead of loading it.** The
  *Walkthrough* and *Sensitivity Analysis* vignettes previously loaded the
  bundled `pretrained_gan.rds` for speed; they now call
  `train_gan_on_real_data()` for a small number of epochs (guarded by
  `check_torch_setup()`), matching the *Texture Model* vignette.

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
