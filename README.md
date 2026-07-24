# ICONIC

**I**nference with **C**ausal **O**bservational **N**egative-control **I**nstruments and **C**ontrols

An R package for causal model selection in multi-omic observational studies with unmeasured confounding, motivated by studies of maternal PFAS exposure and placental transcription. ICONIC provides a modular workflow that diagnoses which estimators are valid for your data, fits all eligible estimators, stress-tests them against assumption violations, and recommends the most robust identified model — backed by a hybrid GAN + copula simulation engine for benchmarking under known ground truth.

---

## Background

Estimating causal effects of prenatal exposures on fetal outcomes is complicated by unmeasured confounding — factors like socioeconomic stress or diet that affect both maternal PFAS levels and placental gene expression. Standard regression adjusts for observed covariates but cannot remove this residual confounding, which does not shrink with sample size.

The ICONIC framework leverages three sources of identification:

- **G** — a polygenic risk score for PFAS metabolism as a genetic instrument (Mendelian randomization)
- **G_m** — a mediator-specific genetic instrument (e.g. an eQTL from fetal genotype) for placental isoform expression
- **W** — a panel of negative-control outcomes (transcripts not on the PFAS causal pathway) that share the same unmeasured confounders as Y

The package provides a toy simulation where the ground truth is known exactly, allowing evaluation of estimators across confounding scenarios. It also supports training a hybrid generative texture model (sample-level GAN + feature-level Gaussian copula) on real data for sensitivity analysis calibrated to the user's own cohort.

---

## Installation

```r
remotes::install_github("karrixxa/iconic")
```

Dependencies: `AER`, `MASS`, `parallel`, `survival`, `torch` (all on CRAN). The `torch` package enables the GAN-based texture generation; the feature-level Gaussian copula learns the mediator panel's full joint distribution. The `survival` package provides Cox proportional hazards and Kaplan-Meier estimators for time-to-event outcomes (v0.9.4).

---

## Model selection workflow

ICONIC provides a modular four-step workflow that uses the simulation framework as an engine for recommending the most appropriate estimator for a given dataset. Rather than requiring the user to choose an estimator a priori, ICONIC diagnoses which estimators are valid for their data, stress-tests them against plausible assumption violations, and recommends the most robust option.

```
iconic_data()  →  iconic_diagnose()  →  iconic_estimate()  →  iconic_recommend()
                                                       ↗
                    iconic_sensitivity()  ───────────────
                                                       ↗
iconic_prospect()  (for users without instruments/NCs yet)
```

### Step 1: Diagnose

`iconic_diagnose()` runs all applicable empirical diagnostics on the user's data and returns an eligibility report for all eight estimators. For each estimator, it checks whether the required inputs are present and whether the identifying assumptions pass empirical screens: instrument strength (first-stage partial F ≥ 10), negative-control exposure-independence (A1: W ⊥ Z | C), instrument-independence (A2: W ⊥ G | C, A2': W ⊥ Gm | C), and proximal completeness (A3: dim(W_valid) ≥ k). Estimators failing any required check are marked ineligible with a reason.

### Step 2: Estimate

`iconic_estimate()` fits all eligible estimators on the user's real data, returning per-feature (and per-mediator) point estimates, standard errors, and p-values. For mediation, the function loops over outcome features × mediators, applying each eligible estimator to each combination.

Three `se_method` options control how the NIE standard error and p-value are computed:

| `se_method` | NIE SE | NIE p-value | Speed |
|---|---|---|---|
| `"delta"` (default) | Delta-method | Sobel / Wald z-test | Fast |
| `"bootstrap"` | Nonparametric bootstrap SD | Wald z-test with bootstrap SE | Slow |
| `"composite"` | Delta-method (unchanged) | Huang (2019) JT-comp composite null test | Fast |

The NIE = α_M × β_M is a product of two coefficients, so H₀: α_M·β_M = 0 is a **composite null** (holds if either coefficient is zero). The Sobel test (delta) approximates the product as normal, which is conservative when both coefficients are zero (sparse signals). The JT-comp test (Huang 2019) accounts for all three null cases via a closed-form p-value based on the normal product distribution, providing higher power when signals are sparse. When `se_method = "composite"`, only `NIE_p` changes; `NDE_p`, `NIE_se`, and point estimates are unchanged. The output includes `var_a` and `var_b` columns for transparency.

**Limitation:** the JT-comp test cannot distinguish the partial-null case (α≠0, β=0) from the alternative when α_M is large, because the product ab can be large even when b ~ N(0,1). This inflates type I error under H₀(2) for estimators with strong stage-1 instruments. The Sobel/delta test handles H₀(2) correctly via the delta-method SE. When α_M is expected to be nonzero and signals are not sparse, `"delta"` may be preferable; `"composite"` is most beneficial when signals are sparse. The JT-comp approximation is valid for Var(a), Var(b) < 1.5 (approximately n < 2000 with sparse signals); ICONIC clamps variances to [1, 1.5].

### Step 3: Stress-test

`iconic_sensitivity()` maps the degradation surface calibrated to the user's data. It generates synthetic data across a grid of instrument exogeneity violations (ρ_G1 × ρ_G2) and negative-control coverage (ω) values, running all eligible estimators at each grid point. The resulting surface shows how much each estimator's bias grows under plausible assumption violations.

### Step 4: Recommend

`iconic_recommend()` combines the diagnostic eligibility, point estimates, and degradation surface to produce a model recommendation. Estimators are ranked by identification tier (identified > negative-control > instrument-based > bias-reference) and within tier by robustness (worst-case bias across the degradation surface).

### Prospective analysis

For users who do not yet have instruments or negative controls, `iconic_prospect()` provides a prospective analysis: it estimates the mediation effect with available methods (UNADJ, DIRECT), maps the confounding sensitivity surface, and simulates how estimates would shift if instruments and negative controls were added at varying quality.

### Quick start

```r
library(iconic)

# 1. Construct a standardized data object from your real data
data <- iconic_data(
  Z = pfas_vec,              # exposure (vector or features x samples matrix)
  Y = expr_matrix,           # outcome (features x samples)
  M = mediator_vec,          # mediator (optional, for mediation)
  G = prs_vec,               # exposure instrument (optional)
  Gm = eqtl_vec,             # mediator instrument (optional)
  W = nc_matrix,             # negative controls (features x samples, optional)
  W1 = nc_zm_matrix,         # path-specific NCs for Z->M path (optional)
  W2 = nc_my_matrix,         # path-specific NCs for M->Y path (optional)
  covariates = meta_df       # sample-level covariates (optional)
)

# 2. Diagnose: which estimators are eligible?
diag <- iconic_diagnose(data)
print(diag)
#   G (exposure instrument):  partial F = 347.1 (ok)
#   Gm (mediator instrument): partial F = 951.4 (ok)
#   NC independence (A2):     10/10 controls valid
#   Completeness:             10 valid NCs vs k=1 -> satisfied
#   Eligible estimators: 8/8

# 3. Estimate: fit all eligible estimators
est <- iconic_estimate(data, diagnosis = diag)
head(est)
#   feature  mediator  method     NDE  NDE_se   NDE_p    NIE  NIE_se   NIE_p
#   gene1    med1      IV2SLS2   0.112  0.031  0.0003   0.136  0.042  0.0012
#   gene1    med1      PGC2Gm    0.130  0.035  0.0002   0.134  0.038  0.0004
#   ...

# 4. Stress-test: map the degradation surface
sens <- iconic_sensitivity(data, diagnosis = diag, n_iter = 50)
print(sens)
#   IV2SLS2: max|NDE bias|=0.090, max|NIE bias|=0.068
#   PGC2Gm:  max|NDE bias|=0.124, max|NIE bias|=0.060

# 5. Recommend: which estimator should you use?
rec <- iconic_recommend(data, diagnosis = diag, estimate = est,
                        sensitivity = sens)
print(rec)
#   Recommended: IV2SLS2 (tier A — identified)
#   Estimated NDE (mean): 0.123
#   Estimated NIE (mean): 0.148
#   Max |bias| across sensitivity surface: 0.090
```

### Estimator eligibility rules

| Estimator | Requires | Tier |
|---|---|---|
| UNADJ | nothing | D (bias reference) |
| DIRECT | G + W | D (bias reference) |
| COCA | W, valid NCs (A2) | B (negative-control) |
| IV2SLS | G, F_G ≥ 10, W | C (instrument-based) |
| PGC | G, F_G ≥ 10, W, completeness | B (negative-control) |
| IV2SLS2 | G + Gm + W, F_G ≥ 10, F_Gm ≥ 10 | A (identified) |
| PGC2 | G + W1/W2, F_G ≥ 10, path completeness | B (negative-control) |
| PGC2Gm | G + Gm + W1/W2, F_G ≥ 10, path completeness | A (identified) |

The A1 screen (W ~ Z | C) flags all NCs that share the confounder U with Z — this is the **intended** negative-control behavior, not a violation. A1 cannot distinguish "W downstream of Z" (true violation) from "W shares a cause with Z" (intended). Therefore **A2 (not A1) is the hard gate** for eligibility.

### Recommendation tiers

| Tier | Estimators | Description |
|---|---|---|
| A (identified) | IV2SLS2, PGC2Gm | Two instruments + path-specific NCs; point-identified under M-O confounding |
| B (negative-control) | PGC, PGC2, COCA | NC-based confounding correction |
| C (instrument-based) | IV2SLS | Single instrument, no NC correction |
| D (bias reference) | UNADJ, DIRECT | Naive; no causal identification |

When sensitivity is supplied, estimators within the same tier are ranked by robustness (lower max |bias| across the violation surface = higher rank).

---

## Causal model

All estimators are judged against a structural causal model (SCM):

```
G1 ~ rho_G1*U_XM + rho_pop*P + sqrt(1-rho_G1^2-rho_pop^2)*N(0,1)  Exposure instrument
G2 ~ rho_G2*U_MY + rho_pop*P + sqrt(1-rho_G2^2-rho_pop^2)*N(0,1)  Mediator instrument (when phi > 0)
U_XM, U_MY ~ N(0,1)                          Path-specific unmeasured confounders
Z = scale(0.6*G1 + delta*U_XM + eps)         Maternal PFAS (scaled)
M = alpha_M*Z + mo_conf*0.5*U_MY + phi*G2 + eps    Mediator
Y_f = beta_M*M + beta_Z*Z + gamma_f*U_MY + pleio*G1 + eps   Placental transcript f
W1_f = omega_1*U_XM + (1-omega_1)*U2 + eps   Negative-control transcript f (Z->M path)
W2_f = omega_2*U_MY + (1-omega_2)*N(0,1) + eps   Negative-control transcript f (M->Y path)
```

When `separate_U = FALSE` (default), `U_XM = U_MY = U1` and `W1 = W2 = W`, reproducing the single-confounder model. When `separate_U = TRUE`, independent confounders drive the Z->M and M->Y paths, and path-specific negative controls W1 and W2 are generated with independent coverage.

**Total effect:** `tau = beta_Z + alpha_M * beta_M`

**Natural direct effect (NDE):** `beta_Z` (Z -> Y not through M)

**Natural indirect effect (NIE):** `alpha_M * beta_M` (Z -> M -> Y)

Key parameters:

- `mo_confounding` — controls whether U also affects M, creating mediator-outcome (M-O) confounding. When > 0, natural effects are not point-identified by a single genetic instrument.
- `phi` — controls the strength of the mediator instrument Gm -> M. When > 0, Gm is generated as an instrument for M, enabling the 2-stage MR estimator to point-identify NDE and NIE even under M-O confounding.
- `pleio` — controls horizontal pleiotropy (a direct G -> Y path), violating the exclusion restriction.
- `rho_G1` — correlation of the exposure instrument G1 with the Z->M confounder U_XM (instrument exogeneity violation). Default 0.
- `rho_G2` — correlation of the mediator instrument G2 with the M->Y confounder U_MY (instrument exogeneity violation). Default 0.
- `rho_pop` — shared population-structure factor P that induces correlation between G1 and G2 (linkage / stratification). Default 0.
- `separate_U` — if TRUE, draw independent confounders U_XM and U_MY for the Z->M and M->Y paths. Default FALSE.
- `omega_1` — coverage of U_XM by the path-specific negative controls W1. NULL uses `w_signal`.
- `omega_2` — coverage of U_MY by the path-specific negative controls W2. NULL uses `w_signal`.
- `feat_cor` — within-module feature correlation for block-diagonal co-expression structure in Y and W. Default 0 (independent features). When > 0, features are divided into `ceiling(sqrt(p))` modules with within-module pairwise correlation = `feat_cor` and between-module correlation = 0.

When all instrument-independence parameters are at their defaults, the DGP is identical to the single-confounder model. When any is non-default, the output additionally includes `G1`, `G2`, `W1`, `W2`, `U_XM`, `U_MY`, and (when `rho_pop > 0`) `P`.

---

## Total-effect estimators

| Name | Function | Approach |
|---|---|---|
| UNADJ | *(internal)* | Unadjusted OLS — bias reference |
| DIRECT | `fit_direct()` | OLS with G, W, and covariates as controls |
| COCA | `fit_coca()` | Negative-control ratio correction (delta method SE) |
| IV/2SLS | `fit_iv2sls()` | Two-stage least squares using G as instrument |
| PGC | `fit_pgc()` | Proxy G-component correction (matrix bridge, 3-step) |

The matrix-bridge PGC regresses the G-residualised exposure Z_resid on the **full W matrix** to construct the confounding proxy W_hat. This requires `ncol(W_valid) >= k` (the proximal completeness condition): if W has fewer valid columns than confounders, the bridge cannot span the confounder subspace and the estimator is under-identified. A scalar-bridge variant (`fit_pgc_scalar()`) is also exported as a standalone fallback; it collapses W to `rowMeans(W)` and is algebraically equivalent to IV/2SLS.

---

## Mediation estimators (NDE / NIE)

| Name | Function | Approach |
|---|---|---|
| UNADJ | `fit_unadj_mediation()` | Naive Baron-Kenny: `M~Z`, `Y~Z+M` — bias reference |
| DIRECT | `fit_direct_mediation()` | OLS with G, W in both mediator and outcome stages |
| COCA | `fit_coca_mediation()` | NC-calibrated both stages: `W~M+Z`, `W~Y+Z+M` |
| IV2SLS | `fit_iv2sls_mediation()` | `Z~G+W->Z_hat`; `M~Z_hat`; `Y~Z_hat+M+W` (IV-cleansed exposure) |
| PGC | `fit_pgc_mediation()` | Matrix bridge: `Z_resid~W->W_hat`; `M~Z+W_hat`; `Y~Z+M+W_hat` |
| IV2SLS2 | `fit_iv2sls_mediation2()` | 2-stage MR: `Z~G+W->Z_hat`; `M~Z_hat+Gm+W->M_hat`; `Y~Z_hat+M_hat+W` |
| PGC2 | `fit_pgc_mediation2()` (gm=NULL) | Two-stage proximal: path-specific bridges `W1->U_XM`, `W2->U_MY` |
| PGC2Gm | `fit_pgc_mediation2()` (gm supplied) | NC-augmented: Gm isolates U_MY before bridging W2 |

Each mediation estimator returns `list(NDE, NDE_se, NDE_p, NIE, NIE_se, NIE_p, alpha_M, alpha_se, beta_M, beta_M_se)`. The `alpha_M`/`alpha_se`/`beta_M`/`beta_M_se` fields (v0.9.3) expose the stage-1 and stage-2 coefficient estimates and their SEs, enabling the composite null test.

**IV2SLS2** (2-stage MR) uses two instruments: G for Z and Gm for M, instrumenting both endogenous variables. This estimator is included only when `dat$Gm` is present (`phi > 0`). It is point-identified under M-O confounding but requires instrument exogeneity — both G and Gm must be independent of the confounders.

**PGC-2** (`fit_pgc_mediation2`) is a two-stage proximal mediation estimator that uses path-specific negative controls: W1 captures the Z->M confounder U_XM, and W2 captures the M->Y confounder U_MY. The three-stage bridge procedure constructs proxies for both confounders and includes them in the outcome regression. Unlike IV2SLS2, PGC-2 does not require instrument exogeneity — the bridge absorbs confounding regardless of instrument–confounder correlation. When `gm` is supplied (PGC2Gm), the mediator instrument helps isolate U_MY before bridging, improving performance when the mediator instrument is correlated with the confounder.

PGC-2 and PGC2Gm are included in the pipeline only when path-specific negative controls are present (`separate_U = TRUE` or any path-specific parameter is non-default). The five base methods are unchanged whether or not these are available, giving eight estimators total when all inputs are present.

A scalar-bridge mediation variant (`fit_pgc_scalar_mediation()`) is also exported as a standalone fallback.

---

## Survival / time-to-event outcomes (v0.9.4)

ICONIC v0.9.4 extends the framework to support **Cox proportional hazards** and **restricted mean survival time (RMST)** outcomes, motivated by case studies where the outcome is overall survival (e.g., smoking → lung tumor expression → lung cancer survival). All five total-effect estimators and all eight mediation estimators have survival counterparts that share the same 2SPS (two-stage predictor substitution) structure: first-stage regressions (Z~G, M~Z_hat+Gm, bridges) remain OLS because Z, M, and W are continuous; only the outcome stage switches to `survival::coxph` (log-HR scale) or OLS on RMST pseudo-observations (time scale).

### Two effect scales

| `effect_scale` | Outcome stage | NDE/NIE decomposition | Collapsible |
|---|---|---|---|
| `"loghr"` (default) | `coxph(Surv(time, event) ~ ...)` | Product α·β on log-HR scale (approximate) | No |
| `"rmst"` | OLS on RMST pseudo-observations | Product α·β on time scale (exact) | Yes |

The Cox log-HR is **non-collapsible**, so the product-of-coefficients NIE = α·β on the log-HR scale is an approximation. The RMST scale (pseudo-observation OLS) provides the exact decomposition because RMST differences are collapsible. RMST pseudo-observations are computed via leave-one-out jackknife of Kaplan-Meier RMST (Graw et al. 2009), with the truncation time τ defaulting to the 90th percentile of follow-up.

### COCA incompatibility

COCA regresses W on Y (W ~ y + Z), placing the outcome on the right-hand side — impossible with a `Surv` object. Both `fit_coca_surv()` and `fit_coca_mediation_surv()` return `NA` with an informative `"reason"` attribute for survival outcomes.

### Usage

```r
library(iconic)

# Create survival data: pass surv_time and surv_event instead of Y
sdat <- iconic_data(
  Z = pack_years, outcome_type = "survival",
  surv_time = os_time, surv_event = os_event,
  M = tumor_expr, G = smoking_prs, Gm = lung_eqtl,
  W = non_lung_gp_expr)

# Total effect on Cox log-HR scale
est <- iconic_estimate(sdat, effect_scale = "loghr")

# Total effect on RMST scale (collapsible, exact NDE/NIE decomposition)
est_rmst <- iconic_estimate(sdat, effect_scale = "rmst", tau = 365)

# Mediation: NDE and NIE under survival
est_med <- iconic_estimate(sdat, effect_scale = "loghr")
# Columns: NDE, NDE_se, NDE_p, NIE, NIE_se, NIE_p, ...
```

### Survival simulation

The toy simulation supports `outcome_type = "survival"` in both `generate_toy_data()` and `run_single_iteration()`. The linear predictor is converted to time-to-event via an exponential proportional hazards model with configurable baseline hazard (`surv_h0`), event fraction (`surv_event_frac`), and censoring rate (`surv_censor_rate`). The true effects (total, NDE, NIE) are on the Cox log-HR scale.

```r
# Generate survival toy data
dat <- generate_toy_data(n = 500, beta_Z = 0.25, alpha_M = 0.5, beta_M = 0.3,
                         conf_str = 0.8, phi = 0.8, outcome_type = "survival",
                         surv_event_frac = 0.6)
# dat$surv_time, dat$surv_event are now available
```

Sensitivity analysis (`iconic_sensitivity()`, `iconic_prospect()`) also threads `outcome_type` and `effect_scale` through the simulation engine. The GAN texture model is skipped for survival outcomes (it requires continuous Y); the default texture is used instead.

---

## Empirical negative-control validity diagnostics

The NC-based estimators (COCA, PGC, PGC-2) require identifying assumptions that hold by construction in the simulation but must be **tested on real data**:

1. **W perp Z | C** — controls are not affected by the exposure.
2. **W perp G | C** — controls are independent of the instrument (no meQTLs).
3. **W perp Gm | C** — controls are independent of the mediator instrument (no eQTL-CpG structure).

| Function | Screen | Method |
|---|---|---|
| `nc_validity_screen()` | W perp Z \| C | Regress each W_f on Z + C; BH-FDR control at 0.10 |
| `nc_independence_check()` | W perp G \| C | Partial correlation of W_f with G after residualizing on C; BH-FDR |
| `nc_independence_check_gm()` | W perp Gm \| C | Partial correlation of W_f with Gm after residualizing on C; BH-FDR |
| `nc_completeness_check()` | dim(W_valid) >= k | Count valid controls vs number of confounders; flag under-identified |

```r
dat <- run_single_iteration(n_features = 10, n_confounders = 1, phi = 0.8, seed = 1)

# A1: Are any controls associated with the exposure?
screen_Z <- nc_validity_screen(dat)

# A2: Are any controls associated with the instrument?
screen_G <- nc_independence_check(dat)

# A2': Are any controls associated with the mediator instrument?
screen_Gm <- nc_independence_check_gm(dat)

# Completeness: enough valid controls for the confounder dimensionality?
comp <- nc_completeness_check(dat)
comp$completeness   # "satisfied", "borderline", or "under-identified"
```

Controls that fail either screen should be dropped before running COCA/PGC. `nc_independence_check_gm()` returns NULL with a message when `dat$Gm` is absent.

---

## Simulation and sensitivity analysis

### Toy simulation

The toy simulation generates synthetic datasets with known ground truth, allowing evaluation of estimator bias, RMSE, and Type I error across parameter grids.

```r
# Total-effect simulation: sweep confounding strength
res <- sweep_param("conf_str", c(0.2, 0.5, 0.8, 1.0), n_iter = 50)
res$summary

# Mediation simulation: sweep a single parameter
med <- sweep_mediation_param("conf_str", c(0.2, 0.5, 0.8), n_iter = 50,
                             mo_confounding = 0.8, phi = 0.8)
med$summary

# Sweep instrument exogeneity parameters
degradation <- sweep_mediation_param("rho_G2", c(0, 0.1, 0.2, 0.3, 0.5),
    n_iter = 50, mo_confounding = 0.8, phi = 0.8,
    rho_G1 = 0.3, separate_U = TRUE, omega_1 = 0.7, omega_2 = 0.7)
degradation$summary

# Type I error under the null (no true effect)
null_res <- run_null_mediation_sim(n_iter = 200, mo_confounding = 0.8, phi = 0.8)
null_res$rates
```

### Hybrid GAN + copula sensitivity analysis on your own data

Beyond the toy simulation, `iconic` can train a hybrid generative texture model on a **real dataset** and use it to generate realistic in-silico data for a sensitivity analysis: *which estimator is preferred*, and *do the negative controls hold under different confounding scenarios?*

The design keeps a **known ground truth** by separating three concerns: a structural causal skeleton (with a tunable, multi-confounder `U`), a **pluggable negative-control mechanism**, and a generative *texture* model. The texture model is a hybrid of two independently trained components: a sample-level `torch` GAN that learns the joint distribution of per-sample exposure, outcome, mediator level, and encoded covariates, and a feature-level Gaussian copula (`train_feature_texture()`) that learns the full joint distribution of the mediator (M) panel — marginal distributions (empirical CDF by default, with parametric fallback via KS test) and cross-feature dependence. The `torch` package is a hard dependency.

```r
# 1. Load real data (features x samples matrices + covariates), or the example
input <- load_real_input_data(Z_matrix, Y_matrix, covariates_df = meta)

# 2. Train the generator (hybrid torch GAN + Gaussian copula)
gan <- train_gan_on_real_data(input$gan_training_data, epochs = 300)

# 3. Sweep confounding scenarios: strength x NC coverage x #confounders
sens <- gan_sensitivity(gan, conf_grid = c(0.2, 0.5, 0.8),
                        coverage_grid = c(0.3, 0.6, 1.0), k_grid = c(1, 2))

# 4. Which estimator wins, robustly across scenarios?
recommend_estimator(sens)$overall

# 5. Do the negative controls hold as confounding grows?
nc_validity_check(gan, k_grid = c(1, 2, 3), n_valid_controls = 1)$verdict

# 6. How robust are the estimators to pleiotropy?
pleio_sens <- gan_pleiotropy_sensitivity(gan,
    pleio_grid = c(0, 0.05, 0.10), conf_grid = c(0.2, 0.5, 0.8),
    n_iter = 50, n_features = 10)

# 7. Mediation sensitivity (set phi > 0 to include the 2-stage MR estimator)
med_sens <- gan_mediation_sensitivity(gan,
    conf_grid = c(0.2, 0.5, 0.8), coverage_grid = c(0.3, 0.7, 1.0),
    mo_confounding = 0.8, phi = 0.8, n_iter = 50)
head(med_sens$summary)
```

**Pluggable negative controls.** A negative-control model is any function `function(U, covariates, params) -> W`. Built-ins: `nc_proxy` (direct proxy, with `mode = "shared"` or `"distinct"`) and `nc_cpg` (CpG-predicted expression). Negative controls are only valid when they span the confounder subspace driving `Y`; the `coverage`/`k` knobs let you sweep from perfect to broken, and `nc_validity_check()` flags the proximal-inference regime where `#confounders > #valid controls`.

A full runnable example lives in `inst/scripts/run_iconic_pipeline.R`.

---

## Parallelization

ICONIC v0.8.3 introduces cross-platform multi-core parallelization via the `n_cores` argument on all computationally intensive functions. The implementation uses `parallel::mclapply` on Unix (Linux, macOS) and PSOCK clusters (`makeCluster`/`parLapply`/`stopCluster`) on Windows, requiring no additional dependencies beyond the base `parallel` package.

### Functions supporting `n_cores`

| Function | What is parallelized |
|---|---|
| `iconic_estimate()` | Per-feature (and per-mediator × feature) estimation |
| `iconic_sensitivity()` | Replicate iterations within each grid cell |
| `iconic_prospect()` | Replicate iterations in Phase 1 and Phase 2 |
| `iconic_diagnose()` | NC validity screens (A1, A2, A2') across control features |
| `infer_confounding()` | NC-coverage loop (omega) and k permutation analysis |
| `sweep_nc_validity()` | Replicate iterations in all four diagnostic panels |
| `sweep_instrument_strength()` | Replicate iterations per instrument-strength grid point |
| `analyze_methods_robust()` | Per-feature estimation |
| `analyze_mediation_robust()` | Per-feature mediation estimation |
| `run_methods()` / `run_mediation_methods()` | Per-feature estimation (internal drivers) |
| `nc_validity_screen()` | Per-control-feature regression |
| `nc_independence_check()` / `nc_independence_check_gm()` | Per-control-feature partial correlation |
| `nc_completeness_check()` | Delegates to `nc_validity_screen` and `nc_independence_check` |

### Usage

```r
# Use 8 cores on an HPC node
est <- iconic_estimate(data, diagnosis = diag, n_cores = 8)
sens <- iconic_sensitivity(data, diagnosis = diag, n_iter = 1000, n_cores = 8)
conf <- infer_confounding(data, diagnosis = diag, estimate = est, n_cores = 8)

# Auto-detect available cores
n_cores <- parallel::detectCores() - 1
diag <- iconic_diagnose(data, n_cores = n_cores)
```

### Progress messages

When `n_cores > 1`, each parallelized function emits `message()`-based progress milestones: a start line announcing the task count and core count, and completion lines as chunks finish. In sequential mode (`n_cores = 1`, the default), progress is reported every ~10% of tasks. These messages go to `stderr` and do not interfere with `capture.output()` or R Markdown rendering.

### Design notes

- **Default `n_cores = 1`** preserves exact backward compatibility with v0.8.2.
- **No nested parallelism**: simulation functions that are already parallelized at the replicate level (e.g., `run_simulation`, `sweep_param`) call `run_methods()` with `n_cores = 1` inside each worker, avoiding nested fork overhead.
- **Windows PSOCK**: on Windows, each worker loads the `iconic` namespace via `clusterEvalQ(cl, library(iconic))` so that internal helper functions are available. On Unix, `mclapply` inherits the parent's environment by default.
- **No new dependencies**: the implementation uses only the `parallel` package, which is already in `Imports`.

---

## Feature-level correlations

Real omics data exhibit cross-feature correlations — e.g., placental transcript expression forms co-expression modules with high within-module and low between-module correlation. ICONIC supports this structure in two places: a `feat_cor` parameter in the toy simulation, and feature-level texture learning in the generative pipeline (Gaussian copula for the mediator panel, residual correlation matrices for Y and W).

### Block-diagonal co-expression in the toy simulation

`generate_toy_data()` and `sweep_mediation_param()` accept a `feat_cor` parameter that injects block-diagonal correlated noise into the outcome (Y) and negative-control (W) panels. Features are divided into `ceiling(sqrt(p))` equal-sized modules; within-module pairwise correlation equals `feat_cor`, between-module correlation is 0. For `p = 10`: 4 modules of sizes 3, 3, 2, 2.

```r
# Sweep feature correlation strength
# rho_G1 = 0 ensures the NC bridge is live (A2 passes, all 8 estimators eligible)
fc_sweep <- sweep_mediation_param("feat_cor", c(0, 0.2, 0.4, 0.6, 0.8),
    n_iter = 50, mo_confounding = 0.8, phi = 0.8,
    rho_G1 = 0, separate_U = TRUE, omega_1 = 0.7, omega_2 = 0.7)
fc_sweep$summary

# Plot the sweep
plot_feature_correlation_sweep(fc_sweep)
```

The mediator M is scalar by design (one mediator per mediation analysis), so `feat_cor` does not apply to M.

### Feature-level texture learning in the generative pipeline

When training a generative texture model on real data, two mechanisms capture feature-level structure:

- **Mediator panel (M)**: the feature-level Gaussian copula (`train_feature_texture()`) learns each feature's marginal distribution (empirical CDF or parametric family) and their cross-feature dependence (copula correlation matrix). During simulation, copula draws replace the parametric Gaussian noise in the mediator structural equations, so synthetic mediator data preserves the marginal shapes and correlation structure of the user's cohort.

- **Outcome (Y) and negative-control (W) panels**: feature-level residual correlation matrices are learned after residualizing on the sample-level exposure and covariates. During simulation, the noise terms in the outcome and negative-control equations are drawn from a multivariate normal with the learned correlation structure, so that conditional on the confounder U the simulated features retain realistic cross-feature correlations. When feature correlation matrices are unavailable (e.g., panels not supplied or feature count mismatch), the noise defaults to independent draws.

The negative-control (W) panel is not learned from data via the copula; it continues to use the pluggable `nc_model` with residual correlation matrices injected as correlated noise.

The example data (`load_real_input_data(example = TRUE)`) includes block-diagonal co-expression modules (within-block r = 0.5) so that the learned correlations are non-trivial by default.

---

## Plotting

The package exports plotting helpers that accept the data frames returned by `sweep_param()` or `run_simulation()`:

| Function | Plot |
|---|---|
| `plot_bias()` | Mean bias vs swept parameter |
| `plot_bias_boxplot()` | Per-iteration bias distribution (box-and-points) |
| `plot_bias_distribution()` | Bias histogram or density |
| `plot_type1_error()` / `plot_type1_vs_conf()` | Type I error rates |
| `plot_power()` | Power curves |
| `plot_sensitivity_heatmap()` | GAN sensitivity heatmap |
| `plot_gan_diagnostics()` | Generator fidelity (real vs synthetic marginals) |
| `plot_feature_correlation_sweep()` | NDE/NIE bias, RMSE, and Type I error vs feature correlation strength |

The canonical method order and colour palette are exported as `iconic_method_order` and `iconic_method_colors`.

---

## Package structure

```
iconic/
├── R/
│   ├── iconic_data.R            # iconic_data() — standardized data container
│   ├── diagnose.R               # iconic_diagnose() — eligibility diagnostics
│   ├── estimate.R               # iconic_estimate() — fit all eligible estimators
│   ├── model_sensitivity.R      # iconic_sensitivity() — degradation surface mapping
│   ├── recommend.R              # iconic_recommend() — model recommendation
│   ├── prospect.R               # iconic_prospect() — prospective analysis without instruments
│   ├── estimators.R             # fit_direct, fit_coca, fit_iv2sls, fit_pgc, fit_pgc_scalar
│   ├── mediation.R              # fit_*_mediation, fit_iv2sls_mediation2, fit_pgc_mediation2, analyze_mediation_robust
│   ├── surv_estimators.R        # fit_*_surv — Cox/RMST total-effect estimators (v0.9.4)
│   ├── surv_mediation.R         # fit_*_mediation_surv — Cox/RMST mediation estimators (v0.9.4)
│   ├── surv_dgp.R               # .linpred_to_surv — exponential PH + censoring (v0.9.4)
│   ├── generate_data.R          # generate_toy_data (internal toy DGP, feat_cor, survival outcome)
│   ├── run_methods.R            # run_methods, analyze_methods_robust/parallel, summarise_results
│   ├── mediation_simulation.R   # run_mediation_sim, sweep_mediation_param, null sims
│   ├── simulation.R             # run_simulation, sweep_param, run_null_sim, sweep_null_by_conf
│   ├── load_data.R              # load_real_input_data (real-data ingest, example data with co-expression modules)
│   ├── gan.R                    # train_gan_on_real_data, GAN architecture, sample_texture, feature correlation learning
│   ├── feature_texture.R        # train_feature_texture, sample_feature_texture (Gaussian copula for mediator panel)
│   ├── nc_models.R              # pluggable negative-control interface: nc_proxy, nc_cpg
│   ├── nc_diagnostics.R         # nc_validity_screen, nc_independence_check, nc_independence_check_gm, nc_completeness_check
│   ├── run_iteration.R          # run_single_iteration (generalized generator)
│   ├── sensitivity.R            # gan_sensitivity, gan_mediation_sensitivity, recommend_estimator, nc_validity_check
│   ├── pleiotropy_sensitivity.R # gan_pleiotropy_sensitivity
│   ├── infer_confounding.R      # infer_confounding() — data-calibrated confounding parameters
│   ├── plots.R                  # plotting helpers, method order, colour palette
│   └── figures.R                # publication figure functions (plot_feature_correlation_sweep, etc.)
├── inst/scripts/                # run_iconic_pipeline.R end-to-end example
├── vignettes/iconic-walkthrough.Rmd  # package vignette
├── tests/testthat/              # test suite including ground truth regression test
├── DESCRIPTION
└── NAMESPACE
```

---

## License

MIT
