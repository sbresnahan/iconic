# ICONIC: Causal Model Selection with Genetic Instruments and Negative Controls

## Introduction

This vignette walks through the complete `iconic` workflow on a
simulated dataset: assembling an `iconic_data` object, diagnosing
estimator eligibility, estimating natural direct and indirect effects,
stress-testing the result, and planning a prospective study.

## Installation

`iconic` is being submitted to Bioconductor. Once accepted, install the
release version with:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("iconic")
```

The development version is available from GitHub:

``` r

BiocManager::install("sbresnahan/iconic")
```

## What this package is for

`iconic` is a workbench for causal inference in observational omics
studies where unmeasured confounding biases ordinary regression. The
motivating application is estimating the causal effect of an exposure
(e.g. gestational diabetes mellitus) on a panel of molecular outcomes
(e.g. placental transcriptome), using genetic instruments and negative
controls to identify the effect.

The package fits eight estimators of the natural direct and indirect
effects (NDE / NIE), diagnoses which are valid for *your* data,
stress-tests them against confounding and pleiotropy violations, and
recommends the one most likely to be unbiased. A generative texture
model (hybrid GAN + copula) lets the sensitivity analysis mirror the
marginal and joint structure of your cohort instead of relying on
generic simulation defaults.

This vignette walks through the core workflow on simulated data. The
statistical theory – the eight estimators, the negative-control
assumptions (A1, A2, A2’), the completeness condition, the sensitivity
surface, and the model-selection criterion – is developed in full in the
accompanying manuscript and technical supplement.

## Quick start

The workflow is six steps:

1.  **Standardize** your data with
    [`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md).
2.  **Diagnose** which estimators are eligible with
    [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md).
3.  **Estimate** the NDE/NIE with all eligible estimators via
    [`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md).
4.  **Stress-test** robustness with
    [`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md).
5.  **Recommend** a single estimator with
    [`iconic_recommend()`](https://seantbresnahan.com/iconic/reference/iconic_recommend.md).
6.  **Prospect** for future studies with
    [`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md).

We illustrate on a simulated mediation panel with a genetic instrument
(`G`), an exposure (`X`), a mediator panel (`M`), an outcome panel
(`Y`), and a panel of negative controls (`W`).

``` r

set.seed(1)
dat <- generate_toy_data(n = 200, n_features = 10, seed = 1)

# Assemble an iconic_data object. Y and M are features x samples matrices;
# X is a length-n vector; G is an n x n_snps matrix; W is an n x n_ncs matrix.
idat <- iconic_data(
  X = dat$X, Y = dat$Y, M = dat$M,
  G = dat$G, W = dat$W,
  covariates = dat$synthetic_data
)
```

## Step 1: Diagnose

[`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
reports instrument strength (partial F for `G` and `Gm`),
negative-control validity screens (A1, A2, A2’), completeness, and an
eligibility table for all eight estimators.

``` r

diag <- iconic_diagnose(idat)
diag$eligibility
```

    #>   estimator eligible a2_required reason_code
    #> 1     UNADJ     TRUE       FALSE          OK
    #> 2    DIRECT     TRUE       FALSE          OK
    #> 3      COCA     TRUE       FALSE          OK
    #> 4    IV2SLS     TRUE        TRUE          OK
    #> 5       PGC     TRUE        TRUE          OK
    #> 6   IV2SLS2    FALSE        TRUE   NEED_DATA
    #> 7      PGC2     TRUE        TRUE          OK
    #> 8    PGC2Gm    FALSE        TRUE   NEED_DATA
    #>                                                                                                             reason
    #> 1                                                                                                  always eligible
    #> 2                                                                                                    G + W present
    #> 3                                                                       W present, valid NCs available (A2 exempt)
    #> 4                                                                            G present (W augmentation), F_G=192.4
    #> 5                                                                 G + W present, F_G=192.4, completeness satisfied
    #> 6                                                        requires G + Gm + F_G>=10 + F_Gm>=10 (F_G=192.4, F_Gm=NA)
    #> 7                                                             G + W1/W2 present, F_G=192.4, completeness satisfied
    #> 8 requires G + Gm + W1/W2 + F_G>=10 + F_Gm>=10 + completeness (have Gm=FALSE, W1/W2=TRUE, completeness: satisfied)

Each estimator is marked eligible only when its required assumptions are
met. For example, `IV2SLS` requires a valid instrument (F \>= 10), while
`PGC` additionally requires negative-control completeness. The
instrument-variable estimators (`IV2SLS`, `IV2SLS2`) are identified by
the instrument(s) alone, so they remain eligible even when no negative
controls are supplied; `W`, when present, is used as an optional
augmentation.

## Step 2: Estimate

[`iconic_estimate()`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md)
fits every eligible estimator and returns per-feature NDE, NIE, standard
error, and p-value.

``` r

est <- iconic_estimate(idat)
head(est)
```

    #>          feature   mediator method        NDE     NDE_se        NDE_p
    #> UNADJ  feature_1 mediator_1  UNADJ  0.6911231 0.39067040 7.844549e-02
    #> DIRECT feature_1 mediator_1 DIRECT  0.6192770 0.26254141 1.938402e-02
    #> COCA   feature_1 mediator_1   COCA  0.4030486 0.35117559 2.510873e-01
    #> IV2SLS feature_1 mediator_1 IV2SLS -0.2849309 0.07031264 7.458041e-05
    #> PGC    feature_1 mediator_1    PGC  0.7204167 0.28261121 1.157215e-02
    #> PGC2   feature_1 mediator_1   PGC2  3.8613917 2.19080836 7.955246e-02
    #>                 NIE    NIE_se        NIE_p        TE      TE_se         TE_p
    #> UNADJ   0.129374321 0.3885271 7.391444e-01 0.8204974 0.55097791 0.1364437548
    #> DIRECT  0.006415187 0.2577313 9.801419e-01 0.6256922 0.36790405 0.0890004830
    #> COCA    0.177040175 0.6676550 7.908807e-01 0.5800888 0.75437884 0.4419157958
    #> IV2SLS  0.615117689 0.0575226 1.091756e-26 0.3301868 0.09084446 0.0002783708
    #> PGC    -0.232467846 0.2799630 4.063398e-01 0.4879488 0.39780447 0.2199712165
    #> PGC2   -3.506045931 2.2024389 1.114089e-01 0.3553458 3.10650585 0.9089305392
    #>           alpha_M    alpha_se       beta_M  beta_M_se NDE_significant
    #> UNADJ   0.9953781 0.007480503  0.129975047 0.39032991           FALSE
    #> DIRECT  0.9813625 0.014004648  0.006537021 0.26262592            TRUE
    #> COCA   -0.8542739 2.878356797 -0.207240524 0.35105126           FALSE
    #> IV2SLS  1.0007388 0.044228755  0.614663581 0.05065558            TRUE
    #> PGC     0.9868719 0.009760905 -0.235560297 0.28367773            TRUE
    #> PGC2    0.9983030 0.044409983 -3.512005721 2.20064391           FALSE
    #>        NIE_significant
    #> UNADJ            FALSE
    #> DIRECT           FALSE
    #> COCA             FALSE
    #> IV2SLS            TRUE
    #> PGC              FALSE
    #> PGC2             FALSE

The `se_method` argument controls standard errors (`"delta"` by default;
`"bootstrap"` for confidence intervals) and the `p_value_method`
argument controls the null-distribution p-value (`"sobel"` or
`"jt-comp"` for higher power under sparse signals).

## Step 3: Train a texture model once

The sensitivity and prospective analyses use a generative texture model
to simulate data that mirrors your cohort’s marginal and joint
structure. Train it **once** and attach it to the `iconic_data` object
so every downstream function reuses it without retraining.

For a runnable example we train a texture model briefly here; in
practice you would train for more epochs on your own cohort with
[`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md)
(see the *Texture Model* vignette).

``` r

# Train a texture model once (few epochs for the vignette).
input <- load_real_input_data(
  X_matrix = matrix(dat$X, nrow = 1),
  Y_matrix = t(dat$Y), M_matrix = t(dat$M),
  W_matrix = t(dat$W), covariates_df = dat$synthetic_data
)

gan <- train_gan_on_real_data(
  input$gan_training_data,
  feature_correlations = input$feature_correlations,
  feature_texture = input$feature_texture,
  epochs = 5, seed = 1, verbose = FALSE
)

# Attach it: iconic_sensitivity() and iconic_prospect() will now reuse it.
idat$trained_gan <- gan
```

The sensitivity, recommendation, and prospective chunks in this vignette
are skipped because torch is not available. Install torch with
`install.packages("torch")` and run
[`torch::install_torch()`](https://torch.mlverse.org/docs/reference/install_torch.html)
to run them.

The texture model requires `torch`; install it with
`install.packages("torch")` before training.

## Step 4: Sensitivity analysis

[`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
sweeps the exogeneity of the instrument (`rho_G1`, `rho_G2`) and reports
how each estimator’s bias degrades. It reuses the attached GAN, so no
retraining happens here.

``` r

sens <- iconic_sensitivity(idat, n_iter = 3)
sens
```

The output reports the maximum absolute bias across the sensitivity
surface for each estimator, the tipping points (first `rho_G2` at
`rho_G1 = 0` where `|bias| > 0.10`), and the texture source. Use
`confounding = "inferred"` to calibrate the held-fixed confounding
parameters from the data:

``` r

conf <- infer_confounding(idat, diag, est)
conf
```

    #> <iconic_confounding>
    #>  conf_strength    0.833 (UNADJ-IV2SLS gap (NDE))
    #>  mo_confounding   default (0.8) -- no mediator instrument
    #>  omega_1          0.544 (sqrt(R^2) of W on Y residualized on X+C)
    #>  warning: composite: coverage x confounder strength, not pure coverage
    #>  omega_2          0.544 (sqrt(R^2) of W on Y residualized on X+C)
    #>  warning: composite: coverage x confounder strength, not pure coverage
    #>  k                1 [CI: 1, 2] (parallel analysis (Horn, 1965))
    #> 
    #>  Unavailable: rho_G1, rho_G2, mo_confounding 
    #> 
    #>  Warnings:
    #>   mo_confounding: no mediator instrument (Gm), using default 0.8.

[`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)
estimates the confounding strength, mediator-outcome confounding,
negative-control coverage (`omega_1`, `omega_2`), and the number of
latent confounders (`k`) from the gap between estimators and the
negative-control residual variance. Pass `confounding = "inferred"` to
[`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
to use these values instead of defaults.

## Step 5: Recommend

[`iconic_recommend()`](https://seantbresnahan.com/iconic/reference/iconic_recommend.md)
combines eligibility, point estimates, and the sensitivity surface into
a single ranked recommendation.

``` r

rec <- iconic_recommend(idat, diag, est, sens)
rec
```

The recommended estimator is the eligible estimator with the best
per-estimand robustness across the sensitivity surface (NDE and NIE are
ranked separately). The full ranking and the ineligible estimators (with
reasons) are printed below the recommendation.

## Step 6: Prospective analysis

[`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md)
answers the design question: *if you don’t yet have a strong instrument,
how strong does it need to be?* It sweeps the instrument strength
(`gamma_G`) and reports how estimates converge, then simulates a
prospective study at a target strength.

``` r

pros <- iconic_prospect(idat, n_iter = 3)
pros
```

The output reports the F-statistics at each strength level and
recommends which estimator would be best if you collected an instrument
of that quality. It also reuses the attached GAN.

## Survival outcomes

`iconic` supports time-to-event outcomes via Cox log-hazard-ratio and
RMST (restricted mean survival time) effect scales. RMST is collapsible
and admits an exact NDE/NIE decomposition; the log-HR scale is reported
for familiarity.

``` r

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
```

    #>    feature method       beta       se    pvalue significant
    #> 1 survival  UNADJ 0.18036720 1.197657 0.2865681       FALSE
    #> 2 survival IV2SLS 0.09989256 1.105052 0.6895867       FALSE
    #> 3 survival DIRECT 0.05365322 1.055119 0.8337282       FALSE
    #> 4 survival   COCA         NA       NA        NA          NA
    #> 5 survival    PGC 0.05374265 1.055213 0.7747792       FALSE

## Putting it together

The full pipeline, end to end:

``` r

# 1. Standardize data (optionally attach a pre-trained GAN).
dat <- generate_toy_data(n = 200, n_features = 10, seed = 1)
idat <- iconic_data(X = dat$X, Y = dat$Y, M = dat$M,
                    G = dat$G, W = dat$W, covariates = dat$synthetic_data)

# 2. Train a texture model once and attach it.
input <- load_real_input_data(
  X_matrix = matrix(dat$X, nrow = 1), Y_matrix = t(dat$Y),
  M_matrix = t(dat$M), W_matrix = t(dat$W),
  covariates_df = dat$synthetic_data
)
idat$trained_gan <- train_gan_on_real_data(
  input$gan_training_data,
  feature_correlations = input$feature_correlations,
  feature_texture = input$feature_texture,
  epochs = 100, seed = 1
)

# 3. Diagnose, estimate, stress-test, recommend.
diag <- iconic_diagnose(idat)
est <- iconic_estimate(idat)
sens <- iconic_sensitivity(idat, confounding = "inferred")
rec <- iconic_recommend(idat, diag, est, sens)

# 4. Prospect for future studies (reuses the attached GAN).
pros <- iconic_prospect(idat, confounding = "inferred")
```

## Further reading

- **Manuscript and technical supplement** – the eight estimators, the
  negative-control assumptions, the completeness condition, the
  sensitivity surface, and the model-selection criterion.
- **Function reference** –
  [`help(package = "iconic")`](https://seantbresnahan.com/iconic/reference)
  for the full API, including
  [`gan_sensitivity()`](https://seantbresnahan.com/iconic/reference/gan_sensitivity.md),
  `pleiotropy_sweep()`, `mediation_sensitivity()`, and the
  negative-control diagnostics
  ([`nc_validity_screen()`](https://seantbresnahan.com/iconic/reference/nc_validity_screen.md),
  [`nc_independence_check()`](https://seantbresnahan.com/iconic/reference/nc_independence_check.md),
  [`nc_completeness_check()`](https://seantbresnahan.com/iconic/reference/nc_completeness_check.md)).

## Session information

``` r

sessionInfo()
```

    #> R version 4.6.1 (2026-06-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] iconic_0.99.0    BiocStyle_2.40.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] Matrix_1.7-5        bit_4.6.0           jsonlite_2.0.0     
    #>  [4] compiler_4.6.1      BiocManager_1.30.27 Rcpp_1.1.2         
    #>  [7] callr_3.8.0         jquerylib_0.1.4     splines_4.6.1      
    #> [10] systemfonts_1.3.2   textshaping_1.0.5   yaml_2.3.12        
    #> [13] fastmap_1.2.0       lattice_0.22-9      R6_2.6.1           
    #> [16] knitr_1.51          bookdown_0.47       desc_1.4.3         
    #> [19] bslib_0.12.0        rlang_1.3.0         cachem_1.1.0       
    #> [22] xfun_0.60           fs_2.1.0            sass_0.4.10        
    #> [25] bit64_4.8.2         otel_0.2.0          cli_3.6.6          
    #> [28] pkgdown_2.2.1       withr_3.0.3         magrittr_2.0.5     
    #> [31] ps_1.9.3            digest_0.6.39       grid_4.6.1         
    #> [34] processx_3.9.0      torch_0.17.0        lifecycle_1.0.5    
    #> [37] coro_1.1.0          evaluate_1.0.5      ragg_1.5.2         
    #> [40] survival_3.8-6      rmarkdown_2.31      tools_4.6.1        
    #> [43] htmltools_0.5.9
