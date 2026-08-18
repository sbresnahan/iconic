# Survival Outcomes

## Introduction

This vignette covers time-to-event outcomes. `iconic` supports Cox
proportional-hazards and restricted-mean-survival-time (RMST) effect
scales for all instrument- and negative-control-based estimators.

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

## Time-to-event outcomes in ICONIC

ICONIC extends all estimators except `COCA` to time-to-event outcomes
via two-stage predictor substitution (2SPS). The key insight is that the
first-stage regressions — instrumenting $`X`$ on $`G_1 + W + C`$,
instrumenting $`M`$ on $`\hat{X} + G_m + W + C`$, and constructing the
bridge proxies $`\hat{W}`$ — involve only continuous variables ($`X`$,
$`M`$, $`W`$) as responses, so they remain ordinary least squares
regardless of the outcome type. Only the final outcome stage changes.
(For the instrument-variable estimators `IV2SLS` / `IV2SLS2`, the
negative controls are an optional path-specific augmentation — `W1` for
the exposure first stage, `W2` for the mediator and outcome stages: when
no negative controls are supplied, the first stages instrument on the
genetic instrument(s) and covariates alone.)

## Two effect scales

ICONIC provides two effect scales for survival outcomes:

1.  **Cox log-hazard-ratio (log-HR)**: fits
    `coxph(Surv(time, event) ~ X_hat + M_hat + W_hat + C)` and reports
    the coefficients as log-HRs. The log-HR is **non-collapsible**: the
    conditional log-HR given covariates differs from the marginal log-HR
    even under no confounding, so the product decomposition
    $`\text{NIE} = \hat{\alpha}_M \cdot \hat{\beta}_M`$ is approximate
    on this scale. The log-HR scale is the **recommended primary
    scale**.

2.  **RMST (restricted mean survival time)**: computes
    pseudo-observations via the leave-one-out jackknife of the
    Kaplan-Meier RMST estimate, then fits OLS on the
    pseudo-observations. RMST is **collapsible**, so the NDE/NIE
    decomposition is exact on this scale. However, the
    pseudo-observations introduce additional variance, producing wider
    confidence intervals and higher RMSE than the Cox partial-likelihood
    estimator.

## COCA incompatibility

`COCA` regresses the negative control on the outcome ($`W \sim Y + X`$)
and recovers $`\hat{\tau} = -\hat{\beta}_Z / \hat{\beta}_Y`$. This
requires $`Y`$ as a continuous response variable. With a censored
time-to-event outcome, $`Y`$ is not available as a continuous variable —
the outcome is represented as a `Surv(time, event)` object — so `COCA`
is structurally incompatible with survival outcomes. ICONIC returns `NA`
for COCA on survival data and excludes it from the eligibility report.

## A minimal example

We simulate a small survival panel and estimate on both scales.

``` r

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
```

### Cox log-HR scale

``` r

est_loghr <- iconic_estimate(sdata, effect_scale = "loghr")
head(est_loghr)
```

    #>    feature   mediator  method        NDE    NDE_se     NDE_p        NIE
    #> 1 survival mediator_1   UNADJ  0.1920140 1.2116875 0.2473450 0.08606224
    #> 2 survival mediator_1  IV2SLS  0.3437374 1.4102082 0.2553210 0.05136780
    #> 3 survival mediator_1  DIRECT -0.4659324 0.6275497 0.1198468 0.05838168
    #> 4 survival mediator_1    COCA         NA        NA        NA         NA
    #> 5 survival mediator_1     PGC -0.0486399 0.9525241 0.8089952 0.09005349
    #> 6 survival mediator_1 IV2SLS2  0.2588500 1.2954394 0.3994431 0.06906628
    #>      NIE_se     NIE_p          TE     TE_se      TE_p   alpha_M   alpha_se
    #> 1 0.6429616 0.8935189  0.27807623 1.3717092 0.8393520 0.5496612 0.08406599
    #> 2 0.6668546 0.9385997  0.39510515 1.5599302 0.8000489 0.6132184 0.10274323
    #> 3 0.4606959 0.8991581 -0.40755072 0.7784981 0.6006201 0.3972325 0.15633256
    #> 4        NA        NA          NA        NA        NA        NA         NA
    #> 5 0.5651781 0.8734038  0.04141358 1.1075777 0.9701731 0.4654260 0.10938650
    #> 6 0.5701464 0.9035819  0.32791626 1.4153552 0.8167828 0.4959191 0.08412030
    #>       beta_M beta_M_se NDE_significant NIE_significant
    #> 1 0.15657324  1.169496           FALSE           FALSE
    #> 2 0.08376754  1.087376           FALSE           FALSE
    #> 3 0.14697105  1.158320           FALSE           FALSE
    #> 4         NA        NA              NA              NA
    #> 5 0.19348616  1.213473           FALSE           FALSE
    #> 6 0.13926926  1.149434           FALSE           FALSE

### RMST scale

``` r

est_rmst <- iconic_estimate(sdata, effect_scale = "rmst")
head(est_rmst)
```

    #>    feature   mediator  method        NDE    NDE_se      NDE_p         NIE
    #> 1 survival mediator_1   UNADJ -0.8974612 0.7343346 0.22467757 -0.28876647
    #> 2 survival mediator_1  IV2SLS -1.5477523 1.1920053 0.19738017  0.06450026
    #> 3 survival mediator_1  DIRECT  1.9759180 1.1161907 0.08003739 -0.08779118
    #> 4 survival mediator_1    COCA         NA        NA         NA          NA
    #> 5 survival mediator_1     PGC  0.2825857 0.8436854 0.73841512 -0.13671973
    #> 6 survival mediator_1 IV2SLS2 -1.2495133 1.2202294 0.30843633 -0.30501049
    #>      NIE_se     NIE_p         TE     TE_se      TE_p   alpha_M   alpha_se
    #> 1 0.4100006 0.4812413 -1.1862277 0.8410397 0.1584127 0.5496612 0.08406599
    #> 2 0.4315876 0.8811995 -1.4832520 1.2677321 0.2419992 0.6132184 0.10274323
    #> 3 0.2879132 0.7604253  1.8881268 1.1527253 0.1014284 0.3972325 0.15633256
    #> 4        NA        NA         NA        NA        NA        NA         NA
    #> 5 0.3390685 0.6867853  0.1458659 0.9092703 0.8725495 0.4654260 0.10938650
    #> 6 0.3663229 0.4050549 -1.5545238 1.2740299 0.2224032 0.4959191 0.08412030
    #>       beta_M beta_M_se NDE_significant NIE_significant
    #> 1 -0.5253535 0.7415750           FALSE           FALSE
    #> 2  0.1051832 0.7035867           FALSE           FALSE
    #> 3 -0.2210070 0.7195599           FALSE           FALSE
    #> 4         NA        NA              NA              NA
    #> 5 -0.2937518 0.7252336           FALSE           FALSE
    #> 6 -0.6150408 0.7312704           FALSE           FALSE

On the log-HR scale, the true total effect is $`\beta_X`$ and the
NDE/NIE decomposition is approximate. On the RMST scale, the truth is
computed empirically at zero confounding (because RMST depends on the
baseline hazard and truncation time), and the decomposition is exact.

## Diagnosis and eligibility

``` r

diag <- iconic_diagnose(sdata)
diag$eligibility
```

    #>   estimator eligible a2_required reason_code
    #> 1     UNADJ     TRUE       FALSE          OK
    #> 2    DIRECT     TRUE       FALSE          OK
    #> 3      COCA     TRUE       FALSE          OK
    #> 4    IV2SLS     TRUE        TRUE          OK
    #> 5       PGC     TRUE        TRUE          OK
    #> 6   IV2SLS2     TRUE        TRUE          OK
    #> 7      PGC2     TRUE        TRUE          OK
    #> 8    PGC2Gm     TRUE        TRUE          OK
    #>                                                                 reason
    #> 1                                                      always eligible
    #> 2                                                        G + W present
    #> 3                           W present, valid NCs available (A2 exempt)
    #> 4                                 G present (W augmentation), F_G=94.9
    #> 5                      G + W present, F_G=94.9, completeness satisfied
    #> 6                G + Gm present (W augmentation), F_G=94.9, F_Gm=752.3
    #> 7                  G + W1/W2 present, F_G=94.9, completeness satisfied
    #> 8 G + Gm + W1/W2 present, F_G=94.9, F_Gm=752.3, completeness satisfied

Note that `COCA` is marked ineligible for survival outcomes. All other
estimators are assessed on the same instrument-strength and NC-validity
criteria as for continuous outcomes.

## Further reading

- **Walkthrough vignette**:
  [`vignette("iconic-walkthrough")`](https://seantbresnahan.com/iconic/articles/iconic-walkthrough.md)
  — the full model-selection workflow, including a survival example.
- **Function reference**:
  [`?iconic_estimate`](https://seantbresnahan.com/iconic/reference/iconic_estimate.md)
  for the `effect_scale` argument and
  [`?generate_toy_data`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md)
  for the survival DGP parameters (`surv_h0`, `surv_event_frac`,
  `surv_censor_rate`).

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
    #> [1] iconic_0.99.2    BiocStyle_2.40.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] cli_3.6.6           knitr_1.51          rlang_1.3.0        
    #>  [4] xfun_0.60           otel_0.2.0          textshaping_1.0.5  
    #>  [7] jsonlite_2.0.0      htmltools_0.5.9     ragg_1.5.2         
    #> [10] sass_0.4.10         rmarkdown_2.31      grid_4.6.1         
    #> [13] evaluate_1.0.5      jquerylib_0.1.4     fastmap_1.2.0      
    #> [16] yaml_2.3.12         lifecycle_1.0.5     bookdown_0.47      
    #> [19] BiocManager_1.30.27 compiler_4.6.1      fs_2.1.0           
    #> [22] lattice_0.22-9      systemfonts_1.3.2   digest_0.6.39      
    #> [25] R6_2.6.1            splines_4.6.1       Matrix_1.7-5       
    #> [28] bslib_0.12.0        withr_3.0.3         tools_4.6.1        
    #> [31] survival_3.8-6      pkgdown_2.2.1       cachem_1.1.0       
    #> [34] desc_1.4.3
