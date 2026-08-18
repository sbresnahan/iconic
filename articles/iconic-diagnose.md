# Diagnosing Estimator Eligibility

## Introduction

This vignette explains how `iconic` decides which causal estimators are
eligible for a given dataset. Eligibility is diagnosed from the data
themselves: instrument strength, negative-control validity and
completeness, and the availability of mediation data each gate a
different subset of the eight estimators.

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

## What diagnosis checks

Before fitting any estimator,
[`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
examines your data and reports which of the eight ICONIC estimators are
*eligible* — that is, which satisfy the assumptions they require to be
valid. These assumptions are not directly testable (the confounder $`U`$
is unobserved); diagnosis checks their empirically observable
projections. Diagnosis checks four things:

1.  **Instrument strength.** The partial $`F`$-statistic for each
    genetic instrument ($`G_1`$ for the exposure, $`G_m`$ for the
    mediator) against the conventional weak-instrument threshold
    ($`F \geq 10`$).
2.  **Negative-control validity.** Three screens — A1 (NCs associate
    with the exposure), A2 (NCs are independent of the instrument given
    covariates), and A2’ (NCs associate with the outcome given the
    exposure) — assessed via partial correlations with BH-FDR control.
3.  **Completeness.** Whether the number of valid negative controls is
    sufficient for the proximal estimators (PGC, PGC2, PGC2Gm), which
    require $`\dim(W_\text{valid}) \geq k`$ (the number of latent
    confounders).
4.  **Eligibility table.** A per-estimator verdict: eligible, ineligible
    (with reason), or not applicable.

## A minimal example

We simulate a small mediation panel with a genetic instrument and
negative controls, then run the diagnosis.

``` r

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
```

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
    #> 6   IV2SLS2     TRUE        TRUE          OK
    #> 7      PGC2     TRUE        TRUE          OK
    #> 8    PGC2Gm     TRUE        TRUE          OK
    #>                                                                 reason
    #> 1                                                      always eligible
    #> 2                                                        G + W present
    #> 3                           W present, valid NCs available (A2 exempt)
    #> 4                                 G present (W augmentation), F_G=98.9
    #> 5                      G + W present, F_G=98.9, completeness satisfied
    #> 6                G + Gm present (W augmentation), F_G=98.9, F_Gm=826.4
    #> 7                  G + W1/W2 present, F_G=98.9, completeness satisfied
    #> 8 G + Gm + W1/W2 present, F_G=98.9, F_Gm=826.4, completeness satisfied

Each row reports whether the estimator is eligible and, if not, why. For
example, `IV2SLS` requires a valid exposure instrument ($`F \geq 10`$),
while `PGC` additionally requires negative-control completeness.

## Instrument strength

``` r

diag$instrument_strength
```

    #> $F_G
    #>       G1 
    #> 98.87353 
    #> 
    #> $F_G_panel
    #>       G1 
    #> 98.87353 
    #> 
    #> $F_G_min
    #> [1] 98.87353
    #> 
    #> $F_G_median
    #> [1] 98.87353
    #> 
    #> $F_G_mean
    #> [1] 98.87353
    #> 
    #> $F_G_p_pass
    #> [1] 1
    #> 
    #> $F_G_n_pass
    #> [1] 1
    #> 
    #> $F_G_n_total
    #> [1] 1
    #> 
    #> $F_Gm
    #> mediator_1 
    #>   826.4032 
    #> 
    #> $F_Gm_panel
    #> mediator_1 
    #>   826.4032 
    #> 
    #> $F_Gm_min
    #> [1] 826.4032
    #> 
    #> $F_Gm_median
    #> [1] 826.4032
    #> 
    #> $F_Gm_mean
    #> [1] 826.4032
    #> 
    #> $F_Gm_p_pass
    #> [1] 1
    #> 
    #> $F_Gm_n_pass
    #> [1] 1
    #> 
    #> $F_Gm_n_total
    #> [1] 1
    #> 
    #> $weak_G
    #> [1] FALSE
    #> 
    #> $weak_Gm
    #> [1] FALSE
    #> 
    #> $min_f
    #> [1] 10
    #> 
    #> $g_threshold
    #> NULL
    #> 
    #> $gm_threshold
    #> NULL

The report gives the partial $`F`$-statistic for $`G_1`$ (exposure
instrument) and $`G_m`$ (mediator instrument). Estimators that rely on
an instrument (`IV2SLS`, `IV2SLS2`, `PGC2Gm`) are eligible only when the
relevant $`F`$-statistic meets the threshold.

## Negative-control validity

``` r

diag$nc_validity
```

    #>   feature      p_value          fdr partial_r relative_effect significant
    #> 1       1 5.582072e-12 1.395518e-11 0.6258002              NA        TRUE
    #> 2       2 2.744321e-10 3.430402e-10 0.5840364              NA        TRUE
    #> 3       3 3.569923e-13 1.784962e-12 0.6519426              NA        TRUE
    #> 4       4 4.220271e-10 4.220271e-10 0.5790310              NA        TRUE
    #> 5       5 1.285262e-11 2.142104e-11 0.6173571              NA        TRUE
    #>                   verdict                   verdict_fdr
    #> 1 drop: associated with X drop: associated with X (FDR)
    #> 2 drop: associated with X drop: associated with X (FDR)
    #> 3 drop: associated with X drop: associated with X (FDR)
    #> 4 drop: associated with X drop: associated with X (FDR)
    #> 5 drop: associated with X drop: associated with X (FDR)
    #>                     verdict_magnitude
    #> 1 drop: associated with X (magnitude)
    #> 2 drop: associated with X (magnitude)
    #> 3 drop: associated with X (magnitude)
    #> 4 drop: associated with X (magnitude)
    #> 5 drop: associated with X (magnitude)

The NC validity screen tests the empirically observable projections of
three identifying assumptions for each negative-control feature:

- **A1**: $`W \perp X \mid C, U`$ (the NC is not affected by the
  exposure once the confounder is accounted for). The screen tests the
  observable projection: association of $`W`$ with $`X`$ given $`C`$.
- **A2**: $`W \perp G_1 \mid C, U`$ (the NC is not an instrument
  itself). The screen tests the observable projection: association of
  $`W`$ with $`G_1`$ given $`C`$.
- **A2’**: $`W \perp G_m \mid C, U`$ (the NC is not a mediator
  instrument). The screen tests the observable projection: association
  of $`W`$ with $`G_m`$ given $`C`$.

Features that fail a screen are flagged, and the completeness count is
reduced accordingly. The `COCA` estimator is exempt from the A2 screen
because it calibrates through the $`W \sim Y + X`$ ratio and does not
require $`W \perp G_1`$.

## Proceeding without instruments or controls

If you supply neither genetic instruments nor negative controls, only
the `UNADJ` estimator (the confounded OLS reference) is eligible. By
default,
[`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
proceeds with a message:

``` r

dat0 <- generate_toy_data(n = 100, n_features = 5, seed = 1)
idat0 <- iconic_data(X = dat0$X, Y = dat0$Y, covariates = dat0$synthetic_data)
diag0 <- iconic_diagnose(idat0)
```

Set `allow_no_proxy = FALSE` to make this an error instead, which is
useful in pipelines where you want to enforce that identification
assumptions are met before proceeding.

## Further reading

- **Walkthrough vignette**:
  [`vignette("iconic-walkthrough")`](https://seantbresnahan.com/iconic/articles/iconic-walkthrough.md)
  — the full model-selection workflow.
- **Function reference**:
  [`?iconic_diagnose`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md)
  for all arguments, including `min_f`, `fdr_level`, and the
  per-estimator NC screens.

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
    #> [10] sass_0.4.10         rmarkdown_2.31      evaluate_1.0.5     
    #> [13] jquerylib_0.1.4     fastmap_1.2.0       yaml_2.3.12        
    #> [16] lifecycle_1.0.5     bookdown_0.47       BiocManager_1.30.27
    #> [19] compiler_4.6.1      fs_2.1.0            systemfonts_1.3.2  
    #> [22] digest_0.6.39       R6_2.6.1            bslib_0.12.0       
    #> [25] withr_3.0.3         tools_4.6.1         pkgdown_2.2.1      
    #> [28] cachem_1.1.0        desc_1.4.3
