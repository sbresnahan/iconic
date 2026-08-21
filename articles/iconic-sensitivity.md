# Sensitivity Analysis and Confounding Inference

## Introduction

This vignette shows how `iconic` stress-tests a causal estimate by
replaying the analysis on synthetic data that preserve the observed
feature texture while sweeping the strength of unmeasured confounding
and instrument exogeneity violations.

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

## What sensitivity analysis does

Even after diagnosis declares an estimator eligible, its validity rests
on the assumption that the genetic instrument is exogenous — that the
instrument-confounder correlations $`\rho_{G1}`$ and $`\rho_{G2}`$ are
zero. These correlations cannot be estimated from data (the confounders
are unobserved), so
[`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
sweeps them over a plausible grid and reports how each estimator’s bias
degrades.

The output reports, for each estimator:

- **Maximum absolute bias** across the sensitivity surface.
- **Tipping points**: the first $`\rho_{G2}`$ at $`\rho_{G1} = 0`$ where
  $`|\text{bias}| > 0.10`$.
- **Texture source**: whether the synthetic data used the generative
  texture model (GAN + copula) or the parametric benchmark generator.

## A minimal example

We simulate a small mediation panel, train a texture model briefly, and
run the sensitivity sweep.

``` r

dat <- generate_toy_data(
  n = 200, n_features = 10, seed = 1,
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

# Train a texture model briefly and attach it to the iconic_data object.
input <- load_real_input_data(
  X_matrix = matrix(dat$X, nrow = 1),
  Y_matrix = t(dat$Y), M_matrix = t(dat$M),
  W_matrix = t(dat$W), covariates_df = dat$synthetic_data
)

idat$trained_gan <- train_gan_on_real_data(
  input$gan_training_data,
  feature_correlations = input$feature_correlations,
  feature_texture = input$feature_texture,
  epochs = 5, seed = 1, verbose = FALSE
)
```

The sensitivity-sweep chunks in this vignette are skipped because torch
is not available. Install torch with `install.packages("torch")` and run
[`torch::install_torch()`](https://torch.mlverse.org/docs/reference/install_torch.html)
to run them.

``` r

sens <- iconic_sensitivity(idat, n_iter = 3)
sens
```

## Inferring confounding from the data

[`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)
estimates the held-fixed confounding parameters from your data,
supplying a `confounding = "inferred"` mode for
[`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
and
[`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md).

``` r

diag <- iconic_diagnose(idat)
est <- iconic_estimate(idat, diagnosis = diag)
conf <- infer_confounding(idat, diag, est)
conf
```

    #> <iconic_confounding>
    #>  conf_strength    0.409 (UNADJ-IV2SLS gap (NDE))
    #>  mo_confounding   0.025 (IV2SLS-IV2SLS2 NIE gap)
    #>  omega_1          0.485 (sqrt(R^2) of W on Y residualized on X+C)
    #>  warning: composite: coverage x confounder strength, not pure coverage
    #>  omega_2          0.484 (sqrt(R^2) of W on Y residualized on X+C)
    #>  warning: composite: coverage x confounder strength, not pure coverage
    #>  k                1 [CI: 1, 2] (parallel analysis (Horn, 1965))
    #> 
    #>  Unavailable: rho_G1, rho_G2

The inferred parameters include:

- **Confounding strength** ($`\delta`$): from the gap between unadjusted
  OLS and IV2SLS estimates.
- **Mediator-outcome confounding** ($`\delta_{mo}`$): from the gap
  between IV2SLS and IV2SLS2 natural indirect effects.
- **Negative-control coverage** ($`\omega_1`$, $`\omega_2`$): from the
  $`R^2`$ of each NC feature regressed on the outcome residual.
- **Number of latent confounders** ($`k`$): via parallel analysis on the
  residualized outcome correlation matrix.

Parameters that cannot be inferred (the instrument-confounder
correlations $`\rho_{G1}`$, $`\rho_{G2}`$) fall back to defaults with a
warning and are always swept.

## Inferred vs. default confounding

``` r

sens_inf <- iconic_sensitivity(idat, n_iter = 3, confounding = "inferred")
sens_inf
```

The `confounding = "default"` mode uses fixed defaults
($`\delta_{mo} = 0.8`$, $`\omega_1 = \omega_2 = 0.7`$, $`\phi = 0.8`$)
and is the recommended starting point. The `confounding = "inferred"`
mode calibrates to your data but uses estimator validity (e.g., that
IV2SLS is unbiased) to set up a benchmark whose purpose is to test
estimator validity — so inferred values should be read as best-case
calibrations under the stated assumptions.

## Further reading

- **Texture model vignette**:
  [`vignette("iconic-texture-model")`](https://seantbresnahan.com/iconic/articles/iconic-texture-model.md)
  — how the GAN + copula texture model is trained and why it matters for
  sensitivity.
- **Function reference**:
  [`?iconic_sensitivity`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md)
  for the full argument list, including the
  $`\rho_{G1} \times \rho_{G2}`$ grid and the `n_iter` control.

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
    #> [1] iconic_0.99.3    BiocStyle_2.40.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] cli_3.6.6           knitr_1.51          rlang_1.3.0        
    #>  [4] xfun_0.60           processx_3.9.0      otel_0.2.0         
    #>  [7] torch_0.17.0        coro_1.1.0          textshaping_1.0.5  
    #> [10] jsonlite_2.0.0      bit_4.6.0           htmltools_0.5.9    
    #> [13] ps_1.9.3            ragg_1.5.2          sass_0.4.10        
    #> [16] rmarkdown_2.31      evaluate_1.0.5      jquerylib_0.1.4    
    #> [19] fastmap_1.2.0       yaml_2.3.12         lifecycle_1.0.5    
    #> [22] bookdown_0.47       BiocManager_1.30.27 compiler_4.6.1     
    #> [25] fs_2.1.0            Rcpp_1.1.2          systemfonts_1.3.2  
    #> [28] digest_0.6.39       R6_2.6.1            magrittr_2.0.5     
    #> [31] callr_3.8.0         bslib_0.12.0        withr_3.0.3        
    #> [34] bit64_4.8.4         tools_4.6.1         pkgdown_2.2.1      
    #> [37] cachem_1.1.0        desc_1.4.3
