# The Generative Texture Model

## Introduction

This vignette describes the generative texture model used to simulate
realistic omics features: a torch GAN for the covariate table combined
with a Gaussian copula for the feature panel, so that synthetic datasets
inherit the marginal and dependence structure of the observed data.

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

## Two simulation modes

ICONIC uses two distinct simulation modes that serve different purposes:

1.  **Benchmark mode** (`generate_toy_data`): a structural
    data-generating process with parametric Gaussian noise and no
    learned texture. The ground truth (NDE, NIE, confounding strength,
    instrument strength) is a closed-form function of the parameters.
    This mode is used to validate the estimators under controlled
    confounding scenarios with exactly known ground truth.

2.  **Data-calibrated mode** (the generative texture pipeline): trains a
    texture model on the user’s own data so that the sensitivity and
    prospective analyses are calibrated to realistic covariate, outcome,
    and mediator distributions. The ground truth is still imposed by the
    structural skeleton, but the nuisance texture (marginals,
    correlations, covariate structure) is learned from the data.

This separation is deliberate: estimator validation requires a
transparent, reproducible data-generating process whose parameters are
known exactly, while the user-facing workflow benefits from realistic
data texture even at the cost of a learned noise model.

## The hybrid GAN + copula architecture

The texture model is a hybrid of two independently trained components,
each chosen to match the dimensionality and structure of its target
block:

- **Sample-level GAN**: an MLP generator and discriminator trained with
  a non-saturating BCE-with-logits objective, Adam optimization
  ($`\beta_1 = 0.5`$), and one-sided label smoothing. This GAN learns
  the joint distribution of the per-sample exposure level, outcome
  level, mediator level, and encoded covariates (one row per sample). A
  GAN is used here because the block contains mixed continuous and
  categorical variables whose nonlinear interactions and
  mutual-exclusivity constraints are more naturally preserved by a
  flexible generator than by a Gaussian copula.

- **Feature-level Gaussian copula**: learns the full joint distribution
  of the mediator ($`M`$) panel. Each feature’s marginal is fit by
  empirical CDF, or by a parametric family (normal, log-normal, gamma,
  or beta, selected by AIC) when it passes a Kolmogorov-Smirnov test at
  $`p > 0.05`$. Cross-feature dependence is captured by the correlation
  matrix of the normal-score-transformed marginals.

Neither component learns the causal effect. The learned texture is
injected into the same structural causal model as the benchmark mode,
but drawn independently of the latent confounders $`U`$ so it cannot
open a hidden backdoor. Because the generative model never models the
joint of $`(X, M, Y)`$ and only learns the nuisance texture, the
ground-truth effect is guaranteed by construction.

## Training the texture model

The texture model requires `torch`. Verify your installation first:

``` r

check_torch_setup()
```

    #> [1] FALSE

If this returns `FALSE`, install torch with `install.packages("torch")`
and run
[`torch::install_torch()`](https://torch.mlverse.org/docs/reference/install_torch.html).

We simulate a small panel and train the texture model briefly:

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

gan
```

The GAN-training chunks in this vignette are skipped because torch is
not available. Install torch with `install.packages("torch")` and run
[`torch::install_torch()`](https://torch.mlverse.org/docs/reference/install_torch.html)
to run them.

In practice, use 100–300 epochs for a production-quality texture model.
The `epochs = 5` setting here is for speed only.

## Attaching the texture model

Once trained, attach the texture model to the `iconic_data` object so
that every downstream function
([`iconic_sensitivity()`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md),
[`iconic_prospect()`](https://seantbresnahan.com/iconic/reference/iconic_prospect.md))
reuses it without retraining:

``` r

idat$trained_gan <- gan
```

You can also pass a pre-trained model via the `trained_gan` argument to
[`iconic_data()`](https://seantbresnahan.com/iconic/reference/iconic_data.md)
at construction time.

## When to use each mode

| Mode | When to use | Ground truth | Texture |
|----|----|----|----|
| Benchmark (`generate_toy_data`) | Estimator validation, benchmarks | Known exactly (closed-form) | Parametric Gaussian |
| Data-calibrated (GAN + copula) | Sensitivity, prospective analysis | Imposed by structural skeleton | Learned from your data |

The sensitivity and prospective analyses benefit from realistic texture
because estimator performance depends on the empirical distribution of
the analysis data. The benchmark mode is for transparent, reproducible
validation where the ground truth must be known exactly.

## Further reading

- **Sensitivity vignette**:
  [`vignette("iconic-sensitivity")`](https://seantbresnahan.com/iconic/articles/iconic-sensitivity.md)
  — how the trained texture model is used in the sensitivity sweep.
- **Function reference**:
  [`?train_gan_on_real_data`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md)
  for training arguments,
  [`?check_torch_setup`](https://seantbresnahan.com/iconic/reference/check_torch_setup.md)
  for the torch verification helper.

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
    #> [34] bit64_4.8.2         tools_4.6.1         pkgdown_2.2.1      
    #> [37] cachem_1.1.0        desc_1.4.3
