# ICONIC

`iconic` provides causal model selection for observational omics studies where
unmeasured confounding biases ordinary regression. It fits eight estimators of
the natural direct and indirect effects (NDE / NIE), diagnoses which are valid
for your data using genetic instruments and negative controls, stress-tests
them against confounding and pleiotropy violations, and recommends the one
most likely to be unbiased.

The motivating application is estimating the causal effect of an exposure
(e.g. gestational diabetes mellitus) on a panel of molecular outcomes
(e.g. placental transcriptome), but the package is applicable to any
observational study with a genetic instrument, negative controls, and a
mediator and/or outcome panel.

## Installation

```r
# From CRAN (once accepted):
install.packages("iconic")

# From GitHub:
# install.packages("remotes")
remotes::install_github("sbresnahan/iconic")

# The generative texture model (GAN) requires torch:
install.packages("torch")
```

## Quick start

```r
library(iconic)
set.seed(1)

# Simulated mediation panel: exposure Z, mediator M, outcome Y,
# genetic instrument G, negative controls W.
dat <- generate_toy_data(n = 200, n_features = 10, seed = 1)
idat <- iconic_data(
  Z = dat$Z, Y = dat$Y, M = dat$M,
  G = dat$G, W = dat$W, covariates = dat$synthetic_data
)

# 1. Diagnose which estimators are eligible.
diag <- iconic_diagnose(idat)

# 2. Estimate NDE/NIE with all eligible estimators.
est <- iconic_estimate(idat)

# 3. Train a texture model once and attach it (reused downstream).
input <- load_real_input_data(
  Z_matrix = matrix(dat$Z, nrow = 1), Y_matrix = t(dat$Y),
  M_matrix = t(dat$M), W_matrix = t(dat$W),
  covariates_df = dat$synthetic_data
)
idat$trained_gan <- train_gan_on_real_data(
  input$gan_training_data,
  feature_correlations = input$feature_correlations,
  feature_texture      = input$feature_texture,
  epochs = 100, seed = 1
)

# 4. Stress-test robustness (reuses the attached GAN).
sens <- iconic_sensitivity(idat, confounding = "inferred")

# 5. Recommend a single estimator.
rec <- iconic_recommend(idat, diag, est, sens)

# 6. Prospect for future studies (reuses the attached GAN).
pros <- iconic_prospect(idat, confounding = "inferred")
```

## Key functions

| Function | Purpose |
|---|---|
| `iconic_data()` | Standardize exposure, outcome, mediator, instruments, and negative controls into one object. |
| `iconic_diagnose()` | Report instrument strength, negative-control validity (A1, A2, A2'), completeness, and estimator eligibility. |
| `iconic_estimate()` | Fit all eligible estimators; return per-feature NDE, NIE, SE, p-value. |
| `train_gan_on_real_data()` | Train a hybrid GAN + copula texture model from your data. |
| `iconic_sensitivity()` | Sweep instrument exogeneity; report bias degradation and tipping points. |
| `infer_confounding()` | Estimate confounding strength, NC coverage, and latent confounder count from the data. |
| `iconic_recommend()` | Combine eligibility, estimates, and robustness into a ranked recommendation. |
| `iconic_prospect()` | Sweep instrument strength to plan future studies. |
| `generate_toy_data()` | Simulate a mediation panel with known ground truth. |

## Further reading

- **Vignette**: `vignette("iconic-walkthrough")` -- the full model-selection workflow.
- **Function reference**: `help(package = "iconic")` -- the complete API, including
  `gan_sensitivity()`, `pleiotropy_sweep()`, `mediation_sensitivity()`, and the
  negative-control diagnostics (`nc_validity_screen()`, `nc_independence_check()`,
  `nc_completeness_check()`).
- **Manuscript and technical supplement** -- the statistical theory: the eight
  estimators, the negative-control assumptions, the completeness condition, the
  sensitivity surface, and the model-selection criterion.

## License

MIT (c) Sean T. Bresnahan and Charis Xiong.
