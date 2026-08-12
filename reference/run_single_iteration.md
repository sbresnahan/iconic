# Generate one synthetic dataset under the generalised SCM

Generate one synthetic dataset under the generalised SCM

## Usage

``` r
run_single_iteration(
  trained_gan = NULL,
  n_synthetic_samples = 500,
  n_features = 20,
  n_confounders = 1,
  n_mediators = 1,
  beta_X = 0.1,
  alpha_M = 0.5,
  beta_M = 0.3,
  effect_size = NULL,
  conf_strength = 0.8,
  coverage = 1,
  captured = NULL,
  nc_model = "proxy",
  nc_params = list(),
  mo_confounding = 0,
  pleio = 0,
  phi = 0,
  gamma_G = 0.6,
  rho_G1 = 0,
  rho_G2 = 0,
  rho_pop = 0,
  lambda_XM = NULL,
  lambda_MY = NULL,
  omega_1 = NULL,
  omega_2 = NULL,
  feat_cor = 0,
  separate_U = NULL,
  u_strength = NULL,
  w_coverage_profile = NULL,
  MMExp = 1,
  MMOut = 1,
  MMCon = 1,
  MMCpG = 1,
  outcome_type = c("continuous", "survival"),
  surv_h0 = 0.1,
  surv_event_frac = 0.6,
  surv_censor_rate = NULL,
  seed = NULL
)
```

## Arguments

- trained_gan:

  Optional `iconic_gan` from
  [`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md),
  supplying realistic covariate/outcome texture. If `NULL`, default
  synthetic covariates are used.

- n_synthetic_samples:

  Sample size. Default 500.

- n_features:

  Number of outcome (and control) features. Default 20.

- n_confounders:

  Number of latent confounders `k`. Default 1 (backward-compatible
  single-confounder model).

- n_mediators:

  Number of independent mediators. When \> 1, each mediator has its own
  genetic instrument Gm and contributes additively to Y. Default 1
  (single mediator, backward compatible).

- beta_X:

  Direct effect of X on Y (true NDE). Default 0.10.

- alpha_M:

  Effect of X on the mediator M. Default 0.50.

- beta_M:

  Effect of M on Y. Default 0.30.

- effect_size:

  Optional shortcut: if non-`NULL`, sets a pure direct total effect
  (`beta_X = effect_size`, `alpha_M = beta_M = 0`). Use `0` for null
  (Type I error) simulations. Default `NULL` (use the mediation
  parameters).

- conf_strength:

  Overall confounding strength (analogue of `conf_str`); scales the
  confounder loadings into X and Y. Default 0.80.

- coverage:

  How well the negative controls span the confounder subspace, in
  `[0, 1]` (passed to `nc_model`). Default 1.

- captured:

  Integer indices of the confounders the negative controls see (passed
  to `nc_model`). Default all `k`.

- nc_model:

  Negative-control model: a function `(U, covariates, params) -> W`, or
  a registered name (`"proxy"`, `"cpg"`). Default \`"proxy".

- nc_params:

  Extra named parameters forwarded to `nc_model`.

- mo_confounding:

  Strength of U1 -\> M (mediator-outcome confounding). 0 = no M-O
  confounding (original DGP). Default 0. When \> 0, the first confounder
  U\[,1\] also affects M, creating the M-O confounding that the
  mediation estimators are benchmarked against.

- pleio:

  Strength of a direct G -\> Y path (horizontal pleiotropy), violating
  the exclusion restriction. 0 = no pleiotropy (valid instrument,
  original DGP). Default 0. When \> 0, G affects Y directly in addition
  to through X, allowing benchmarking of IV/2SLS under instrument
  invalidity.

- phi:

  Strength of the mediator instrument Gm -\> M. 0 = no mediator
  instrument (original DGP, no Gm generated). Default 0. When \> 0, a
  valid instrument for M is generated (independent of U and G, no direct
  path to Y), enabling point identification of NDE/NIE via
  [`fit_iv2sls_mediation2()`](https://seantbresnahan.com/iconic/reference/fit_iv2sls_mediation2.md)
  even under M-O confounding.

- gamma_G:

  Strength of the exposure instrument G -\> X. Default 0.6.

- rho_G1:

  Correlation of G1 with conf_XM (instrument exogeneity violation).
  Default 0.

- rho_G2:

  Correlation of G2 with conf_MY (instrument exogeneity violation).
  Default 0.

- rho_pop:

  Shared population structure inducing G1-G2 correlation Default 0.

- lambda_XM:

  Optional length-k loading vector giving each confounder's weight on
  the X-\>M backdoor path. NULL (default) = shared loadings.

- lambda_MY:

  Optional length-k loading vector giving each confounder's weight on
  the M-\>Y backdoor path. NULL (default) = shared loadings.

- omega_1:

  Coverage of conf_XM by W1. NULL = use `coverage`.

- omega_2:

  Coverage of conf_MY by W2. NULL = use `coverage`.

- feat_cor:

  Within-module feature correlation. When \> 0 and the GAN does not
  provide feature_correlations, a block-diagonal correlation matrix is
  used for the NC and outcome noise. GAN-learned correlations take
  precedence. Default 0.

- separate_U:

  Defunct. Passing a value errors with a message pointing to the
  replacement per-path loading vectors `lambda_XM` / `lambda_MY`.
  Retained in the signature only to catch and redirect old calls.

- u_strength:

  Numeric vector: per-confounder strength profile (length k). Default
  NULL → rep(1, k) (equal strength, backward compatible). Recycled to
  length k and normalized so the total confounding budget is unchanged.

- w_coverage_profile:

  A list with optional `w1` and `w2` numeric vectors: per-control
  coverage of conf_XM / conf_MY (length n_features). Default NULL →
  scalar omega applied uniformly.

- MMExp, MMOut, MMCon, MMCpG:

  Per-pathway confounding multipliers (exposure, outcome, controls,
  methylation). Default 1.

- outcome_type:

  `"continuous"` (default) or `"survival"` When `"survival"`, the linear
  predictor is converted to `surv_time` and `surv_event` via an
  exponential PH model.

- surv_h0:

  Baseline hazard for the survival DGP. Default 0.1.

- surv_event_frac:

  Target fraction of observed events. Default 0.6.

- surv_censor_rate:

  Explicit censoring rate. Default NULL.

- seed:

  Optional RNG seed.

## Value

A named list matching
[`generate_toy_data()`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md)
— `X`, `G` (n x n_features), `Y`, `W`, `U1`, `M`, `synthetic_data`,
`true_total`, `true_NDE`, `true_NIE` — plus `U` (full n x k confounder
matrix), `genetic_instrument`, `successful_features`, `failed_features`,
and `params`. When `phi > 0`, also includes `Gm` (numeric vector, length
n, or `n_mediators x n` matrix when `n_mediators > 1`). When any
parameter is non-default, also includes `G1`, `G2`, `W1`, `W2`,
`conf_XM`, `conf_MY`, and (when `rho_pop > 0`) `P`. When
`n_mediators > 1`, `M` is an `n_mediators x n` matrix.

## Examples

``` r
dat <- run_single_iteration(NULL, n_synthetic_samples = 100,
  n_features = 5, n_confounders = 1, seed = 1)
analyze_methods_robust(dat)
#>         feature method       beta         se       pvalue significant
#> UNADJ         1  UNADJ  0.5261014 0.05357861 2.987041e-16        TRUE
#> DIRECT        1 DIRECT  0.3610414 0.10312635 7.234765e-04        TRUE
#> COCA          1   COCA -0.5410075 0.22242994 1.500493e-02        TRUE
#> IV2SLS        1 IV2SLS  0.2176178 0.09301922 2.149797e-02        TRUE
#> PGC           1    PGC  0.2755532 0.06758716 9.474808e-05        TRUE
#> UNADJ1        2  UNADJ  0.5302118 0.05187722 4.022556e-17        TRUE
#> DIRECT1       2 DIRECT  0.3705756 0.09715792 2.501072e-04        TRUE
#> COCA1         2   COCA -0.4275786 0.18451456 2.048649e-02        TRUE
#> IV2SLS1       2 IV2SLS  0.2108103 0.08800189 1.864500e-02        TRUE
#> PGC1          2    PGC  0.2764705 0.06419232 4.030979e-05        TRUE
#> UNADJ2        3  UNADJ  0.4984669 0.05121875 4.617412e-16        TRUE
#> DIRECT2       3 DIRECT  0.4048192 0.09989213 1.073455e-04        TRUE
#> COCA2         3   COCA -0.5739017 0.23948376 1.655648e-02        TRUE
#> IV2SLS2       3 IV2SLS  0.1931068 0.09134674 3.725080e-02        TRUE
#> PGC2          3    PGC  0.2821838 0.06580514 4.326252e-05        TRUE
#> UNADJ3        4  UNADJ  0.4327221 0.04984945 8.773302e-14        TRUE
#> DIRECT3       4 DIRECT  0.3741117 0.10202212 4.149712e-04        TRUE
#> COCA3         4   COCA -0.9309171 0.39474764 1.836088e-02        TRUE
#> IV2SLS3       4 IV2SLS  0.1793920 0.09287439 5.652621e-02       FALSE
#> PGC3          4    PGC  0.2642272 0.06745048 1.686514e-04        TRUE
#> UNADJ4        5  UNADJ  0.5350279 0.05522612 5.758312e-16        TRUE
#> DIRECT4       5 DIRECT  0.4092637 0.10640144 2.234730e-04        TRUE
#> COCA4         5   COCA -0.5686515 0.23588898 1.592305e-02        TRUE
#> IV2SLS4       5 IV2SLS  0.2056640 0.09687073 3.646193e-02        TRUE
#> PGC4          5    PGC  0.2938509 0.06987364 5.897800e-05        TRUE
```
