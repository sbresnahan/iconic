# Sensitivity (degradation) surface and effect-decomposition bias sweep

Sweeps a 2D grid of instrument-independence violations (`rho_G1` x
`rho_G2`) and reports how each estimator's NDE/NIE bias degrades as the
genetic instruments become correlated with the unmeasured confounders.

## Usage

``` r
iconic_sensitivity(
  data,
  diagnosis = NULL,
  trained_gan = NULL,
  confounding = c("default", "inferred", "manual"),
  gan_epochs = 100,
  rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  n_iter = 30,
  n_samples = NULL,
  n_features = 10,
  phi = NULL,
  mo_confounding = 0.8,
  lambda_XM = NULL,
  lambda_MY = NULL,
  omega_1 = c(0.3, 0.7, 1),
  omega_2 = c(0.3, 0.7, 1),
  bias_threshold = 0.1,
  base_seed = 700,
  n_cores = 1,
  verbose = FALSE,
  outcome_type = NULL,
  effect_scale = c("loghr", "rmst", "logor", "riskdiff"),
  surv_h0 = 0.1,
  surv_event_frac = 0.6,
  surv_censor_rate = NULL,
  bin_prev = 0.5
)
```

## Arguments

- data:

  An `iconic_data` object.

- diagnosis:

  Optional `iconic_diagnosis` from
  [`iconic_diagnose()`](https://seantbresnahan.com/iconic/reference/iconic_diagnose.md).
  Used to infer `phi` and `mo_confounding` when not explicitly supplied.

- trained_gan:

  Optional `iconic_gan` from
  [`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md).
  If `NULL`, a texture model is auto-trained from `data`.

- confounding:

  Confounding parameter source: `"default"` (fixed defaults),
  `"inferred"` (data-calibrated via
  [`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)),
  or `"manual"` (use explicitly supplied arguments). Default
  `"default"`. Alternatively, supply a precomputed `iconic_confounding`
  object from a prior
  [`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)
  call to reuse it as-is (equivalent to `"inferred"` but skips
  recomputation). When `"inferred"` runs internally, the gap-based
  calibration uses
  [`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)'s
  documented random subset of mediators/features (`max_infer_tasks`),
  not the full panel. Inferred scalars only fill in `mo_confounding` /
  `omega_1` / `omega_2` when those arguments are left at their defaults;
  explicitly supplied values (e.g. an omega sweep) always take
  precedence.

- gan_epochs:

  Epochs for auto-trained GAN. Default 100.

- rho_G1_grid:

  Values of rho_G1 (G correlation with conf_XM). Default
  `c(0, 0.1, 0.2, 0.3, 0.5)`.

- rho_G2_grid:

  Values of rho_G2 (Gm correlation with conf_MY). Default
  `c(0, 0.1, 0.2, 0.3, 0.5)`.

- n_iter:

  Replicates per grid cell. Default 30.

- n_samples:

  Samples per replicate. If NULL, uses `data$n`.

- n_features:

  Features per replicate. Default 10.

- phi:

  Mediator instrument strength. If NULL, inferred from diagnosis (F_Gm)
  or defaults to 0.8.

- mo_confounding:

  M-O confounding strength. Default 0.8. Used when
  `confounding = "default"` or `"manual"`.

- lambda_XM:

  Optional per-path confounder loading vector (X-\>M path). NULL
  (default) = shared loadings.

- lambda_MY:

  Optional per-path confounder loading vector (M-\>Y path). NULL
  (default) = shared loadings.

- omega_1, omega_2:

  NC coverage of each path's confounder composite (conf_XM / conf_MY).
  Default `c(0.3, 0.7, 1.0)`, swept jointly with the rho grid. When the
  two vectors are identical (the default), the sweep is taken on the
  diagonal (`omega_1 == omega_2`); supply distinct vectors to cross the
  full `omega_1 x omega_2` grid. Used when `confounding = "default"` or
  `"manual"`.

- bias_threshold:

  Absolute bias threshold for tipping-point annotation. Default 0.10.

- base_seed:

  Base RNG seed. Default 700.

- n_cores:

  Number of parallel workers for simulation replicates. Default 1
  (sequential). Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) on
  Unix and a PSOCK cluster on Windows.

- verbose:

  Logical: print progress messages during the sweep. Default `FALSE`
  (quiet).

- outcome_type:

  `NULL` (inherit from `data`, default) or `"continuous"` / `"survival"`
  / `"binary"`. When survival, the sensitivity sweep uses the Cox / RMST
  survival mediation drivers; when binary, the logistic /
  linear-probability-model binary mediation drivers.

- effect_scale:

  `"loghr"` (default), `"rmst"`, `"logor"`, or `"riskdiff"`.
  `"loghr"`/`"rmst"` apply to survival outcomes; `"logor"`/`"riskdiff"`
  apply to binary outcomes (an inert default is remapped with a
  message).

- surv_h0:

  Baseline hazard for survival DGP. Default 0.1.

- surv_event_frac:

  Target fraction of observed events. Default 0.6.

- surv_censor_rate:

  Explicit censoring rate. Default NULL.

- bin_prev:

  Target prevalence for the binary DGP. Default 0.5.

## Value

An `iconic_sensitivity` S3 object: a named list with `$surface` (data
frame: rho_G1, rho_G2, method, NDE_bias, NIE_bias, NDE_rmse, NIE_rmse,
NDE_type1, NIE_type1, tipped), `$grid`, `$tipping_points`, `$summary`,
`$texture_source` (how the GAN was obtained), and
`$inferred_confounding` (when `confounding = "inferred"`).

## Details

The grid is calibrated to the user's data: sample size, mediator
instrument strength (`phi`), confounding level, and covariate / outcome
texture are inferred from the `iconic_data` object and diagnosis when
available.

## When to use

Use this as an **effect-decomposition bias sweep** when you need to know
how robust a mediation estimate (NDE/NIE) is to violations of the
genetic instruments' independence from the confounders. The degradation
surface shows, for each estimator, the point at which bias becomes
material — the "tipping point" — so you can report not just a point
estimate but the range of violations it tolerates. This is the mediation
analogue of a total-effect sensitivity sweep and is most informative
when comparing estimators that make different independence assumptions
(e.g. IV2SLS2 vs PGC2Gm).

## Texture model

When `trained_gan` is `NULL` and no GAN is attached to `data`, a texture
model is auto-trained from the user's data (GAN via `torch` if
available, otherwise a multivariate-normal fallback). This ensures the
synthetic covariate and outcome texture matches the user's cohort.
Supply `trained_gan` explicitly or attach one via
[`iconic_data`](https://seantbresnahan.com/iconic/reference/iconic_data.md)`(trained_gan = ...)`
to reuse a pre-trained model and avoid retraining.

## Confounding calibration

The `confounding` argument controls how the held-fixed confounding
parameters (`mo_confounding`, `omega_1`, `omega_2`) are set:

- `"default"`: fixed defaults (0.8, 0.7, 0.7).

- `"inferred"`: calls
  [`infer_confounding()`](https://seantbresnahan.com/iconic/reference/infer_confounding.md)
  to estimate parameters from the data. Parameters that cannot be
  inferred fall back to defaults with warnings.

- `"manual"`: use the explicitly supplied arguments.

## Defaults

|  |  |  |
|----|----|----|
| **Parameter** | **Default** | **Source** |
| `confounding` | "default" | Use DGP defaults below |
| `gan_epochs` | 100 | Texture-model training budget |
| `rho_G1_grid` | c(0,0.1,0.2,0.3,0.5) | Independence-violation sweep |
| `rho_G2_grid` | c(0,0.1,0.2,0.3,0.5) | Independence-violation sweep |
| `n_iter` | 30 | Replicates per grid cell |
| `n_features` | 10 | Features per replicate |
| `mo_confounding` | 0.8 | Simulation calibration (delta_mo) |
| `lambda_XM`, `lambda_MY` | shared | Per-path confounder loadings |
| `omega_1, omega_2` | c(0.3,0.7,1.0) | NC coverage (swept on the diagonal) |
| `bias_threshold` | 0.10 | Tipping-point threshold |

## Examples

``` r
if (check_torch_setup()) {
  data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
    W = matrix(rnorm(100 * 10), 10, 100))
  sens <- iconic_sensitivity(data, n_iter = 2, gan_epochs = 5,
    rho_G1_grid = c(0, 0.2), rho_G2_grid = c(0, 0.2))
  print(sens)
}
```
