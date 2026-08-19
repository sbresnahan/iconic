# Prospective bias-reduction analysis for data without instruments or negative controls

For users who have exposure (X), mediator (M), and outcome (Y) but no
genetic instruments (G, Gm) or negative controls (W), this function
simulates what estimates they could expect if they were to collect such
data.

## Usage

``` r
iconic_prospect(
  data,
  trained_gan = NULL,
  confounding = c("default", "inferred", "manual"),
  gan_epochs = 100,
  gamma_G_grid = c(0.2, 0.4, 0.6, 0.8, 1),
  target_gamma_G = 0.6,
  n_iter = 30,
  n_features = 10,
  mo_confounding = 0.8,
  phi = 0.8,
  lambda_XM = NULL,
  lambda_MY = NULL,
  omega_1 = 0.7,
  omega_2 = 0.7,
  rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
  omega_grid_rho = c(0.3, 0.7, 1),
  run_rho_sweep = TRUE,
  bias_threshold = 0.1,
  base_seed = 500,
  verbose = FALSE,
  allow_no_proxy = TRUE,
  outcome_type = c("continuous", "survival", "binary"),
  effect_scale = c("loghr", "rmst", "logor", "riskdiff"),
  surv_h0 = 0.1,
  surv_event_frac = 0.6,
  surv_censor_rate = NULL,
  bin_prev = 0.5
)
```

## Arguments

- data:

  An `iconic_data` object (must have X and Y; M is required for
  mediation prospect).

- trained_gan:

  Optional `iconic_gan` from
  [`train_gan_on_real_data()`](https://seantbresnahan.com/iconic/reference/train_gan_on_real_data.md).
  If `NULL`, a texture model is auto-trained from `data`.

- confounding:

  Confounding parameter source: `"default"`, `"inferred"`, or
  `"manual"`. Default `"default"`.

- gan_epochs:

  Epochs for auto-trained GAN. Default 100.

- gamma_G_grid:

  Instrument strength values to sweep (Phase 1). Default
  `c(0.2, 0.4, 0.6, 0.8, 1.0)`.

- target_gamma_G:

  Target instrument strength for Phase 2. Default 0.6 (matching the DGP
  default).

- n_iter:

  Replicates per grid cell. Default 30.

- n_features:

  Features per replicate. Default 10.

- mo_confounding:

  Assumed M-O confounding strength. Default 0.8.

- phi:

  Assumed mediator instrument strength. Default 0.8.

- lambda_XM:

  Optional per-path confounder loading vector (X-\>M path).

- lambda_MY:

  Optional per-path confounder loading vector (M-\>Y path). Default
  TRUE.

- omega_1, omega_2:

  Assumed NC coverage. Default 0.7.

- rho_G1_grid, rho_G2_grid:

  Instrument-exogeneity violation values (correlation of each instrument
  with its path's confounder composite) to sweep in Phase 3 at the
  target instrument strength. Default `c(0, 0.1, 0.2, 0.3, 0.5)`,
  matching
  [`iconic_sensitivity`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md).

- omega_grid_rho:

  Negative-control coverage values swept jointly with the rho grid in
  Phase 3, on the diagonal (`omega_1 == omega_2`). Default
  `c(0.3, 0.7, 1.0)`, matching
  [`iconic_sensitivity`](https://seantbresnahan.com/iconic/reference/iconic_sensitivity.md).
  This is independent of the Phase 1 `omega_1`/`omega_2` sweep.

- run_rho_sweep:

  Logical: run the Phase 3 robustness sweep (default `TRUE`). The sweep
  crosses instrument-exogeneity violations (rho) with negative-control
  coverage (omega) and feeds the resulting degradation surface into
  [`iconic_recommend()`](https://seantbresnahan.com/iconic/reference/iconic_recommend.md)
  so the recommended estimator is chosen by robustness to both imperfect
  instruments and weakening controls, rather than by eligibility alone.
  Set `FALSE` to skip (faster, but the recommendation then falls back to
  a single-point / eligibility ranking).

- bias_threshold:

  Tipping-point threshold. Default 0.10.

- base_seed:

  Base RNG seed. Default 500.

- verbose:

  Logical: print progress messages during the sweep. Default `FALSE`
  (quiet).

- allow_no_proxy:

  Logical: when `TRUE` (default), proceed with the prospective sweep
  even if the data already has instruments/NCs (with a message). When
  `FALSE`, error if the data already has IV+NC (use iconic_estimate
  instead).

- outcome_type:

  `"continuous"` (default), `"survival"`, or `"binary"`. Threads through
  to the simulation DGP.

- effect_scale:

  `"loghr"` (default), `"rmst"`, `"logor"`, or `"riskdiff"`.
  `"loghr"`/`"rmst"` apply to survival outcomes; `"logor"`/`"riskdiff"`
  apply to binary outcomes (an inert default is remapped with a
  message).

- surv_h0, surv_event_frac, surv_censor_rate:

  Survival DGP parameters See
  [`generate_toy_data`](https://seantbresnahan.com/iconic/reference/generate_toy_data.md).

- bin_prev:

  Target prevalence for the binary DGP. Default 0.5.

## Value

An `iconic_prospect` S3 object: a named list with `$strength_surface`
(Phase 1: gamma_G x method estimates), `$prospective` (Phase 2: full
simulation at target strength), `$rho_surface` (Phase 3: rho_G1 x rho_G2
exogeneity-robustness surface at the target strength; `NULL` when
`run_rho_sweep = FALSE`), `$summary`, `$recommendation`,
`$texture_source`, and `$inferred_confounding` (when
`confounding = "inferred"`).

## Details

Phase 1 sweeps instrument strength (`gamma_G`) to show how estimates
converge as the instrument strengthens. Phase 2 runs a full prospective
simulation at a target strength, generating synthetic instruments and
NCs calibrated to the user's sample size and confounding level.

## When to use

Use this as a **bias-reduction prospective** when you have an
observational exposure-mediator-outcome triplet but lack the genetic
instruments and negative controls that the core estimators require. It
quantifies the **relative bias improvement** you could expect by
collecting such data: the sweep shows how much of the naive confounding
bias is removed as instrument strength and NC coverage increase, letting
you decide whether the marginal gain justifies the cost of genotyping /
profiling the additional assays. It is a planning tool, not an estimator
– it does not produce a causal estimate from your current data, but
tells you what a future instrumented study would yield.

## Texture model

When `trained_gan` is `NULL` and no GAN is attached to `data`, a texture
model is auto-trained from the user's data.

## Confounding calibration

The `confounding` argument controls how the held-fixed confounding
parameters are set. In the prospective setting (no instruments or NCs),
most parameters cannot be inferred and will fall back to defaults with
warnings – this is an honest limitation, not a silent failure.

## Defaults

|  |  |  |
|----|----|----|
| **Parameter** | **Default** | **Source** |
| `confounding` | "default" | Use DGP defaults below |
| `gan_epochs` | 100 | Texture-model training budget |
| `gamma_G_grid` | c(0.2,0.4,0.6,0.8,1.0) | Instrument-strength sweep |
| `target_gamma_G` | 0.6 | DGP default (gamma_G) |
| `n_iter` | 30 | Replicates per grid cell |
| `n_features` | 10 | Features per replicate |
| `mo_confounding` | 0.8 | Simulation calibration (delta_mo) |
| `phi` | 0.8 | Strong mediator instrument assumption |
| `lambda_XM`, `lambda_MY` | shared | Per-path confounder loadings |
| `omega_1, omega_2` | 0.7 | NC coverage (simulation calibration) |
| `rho_G1_grid`, `rho_G2_grid` | c(0,0.1,0.2,0.3,0.5) | Phase 3 exogeneity sweep |
| `omega_grid_rho` | c(0.3,0.7,1.0) | Phase 3 NC-coverage sweep (diagonal) |
| `run_rho_sweep` | TRUE | Run Phase 3 robustness sweep |
| `bias_threshold` | 0.10 | Tipping-point threshold |
| `allow_no_proxy` | TRUE | Proceed in prospective setting |

## Examples

``` r
if (check_torch_setup()) {
  data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
    M = rnorm(100))
  result <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
    gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
  print(result)
}
```
