# Composite data-driven robustness score (internal)

Combines per-estimand robustness into a single score and discounts each
estimator by the graded strength of the diagnostic assumptions it
depends on. Three steps:

1.  **Worst-estimand composite**:
    `pmin(robustness_NDE, robustness_NIE)`. An estimator is only as
    trustworthy as its weakest estimand, so a strong NDE cannot
    compensate for a weak NIE (the usual object of a mediation
    analysis).

2.  **Confidence multiplier**: each estimator's composite is multiplied
    by a factor in \\(0, 1\]\\ read from the graded completeness verdict
    (`diagnosis$completeness$completeness`). Only estimators whose
    identification routes through the negative-control bridge (DIRECT,
    COCA, PGC, PGC2, PGC2Gm) are discounted; instrument-only estimators
    (IV2SLS, IV2SLS2) and UNADJ keep their full score. When no
    completeness object exists (no NC panel), bridge-dependent
    estimators are already ineligible, so the multiplier is never
    misapplied.

3.  **Final score**: `composite * confidence_mult`.

## Usage

``` r
.composite_robustness(ranking, diagnosis, completeness_penalty)
```

## Arguments

- ranking:

  Data frame with `estimator`, `robustness_NDE`, `robustness_NIE`
  columns.

- diagnosis:

  An `iconic_diagnosis` object (may carry a NULL `$completeness`).

- completeness_penalty:

  Named numeric vector mapping completeness verdicts to multipliers in
  \\\[0, 1\]\\.

## Value

Data frame with `composite`, `confidence_mult`, and `final_score`
aligned to `ranking$estimator`.
