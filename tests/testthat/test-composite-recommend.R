# Tests for the v0.9.9.2 data-driven composite recommendation:
#  1. .composite_robustness() computes the worst-estimand (min) composite.
#  2. The confidence multiplier discounts bridge-dependent estimators by the
#     graded completeness verdict; instrument-only estimators are not discounted.
#  3. Structurally-naive estimators (UNADJ, DIRECT) are demoted below eligible
#     IV/NC estimators in the ranking.
#  4. iconic_recommend() headline = top eligible non-naive estimator by final
#     score (recovers IV2SLS2 on the GUSTO-like scenario).
#  5. completeness_penalty is user-overridable.
#  6. Requirement labels are threshold-aware (interpolate actual min_f).
#  7. min_f pass-through: auto-diagnose honours a non-default min_f.

# ── Helper: build a minimal iconic_diagnosis stub ──
.make_diag <- function(completeness = "weak-capture", min_f = 5,
                       eligible = rep(TRUE, 8)) {
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")
  elig <- data.frame(
    estimator = methods,
    eligible = eligible,
    a2_required = methods %in% c("IV2SLS", "PGC", "IV2SLS2", "PGC2", "PGC2Gm"),
    reason_code = ifelse(eligible, "OK", "NEED_DATA"),
    reason = ifelse(eligible, "ok", "missing data"),
    stringsAsFactors = FALSE
  )
  comp <- if (is.null(completeness)) NULL else
    list(completeness = completeness, n_valid_controls = 20, k = 1,
         dim_W = 164, n_valid_coca = 20)
  structure(list(
    instrument_strength = list(weak_G = FALSE, weak_Gm = FALSE),
    completeness = comp,
    eligibility = elig,
    min_f = min_f
  ), class = c("iconic_diagnosis", "list"))
}

# ── Helper: build a minimal iconic_sensitivity stub with chosen scores ──
# nde_bias / nie_bias are max|bias| per method; combined criterion with no
# coverage columns yields score = 1/(1 + norm(bias)). We instead inject
# coverage-free surfaces and check the composite on the resulting scores.
.make_sens <- function(methods, nde_bias, nie_bias) {
  surface <- data.frame(
    rho_G1 = 0, rho_G2 = 0, omega_1 = 0.7, omega_2 = 0.7,
    method = methods,
    NDE_bias = nde_bias, NIE_bias = nie_bias,
    stringsAsFactors = FALSE
  )
  structure(list(surface = surface), class = c("iconic_sensitivity", "list"))
}

# Minimal iconic_data stub (only fields recommend touches).
.make_data <- function() {
  structure(list(is_mediation = TRUE), class = c("iconic_data", "list"))
}

test_that(".composite_robustness computes worst-estimand composite", {
  ranking <- data.frame(
    estimator = c("IV2SLS", "IV2SLS2"),
    robustness_NDE = c(0.994, 0.892),
    robustness_NIE = c(0.499, 0.739),
    stringsAsFactors = FALSE
  )
  comp <- iconic:::.composite_robustness(ranking, .make_diag("satisfied"),
                                         NULL)
  # worst-estimand = min
  expect_equal(comp$composite, c(0.499, 0.739), tolerance = 1e-8)
  # satisfied verdict -> multiplier 1 for everyone
  expect_equal(comp$confidence_mult, c(1, 1))
  expect_equal(comp$final_score, comp$composite, tolerance = 1e-8)
})

test_that("confidence multiplier discounts bridge-dependent estimators only", {
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")
  ranking <- data.frame(
    estimator = methods,
    robustness_NDE = rep(0.8, 8),
    robustness_NIE = rep(0.8, 8),
    stringsAsFactors = FALSE
  )
  comp <- iconic:::.composite_robustness(ranking, .make_diag("weak-capture"),
                                         NULL)
  # bridge-dependent: DIRECT, COCA, PGC, PGC2, PGC2Gm -> 0.5
  # instrument-only / naive: UNADJ, IV2SLS, IV2SLS2 -> 1.0
  expect_equal(comp$confidence_mult,
               c(1.0, 0.5, 0.5, 1.0, 0.5, 1.0, 0.5, 0.5))
  expect_equal(comp$final_score, 0.8 * comp$confidence_mult, tolerance = 1e-8)
})

test_that("completeness verdicts map to documented penalties", {
  methods <- c("PGC", "IV2SLS2")
  ranking <- data.frame(estimator = methods, robustness_NDE = 0.8,
                        robustness_NIE = 0.8, stringsAsFactors = FALSE)
  f <- function(v) iconic:::.composite_robustness(ranking, .make_diag(v), NULL)$confidence_mult[1]
  expect_equal(f("satisfied"), 1.0)
  expect_equal(f("borderline"), 0.7)
  expect_equal(f("weak-capture"), 0.5)
  expect_equal(f("under-identified"), 0)
  # IV2SLS2 (instrument-only) is never discounted
  g <- function(v) iconic:::.composite_robustness(ranking, .make_diag(v), NULL)$confidence_mult[2]
  expect_equal(g("weak-capture"), 1.0)
})

test_that("completeness_penalty is user-overridable", {
  ranking <- data.frame(estimator = c("PGC", "IV2SLS2"),
                        robustness_NDE = 0.8, robustness_NIE = 0.8,
                        stringsAsFactors = FALSE)
  comp <- iconic:::.composite_robustness(
    ranking, .make_diag("weak-capture"),
    completeness_penalty = c("weak-capture" = 0.9))
  expect_equal(comp$confidence_mult[1], 0.9)  # PGC overridden
  expect_equal(comp$confidence_mult[2], 1.0)  # IV2SLS2 untouched
})

test_that("NULL completeness leaves multiplier at 1 (bridge already ineligible)", {
  ranking <- data.frame(estimator = c("PGC", "IV2SLS2"),
                        robustness_NDE = 0.8, robustness_NIE = 0.8,
                        stringsAsFactors = FALSE)
  comp <- iconic:::.composite_robustness(ranking, .make_diag(NULL), NULL)
  expect_equal(comp$confidence_mult, c(1.0, 1.0))
})

test_that("non-finite bias ranks last (survival COCA case)", {
  ranking <- data.frame(estimator = c("COCA", "IV2SLS2"),
                        robustness_NDE = c(-Inf, 0.7),
                        robustness_NIE = c(0.6, 0.7),
                        stringsAsFactors = FALSE)
  comp <- iconic:::.composite_robustness(ranking, .make_diag("satisfied"), NULL)
  expect_equal(comp$final_score[1], -Inf)
  expect_true(comp$final_score[2] > comp$final_score[1])
})

test_that("headline recovers IV2SLS2 and demotes naive estimators", {
  # GUSTO-like robustness profile
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")
  # max|bias| chosen so combined score (no coverage) = 1/(1+norm) reproduces
  # the qualitative ordering; we just need IV2SLS high-NDE/low-NIE, IV2SLS2
  # balanced, PGC high-NIE but bridge-dependent.
  sens <- .make_sens(methods,
                     nde_bias = c(0.30, 0.10, 0.30, 0.001, 0.10, 0.05, 0.02, 0.01),
                     nie_bias = c(0.30, 0.02, 0.20, 0.30, 0.02, 0.10, 0.20, 0.20))
  diag <- .make_diag("weak-capture", min_f = 5,
                     eligible = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  rec <- iconic_recommend(.make_data(), diagnosis = diag,
                          sensitivity = sens, auto_sensitivity = FALSE)
  # Headline is the composite winner among eligible non-naive estimators.
  expect_equal(rec$recommended, "IV2SLS2")
  # Naive estimators never out-rank an eligible IV/NC estimator.
  rk <- rec$ranking
  elig <- rk[rk$eligible, ]
  naive_pos <- which(elig$estimator %in% c("UNADJ", "DIRECT"))
  ivnc_pos <- which(elig$estimator %in% c("IV2SLS", "IV2SLS2", "PGC", "COCA"))
  expect_true(all(ivnc_pos < min(naive_pos)))
  # ranking carries the composite columns
  expect_true(all(c("composite", "confidence_mult", "final_score") %in% names(rk)))
})

test_that("per-estimand recommendations are retained", {
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")
  sens <- .make_sens(methods,
                     nde_bias = c(0.30, 0.10, 0.30, 0.001, 0.10, 0.05, 0.02, 0.01),
                     nie_bias = c(0.30, 0.02, 0.20, 0.30, 0.02, 0.10, 0.20, 0.20))
  diag <- .make_diag("weak-capture", min_f = 5,
                     eligible = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  rec <- iconic_recommend(.make_data(), diagnosis = diag,
                          sensitivity = sens, auto_sensitivity = FALSE)
  expect_true(!is.null(rec$recommended_NDE))
  expect_true(!is.null(rec$recommended_NIE))
  # NDE-best should be IV2SLS (lowest NDE bias), NIE-best PGC or DIRECT
  expect_equal(rec$recommended_NDE, "IV2SLS")
})

test_that("requirement labels are threshold-aware", {
  reqs5 <- iconic:::.estimator_requirements(5)
  reqs10 <- iconic:::.estimator_requirements(10)
  expect_match(reqs5[["IV2SLS"]], "F>=5", fixed = TRUE)
  expect_match(reqs10[["IV2SLS"]], "F>=10", fixed = TRUE)
  expect_match(reqs5[["IV2SLS2"]], "F>=5", fixed = TRUE)
  # UNADJ has no F label
  expect_false(grepl("F>=", reqs5[["UNADJ"]]))
})

test_that("rationale uses the diagnosis min_f, not a hardcoded value", {
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")
  sens <- .make_sens(methods,
                     nde_bias = rep(0.1, 8), nie_bias = rep(0.1, 8))
  diag <- .make_diag("satisfied", min_f = 5)
  rec <- iconic_recommend(.make_data(), diagnosis = diag,
                          sensitivity = sens, auto_sensitivity = FALSE)
  iv2sls_rat <- rec$ranking$rationale[rec$ranking$estimator == "IV2SLS"]
  expect_match(iv2sls_rat, "F>=5", fixed = TRUE)
  expect_false(grepl("F>=10", iv2sls_rat, fixed = TRUE))
})
