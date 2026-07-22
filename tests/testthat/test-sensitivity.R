test_that("gan_sensitivity returns a scenario x method summary with five methods", {
  sens <- gan_sensitivity(NULL, conf_grid = c(0.3, 0.8), coverage_grid = c(0.5, 1),
                          k_grid = 1, n_iter = 5, n_samples = 250, n_features = 5,
                          base_seed = 11)
  expect_named(sens, c("summary", "grid"))
  expect_true(all(c("conf_strength", "coverage", "k", "true_total", "method",
                    "bias", "rmse", "power") %in% names(sens$summary)))
  # 2 conf x 2 coverage x 1 k x 5 methods
  expect_equal(nrow(sens$summary), 2 * 2 * 1 * 5)
})

test_that("recommend_estimator excludes UNADJ, picks a robust method", {
  sens <- gan_sensitivity(NULL, conf_grid = c(0.3, 0.8), coverage_grid = c(0.5, 1),
                          k_grid = 1, n_iter = 6, n_samples = 300, n_features = 5,
                          base_seed = 12)
  rec <- recommend_estimator(sens)
  expect_named(rec, c("per_scenario", "worst_case", "overall"))
  expect_false("UNADJ" %in% rec$worst_case$method)
  expect_true(rec$overall %in% c("DIRECT", "COCA", "IV2SLS", "PGC"))
  # IV2SLS is designed to be the robust winner here
  expect_equal(rec$overall, "IV2SLS")
})

test_that("nc_validity_check flags under-identification when k exceeds valid controls", {
  nv <- suppressWarnings(
    nc_validity_check(NULL, coverage_grid = c(0.5, 1), k_grid = c(1, 3),
                      conf_strength = 0.8, n_valid_controls = 1,
                      n_iter = 5, n_samples = 250, n_features = 5, base_seed = 13))
  expect_named(nv, c("summary", "verdict"))
  # k = 1 identified, k = 3 not (only 1 valid control)
  expect_true(all(nv$verdict$identified[nv$verdict$k == 1]))
  expect_false(any(nv$verdict$identified[nv$verdict$k == 3]))
  expect_true(any(grepl("under-identified", nv$verdict$diagnosis)))
  # PGC and COCA should be in the summary (NC-dependent estimators)
  expect_true("PGC" %in% nv$summary$method)
  expect_true("COCA" %in% nv$summary$method)
})
