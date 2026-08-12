# Tests for the as_iconic_data() S3 generic and its methods

test_that("as_iconic_data passes iconic_data objects through unchanged", {
  d <- iconic_data(X = rnorm(50), Y = matrix(rnorm(50 * 5), 5, 50),
                   G = rnorm(50), W = matrix(rnorm(50 * 5), 5, 50))
  expect_identical(as_iconic_data(d), d)
})

test_that("as_iconic_data.default accepts an exposure vector + named args", {
  d <- as_iconic_data(rnorm(50), Y = matrix(rnorm(50 * 5), 5, 50),
                      G = rnorm(50), W = matrix(rnorm(50 * 5), 5, 50))
  expect_s3_class(d, "iconic_data")
  expect_equal(d$n, 50)
})

test_that("as_iconic_data.default accepts a load_real_input_data() result", {
  raw <- load_real_input_data(example = TRUE)
  # W panel is smaller than Y -> iconic_data() recycles (expected)
  d <- suppressWarnings(as_iconic_data(raw))
  expect_s3_class(d, "iconic_data")
  expect_equal(d$n, raw$n_samples)
})

test_that("as_iconic_data.SummarizedExperiment maps assays and colData", {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(31)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = matrix(rnorm(20 * 60), 20, 60,
                                dimnames = list(paste0("gene", 1:20),
                                                paste0("S", 1:60)))),
    colData = S4Vectors::DataFrame(
      bmi = rnorm(60), prs = rnorm(60),
      nc1 = rnorm(60), nc2 = rnorm(60), age = rnorm(60))
  )
  # 2 negative controls vs 20 outcomes -> iconic_data() recycles (expected)
  d <- suppressWarnings(
    as_iconic_data(se, assay = "expr", exposure = "bmi",
                   instrument = "prs",
                   negative_controls = c("nc1", "nc2"),
                   covariates = "age")
  )
  expect_s3_class(d, "iconic_data")
  expect_equal(d$n, 60)
  expect_equal(ncol(d$Y), 60)
  expect_equal(nrow(d$Y), 20)
  expect_length(d$X, 60)
  expect_true(!is.null(d$G))
  expect_true(!is.null(d$W))
  expect_true("age" %in% names(d$covariates))
})

test_that("as_iconic_data.SummarizedExperiment supports survival outcomes", {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(32)
  se <- SummarizedExperiment::SummarizedExperiment(
    colData = S4Vectors::DataFrame(
      bmi = rnorm(50), prs = rnorm(50),
      time = rexp(50, 0.1), event = rbinom(50, 1, 0.6))
  )
  d <- as_iconic_data(se, assay = NULL, exposure = "bmi",
                      instrument = "prs",
                      surv_time = "time", surv_event = "event",
                      outcome_type = "survival")
  expect_s3_class(d, "iconic_data")
  expect_equal(d$outcome_type, "survival")
})

test_that("as_iconic_data.SummarizedExperiment errors on missing columns", {
  skip_if_not_installed("SummarizedExperiment")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = matrix(rnorm(10 * 20), 10, 20)),
    colData = S4Vectors::DataFrame(bmi = rnorm(20))
  )
  expect_error(as_iconic_data(se, assay = "expr", exposure = "missing_col"))
})
