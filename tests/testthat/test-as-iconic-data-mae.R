# Tests for the MultiAssayExperiment interface:
# as_iconic_data.MultiAssayExperiment() and as_MultiAssayExperiment().

.mk_mat <- function(p, pre, ids) {
    matrix(rnorm(p * length(ids)), p, length(ids),
        dimnames = list(paste0(pre, seq_len(p)), ids))
}
.se <- function(m, name = "value") {
    SummarizedExperiment::SummarizedExperiment(
        assays = stats::setNames(list(m), name))
}

.toy_mae <- function(n = 80, seed = 11) {
    set.seed(seed)
    ids <- paste0("S", seq_len(n))
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(
            rna    = .se(.mk_mat(8, "gene", ids), "expr"),
            methyl = .se(.mk_mat(1, "cg", ids), "m"),
            nc     = .se(.mk_mat(8, "PC", ids), "pcs"),
            nc_xm  = .se(.mk_mat(3, "xm", ids)),
            nc_my  = .se(.mk_mat(3, "my", ids))),
        colData = S4Vectors::DataFrame(
            bmi = rnorm(n), prs = rnorm(n), cis_score = rnorm(n),
            age = rnorm(n), time = rexp(n, 0.1), event = rbinom(n, 1, 0.6),
            case = rbinom(n, 1, 0.4), row.names = ids))
}

test_that("MAE method matches a direct iconic_data() call", {
    mae <- .toy_mae()
    cd <- MultiAssayExperiment::colData(mae)
    ex <- MultiAssayExperiment::experiments(mae)
    d_mae <- as_iconic_data(mae, exposure = "bmi", outcome = "rna",
        mediator = "methyl", instrument = "prs",
        mediator_instrument = "cis_score",
        negative_controls = "nc", covariates = "age")
    d_dir <- iconic_data(
        X = cd$bmi, Y = SummarizedExperiment::assay(ex$rna),
        M = SummarizedExperiment::assay(ex$methyl), G = cd$prs,
        Gm = cd$cis_score, W = SummarizedExperiment::assay(ex$nc),
        covariates = data.frame(age = cd$age))
    expect_s3_class(d_mae, "iconic_data")
    for (f in c("X", "Y", "M", "G", "Gm", "W", "W1", "W2", "n",
        "n_features", "n_mediators", "has_path_nc",
        "feature_names", "mediator_names"))
        expect_equal(d_mae[[f]], d_dir[[f]], ignore_attr = TRUE, info = f)
    expect_equal(d_mae$sample_names, rownames(cd))

    e_mae <- suppressMessages(suppressWarnings(
        iconic_estimate(d_mae, methods = c("UNADJ", "IV2SLS2"))))
    e_dir <- suppressMessages(suppressWarnings(
        iconic_estimate(d_dir, methods = c("UNADJ", "IV2SLS2"))))
    expect_equal(e_mae$NDE, e_dir$NDE)
    expect_equal(e_mae$NIE, e_dir$NIE)
})

test_that("MAE method aligns incomplete and reordered experiments", {
    set.seed(12)
    n <- 60
    ids <- paste0("P", seq_len(n))
    y <- .mk_mat(4, "gene", ids)
    m <- .mk_mat(1, "cg", ids)
    # mediator assay: samples 1-5 missing, columns reversed, own colnames
    m_sub <- m[, rev(ids[-(1:5)]), drop = FALSE]
    colnames(m_sub) <- paste0("meth_", colnames(m_sub))
    smap <- rbind(
        data.frame(assay = "rna", primary = ids, colname = ids),
        data.frame(assay = "methyl", primary = rev(ids[-(1:5)]),
            colname = colnames(m_sub)))
    smap$assay <- factor(smap$assay, levels = c("rna", "methyl"))
    cdat <- S4Vectors::DataFrame(bmi = rnorm(n), prs = rnorm(n),
        row.names = ids)
    cdat$prs[10] <- NA # dropped: incomplete required colData
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(rna = .se(y), methyl = .se(m_sub)),
        colData = cdat, sampleMap = S4Vectors::DataFrame(smap))
    expect_message(
        d <- as_iconic_data(mae, exposure = "bmi", outcome = "rna",
            mediator = "methyl", instrument = "prs",
            scale = FALSE),
        "kept 54 of 60")
    keep <- setdiff(ids[-(1:5)], "P10")
    expect_equal(d$sample_names, keep)
    expect_equal(unname(d$Y), unname(y[, keep]))
    expect_equal(as.numeric(d$M), as.numeric(m[, keep]))
    expect_equal(d$X, as.numeric(cdat[keep, "bmi"]))
})

test_that("MAE method supports path-specific W1/W2 and lone panels", {
    mae <- .toy_mae()
    # 3-column panels vs 8 outcomes -> iconic_data() recycles (expected)
    d <- suppressWarnings(
        as_iconic_data(mae, exposure = "bmi", outcome = "rna",
            mediator = "methyl", instrument = "prs",
            mediator_instrument = "cis_score",
            W1 = "nc_xm", W2 = "nc_my"))
    expect_true(d$has_path_nc)
    expect_equal(nrow(d$W1), 3)
    expect_false(isTRUE(all.equal(d$W1, d$W2)))
    d_lone <- suppressWarnings(
        as_iconic_data(mae, exposure = "bmi", outcome = "rna",
            mediator = "methyl", instrument = "prs",
            W2 = "nc_my"))
    expect_false(d_lone$has_path_nc)
    expect_null(d_lone$W1)
    expect_false(is.null(d_lone$W2))
})

test_that("MAE method supports survival and binary outcomes", {
    mae <- .toy_mae()
    ds <- as_iconic_data(mae, exposure = "bmi", mediator = "methyl",
        instrument = "prs", outcome_type = "survival",
        surv_time = "time", surv_event = "event")
    expect_identical(ds$outcome_type, "survival")
    expect_length(ds$surv_time, 80)
    db <- as_iconic_data(mae, exposure = "bmi", mediator = "methyl",
        instrument = "prs", outcome_type = "binary",
        bin_outcome = "case")
    expect_identical(db$outcome_type, "binary")
    expect_true(all(db$Y_bin %in% c(0, 1)))
})

test_that("MAE method gives informative errors", {
    mae <- .toy_mae()
    expect_error(as_iconic_data(mae, exposure = "bmi", outcome = "nope"),
        "not an experiment")
    expect_error(as_iconic_data(mae, exposure = "bmi"),
        "`outcome`")
    expect_error(as_iconic_data(mae, exposure = "bmi", outcome = "rna",
        instrument = "missing_col"),
    "neither an experiment nor a colData column")
    expect_error(as_iconic_data(mae, exposure = "bmi", outcome = "rna",
        instrument = "rna"),
    "single variable")
    expect_error(as_iconic_data(mae, outcome = "rna"), "`exposure` is required")
})

test_that("as_MultiAssayExperiment round-trips iconic_data objects", {
    set.seed(13)
    n <- 70
    d <- iconic_data(X = rnorm(n), Y = matrix(rnorm(n * 5), 5, n),
        M = rnorm(n), G = rnorm(n), Gm = rnorm(n),
        W = matrix(rnorm(n * 5), 5, n),
        covariates = data.frame(age = rnorm(n)))
    mae <- as_MultiAssayExperiment(d)
    expect_s4_class(mae, "MultiAssayExperiment")
    expect_setequal(names(mae), c("Y", "M", "Gm", "W"))
    expect_true(all(c("X", "G", "age") %in%
        colnames(MultiAssayExperiment::colData(mae))))
    d2 <- as_iconic_data(mae)
    for (f in c("X", "Y", "M", "G", "Gm", "W", "W1", "W2", "covariates",
        "has_path_nc", "is_mediation", "scaling", "feature_names"))
        expect_equal(d2[[f]], d[[f]], ignore_attr = TRUE, info = f)

    # path-specific panels
    dp <- iconic_data(X = rnorm(n), Y = matrix(rnorm(n * 3), 3, n),
        M = rnorm(n), G = rnorm(n), Gm = rnorm(n),
        W1 = matrix(rnorm(n * 2), 2, n),
        W2 = matrix(rnorm(n * 2), 2, n))
    mp <- as_MultiAssayExperiment(dp)
    expect_setequal(names(mp), c("Y", "M", "Gm", "W1", "W2"))
    dp2 <- as_iconic_data(mp)
    for (f in c("W", "W1", "W2", "has_path_nc"))
        expect_equal(dp2[[f]], dp[[f]], ignore_attr = TRUE, info = f)

    # lone panel keeps has_path_nc = FALSE
    dl <- suppressWarnings(
        iconic_data(X = rnorm(n), Y = matrix(rnorm(n * 3), 3, n),
            M = rnorm(n), G = rnorm(n), Gm = rnorm(n),
            W2 = matrix(rnorm(n * 2), 2, n)))
    dl2 <- suppressWarnings(as_iconic_data(as_MultiAssayExperiment(dl)))
    expect_false(dl2$has_path_nc)
    expect_equal(dl2$W2, dl$W2, ignore_attr = TRUE)
})

test_that("as_MultiAssayExperiment handles survival and binary outcomes", {
    set.seed(14)
    n <- 60
    ds <- iconic_data(X = rnorm(n), outcome_type = "survival",
        surv_time = rexp(n), surv_event = rbinom(n, 1, 0.6),
        M = rnorm(n), G = rnorm(n))
    ms <- as_MultiAssayExperiment(ds)
    expect_false("Y" %in% names(ms))
    ds2 <- as_iconic_data(ms)
    expect_equal(ds2$surv_time, ds$surv_time)
    expect_equal(ds2$surv_event, ds$surv_event)
    db <- iconic_data(X = rnorm(n), Y = rbinom(n, 1, 0.4),
        outcome_type = "binary", M = rnorm(n), G = rnorm(n))
    db2 <- as_iconic_data(as_MultiAssayExperiment(db))
    expect_equal(db2$Y_bin, db$Y_bin)
    expect_identical(db2$outcome_type, "binary")
})

test_that("as_MultiAssayExperiment rejects non-iconic_data input", {
    expect_error(as_MultiAssayExperiment(list(a = 1)), "iconic_data")
})
