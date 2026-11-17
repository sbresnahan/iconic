# ============================================================
# mae.R: MultiAssayExperiment interface.
#
# as_iconic_data.MultiAssayExperiment() – MAE -> iconic_data
# as_MultiAssayExperiment()             – iconic_data -> MAE
#
# MultiAssayExperiment is the user-facing Bioconductor container:
# omics panels (outcome Y, mediator M, negative controls W / W1 / W2,
# mediator instruments Gm) are experiments, and sample-level
# variables (exposure X, instrument G, covariates, survival / binary
# outcome) are colData columns. Samples are matched across
# experiments through the sampleMap. iconic_data is the derived,
# standardised (scaled, aligned, validated) analysis view consumed by
# the estimators.
# ============================================================


#' @rdname as_iconic_data
#' @param outcome For the MultiAssayExperiment method: name of the
#' experiment holding the outcome panel Y (features x samples).
#' Required when `outcome_type = "continuous"`.
#' @param mediator Optional: name of the experiment (or character
#' vector of `colData` columns) holding the mediator panel M.
#' @param W1,W2 Optional: names of experiments (or character vectors of
#' `colData` columns) holding the path-specific negative-control panels
#' for the X->M (`W1`) and M->Y (`W2`) paths.
#' @param verbose Logical: report samples dropped during alignment.
#' Default `TRUE`.
#'
#' @details
#' **MultiAssayExperiment method.** Each role (`outcome`, `mediator`,
#' `negative_controls`, `W1`, `W2`, `mediator_instrument`, and
#' optionally `exposure` / `instrument`) names either one experiment of
#' `input` or one or more `colData` columns; experiment names take
#' precedence. For SummarizedExperiment experiments the assay selected
#' by `assay` is used; matrix experiments are used directly. Samples are
#' matched through `sampleMap(input)` to primary IDs (the rownames of
#' `colData`). Only primaries present in every experiment used, and
#' with non-missing values in the sample-level `colData` columns used,
#' are kept. Experiments with more than one column per primary (e.g.
#' technical replicates) must be collapsed first, for example with
#' [MultiAssayExperiment::mergeReplicates()]. The aligned components are
#' passed to [iconic_data()], so validation and scaling are identical to
#' calling it directly.
#'
#' @method as_iconic_data MultiAssayExperiment
#' @export
#'
#' @examples
#' # From a MultiAssayExperiment: omics panels are experiments,
#' # sample-level variables are colData columns
#' set.seed(1)
#' ids <- paste0("S", 1:80)
#' mk <- function(p, pre) matrix(rnorm(p * 80), p, 80,
#'     dimnames = list(paste0(pre, seq_len(p)), ids))
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'     experiments = list(
#'         rna = SummarizedExperiment::SummarizedExperiment(
#'             assays = list(expr = mk(10, "gene"))),
#'         methyl = SummarizedExperiment::SummarizedExperiment(
#'             assays = list(m = mk(1, "cg"))),
#'         nc = SummarizedExperiment::SummarizedExperiment(
#'             assays = list(pcs = mk(10, "PC")))),
#'     colData = S4Vectors::DataFrame(
#'         bmi = rnorm(80), prs = rnorm(80), cis_score = rnorm(80),
#'         age = rnorm(80), row.names = ids))
#' data <- as_iconic_data(mae, exposure = "bmi", outcome = "rna",
#'     mediator = "methyl", instrument = "prs",
#'     mediator_instrument = "cis_score", negative_controls = "nc",
#'     covariates = "age")
#' print(data)
as_iconic_data.MultiAssayExperiment <- function(input, exposure,
                                                outcome = NULL,
                                                mediator = NULL,
                                                instrument = NULL,
                                                mediator_instrument = NULL,
                                                negative_controls = NULL,
                                                W1 = NULL, W2 = NULL,
                                                covariates = NULL,
                                                surv_time = NULL,
                                                surv_event = NULL,
                                                bin_outcome = NULL,
                                                outcome_type = c("continuous",
                                                    "survival",
                                                    "binary"),
                                                assay = 1L,
                                                verbose = TRUE, ...) {
    ## Roles recorded by as_MultiAssayExperiment() are used when the
    ## caller does not specify them (round trip).
    meta <- S4Vectors::metadata(input)$iconic
    restore <- FALSE
    if (missing(exposure)) {
        if (is.null(meta$roles))
            stop("`exposure` is required (name of a colData column or ",
                "experiment holding the exposure).")
        r <- meta$roles
        exposure <- r$exposure; outcome <- r$outcome; mediator <- r$mediator
        instrument <- r$instrument; mediator_instrument <- r$mediator_instrument
        negative_controls <- r$negative_controls; W1 <- r$W1; W2 <- r$W2
        covariates <- r$covariates; surv_time <- r$surv_time
        surv_event <- r$surv_event; bin_outcome <- r$bin_outcome
        outcome_type <- meta$outcome_type
        restore <- TRUE
    }
    outcome_type <- match.arg(outcome_type, c("continuous", "survival",
        "binary"))

    exps <- MultiAssayExperiment::experiments(input)
    exp_names <- names(exps)
    cd <- MultiAssayExperiment::colData(input)
    cd_names <- colnames(cd)
    sm <- MultiAssayExperiment::sampleMap(input)

    ## Classify each role as an experiment or a set of colData columns.
    classify <- function(spec, label) {
        if (is.null(spec)) return(NULL)
        spec <- as.character(spec)
        if (length(spec) == 1L && spec %in% exp_names)
            return(list(type = "experiment", name = spec))
        miss <- setdiff(spec, cd_names)
        if (length(miss))
            stop(label, ": '", paste(miss, collapse = "', '"), "' is neither an ",
                "experiment nor a colData column of the MultiAssayExperiment.\n",
                "  Experiments: ", paste(exp_names, collapse = ", "), "\n",
                "  colData columns: ", paste(cd_names, collapse = ", "),
                call. = FALSE)
        list(type = "coldata", cols = spec)
    }

    if (outcome_type == "continuous" && is.null(outcome))
        stop("`outcome` (name of the experiment holding the outcome panel) ",
            "is required when outcome_type = 'continuous'.", call. = FALSE)
    if (outcome_type == "continuous" && !is.null(outcome) &&
        !(length(outcome) == 1L && outcome %in% exp_names))
        stop("outcome: '", paste(outcome, collapse = "', '"), "' is not an ",
            "experiment. Available: ", paste(exp_names, collapse = ", "),
            call. = FALSE)
    if (outcome_type == "survival" && (is.null(surv_time) || is.null(surv_event)))
        stop("When outcome_type = 'survival', supply surv_time and surv_event ",
            "(colData column names).", call. = FALSE)
    if (outcome_type == "binary" && is.null(bin_outcome))
        stop("When outcome_type = 'binary', supply bin_outcome = '<colData ",
            "column>' naming the 0/1 outcome column.", call. = FALSE)

    roles <- list(
        X  = classify(exposure, "exposure"),
        Y  = if (outcome_type == "continuous") classify(outcome, "outcome"),
        M  = classify(mediator, "mediator"),
        G  = classify(instrument, "instrument"),
        Gm = classify(mediator_instrument, "mediator_instrument"),
        W  = classify(negative_controls, "negative_controls"),
        W1 = classify(W1, "W1"),
        W2 = classify(W2, "W2")
    )
    roles <- Filter(Negate(is.null), roles)
    sample_cols <- c(
        unlist(lapply(roles, function(r) if (r$type == "coldata") r$cols)),
        if (outcome_type == "survival") c(surv_time, surv_event),
        if (outcome_type == "binary") bin_outcome,
        covariates)
    miss <- setdiff(sample_cols, cd_names)
    if (length(miss))
        stop("colData column(s) not found: ", paste(miss, collapse = ", "),
            ". Available: ", paste(cd_names, collapse = ", "), call. = FALSE)

    ## Align samples via the sampleMap.
    keep <- rownames(cd)
    used_exps <- unique(unlist(lapply(roles, function(r)
        if (r$type == "experiment") r$name)))
    for (nm in used_exps) {
        smi <- sm[sm$assay == nm, , drop = FALSE]
        if (anyDuplicated(smi$primary))
            stop("Experiment '", nm, "' has more than one column for some ",
                "primary samples. Collapse replicates first, e.g. with ",
                "MultiAssayExperiment::mergeReplicates().", call. = FALSE)
        keep <- keep[keep %in% smi$primary]
    }
    if (length(sample_cols)) {
        cdk <- as.data.frame(cd[keep, unique(sample_cols), drop = FALSE])
        keep <- keep[stats::complete.cases(cdk)]
    }
    n_dropped <- nrow(cd) - length(keep)
    if (n_dropped > 0 && isTRUE(verbose))
        message("as_iconic_data(): kept ", length(keep), " of ", nrow(cd),
            " samples present in all experiments used and complete for ",
            "the colData columns used (", n_dropped, " dropped).")
    if (length(keep) == 0)
        stop("No samples are shared by all experiments used.", call. = FALSE)

    ## Extract one role as a features x samples matrix over `keep`.
    get_matrix <- function(r) {
        if (r$type == "experiment") {
            e <- exps[[r$name]]
            m <- if (methods::is(e, "SummarizedExperiment"))
                SummarizedExperiment::assay(e, assay, withDimnames = TRUE)
            else as.matrix(e)
            smi <- sm[sm$assay == r$name, , drop = FALSE]
            m <- as.matrix(m[, smi$colname[match(keep, smi$primary)], drop = FALSE])
            colnames(m) <- keep
            storage.mode(m) <- "double"
            m
        } else {
            v <- vapply(r$cols, function(cl) {
                x <- cd[keep, cl]
                if (!is.numeric(x))
                    stop("colData column '", cl, "' must be numeric.", call. = FALSE)
                as.numeric(x)
            }, numeric(length(keep)))
            m <- t(matrix(v, nrow = length(keep), dimnames = list(keep, r$cols)))
            m
        }
    }
    get_vector <- function(r, label) {
        m <- get_matrix(r)
        if (nrow(m) != 1L)
            stop(label, " must be a single variable (one colData column or a ",
                "one-row experiment); got ", nrow(m), ".", call. = FALSE)
        stats::setNames(as.numeric(m), keep)
    }
    num_col <- function(cl) as.numeric(cd[keep, cl])

    X <- if (roles[["X"]]$type == "experiment") get_matrix(roles[["X"]])
    else get_vector(roles[["X"]], "exposure")
    Y <- NULL
    if (outcome_type == "continuous") Y <- get_matrix(roles[["Y"]])
    if (outcome_type == "binary") Y <- num_col(bin_outcome)
    M <- if (!is.null(roles[["M"]])) get_matrix(roles[["M"]])
    G <- if (!is.null(roles[["G"]])) get_vector(roles[["G"]], "instrument")
    Gm <- if (!is.null(roles[["Gm"]])) get_matrix(roles[["Gm"]])
    W <- if (!is.null(roles[["W"]])) get_matrix(roles[["W"]])
    Wp1 <- if (!is.null(roles[["W1"]])) get_matrix(roles[["W1"]])
    Wp2 <- if (!is.null(roles[["W2"]])) get_matrix(roles[["W2"]])
    cv <- if (!is.null(covariates))
        as.data.frame(cd[keep, covariates, drop = FALSE])

    out <- iconic_data(
        X = X, Y = Y, M = M, G = G, Gm = Gm, W = W, W1 = Wp1, W2 = Wp2,
        covariates = cv,
        outcome_type = outcome_type,
        surv_time = if (outcome_type == "survival") num_col(surv_time),
        surv_event = if (outcome_type == "survival") num_col(surv_event),
        ...)
    out$sample_names <- keep
    if (restore) {
        if (!is.null(meta$scaling)) out$scaling <- meta$scaling
        if (!is.null(meta$feature_names) && outcome_type != "continuous")
            out$feature_names <- meta$feature_names
    }
    out
}


#' Convert an iconic_data object to a MultiAssayExperiment
#'
#' Exports the standardised analysis data held in an `iconic_data`
#' object as a [MultiAssayExperiment::MultiAssayExperiment], the
#' Bioconductor container for multi-omic data on shared samples. Use it
#' to store, share, subset, or inspect the exact data an ICONIC analysis
#' was run on with standard Bioconductor tools. The result can be passed
#' straight back to [as_iconic_data()] (the role mapping is recorded in
#' the metadata).
#'
#' @param data An `iconic_data` object (from [iconic_data()] or
#' [as_iconic_data()]).
#'
#' @return A `MultiAssayExperiment` with
#' \describe{
#'   \item{experiments}{`SummarizedExperiment`s (assay `value`, features
#'   x samples) named `Y` (outcome panel; continuous outcomes only), `M`
#'   (mediators), `Gm` (mediator instruments), `W`, `W1`, `W2`
#'   (negative-control panels), each present only when supplied.}
#'   \item{colData}{`X` (exposure), `G` (exposure instrument),
#'   `surv_time` / `surv_event` (survival), `Y_bin` (binary outcome),
#'   and the encoded covariates (prefixed `cov_` if a name collides).}
#'   \item{metadata}{`iconic`: the outcome type, scaling parameters,
#'   feature / mediator names, availability flags, the role mapping used
#'   by [as_iconic_data()], and the attached GAN (if any).}
#' }
#' Values are those stored in `data`, i.e. centred and scaled when
#' `iconic_data(scale = TRUE)` (the default); original-scale parameters
#' are kept in `metadata(mae)$iconic$scaling`.
#' @export
#'
#' @examples
#' d <- iconic_data(X = rnorm(60), Y = matrix(rnorm(60 * 5), 5, 60),
#'     M = rnorm(60), G = rnorm(60), Gm = rnorm(60),
#'     W = matrix(rnorm(60 * 5), 5, 60))
#' mae <- as_MultiAssayExperiment(d)
#' mae
#' d2 <- as_iconic_data(mae) # round trip
#' all.equal(d$Y, d2$Y, check.attributes = FALSE)
as_MultiAssayExperiment <- function(data) {
    if (!inherits(data, "iconic_data"))
        stop("data must be an iconic_data object.", call. = FALSE)
    n <- data$n
    sn <- data$sample_names
    if (is.null(sn) || length(sn) != n) sn <- paste0("S", seq_len(n))

    as_se <- function(m, names_vec, prefix) {
        m <- as.matrix(m)
        rn <- if (!is.null(names_vec) && length(names_vec) == nrow(m)) names_vec
        else if (!is.null(rownames(m))) rownames(m)
        else paste0(prefix, seq_len(nrow(m)))
        dimnames(m) <- list(as.character(rn), sn)
        SummarizedExperiment::SummarizedExperiment(assays = list(value = m))
    }

    ## Which NC panels to export so that as_iconic_data() rebuilds the
    ## same configuration (pooled W, path-specific W1/W2, or a lone panel).
    same <- function(a, b) !is.null(a) && !is.null(b) &&
        isTRUE(all.equal(unname(a), unname(b)))
    W <- data$W; W1 <- data$W1; W2 <- data$W2
    exp_W <- exp_W1 <- exp_W2 <- FALSE
    if (!is.null(W1) && !is.null(W2) && data$has_path_nc) {
        if (same(W1, W) && same(W2, W)) {
            exp_W <- TRUE # pooled single panel
        } else {
            exp_W1 <- exp_W2 <- TRUE
            exp_W <- !is.null(W) && !same(W, (W1 + W2) / 2)
        }
    } else if (!is.null(W1) || !is.null(W2)) {
        exp_W1 <- !is.null(W1); exp_W2 <- !is.null(W2) # lone panel
        lone <- if (!is.null(W2)) W2 else W1
        exp_W <- !is.null(W) && !same(W, lone)
    } else if (!is.null(W)) {
        exp_W <- TRUE
    }

    ex <- list()
    if (!is.null(data$Y)) ex$Y <- as_se(data$Y, data$feature_names, "feature_")
    if (!is.null(data$M)) ex$M <- as_se(data$M, data$mediator_names, "mediator_")
    if (!is.null(data$Gm)) ex$Gm <- as_se(data$Gm, NULL, "Gm_")
    if (exp_W) ex$W <- as_se(W, NULL, "W_")
    if (exp_W1) ex$W1 <- as_se(W1, NULL, "W1_")
    if (exp_W2) ex$W2 <- as_se(W2, NULL, "W2_")

    cdf <- S4Vectors::DataFrame(X = data$X, row.names = sn)
    if (!is.null(data$G)) cdf$G <- data$G
    if (identical(data$outcome_type, "survival")) {
        cdf$surv_time <- data$surv_time
        cdf$surv_event <- data$surv_event
    }
    if (identical(data$outcome_type, "binary")) cdf$Y_bin <- data$Y_bin
    cov_names <- NULL
    if (!is.null(data$covariates) && ncol(data$covariates) > 0) {
        cv <- data$covariates
        nm <- names(cv)
        nm[nm %in% colnames(cdf)] <- paste0("cov_", nm[nm %in% colnames(cdf)])
        for (j in seq_along(nm)) cdf[[nm[j]]] <- cv[[j]]
        cov_names <- nm
    }

    roles <- list(
        exposure = "X",
        outcome = if (!is.null(ex$Y)) "Y",
        mediator = if (!is.null(ex$M)) "M",
        instrument = if (!is.null(data$G)) "G",
        mediator_instrument = if (!is.null(ex$Gm)) "Gm",
        negative_controls = if (exp_W) "W",
        W1 = if (exp_W1) "W1",
        W2 = if (exp_W2) "W2",
        covariates = cov_names,
        surv_time = if (identical(data$outcome_type, "survival")) "surv_time",
        surv_event = if (identical(data$outcome_type, "survival")) "surv_event",
        bin_outcome = if (identical(data$outcome_type, "binary")) "Y_bin")

    meta <- list(iconic = list(
        outcome_type = data$outcome_type,
        scaling = data$scaling,
        feature_names = data$feature_names,
        mediator_names = data$mediator_names,
        flags = data[c("has_instrument", "has_mediator_instrument", "has_nc",
            "has_path_nc", "recycled_lone_panel", "is_mediation")],
        roles = roles,
        trained_gan = data$trained_gan,
        iconic_version = as.character(utils::packageVersion("iconic"))))

    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(ex),
        colData = cdf, metadata = meta)
}
