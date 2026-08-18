# ============================================================
# Figure plotting functions (ggplot2 + patchwork)
#
# These functions produce publication-quality multi-panel figures
# from the simulation results returned by iconic's computational
# functions. They complement the base-R diagnostic plots in
# plots.R.
#
# Dependencies: ggplot2, patchwork, scales (Suggests)
# ============================================================


# --- Dependency check helper -----------------------------------------------

.figures_check_deps <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for this function. ",
         "Install it with install.packages('ggplot2').")
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for this function. ",
         "Install it with install.packages('patchwork').")
}

.figures_check_scales <- function() {
  if (!requireNamespace("scales", quietly = TRUE))
    stop("Package 'scales' is required for this function. ",
         "Install it with install.packages('scales').")
}

# --- Shared theme ----------------------------------------------------------

.figure_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_rect(colour = "grey35", fill = NA,
                                                linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "grey25", size = 8),
      axis.title = ggplot2::element_text(colour = "grey15", size = 9.5),
      axis.ticks = ggplot2::element_line(colour = "grey35", linewidth = 0.3),
      axis.ticks.length = ggplot2::unit(2, "pt"),
      plot.title = ggplot2::element_text(colour = "grey5", face = "bold",
                                                size = 11),
      plot.tag = ggplot2::element_text(colour = "grey10", face = "bold",
                                                size = 13),
      plot.tag.position = c(0.02, 0.97),
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(9, "pt"),
      legend.text = ggplot2::element_text(colour = "grey25", size = 7.5),
      legend.margin = ggplot2::margin(1, 1, 1, 1)
    )
}

# --- Shared colour scale for method colours --------------------------------

.method_colour_scale <- function(show_legend = FALSE) {
  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")
  ggplot2::scale_colour_manual(
    values = iconic_method_colors,
    limits = methods,
    name = NULL,
    guide = if (show_legend)
      ggplot2::guide_legend(override.aes = list(shape = 19, size = 2.4, alpha = 1))
    else "none"
  )
}

# --- Save helper -----------------------------------------------------------

.save_figure <- function(fig, file, width, height, dpi = 300) {
  if (!is.null(file)) {
    ggplot2::ggsave(file, fig, width = width, height = height,
                    units = "in", bg = "white", dpi = dpi)
    message("Saved ", file)
  }
}

# ============================================================
# Estimator benchmark (4-panel simulation figure)
# ============================================================

# Internal: Panel A — total-effect bias vs confounding strength
.plot_benchmark_panel_a <- function(iter_df, conf_grid) {
  te_methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC")
  df <- iter_df[iter_df$method %in% te_methods & iter_df$pval %in% conf_grid, ]
  df$method <- factor(df$method, levels = te_methods)
  df$xg <- factor(df$pval, levels = conf_grid)

  fill_vals <- setNames(adjustcolor(iconic_method_colors[te_methods], alpha.f = 0.25),
                        te_methods)
  ylim <- range(df$bias, na.rm = TRUE) + c(-0.05, 0.08)

  ggplot2::ggplot(df, ggplot2::aes(x = xg, y = bias, colour = method, fill = method)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, colour = "grey40") +
    ggplot2::geom_boxplot(width = 0.6, position = ggplot2::position_dodge(width = 0.7),
                          outlier.shape = NA, linewidth = 0.35, fatten = 1.6) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.15,
                          dodge.width = 0.7),
                        size = 0.5, alpha = 0.5, show.legend = FALSE) +
    .method_colour_scale() +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(x = "Confounding strength (delta)", y = "Bias", tag = "A") +
    .figure_theme()
}

# Internal: Mediation facet — NDE or NIE bias
.plot_benchmark_mediation_facet <- function(iter_df, xvals, bias_col, xlab,
                                             estimand_lab, panel_lab = NULL) {
  med_methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC", "IV2SLS2",
                   "PGC2", "PGC2Gm")
  methods <- intersect(med_methods, unique(iter_df$method))
  df <- iter_df[iter_df$method %in% methods & iter_df$pval %in% xvals, ]
  df$method <- factor(df$method, levels = methods)
  df$xg <- factor(df$pval, levels = xvals)
  df <- df[!is.na(df$xg), ]
  df$xg <- droplevels(df$xg)

  fill_vals <- setNames(adjustcolor(iconic_method_colors[methods], alpha.f = 0.25),
                        methods)

  keep <- ave(df[[bias_col]], df$method, df$xg, FUN = function(v) {
    q <- quantile(v, c(0.25, 0.75), na.rm = TRUE)
    iqr <- diff(q)
    v >= q[1] - 1.5 * iqr & v <= q[2] + 1.5 * iqr
  })
  df_pts <- df[as.logical(keep), ]

  ylim <- range(df_pts[[bias_col]], na.rm = TRUE)
  if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)
  ylim <- ylim + c(-0.06, 0.06) * diff(ylim)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = xg, y = .data[[bias_col]],
                          colour = method, fill = method)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, colour = "grey40") +
    ggplot2::geom_boxplot(width = 0.62, position = ggplot2::position_dodge(width = 0.72),
                          outlier.shape = NA, linewidth = 0.3, fatten = 1.4) +
    ggplot2::geom_point(data = df_pts,
                        position = ggplot2::position_jitterdodge(jitter.width = 0.15,
                          dodge.width = 0.72),
                        size = 0.45, alpha = 0.5) +
    .method_colour_scale() +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::scale_x_discrete(na.translate = FALSE) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = estimand_lab,
                      hjust = 1.15, vjust = 1.4, size = 2.6, fontface = 2,
                      colour = "grey25") +
    ggplot2::coord_cartesian(ylim = ylim, clip = "off") +
    ggplot2::labs(x = if (nzchar(xlab)) xlab else NULL,
                  y = paste0(estimand_lab, " bias")) +
    .figure_theme()
  if (!is.null(panel_lab)) p <- p + ggplot2::labs(tag = panel_lab)
  p
}

# Internal: Panel D — NIE Type I error vs confounding strength
.plot_benchmark_panel_d <- function(t1e_df, conf_grid) {
  med_methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC", "IV2SLS2",
                   "PGC2", "PGC2Gm")
  methods <- intersect(med_methods, unique(t1e_df$method))
  df <- t1e_df[t1e_df$method %in% methods, ]
  df$method <- factor(df$method, levels = methods)
  df <- df[order(df$method, df$conf_str), ]

  ggplot2::ggplot(df, ggplot2::aes(x = conf_str, y = NIE_type1, colour = method)) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = 2, linewidth = 0.7,
                        colour = "#c0392b") +
    ggplot2::geom_line(linewidth = 0.8, key_glyph = "point") +
    ggplot2::geom_point(size = 1.6, show.legend = FALSE) +
    .method_colour_scale(show_legend = TRUE) +
    ggplot2::coord_cartesian(xlim = range(conf_grid) + c(-0.03, 0.03)) +
    ggplot2::labs(x = "Confounding strength (delta)",
                  y = "NIE Type I error (prop. p < 0.05)",
                  title = "Red dashed line = alpha = 0.05", tag = "D") +
    .figure_theme() +
    ggplot2::theme(plot.title = ggplot2::element_text(colour = "#c0392b",
                     face = "plain", size = 8, hjust = 0.5))
}

#' Estimator benchmark figure
#'
#' Produces a 4-panel publication figure benchmarking the eight ICONIC
#' estimators under unmeasured confounding and mediator-outcome
#' confounding. Panels: (A) total-effect bias vs confounding strength,
#' (B) NDE/NIE bias vs confounding strength, (C) NDE/NIE bias vs sample
#' size, (D) NIE Type I error vs confounding strength.
#'
#' @param panel_a Result of \code{sweep_param("conf_str", ...)} for the
#' total-effect panel.
#' @param panel_b Result of \code{sweep_mediation_param("conf_str", ...)}
#' for the mediation confounding-strength panel.
#' @param panel_c Result of \code{sweep_mediation_param("n_samples", ...)}
#' for the mediation sample-size panel.
#' @param panel_d Result of \code{sweep_mediation_null_by_conf(...)} for
#' the Type I error panel.
#' @param conf_grid Numeric vector of confounding-strength values.
#' @param n_grid Numeric vector of sample sizes.
#' @param file Optional file path to save the figure (PDF or PNG).
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return A \code{patchwork} ggplot object.
#' @export
#' @examples
#' # Toy inputs standing in for sweep_param()/sweep_mediation_param() output
#' methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
#'   "IV2SLS2", "PGC2", "PGC2Gm")
#' it <- expand.grid(method = methods, pval = c(0.5, 0.8), rep = 1:3)
#' it$bias <- rnorm(nrow(it), 0, 0.05)
#' it$NDE_bias <- rnorm(nrow(it), 0, 0.05)
#' it$NIE_bias <- rnorm(nrow(it), 0, 0.05)
#' panel <- list(iter_bias = it)
#' panel_d <- expand.grid(method = methods, conf_str = c(0.5, 0.8))
#' panel_d$NIE_type1 <- runif(nrow(panel_d), 0, 0.1)
#' plot_estimator_benchmark(panel, panel, panel, panel_d,
#'   conf_grid = c(0.5, 0.8), n_grid = c(0.5, 0.8))
plot_estimator_benchmark <- function(panel_a, panel_b, panel_c, panel_d,
                                      conf_grid = c(0.2, 0.4, 0.6, 0.8, 1.0),
                                      n_grid = c(100, 200, 500, 1000),
                                      file = NULL, width = 8, height = 6) {
  .figures_check_deps()

  pA <- .plot_benchmark_panel_a(panel_a$iter_bias, conf_grid)
  pBnde <- .plot_benchmark_mediation_facet(panel_b$iter_bias, conf_grid,
                                            "NDE_bias", "", "NDE", "B")
  pBnie <- .plot_benchmark_mediation_facet(panel_b$iter_bias, conf_grid,
                                            "NIE_bias", "Confounding strength (delta)", "NIE")
  pCnde <- .plot_benchmark_mediation_facet(panel_c$iter_bias, n_grid,
                                            "NDE_bias", "", "NDE", "C")
  pCnie <- .plot_benchmark_mediation_facet(panel_c$iter_bias, n_grid,
                                            "NIE_bias", "Sample size (n)", "NIE")
  pD <- .plot_benchmark_panel_d(panel_d, conf_grid)

  design <- "
    AB
    AC
    DF
    EF
  "

  fig <- patchwork::wrap_plots(A = pA, B = pBnde, C = pBnie, D = pCnde,
                                E = pCnie, F = pD, design = design) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "ICONIC simulation: estimator bias and Type I error under unmeasured confounding",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(
        face = "bold", size = 11, colour = "grey10", hjust = 0.5))) &
    ggplot2::theme(legend.position = "bottom")

  .save_figure(fig, file, width, height)
  fig
}


# ============================================================
# Feature correlation sweep figure
# ============================================================

#' Feature correlation sweep figure
#'
#' Produces a 3-panel figure showing how estimator performance changes as
#' within-module feature correlation increases, modelling co-expression
#' modules in a correlated omics panel.
#'
#' Panel A: NDE and NIE bias vs feature correlation strength, faceted by
#' estimand. Panel B: NDE and NIE RMSE vs feature correlation strength,
#' faceted by estimand. Panel C: NIE Type I error vs feature correlation
#' strength.
#'
#' @param sweep_results Results from
#' \code{sweep_mediation_param("feat_cor", ...)}: a list with
#' \code{$summary} (data frame with columns \code{param_value},
#' \code{method}, \code{NDE_bias}, \code{NIE_bias}, \code{NDE_rmse},
#' \code{NIE_rmse}).
#' @param null_results Optional results from
#' \code{sweep_mediation_null_by_conf(..., feat_cor = ...)} run with
#' \code{feat_cor} swept: a data frame with columns \code{feat_cor},
#' \code{method}, \code{NIE_type1}. If \code{NULL}, Panel C is omitted
#' and the figure has 2 panels.
#' @param file Optional file path to save the figure.
#' @param width Figure width in inches. Default 10.
#' @param height Figure height in inches. Default 12.
#' @return A \code{patchwork} ggplot object.
#' @export
#'
#' @examples
#' # Toy inputs standing in for sweep_mediation_param("feat_cor", ...) output
#' sweep <- list(summary = expand.grid(param_value = c(0, 0.5),
#'   method = c("UNADJ", "IV2SLS2", "PGC2Gm")))
#' sweep$summary$NDE_bias <- rnorm(nrow(sweep$summary), 0, 0.05)
#' sweep$summary$NIE_bias <- rnorm(nrow(sweep$summary), 0, 0.05)
#' sweep$summary$NDE_rmse <- runif(nrow(sweep$summary), 0.05, 0.15)
#' sweep$summary$NIE_rmse <- runif(nrow(sweep$summary), 0.05, 0.15)
#' null <- expand.grid(feat_cor = c(0, 0.5),
#'   method = c("UNADJ", "IV2SLS2", "PGC2Gm"))
#' null$NIE_type1 <- runif(nrow(null), 0, 0.1)
#' plot_feature_correlation_sweep(sweep, null)
plot_feature_correlation_sweep <- function(sweep_results,
                                           null_results = NULL,
                                           file = NULL,
                                           width = 10,
                                           height = 12) {
  .figures_check_deps()

  s <- sweep_results$summary
  if (is.null(s)) stop("sweep_results must be a list with $summary (from sweep_mediation_param).")

  methods <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
               "IV2SLS2", "PGC2", "PGC2Gm")
  methods <- intersect(methods, unique(s$method))

  # --- Panel A: NDE/NIE bias vs feat_cor ---
  bias_df <- rbind(
    data.frame(feat_cor = s$param_value, method = s$method,
               estimand = "NDE", value = s$NDE_bias),
    data.frame(feat_cor = s$param_value, method = s$method,
               estimand = "NIE", value = s$NIE_bias)
  )
  bias_df <- bias_df[bias_df$method %in% methods, ]
  bias_df$method <- factor(bias_df$method, levels = methods)
  bias_df$estimand <- factor(bias_df$estimand, levels = c("NDE", "NIE"))

  pA <- ggplot2::ggplot(bias_df,
                        ggplot2::aes(x = feat_cor, y = value,
                                     color = method, group = method)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey50", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ estimand, ncol = 2, scales = "free_y") +
    .method_colour_scale(show_legend = TRUE) +
    ggplot2::labs(x = "Within-module feature correlation (rho)",
                  y = "Bias", tag = "A",
                  title = "NDE and NIE bias vs feature correlation") +
    .figure_theme()

  # --- Panel B: NDE/NIE RMSE vs feat_cor ---
  rmse_df <- rbind(
    data.frame(feat_cor = s$param_value, method = s$method,
               estimand = "NDE", value = s$NDE_rmse),
    data.frame(feat_cor = s$param_value, method = s$method,
               estimand = "NIE", value = s$NIE_rmse)
  )
  rmse_df <- rmse_df[rmse_df$method %in% methods, ]
  rmse_df$method <- factor(rmse_df$method, levels = methods)
  rmse_df$estimand <- factor(rmse_df$estimand, levels = c("NDE", "NIE"))

  pB <- ggplot2::ggplot(rmse_df,
                        ggplot2::aes(x = feat_cor, y = value,
                                     color = method, group = method)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ estimand, ncol = 2, scales = "free_y") +
    .method_colour_scale() +
    ggplot2::labs(x = "Within-module feature correlation (rho)",
                  y = "RMSE", tag = "B",
                  title = "NDE and NIE RMSE vs feature correlation") +
    .figure_theme()

  # --- Assemble ---
  if (!is.null(null_results)) {
    # Panel C: NIE Type I error vs feat_cor
    # null_results is expected to have a feat_cor column when swept
    t1e <- null_results
    t1e <- t1e[t1e$method %in% methods, ]
    t1e$method <- factor(t1e$method, levels = methods)

    # Determine which column holds the feat_cor values
    if ("feat_cor" %in% names(t1e)) {
      t1e$x <- t1e$feat_cor
    } else {
      # If null_results was run at a single feat_cor, use the sweep grid
      # from sweep_results instead (requires matching by method)
      t1e$x <- NA_real_
    }

    pC <- ggplot2::ggplot(t1e,
                          ggplot2::aes(x = x, y = NIE_type1,
                                       color = method, group = method)) +
      ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed",
                          color = "#c0392b", linewidth = 0.6) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 2) +
      .method_colour_scale() +
      ggplot2::labs(x = "Within-module feature correlation (rho)",
                    y = "NIE Type I error (prop. p < 0.05)", tag = "C",
                    title = "NIE Type I error vs feature correlation") +
      .figure_theme() +
      ggplot2::theme(legend.position = "bottom")

    fig <- (pA / pB / pC) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = "Impact of feature-level correlations on estimator performance",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 12,
                                             colour = "grey10", hjust = 0.5)))
  } else {
    fig <- (pA / pB) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = "Impact of feature-level correlations on estimator performance",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 12,
                                             colour = "grey10", hjust = 0.5)))
  }

  .save_figure(fig, file, width, height)
  fig
}

# ============================================================
# Degradation surface (3-panel heatmap figure)
# ============================================================

#' Degradation surface figure
#'
#' Produces a 3-panel figure showing estimator bias under simultaneous
#' instrument exogeneity violations. Panels: (A) IV2SLS2 NDE bias
#' surface, (B) PGC2Gm NDE bias surface, (C) crossover map showing
#' which estimator has lower |NDE bias| at each grid cell.
#'
#' @param results Data frame with columns \code{rho_G1}, \code{rho_G2},
#' \code{method}, \code{NDE_bias}, \code{NIE_bias}, and derived
#' \code{NDE_abs}, \code{NIE_abs}. Typically built by sweeping
#' \code{sweep_mediation_param()} across a rho_G1 x rho_G2 grid.
#' @param rho_G1_grid Numeric vector of exposure-instrument violation
#' values.
#' @param rho_G2_grid Numeric vector of mediator-instrument violation
#' values.
#' @param omega_facet Logical: when the sensitivity surface swept
#' \code{omega_1} / \code{omega_2} (more than one distinct value), facet the
#' rho_G1 x rho_G2 heatmaps by the omega cell. Default \code{TRUE}. Ignored
#' when omega was not swept.
#' @param file Optional file path to save the figure.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return A \code{patchwork} ggplot object.
#' @export
#' @examples
#' results <- expand.grid(rho_G1 = c(0, 0.2), rho_G2 = c(0, 0.2),
#'   method = c("IV2SLS2", "PGC2Gm"))
#' results$NDE_bias <- rnorm(nrow(results), 0, 0.05)
#' results$NIE_bias <- rnorm(nrow(results), 0, 0.05)
#' plot_degradation_surface(results, rho_G1_grid = c(0, 0.2),
#'   rho_G2_grid = c(0, 0.2))
plot_degradation_surface <- function(results,
                                      rho_G1_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                                      rho_G2_grid = c(0, 0.1, 0.2, 0.3, 0.5),
                                      omega_facet = TRUE,
                                      file = NULL, width = 15, height = 5.5) {
  .figures_check_deps()

  # When omega is swept, facet on the omega cell. Use a single composite
  # omega label (omega_1 x omega_2) as the facet dimension.
  has_omega <- "omega_1" %in% names(results) && "omega_2" %in% names(results)
  omega_swept <- has_omega &&
    (length(unique(results$omega_1)) > 1 || length(unique(results$omega_2)) > 1)
  do_facet <- omega_facet && omega_swept

  if (do_facet) {
    results$omega_cell <- sprintf("omega[1]==%.2g ~ omega[2]==%.2g",
                                  results$omega_1, results$omega_2)
    # order facets by omega_1 then omega_2
    oc <- unique(results[, c("omega_1", "omega_2", "omega_cell")])
    oc <- oc[order(oc$omega_1, oc$omega_2), ]
    results$omega_cell <- factor(results$omega_cell, levels = oc$omega_cell)
  }

  iv_df <- results[results$method == "IV2SLS2", ]
  pg_df <- results[results$method == "PGC2Gm", ]
  iv_df$rho_G1 <- factor(iv_df$rho_G1, levels = rho_G1_grid)
  iv_df$rho_G2 <- factor(iv_df$rho_G2, levels = rho_G2_grid)
  pg_df$rho_G1 <- factor(pg_df$rho_G1, levels = rho_G1_grid)
  pg_df$rho_G2 <- factor(pg_df$rho_G2, levels = rho_G2_grid)

  # Per-estimator robust colour scales. The two estimators can have very
  # different bias magnitudes (IV2SLS2 degrades under instrument violation
  # while PGC2Gm stays near zero), so a single shared scale washes out the
  # better-behaved estimator. Give each panel its own diverging scale capped
  # at that estimator's 95% |bias| quantile (floored at 0.10 so a near-perfect
  # estimator still shows its structure), squishing out-of-bounds cells to the
  # extreme colour. The crossover map (panel C) compares |bias| directly, so
  # the cross-estimator comparison is preserved there.
  squish <- if (requireNamespace("scales", quietly = TRUE)) scales::squish else NULL
  cap_for <- function(x) {
    v <- abs(x[is.finite(x)])
    if (length(v) == 0) return(0.25)
    max(0.10, as.numeric(stats::quantile(v, 0.95, na.rm = TRUE)))
  }
  bias_lim_A <- cap_for(iv_df$NDE_bias)
  bias_lim_B <- cap_for(pg_df$NDE_bias)

  facet_layer <- if (do_facet)
    ggplot2::facet_wrap(~ omega_cell, labeller = ggplot2::label_parsed, nrow = 1) else NULL

  pA <- ggplot2::ggplot(iv_df, ggplot2::aes(x = rho_G2, y = rho_G1, fill = NDE_bias)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", NDE_bias)),
                       size = 2.6, fontface = "bold", colour = "grey20") +
    ggplot2::scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#FF9400",
                                  midpoint = 0, limits = c(-bias_lim_A, bias_lim_A),
                                  oob = squish,
                                  name = "NDE\nbias") +
    facet_layer +
    ggplot2::labs(x = expression(rho[G2]~" (mediator-instrument violation)"),
                  y = expression(rho[G1]~" (exposure-instrument violation)"),
                  title = "IV2SLS2 NDE bias", tag = "A") +
    .figure_theme() +
    ggplot2::theme(legend.position = "right")

  pB <- ggplot2::ggplot(pg_df, ggplot2::aes(x = rho_G2, y = rho_G1, fill = NDE_bias)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", NDE_bias)),
                       size = 2.6, fontface = "bold", colour = "grey20") +
    ggplot2::scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#FF9400",
                                  midpoint = 0, limits = c(-bias_lim_B, bias_lim_B),
                                  oob = squish,
                                  name = "NDE\nbias") +
    facet_layer +
    ggplot2::labs(x = expression(rho[G2]~" (mediator-instrument violation)"),
                  y = expression(rho[G1]~" (exposure-instrument violation)"),
                  title = "PGC2Gm NDE bias", tag = "B") +
    .figure_theme() +
    ggplot2::theme(legend.position = "right")

  # Crossover map
  crossover_df <- data.frame(
    rho_G1 = iv_df$rho_G1, rho_G2 = iv_df$rho_G2,
    iv_abs = abs(iv_df$NDE_bias), pg_abs = abs(pg_df$NDE_bias),
    iv_nde = iv_df$NDE_bias, pg_nde = pg_df$NDE_bias
  )
  if (do_facet) crossover_df$omega_cell <- iv_df$omega_cell
  crossover_df$winner <- ifelse(crossover_df$iv_abs <= crossover_df$pg_abs,
                                "IV2SLS2", "PGC2Gm")
  crossover_df$iv_exceeds <- crossover_df$iv_abs > 0.10
  crossover_df$pg_exceeds <- crossover_df$pg_abs > 0.10
  crossover_df$label <- ifelse(crossover_df$iv_exceeds | crossover_df$pg_exceeds, "!", "")

  winner_cols <- c("IV2SLS2" = "#0279EE", "PGC2Gm" = "#FD9BED")

  pC <- ggplot2::ggplot(crossover_df, ggplot2::aes(x = rho_G2, y = rho_G1, fill = winner)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       size = 3.4, fontface = "bold", colour = "black") +
    ggplot2::scale_fill_manual(values = winner_cols, name = "Preferred\nestimator",
                               limits = c("IV2SLS2", "PGC2Gm")) +
    facet_layer +
    ggplot2::labs(x = expression(rho[G2]~" (mediator-instrument violation)"),
                  y = expression(rho[G1]~" (exposure-instrument violation)"),
                  title = "Crossover: preferred estimator by |NDE bias|\n(! = preferred estimator exceeds |bias| = 0.10)",
                  tag = "C") +
    .figure_theme() +
    ggplot2::theme(legend.position = "right")

  fig <- (pA | pB | pC) +
    patchwork::plot_annotation(
      title = "Degradation surface: estimator bias under simultaneous instrument exogeneity violations",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(
        face = "bold", size = 11, colour = "grey10", hjust = 0.5)))

  # widen when faceting on omega
  if (do_facet) {
    n_facets <- length(unique(results$omega_cell))
    width <- max(width, 5 * n_facets)
    height <- max(height, 6.5)
  }
  .save_figure(fig, file, width, height)
  fig
}

# ============================================================
# NC validity diagnostics (5-panel sweep figure)
# ============================================================

#' Sweep negative-control validity diagnostics
#'
#' Runs simulation sweeps for the empirical NC validity diagnostics
#' and returns data frames ready for plotting with
#' \code{plot_nc_validity_diagnostics()}. Panels: (A) A1 W perp X | C,
#' (B) A2 W perp G | C, (C) A2' W perp Gm | C, (D) A3
#' covariance-capture versus true coverage omega, (E) A3 support
#' R2(Utilde | W) versus omega.
#'
#' @param n_samples Number of synthetic samples per iteration.
#' @param n_iter Number of replications per sweep point.
#' @param phi_val Mediator-instrument strength (for the A2' panel and
#' the capture/support sweeps).
#' @param contam_grid X->W contamination strength grid (A1).
#' @param meqtl_grid G->W (meQTL) strength grid (A2).
#' @param eqtl_grid Gm->W (eQTL) strength grid (A2').
#' @param omega_grid True negative-control coverage grid for the
#' capture (D) and support (E) sweeps. Includes 0 by default so the
#' permutation-null calibration is visible.
#' @param n_perm Permutations per replicate for the capture test.
#' Default 200, matching \code{iconic_diagnose()}.
#' @param cs_confounders Number of latent confounders for the
#' capture/support sweeps (distinct loadings). Default 2.
#' @param include_a3_grid Logical; also run the legacy A3 dimensional
#' grid (PGC bias over n_valid x k). Default FALSE (no longer plotted).
#' @param k_grid,n_valid_grid Grids for the legacy A3 dimensional sweep.
#' @param n_cores Number of parallel workers for the replicate loops.
#' Default 1 (sequential). Uses \code{parallel::mclapply} on Unix
#' and a PSOCK cluster on Windows.
#' @return A list with elements \code{panel_a}, \code{panel_b},
#' \code{panel_d}, \code{panel_capture_support} (per-replicate raw),
#' and optionally \code{panel_c} when \code{include_a3_grid = TRUE}.
#' @export
#' @examples
#' \donttest{
#' panels <- sweep_nc_validity(n_samples = 100, n_iter = 2,
#'   contam_grid = c(0, 0.1), meqtl_grid = c(0, 0.1),
#'   eqtl_grid = c(0, 0.1), omega_grid = c(0, 0.5), n_perm = 20)
#' panels$panel_a
#' panels$panel_capture_support
#' }
sweep_nc_validity <- function(n_samples = 500, n_iter = 50, phi_val = 0.8,
                              contam_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
                              meqtl_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
                              eqtl_grid = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5),
                              omega_grid = c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0),
                              n_perm = 200, cs_confounders = 2,
                              include_a3_grid = FALSE,
                              k_grid = c(1, 2, 3), n_valid_grid = c(1, 2, 3),
                              n_cores = 1) {

  message("Running NC validity sweeps (n_iter = ", n_iter, ", n = ", n_samples, ") ...")
  t0 <- Sys.time()

  ## Panel A: A1 (W perp X | C) — X->W contamination  [unchanged]
  message(" Panel A: nc_validity_screen (A1)")
  panel_a <- do.call(rbind, lapply(contam_grid, function(cs) {
    res_a <- .parallel_lapply(seq_len(n_iter), function(i) {
      dat <- run_single_iteration(NULL, n_synthetic_samples = n_samples,
                                  n_features = 10, n_confounders = 1,
                                  coverage = 0.7, seed = i)
      dat$W[, 6:10] <- matrix(rnorm(n_samples * 5), n_samples, 5)
      dat$W[, 6:10] <- dat$W[, 6:10] + cs * dat$X
      s <- nc_validity_screen(dat)
      c(sum(s$significant[seq.int(6, 10)]) / 5, sum(s$significant[seq_len(5)]) / 5)
    }, n_cores = n_cores, progress = " Panel A replicates")
    det_violated <- vapply(res_a, function(x) x[1], numeric(1))
    det_confounding <- vapply(res_a, function(x) x[2], numeric(1))
    data.frame(contamination = cs,
               violated_mean = mean(det_violated), violated_sd = sd(det_violated),
               confounding_mean = mean(det_confounding), confounding_sd = sd(det_confounding),
               stringsAsFactors = FALSE)
  }))

  ## Panel B: A2 (W perp G | C) — G->W (meQTL)  [unchanged]
  message(" Panel B: nc_independence_check (A2)")
  panel_b <- do.call(rbind, lapply(meqtl_grid, function(ms) {
    res_b <- .parallel_lapply(seq_len(n_iter), function(i) {
      dat <- run_single_iteration(NULL, n_synthetic_samples = n_samples,
                                  n_features = 10, n_confounders = 1, seed = i)
      dat$W[, seq_len(5)] <- dat$W[, seq_len(5)] + ms * dat$G[, 1]
      s <- nc_independence_check(dat)
      c(sum(s$significant[seq_len(5)]) / 5, sum(s$significant[seq.int(6, 10)]) / 5)
    }, n_cores = n_cores, progress = " Panel B replicates")
    det_violated <- vapply(res_b, function(x) x[1], numeric(1))
    det_clean <- vapply(res_b, function(x) x[2], numeric(1))
    data.frame(meqtl = ms,
               violated_mean = mean(det_violated), violated_sd = sd(det_violated),
               clean_mean = mean(det_clean), clean_sd = sd(det_clean),
               stringsAsFactors = FALSE)
  }))

  ## Panel C: A2' (W perp Gm | C) — Gm->W (eQTL)  [unchanged]
  message(" Panel C: nc_independence_check_gm (A2')")
  panel_d <- do.call(rbind, lapply(eqtl_grid, function(es) {
    res_d <- .parallel_lapply(seq_len(n_iter), function(i) {
      dat <- run_single_iteration(NULL, n_synthetic_samples = n_samples,
                                  n_features = 10, n_confounders = 1,
                                  phi = phi_val, seed = i)
      dat$W[, seq_len(5)] <- dat$W[, seq_len(5)] + es * dat$Gm
      s <- nc_independence_check_gm(dat)
      c(sum(s$significant[seq_len(5)]) / 5, sum(s$significant[seq.int(6, 10)]) / 5)
    }, n_cores = n_cores, progress = " Panel C replicates")
    det_violated <- vapply(res_d, function(x) x[1], numeric(1))
    det_clean <- vapply(res_d, function(x) x[2], numeric(1))
    data.frame(eqtl = es,
               violated_mean = mean(det_violated), violated_sd = sd(det_violated),
               clean_mean = mean(det_clean), clean_sd = sd(det_clean),
               stringsAsFactors = FALSE)
  }))

  ## Panels D-E (new): capture + support versus true coverage omega.
  ## Each replicate draws one dataset (k = cs_confounders, distinct NC
  ## loadings) and runs nc_completeness_capture() and nc_support_check()
  ## on it, so both diagnostics share the same realized data.
  message(" Panels D-E: nc_completeness_capture + nc_support_check vs omega")
  panel_capture_support <- do.call(rbind, lapply(seq_along(omega_grid), function(oi) {
    om <- omega_grid[oi]
    res <- .parallel_lapply(seq_len(n_iter), function(i) {
      dat <- run_single_iteration(NULL, n_synthetic_samples = n_samples,
                                  n_features = 10,
                                  n_confounders = cs_confounders,
                                  coverage = om, phi = phi_val,
                                  seed = 3000L + oi * 1000L + i,
                                  nc_params = list(mode = "distinct"))
      cap <- tryCatch(
        nc_completeness_capture(dat, outcome = "Y", n_perm = n_perm, n_cores = 1),
        error = function(e) NULL)
      sup <- tryCatch(nc_support_check(dat), error = function(e) NULL)
      data.frame(
        omega = om, rep = i,
        capture_R2 = if (!is.null(cap)) cap$capture_R2 else NA_real_,
        capture_p = if (!is.null(cap)) cap$capture_pvalue else NA_real_,
        null_R2 = if (!is.null(cap)) mean(cap$null_distribution, na.rm = TRUE) else NA_real_,
        support_R2 = if (!is.null(sup)) sup$R2_utilde_given_W else NA_real_,
        frac_adds_coverage = if (!is.null(sup)) mean(sup$support$adds_coverage, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE)
    }, n_cores = n_cores, progress = paste0("  omega = ", om))
    do.call(rbind, res)
  }))

  out <- list(panel_a = panel_a, panel_b = panel_b, panel_d = panel_d,
              panel_capture_support = panel_capture_support)

  ## Legacy A3 dimensional grid (PGC bias over n_valid x k) — optional.
  if (isTRUE(include_a3_grid)) {
    message(" Legacy A3 grid: nc_completeness_check + fit_pgc")
    grid <- expand.grid(n_valid = n_valid_grid, k = k_grid, KEEP.OUT.ATTRS = FALSE)
    panel_c <- do.call(rbind, lapply(seq_len(nrow(grid)), function(gi) {
      nv <- grid$n_valid[gi]; kk <- grid$k[gi]
      biases <- unlist(.parallel_lapply(seq_len(n_iter), function(i) {
        dat <- run_single_iteration(NULL, n_synthetic_samples = n_samples,
                                    n_features = 10, n_confounders = kk, seed = i,
                                    nc_params = list(mode = "distinct"))
        W_valid <- dat$W[, seq_len(nv), drop = FALSE]
        res <- fit_pgc(dat$Y[, 1], dat$X, dat$G[, 1], W_valid)
        res$beta - dat$true_total
      }, n_cores = n_cores, progress = " A3 grid replicates"))
      dat_comp <- run_single_iteration(NULL, n_synthetic_samples = n_samples,
                                       n_features = 10, n_confounders = kk, seed = 1,
                                       nc_params = list(mode = "distinct"))
      comp <- nc_completeness_check(dat_comp, n_valid_controls = nv)
      data.frame(n_valid = nv, k = kk,
                 pgc_bias = mean(biases), pgc_bias_sd = sd(biases),
                 completeness = comp$completeness, stringsAsFactors = FALSE)
    }))
    out$panel_c <- panel_c
  }

  message(sprintf(" Done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out
}

#' NC validity diagnostics figure (5 panels)
#'
#' Produces a 5-panel figure sweeping the empirical NC validity
#' diagnostics. Panels: (A) A1 W perp X | C, (B) A2 W perp G | C,
#' (C) A2' W perp Gm | C, (D) A3 covariance-capture versus true
#' coverage omega (with permutation-null mean and the fraction of
#' replicates with permutation p < 0.05 on a secondary axis), and
#' (E) A3 support R2(Utilde | W) versus omega.
#'
#' @param panels List returned by \code{sweep_nc_validity()}.
#' @param file Optional file path to save the figure (cairo_pdf, so the
#' omega and Utilde glyphs render).
#' @param width Figure width in inches. Default 10.
#' @param height Figure height in inches. Default 6.5.
#' @return A \code{patchwork} ggplot object.
#' @export
#' @examples
#' # Toy panels standing in for sweep_nc_validity() output
#' panels <- list(
#'   panel_a = data.frame(contamination = c(0, 0.1),
#'     violated_mean = c(0.1, 0.3), violated_sd = 0.05,
#'     confounding_mean = c(0.1, 0.12), confounding_sd = 0.05),
#'   panel_b = data.frame(meqtl = c(0, 0.1),
#'     violated_mean = c(0.1, 0.3), violated_sd = 0.05,
#'     clean_mean = c(0.1, 0.12), clean_sd = 0.05),
#'   panel_d = data.frame(eqtl = c(0, 0.1),
#'     violated_mean = c(0.1, 0.3), violated_sd = 0.05,
#'     clean_mean = c(0.1, 0.12), clean_sd = 0.05),
#'   panel_capture_support = data.frame(
#'     omega = rep(c(0, 0.5), each = 2), rep = 1:2,
#'     capture_R2 = c(0.02, 0.03, 0.25, 0.30),
#'     capture_p = c(0.40, 0.60, 0.01, 0.02),
#'     null_R2 = 0.02,
#'     support_R2 = c(0.02, 0.03, 0.35, 0.40),
#'     frac_adds_coverage = 1))
#' plot_nc_validity_diagnostics(panels)
plot_nc_validity_diagnostics <- function(panels, file = NULL,
                                         width = 10, height = 6.5) {
  .figures_check_deps()
  .figures_check_scales()

  col_violated <- "#FF9400"; col_confounding <- "#888888"
  col_clean <- "#27A062"; col_capture <- "#0279EE"; col_support <- "#75A025"

  build_screen_panel <- function(df, xvar, xlabel, title_txt, groups) {
    long <- do.call(rbind, lapply(names(groups), function(g) {
      spec <- groups[[g]]
      data.frame(x = df[[xvar]], rate = df[[spec$mean]], sd = df[[spec$sd]],
                 group = g, stringsAsFactors = FALSE)
    }))
    cols <- vapply(names(groups), function(g) groups[[g]]$col, character(1))
    ggplot2::ggplot(long, ggplot2::aes(x = x, y = rate, color = group, fill = group)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(rate - sd, 0),
                                        ymax = pmin(rate + sd, 1)),
                           alpha = 0.15, color = NA) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2, shape = 16) +
      ggplot2::geom_hline(yintercept = 0.10, linetype = "dotted", color = "grey50",
                          linewidth = 0.4) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::scale_y_continuous(limits = c(0, 1.05),
                                  labels = scales::percent_format(accuracy = 1)) +
      ggplot2::labs(x = xlabel, y = "Proportion flagged \"drop\"",
                    title = title_txt, color = NULL, fill = NULL) +
      .figure_theme() +
      ggplot2::theme(legend.position = "bottom",
                     legend.text = ggplot2::element_text(size = 7.5)) +
      ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
                      fill = "none")
  }

  pA <- build_screen_panel(
    panels$panel_a, "contamination",
    "X\u2192W contamination strength", "A1: W perp X | C",
    groups = list(
      "Confounder-sharing (expected)" = list(mean = "confounding_mean", sd = "confounding_sd", col = col_confounding),
      "Violated (injected X\u2192W)" = list(mean = "violated_mean", sd = "violated_sd", col = col_violated)))
  pB <- build_screen_panel(
    panels$panel_b, "meqtl",
    "G\u2192W (meQTL) strength", "A2: W perp G | C",
    groups = list(
      "Violated (injected G\u2192W)" = list(mean = "violated_mean", sd = "violated_sd", col = col_violated),
      "Clean (false positive)" = list(mean = "clean_mean", sd = "clean_sd", col = col_clean)))
  pC <- build_screen_panel(
    panels$panel_d, "eqtl",
    "Gm\u2192W strength", "A2': W perp Gm | C",
    groups = list(
      "Violated (injected Gm\u2192W)" = list(mean = "violated_mean", sd = "violated_sd", col = col_violated),
      "Clean (false positive)" = list(mean = "clean_mean", sd = "clean_sd", col = col_clean)))

  ## Summarize the capture/support sweep
  raw <- panels$panel_capture_support
  summ <- do.call(rbind, lapply(split(raw, raw$omega), function(d) {
    data.frame(
      omega = d$omega[1],
      capture_R2_mean = mean(d$capture_R2, na.rm = TRUE),
      capture_R2_sd = stats::sd(d$capture_R2, na.rm = TRUE),
      null_R2_mean = mean(d$null_R2, na.rm = TRUE),
      frac_sig = mean(d$capture_p < 0.05, na.rm = TRUE),
      support_R2_mean = mean(d$support_R2, na.rm = TRUE),
      support_R2_sd = stats::sd(d$support_R2, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))

  ## Panel D: capture vs omega
  sec_scale <- 0.6
  pD <- ggplot2::ggplot(summ, ggplot2::aes(x = omega)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(capture_R2_mean - capture_R2_sd, 0),
                                      ymax = capture_R2_mean + capture_R2_sd),
                         alpha = 0.15, fill = col_capture, color = NA) +
    ggplot2::geom_line(ggplot2::aes(y = capture_R2_mean, color = "Observed capture R\u00b2"),
                       linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = capture_R2_mean, color = "Observed capture R\u00b2"),
                        size = 2, shape = 16) +
    ggplot2::geom_line(ggplot2::aes(y = null_R2_mean, color = "Permutation null (mean)"),
                       linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = null_R2_mean, color = "Permutation null (mean)"),
                        size = 1.6, shape = 17) +
    ggplot2::geom_line(ggplot2::aes(y = frac_sig * sec_scale,
                                    color = "Fraction permutation p<0.05"),
                       linewidth = 0.6, linetype = "twodash") +
    ggplot2::geom_hline(yintercept = 0.3, linetype = "dashed", colour = "grey40",
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.1, linetype = "dotted", colour = "grey60",
                        linewidth = 0.4) +
    ggplot2::scale_color_manual(values = c(
      "Observed capture R\u00b2" = col_capture,
      "Permutation null (mean)" = "grey50",
      "Fraction permutation p<0.05" = col_violated)) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      sec.axis = ggplot2::sec_axis(~ . / sec_scale, name = "Fraction p<0.05",
                                   labels = scales::percent_format(accuracy = 1))) +
    ggplot2::labs(x = "True NC coverage (\u03c9)",
                  y = "Incremental R\u00b2 of W for Y | C",
                  title = "A3 capture: R\u00b2(W | C) vs coverage", color = NULL) +
    .figure_theme() +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 7.5)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  ## Panel E: support vs omega
  pE <- ggplot2::ggplot(summ, ggplot2::aes(x = omega)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(support_R2_mean - support_R2_sd, 0),
                                      ymax = pmin(support_R2_mean + support_R2_sd, 1)),
                         alpha = 0.15, fill = col_support, color = NA) +
    ggplot2::geom_line(ggplot2::aes(y = support_R2_mean), color = col_support,
                       linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = support_R2_mean), color = col_support,
                        size = 2, shape = 16) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40",
                        linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.2, linetype = "dotted", colour = "grey60",
                        linewidth = 0.4) +
    ggplot2::scale_y_continuous(limits = c(0, 1.0)) +
    ggplot2::labs(x = "True NC coverage (\u03c9)",
                  y = expression(R^2 * (tilde(U) ~ "|" ~ W)),
                  title = "A3 support: R\u00b2(\u0168 | W) vs coverage") +
    .figure_theme()

  fig <- ((pA | pB | pC) / (pD | pE | patchwork::plot_spacer())) +
    patchwork::plot_annotation(
      title = "ICONIC negative-control validity diagnostics: simulation sweeps",
      tag_levels = "A",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(
        size = 11, face = "bold", hjust = 0.5))) &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))

  if (!is.null(file)) {
    ggplot2::ggsave(file, fig, width = width, height = height, units = "in",
                    bg = "white", dpi = 300, device = grDevices::cairo_pdf)
    message("Saved ", file)
  }
  fig
}

# ============================================================
# Model selection workflow (3-panel figure)
# ============================================================

#' Model selection workflow figure
#'
#' Produces a 3-panel figure showing the ICONIC model selection
#' workflow on example data. Panels: (A) eligibility table from
#' \code{iconic_diagnose()}, (B) NDE/NIE forest plot from
#' \code{iconic_estimate()}, (C) degradation surface heatmap from
#' \code{iconic_sensitivity()}.
#'
#' @param diagnosis Result of \code{iconic_diagnose()}.
#' @param estimate Result of \code{iconic_estimate()}.
#' @param sensitivity Result of \code{iconic_sensitivity()}.
#' @param recommendation Result of \code{iconic_recommend()}.
#' @param file Optional file path to save the figure.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return A \code{patchwork} ggplot object.
#' @export
#' @examples
#' if (check_torch_setup()) {
#'   data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#'     M = rnorm(100), G = rnorm(100), Gm = rnorm(100),
#'     W = matrix(rnorm(100 * 10), 10, 100))
#'   diag <- iconic_diagnose(data)
#'   est <- iconic_estimate(data, diagnosis = diag)
#'   sens <- iconic_sensitivity(data, n_iter = 2, gan_epochs = 5,
#'     rho_G1_grid = c(0, 0.2), rho_G2_grid = c(0, 0.2))
#'   rec <- iconic_recommend(data, diagnosis = diag, estimate = est,
#'     sensitivity = sens, auto_sensitivity = FALSE)
#'   plot_model_selection(diag, est, sens, rec)
#' }
plot_model_selection <- function(diagnosis, estimate, sensitivity, recommendation,
                                  file = NULL, width = 12, height = 14) {
  .figures_check_deps()

  method_order <- c("UNADJ", "DIRECT", "COCA", "IV2SLS", "PGC",
                    "IV2SLS2", "PGC2", "PGC2Gm")

  # Panel A: Eligibility table
  elig <- diagnosis$eligibility
  elig$label <- ifelse(elig$eligible, "Yes", "No")
  elig$estimator <- factor(elig$estimator, levels = method_order)

  pA <- ggplot2::ggplot(elig, ggplot2::aes(x = estimator, y = 1, fill = eligible)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label, colour = eligible),
                       size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#75A025", "FALSE" = "#E07B00"),
                               labels = c("Eligible", "Ineligible"), name = "") +
    ggplot2::scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "white"),
                                 guide = "none") +
    ggplot2::labs(x = "Estimator", y = NULL,
                  title = "Eligibility (iconic_diagnose)", tag = "A") +
    .figure_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) +
    ggplot2::coord_fixed(ratio = 3)

  # Panel B: NDE/NIE forest plot
  est_agg <- aggregate(cbind(NDE, NIE) ~ method, data = estimate,
                       FUN = mean, na.rm = TRUE)
  est_agg$NDE_se <- aggregate(NDE_se ~ method, data = estimate,
                              FUN = mean, na.rm = TRUE)$NDE_se
  est_agg$NIE_se <- aggregate(NIE_se ~ method, data = estimate,
                              FUN = mean, na.rm = TRUE)$NIE_se
  est_agg$recommended <- est_agg$method == recommendation$recommended
  est_agg$method <- factor(est_agg$method, levels = method_order)

  forest_data <- rbind(
    data.frame(method = est_agg$method, effect = "NDE",
               estimate = est_agg$NDE, se = est_agg$NDE_se,
               recommended = est_agg$recommended),
    data.frame(method = est_agg$method, effect = "NIE",
               estimate = est_agg$NIE, se = est_agg$NIE_se,
               recommended = est_agg$recommended))
  forest_data <- forest_data[!is.na(forest_data$estimate), ]
  forest_data$effect <- factor(forest_data$effect, levels = c("NDE", "NIE"))

  pB <- ggplot2::ggplot(forest_data,
                        ggplot2::aes(x = estimate, y = method, colour = method)) +
    ggplot2::geom_vline(xintercept = 0.10, linetype = "dashed",
                        colour = "grey50", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0.15, linetype = "dashed",
                        colour = "grey50", linewidth = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = estimate - 1.96 * se,
                          xmax = estimate + 1.96 * se),
                          height = 0.2, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(size = recommended, shape = recommended)) +
    ggplot2::scale_colour_manual(values = iconic_method_colors, guide = "none") +
    ggplot2::scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
    ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    ggplot2::facet_wrap(~ effect, ncol = 2, scales = "free_x") +
    ggplot2::labs(x = "Estimate (mean across features)", y = "Method",
                  title = "Estimates with 95% CI (iconic_estimate)", tag = "B") +
    .figure_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(
                     colour = "grey92", linewidth = 0.3))

  # Panel C: Degradation surface heatmap for the recommended estimator,
  # faceted by negative-control coverage (omega) when omega was swept. This
  # mirrors the Figure 3 layout: a rho_G1 x rho_G2 plane per omega cell.
  surf <- sensitivity$surface
  rec_method <- recommendation$recommended
  if (!rec_method %in% surf$method) rec_method <- "PGC2Gm"
  surf_sub <- surf[surf$method == rec_method, ]

  # Build the omega facet label when omega was swept (more than one distinct
  # omega value). Use a single composite omega_1 x omega_2 label, ordered by
  # omega_1 then omega_2, as in plot_degradation_surface().
  has_omega <- "omega_1" %in% names(surf_sub) && "omega_2" %in% names(surf_sub)
  omega_swept <- has_omega &&
    (length(unique(surf_sub$omega_1)) > 1 || length(unique(surf_sub$omega_2)) > 1)
  if (omega_swept) {
    surf_sub$omega_cell <- sprintf("omega[1]==%.2g ~ omega[2]==%.2g",
                                   surf_sub$omega_1, surf_sub$omega_2)
    oc <- unique(surf_sub[, c("omega_1", "omega_2", "omega_cell")])
    oc <- oc[order(oc$omega_1, oc$omega_2), ]
    surf_sub$omega_cell <- factor(surf_sub$omega_cell, levels = oc$omega_cell)
  }

  surf_sub$rho_G1 <- factor(surf_sub$rho_G1, levels = c(0, 0.1, 0.2, 0.3, 0.5))
  surf_sub$rho_G2 <- factor(surf_sub$rho_G2, levels = c(0, 0.1, 0.2, 0.3, 0.5))
  bias_lim <- max(abs(surf_sub$NDE_bias), na.rm = TRUE)

  facet_layer <- if (omega_swept)
    ggplot2::facet_wrap(~ omega_cell, labeller = ggplot2::label_parsed, nrow = 1) else NULL

  pC <- ggplot2::ggplot(surf_sub, ggplot2::aes(x = rho_G2, y = rho_G1, fill = NDE_bias)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f", NDE_bias)),
                       size = 2.4, fontface = "bold", colour = "grey20") +
    ggplot2::scale_fill_gradient2(low = "#0279EE", mid = "white", high = "#FF9400",
                                  midpoint = 0, limits = c(-bias_lim, bias_lim),
                                  name = "NDE\nbias") +
    facet_layer +
    ggplot2::labs(x = expression(rho[G2]~" (Gm violation)"),
                  y = expression(rho[G1]~" (G violation)"),
                  title = sprintf("Degradation surface (iconic_sensitivity): %s",
                                  rec_method),
                  tag = "C") +
    .figure_theme() +
    ggplot2::theme(legend.position = "right")

  fig <- (pA / pB / pC) +
    patchwork::plot_annotation(
      title = "ICONIC model selection workflow on example data",
      subtitle = paste0("Recommended: ", recommendation$recommended),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13,
                                           colour = "grey10", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, colour = "grey30",
                                              hjust = 0.5)))

  .save_figure(fig, file, width, height)
  fig
}

# ============================================================
# Prospective analysis (2-panel figure)
# ============================================================

#' Prospective analysis figure
#'
#' Produces a 2-panel figure showing expected gains from adding
#' genetic instruments. Panel A: NDE bias vs instrument strength.
#' Panel B: prospective NDE/NIE estimates at target instrument
#' strength.
#'
#' @param prospect Result of \code{iconic_prospect()}.
#' @param file Optional file path to save the figure.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return A \code{patchwork} ggplot object.
#' @export
#' @examples
#' if (check_torch_setup()) {
#'   data <- iconic_data(X = rnorm(100), Y = matrix(rnorm(100 * 10), 10, 100),
#'     M = rnorm(100))
#'   pros <- iconic_prospect(data, n_iter = 2, gan_epochs = 5,
#'     gamma_G_grid = c(0.4, 0.8), run_rho_sweep = FALSE)
#'   plot_prospective_analysis(pros)
#' }
plot_prospective_analysis <- function(prospect, file = NULL,
                                       width = 12, height = 8) {
  .figures_check_deps()

  method_cols <- iconic_method_colors
  varying_methods <- c("IV2SLS", "IV2SLS2", "PGC2Gm")

  # Panel A: NDE bias vs instrument strength
  surf <- prospect$strength_surface
  surf_varying <- surf[surf$method %in% varying_methods, ]
  surf_varying$method <- factor(surf_varying$method, levels = varying_methods)
  unadj_ref_bias <- surf$NDE_bias[surf$method == "UNADJ" & surf$gamma_G == min(surf$gamma_G)]

  pA <- ggplot2::ggplot(surf_varying, ggplot2::aes(x = gamma_G, y = NDE_bias,
                        colour = method)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed",
                        linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = unadj_ref_bias, colour = method_cols["UNADJ"],
                        linetype = "dotted", linewidth = 0.8) +
    ggplot2::annotate("text", x = 0.25, y = unadj_ref_bias + 0.02,
                      label = "Your current estimate\n(no instrument, confounded)",
                      size = 2.8, colour = method_cols["UNADJ"], hjust = 0,
                      fontface = "bold") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(
      values = method_cols[varying_methods],
      labels = c("IV2SLS" = "2SLS\n(single instrument)",
                 "IV2SLS2" = "2-stage MR\n(two instruments)",
                 "PGC2Gm" = "NC-augmented\n(instruments + NCs)"),
      name = "Method") +
    ggplot2::labs(x = "Instrument strength (gamma_G)",
                  y = "NDE bias (estimate - true)",
                  title = "How estimates converge as the instrument strengthens",
                  tag = "A") +
    .figure_theme()

  # Panel B: Prospective estimates at target strength
  prosp <- prospect$prospective
  unadj_nde <- prosp$NDE_mean[prosp$method == "UNADJ"]
  unadj_nie <- prosp$NIE_mean[prosp$method == "UNADJ"]
  prosp_sub <- prosp[prosp$method %in% varying_methods, ]
  prosp_sub$method <- factor(prosp_sub$method, levels = varying_methods)

  prosp_forest <- rbind(
    data.frame(method = prosp_sub$method, effect = "NDE",
               estimate = prosp_sub$NDE_mean),
    data.frame(method = prosp_sub$method, effect = "NIE",
               estimate = prosp_sub$NIE_mean))
  prosp_forest$effect <- factor(prosp_forest$effect, levels = c("NDE", "NIE"))

  true_refs <- data.frame(effect = factor(c("NDE", "NIE"), levels = c("NDE", "NIE")),
                          xint = c(0.10, 0.15))
  unadj_refs <- data.frame(effect = factor(c("NDE", "NIE"), levels = c("NDE", "NIE")),
                           xint = c(unadj_nde, unadj_nie))

  pB <- ggplot2::ggplot(prosp_forest, ggplot2::aes(x = estimate, y = method,
                        colour = method)) +
    ggplot2::geom_vline(data = true_refs, ggplot2::aes(xintercept = xint),
                        linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    ggplot2::geom_vline(data = unadj_refs, ggplot2::aes(xintercept = xint),
                        linetype = "dotted", colour = method_cols["UNADJ"],
                        linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(
      values = method_cols[c("IV2SLS", "IV2SLS2", "PGC2Gm")],
      labels = c("IV2SLS" = "2SLS", "IV2SLS2" = "2-stage MR", "PGC2Gm" = "NC-augmented"),
      name = "Method") +
    ggplot2::facet_wrap(~ effect, ncol = 2, scales = "free_x") +
    ggplot2::labs(
      x = "Estimate (dashed = true effect, dotted = your current confounded estimate)",
      y = "Method",
      title = "With strong instruments (F = 189 for exposure, F = 866 for mediator)",
      tag = "B") +
    .figure_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_line(
                     colour = "grey92", linewidth = 0.3))

  # Panel C: CI coverage vs instrument strength (NDE and NIE)
  has_cov <- all(c("NDE_coverage", "NIE_coverage") %in% names(surf))
  pC <- NULL
  if (has_cov) {
    cov_df <- rbind(
      data.frame(gamma_G = surf_varying$gamma_G, method = surf_varying$method,
                 estimand = "NDE", coverage = surf_varying$NDE_coverage),
      data.frame(gamma_G = surf_varying$gamma_G, method = surf_varying$method,
                 estimand = "NIE", coverage = surf_varying$NIE_coverage))
    cov_df$estimand <- factor(cov_df$estimand, levels = c("NDE", "NIE"))
    pC <- ggplot2::ggplot(cov_df, ggplot2::aes(x = gamma_G, y = coverage,
                          colour = method, linetype = estimand)) +
      ggplot2::geom_hline(yintercept = 0.95, colour = "grey50", linetype = "dashed",
                          linewidth = 0.5) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.2) +
      ggplot2::scale_colour_manual(
        values = method_cols[varying_methods],
        labels = c("IV2SLS" = "2SLS", "IV2SLS2" = "2-stage MR",
                   "PGC2Gm" = "NC-augmented"),
        name = "Method") +
      ggplot2::scale_linetype_discrete(name = "Estimand") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "Instrument strength (gamma_G)",
                    y = "95% CI coverage",
                    title = "Confidence-interval coverage vs instrument strength",
                    tag = "C") +
      .figure_theme()
  }

  fig <- if (!is.null(pC)) (pA / pB / pC) else (pA / pB)
  fig <- fig +
    patchwork::plot_annotation(
      title = "Prospective analysis: what adding instruments would do to your estimates",
      subtitle = paste0("Current data has X, M, Y only. The unadjusted estimate is confounded;",
                        "the identified estimators converge to the true effect as the instrument strengthens."),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12,
                                           colour = "grey10", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30", hjust = 0.5)))

  if (!is.null(pC)) height <- max(height, 11)
  .save_figure(fig, file, width, height)
  fig
}

# ============================================================
# Pleiotropy sweep
# ============================================================

#' Pleiotropy sweep figure
#'
#' Produces a faceted figure showing estimator bias and Type I error
#' under horizontal pleiotropy, across three confounding strengths.
#'
#' @param sensitivity Result of \code{gan_pleiotropy_sensitivity()}.
#' @param file Optional file path to save the figure.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' sens <- gan_pleiotropy_sensitivity(NULL, pleio_grid = c(0, 0.1),
#'   conf_grid = 0.8, n_iter = 2, n_samples = 100, n_features = 5)
#' plot_pleiotropy_sweep(sens)
plot_pleiotropy_sweep <- function(sensitivity, file = NULL,
                                   width = 10, height = 5.5) {
  .figures_check_deps()

  s <- sensitivity$summary

  bias_df <- s[s$arm == "alt", c("pleio", "conf_strength", "method", "bias")]
  bias_df$metric <- "Bias"
  bias_df$value <- bias_df$bias

  t1e_df <- s[s$arm == "null", c("pleio", "conf_strength", "method", "power")]
  t1e_df$metric <- "Type I Error"
  t1e_df$value <- t1e_df$power

  plot_df <- rbind(
    bias_df[, c("pleio", "conf_strength", "method", "metric", "value")],
    t1e_df[, c("pleio", "conf_strength", "method", "metric", "value")])

  plot_df$metric <- factor(plot_df$metric, levels = c("Bias", "Type I Error"))
  conf_levels <- sort(unique(plot_df$conf_strength))
  conf_labels <- paste0("Conf = ", conf_levels)
  plot_df$conf_strength <- factor(plot_df$conf_strength,
                                  levels = conf_levels, labels = conf_labels)

  method_colors <- c(UNADJ = "#888888", DIRECT = "#E07B00", COCA = "#3A9EC2",
                     IV2SLS = "#27A062", PGC = "#C455A8")
  method_colors <- method_colors[names(method_colors) %in% unique(plot_df$method)]
  method_order <- c("IV2SLS", "PGC", "COCA", "DIRECT", "UNADJ")
  method_order <- method_order[method_order %in% unique(plot_df$method)]
  plot_df$method <- factor(plot_df$method, levels = method_order)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pleio, y = value,
                       color = method, group = method)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(metric ~ conf_strength, scales = "free_y") +
    ggplot2::scale_color_manual(values = method_colors, name = "Estimator") +
    ggplot2::geom_hline(data = data.frame(metric = "Type I Error", yint = 0.05),
                        ggplot2::aes(yintercept = yint), linetype = "dashed",
                        color = "grey50", inherit.aes = FALSE) +
    ggplot2::geom_hline(data = data.frame(metric = "Bias", yint = 0),
                        ggplot2::aes(yintercept = yint), linetype = "dashed",
                        color = "grey50", inherit.aes = FALSE) +
    ggplot2::labs(x = "Pleiotropy strength (direct G\u2192Y coefficient)",
                  y = NULL) +
    ggplot2::scale_x_continuous(breaks = c(0, 0.05, 0.10)) +
    .figure_theme() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 10, face = "bold"),
                   legend.position = "bottom",
                   panel.border = ggplot2::element_rect(color = "black", fill = NA,
                                                        linewidth = 0.5),
                   panel.spacing.y = ggplot2::unit(0.5, "lines"))

  .save_figure(p, file, width, height)
  p
}

# ============================================================
# Instrument strength sweep
# ============================================================

#' Sweep instrument strength
#'
#' Runs a custom DGP that varies the G->X coefficient (pi_GX) and
#' records IV2SLS performance binned by first-stage partial F.
#' Uses \code{fit_iv2sls()} with \code{min_f = 0} (no weak-IV guard)
#' so bias is visible across all instrument strengths.
#'
#' @param pi_GX_grid Numeric vector of G->X coefficients to sweep.
#' @param n_iter Number of replications per grid point.
#' @param n Sample size.
#' @param k Number of confounders.
#' @param conf_str Confounding strength.
#' @param tau True total effect.
#' @param coverage Negative-control coverage.
#' @param n_cores Number of parallel workers for the replicate loops.
#' Default 1 (sequential). Uses \code{parallel::mclapply} on Unix
#' and a PSOCK cluster on Windows.
#' @return A data frame with columns \code{pi_GX}, \code{arm},
#' \code{iter}, \code{partial_F}, \code{beta}, \code{se},
#' \code{pvalue}, \code{rejected}.
#' @export
#' @examples
#' res <- sweep_instrument_strength(pi_GX_grid = c(0.05, 0.2), n_iter = 2)
#' head(res)
sweep_instrument_strength <- function(pi_GX_grid = c(0.02, 0.05, 0.10, 0.15,
                                        0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80),
                                       n_iter = 50, n = 200, k = 1,
                                       conf_str = 0.8, tau = 0.25, coverage = 0.7,
                                       n_cores = 1) {

  message("Running instrument-strength sweep...")
  results <- data.frame()

  run_one <- function(pi_GX, arm, seed) {
    withr::local_seed(seed)
    U <- matrix(rnorm(n * k), n, k)
    gi <- simulate_single_genetic_instrument(n, seed = seed)
    G <- gi$G
    aX <- conf_str * rep(1, k) / sqrt(k)
    X <- as.numeric(scale(pi_GX * G + as.numeric(U %*% aX) + rnorm(n, 0, 0.5)))
    gY <- conf_str * runif(k, 0.4, 0.8) / sqrt(k)
    exo_Y <- rnorm(n)
    effect <- if (arm == "alt") tau else 0
    Y <- effect * X + as.numeric(U %*% gY) + 0.4 * exo_Y + rnorm(n, 0, 0.2)
    W <- coverage * U[, 1] + (1 - coverage) * rnorm(n) + rnorm(n, 0, 0.3)
    W <- as.numeric(scale(W))
    list(Y = Y, X = X, G = G, W = W)
  }

  for (pi_GX in pi_GX_grid) {
    arms <- c("alt", "null")
    task_grid <- expand.grid(arm = arms, iter = seq_len(n_iter),
                             KEEP.OUT.ATTRS = FALSE)
    rows <- .parallel_lapply(seq_len(nrow(task_grid)), function(ti) {
      arm <- task_grid$arm[ti]
      iter <- task_grid$iter[ti]
      seed <- as.integer(10000 * pi_GX + 100 * (arm == "null") + iter)
      d <- run_one(pi_GX, arm, seed)
      fs <- lm(d$X ~ d$G + d$W)
      sm <- summary(fs)$coefficients
      partial_F <- as.numeric(sm["d$G", "t value"]^2)
      iv <- fit_iv2sls(d$Y, d$X, d$G, d$W, min_f = 0)
      data.frame(
        pi_GX = pi_GX, arm = arm, iter = iter,
        partial_F = partial_F, beta = iv$beta, se = iv$se,
        pvalue = iv$pvalue,
        rejected = !is.na(iv$pvalue) & iv$pvalue < 0.05)
    }, n_cores = n_cores, progress = paste0(" pi_GX=", pi_GX))
    results <- rbind(results, do.call(rbind, rows))
    message(sprintf(" pi_GX = %.2f done (mean F = %.1f)", pi_GX,
                    mean(results$partial_F[results$pi_GX == pi_GX])))
  }
  results
}

#' Instrument strength sweep figure
#'
#' Produces a 2-panel figure showing IV2SLS bias and Type I error
#' as a function of first-stage partial F-statistic, with the
#' conventional weak-instrument threshold (F = 10) marked.
#'
#' @param results Data frame from \code{sweep_instrument_strength()}.
#' @param tau True total effect used in the sweep.
#' @param file Optional file path to save the figure.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return A \code{patchwork} ggplot object.
#' @export
#' @examples
#' res <- sweep_instrument_strength(pi_GX_grid = c(0.05, 0.2), n_iter = 2)
#' plot_instrument_strength_sweep(res)
plot_instrument_strength_sweep <- function(results, tau = 0.25,
                                            file = NULL, width = 10, height = 4.5) {
  .figures_check_deps()

  alt <- results[results$arm == "alt", ]
  nul <- results[results$arm == "null", ]
  lv <- sort(unique(results$pi_GX))

  plot_df <- data.frame(
    pi_GX = lv,
    mean_F = vapply(lv, function(v) mean(results$partial_F[results$pi_GX == v], na.rm = TRUE), numeric(1)),
    bias = vapply(lv, function(v) mean(alt$beta[alt$pi_GX == v], na.rm = TRUE) - tau, numeric(1)),
    bias_sd = vapply(lv, function(v) sd(alt$beta[alt$pi_GX == v]), numeric(1)),
    t1e = vapply(lv, function(v) mean(nul$rejected[nul$pi_GX == v], na.rm = TRUE), numeric(1)))

  pts_bias <- results[results$arm == "alt" & !is.na(results$beta), ]
  pts_bias$bias <- pts_bias$beta - tau
  pts_t1e <- results[results$arm == "null" & !is.na(results$beta), ]

  pA <- ggplot2::ggplot() +
    ggplot2::geom_point(data = pts_bias,
                        ggplot2::aes(x = partial_F, y = bias),
                        alpha = 0.15, size = 0.8, color = "#27A062") +
    ggplot2::geom_point(data = plot_df,
                        ggplot2::aes(x = mean_F, y = bias), size = 2.5, color = "#27A062") +
    ggplot2::geom_line(data = plot_df,
                       ggplot2::aes(x = mean_F, y = bias), color = "#27A062", linewidth = 0.8) +
    ggplot2::geom_errorbar(data = plot_df,
                           ggplot2::aes(x = mean_F, ymin = bias - bias_sd, ymax = bias + bias_sd),
                           width = 0.05, color = "#27A062", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = 10, linetype = "dotted", color = "#FF9400", linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::ylim(-50, 50) +
    ggplot2::annotate("text", x = 11, y = 50, label = "F = 10",
                      hjust = 0, size = 3, color = "#FF9400") +
    ggplot2::labs(x = "First-stage partial F (log scale)",
                  y = "Bias (estimate \u2212 true effect)", title = "A") +
    .figure_theme()

  pB <- ggplot2::ggplot() +
    ggplot2::geom_point(data = pts_t1e,
                        ggplot2::aes(x = partial_F, y = as.numeric(rejected)),
                        alpha = 0.05, size = 0.8, color = "#27A062") +
    ggplot2::geom_point(data = plot_df,
                        ggplot2::aes(x = mean_F, y = t1e), size = 2.5, color = "#27A062") +
    ggplot2::geom_line(data = plot_df,
                       ggplot2::aes(x = mean_F, y = t1e), color = "#27A062", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = 10, linetype = "dotted", color = "#FF9400", linewidth = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(limits = c(-0.25, 0.25)) +
    ggplot2::annotate("text", x = 11, y = 0.2, label = "F = 10",
                      hjust = 0, size = 3, color = "#FF9400") +
    ggplot2::labs(x = "First-stage partial F (log scale)",
                  y = "Type I error rate", title = "B") +
    .figure_theme()

  fig <- pA + pB + patchwork::plot_layout(nrow = 1)

  .save_figure(fig, file, width, height)
  fig
}

# ============================================================
# NC coverage comparison
# (bias degradation as negative-control coverage drops)
# ============================================================

#' NC coverage comparison figure
#'
#' Produces a 4-panel figure comparing PGC2, PGC2Gm, and IV2SLS2 as
#' negative-control coverage of each path's confounder composite
#' (omega_1, omega_2) drops. This is the coverage-axis analogue of the
#' instrument-violation degradation surface: it shows how bias
#' accumulates as the proxy panel captures less of the confounder.
#'
#' Panels: (A) NDE bias vs omega_1, (B) NIE bias vs omega_2,
#' (C) NDE/NIE bias vs omega (both paths), (D) NIE Type I error vs
#' omega_2.
#'
#' @param omega1_sweep Summary from
#' \code{sweep_mediation_param("omega_1", ...)}: a list with
#' \code{$summary} (columns \code{param_value}, \code{method},
#' \code{NDE_bias}, \code{NIE_bias}).
#' @param omega2_sweep Same, sweeping \code{omega_2}.
#' @param t1e_omega2 Result of
#' \code{sweep_mediation_null_by_conf(...)} run across an \code{omega_2}
#' grid: a data frame with columns \code{omega_2}, \code{method},
#' \code{NIE_type1}. If \code{NULL}, Panel D is omitted.
#' @param file Optional file path to save the figure.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @return A \code{patchwork} ggplot object.
#' @export
#' @examples
#' # Toy sweeps standing in for sweep_mediation_param("omega_1"/"omega_2", ...)
#' mk <- function() {
#'   s <- expand.grid(param_value = c(0.3, 0.7),
#'     method = c("IV2SLS2", "PGC2", "PGC2Gm"))
#'   s$NDE_bias <- rnorm(nrow(s), 0, 0.05)
#'   s$NIE_bias <- rnorm(nrow(s), 0, 0.05)
#'   list(summary = s)
#' }
#' plot_nc_coverage_comparison(mk(), mk())
plot_nc_coverage_comparison <- function(omega1_sweep, omega2_sweep,
                                        t1e_omega2 = NULL,
                                        file = NULL, width = 10, height = 12) {
  .figures_check_deps()

  focus_methods <- c("IV2SLS2", "PGC2", "PGC2Gm")
  method_colors <- c(IV2SLS2 = "#0279EE", PGC2 = "#75A025", PGC2Gm = "#FD9BED")

  base_theme <- .figure_theme() +
    ggplot2::theme(legend.position = "bottom")

  make_panel <- function(df, xlab, ylab, hline_val = NULL, yrange = NULL) {
    df$method <- factor(df$method, levels = focus_methods)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = xvar, y = value,
                         color = method, group = method)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = method_colors, name = "Estimator") +
      ggplot2::labs(x = xlab, y = ylab) +
      base_theme
    if (!is.null(hline_val))
      p <- p + ggplot2::geom_hline(yintercept = hline_val, linetype = "dashed",
                                   color = "grey60")
    if (!is.null(yrange))
      p <- p + ggplot2::coord_cartesian(ylim = yrange)
    p
  }

  s1 <- omega1_sweep$summary
  s1 <- s1[s1$method %in% focus_methods, ]
  s2 <- omega2_sweep$summary
  s2 <- s2[s2$method %in% focus_methods, ]

  # Panel A: NDE bias vs omega_1
  pA <- make_panel(
    data.frame(method = s1$method, value = s1$NDE_bias, xvar = s1$param_value),
    xlab = "Coverage of X->M confounder (omega_1)", ylab = "NDE bias",
    hline_val = 0)

  # Panel B: NIE bias vs omega_2
  pB <- make_panel(
    data.frame(method = s2$method, value = s2$NIE_bias, xvar = s2$param_value),
    xlab = "Coverage of M->Y confounder (omega_2)", ylab = "NIE bias",
    hline_val = 0)

  # Panel C: NDE & NIE bias vs omega (both paths, from the omega_1 sweep
  # when omega_1 = omega_2; otherwise uses each path's own sweep)
  nde_c <- data.frame(method = s1$method, metric = "NDE bias",
                      value = s1$NDE_bias, xvar = s1$param_value)
  nie_c <- data.frame(method = s2$method, metric = "NIE bias",
                      value = s2$NIE_bias, xvar = s2$param_value)
  bias_c <- rbind(nde_c, nie_c)
  bias_c$method <- factor(bias_c$method, levels = focus_methods)
  bias_c$metric <- factor(bias_c$metric, levels = c("NDE bias", "NIE bias"))
  pC <- ggplot2::ggplot(bias_c, ggplot2::aes(x = xvar, y = value,
                       color = method, group = method)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ metric, ncol = 2, scales = "free_y") +
    ggplot2::scale_color_manual(values = method_colors, name = "Estimator") +
    ggplot2::labs(x = "Negative-control coverage (omega)", y = "Bias") +
    base_theme

  panels <- list(pA, pB, pC)

  # Panel D: NIE Type I error vs omega_2
  if (!is.null(t1e_omega2)) {
    t1e_df <- t1e_omega2[t1e_omega2$method %in% focus_methods, ]
    xcol <- if ("omega_2" %in% names(t1e_df)) t1e_df$omega_2 else t1e_df$param_value
    pD <- make_panel(
      data.frame(method = t1e_df$method, value = t1e_df$NIE_type1, xvar = xcol),
      xlab = "Coverage of M->Y confounder (omega_2)", ylab = "NIE Type I error",
      hline_val = 0.05, yrange = c(0, 1))
    panels <- c(panels, list(pD))
  }

  if (length(panels) == 4) {
    fig <- (panels[[1]] / panels[[2]] / panels[[3]] / panels[[4]])
  } else {
    fig <- (panels[[1]] / panels[[2]] / panels[[3]])
  }
  fig <- fig + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")

  .save_figure(fig, file, width, height)
  fig
}
