# =============================================================================
# PlantVelocity: diagnostics.R
# Four-panel diagnostics for IR-excluded two-state dynamics
# =============================================================================


#' @title Plot Diagnostics Dashboard
#' @description
#' Plot four diagnostics panels: fitted likelihood distribution, gene-filter
#' counts, fitted-time versus latent-time agreement, and top fitted genes.
#'
#' @param pv A \code{plantvelo} object with recovered dynamics.
#' @param n_top Integer. Number of top-likelihood genes shown.
#' @param group_by Character scalar or \code{NULL}. Metadata column used to
#'   colour the fitted-time panel.
#' @param input.color Character vector or \code{NULL}. Optional colour palette.
#'
#' @return A \code{patchwork} object, or a named list with \code{likelihood},
#'   \code{gene_filter}, \code{fit_t_scatter}, and \code{top_genes} when
#'   \pkg{patchwork} is unavailable.
#'
#' @examples
#' \dontrun{
#' plot_diagnostics(pv)
#' plot_diagnostics(pv, group_by = "cell_type", n_top = 30)
#' }
#' @export
plot_diagnostics <- function(pv,
                             n_top       = 10L,
                             group_by    = NULL,
                             input.color = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  # ----------------extract shared data----------------
  params     <- pv@kinetics[["params"]]
  vgenes     <- pv@velocity[["genes"]]
  all_genes  <- rownames(pv@layers[["spliced"]])
  n_total    <- length(all_genes)

  # ==========================================================================
  # Panel 1 — Likelihood histogram
  # ==========================================================================

  p1 <- .diag_likelihood_hist(params)

  # ==========================================================================
  # Panel 2 — Gene filter summary
  # ==========================================================================

  p2 <- .diag_gene_filter(params, vgenes, n_total)

  # ==========================================================================
  # Panel 3 — mean(fit_t) vs latent_time scatter
  # ==========================================================================

  p3 <- .diag_fit_t_scatter(pv, vgenes, group_by, input.color)

  # ==========================================================================
  # Panel 4: Top-gene likelihood bar chart
  # ==========================================================================

  p4 <- .diag_top_genes(params, n_top)

  # ==========================================================================
  # Assemble panels
  # ==========================================================================

  panels <- list(
    likelihood = p1,
    gene_filter = p2,
    fit_t_scatter = p3,
    top_genes = p4
  )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("Install 'patchwork' for a combined layout. Returning list of panels.")
    return(panels)
  }

  out <- patchwork::wrap_plots(panels, ncol = 2L)

  out
}


# =============================================================================
# Internal panel builders
# =============================================================================

#' Shared ggplot2 theme for diagnostics panels
#'
#' Returns a list of ggplot2 theme components shared across all panels in
#' \code{diagnostics.R}.
#'
#' @return A list of ggplot2 theme objects.
#' @keywords internal
.theme_diag <- function() {
  list(
    ggplot2::theme_classic(base_size = 11),
    ggplot2::theme(
      panel.border    = ggplot2::element_rect(fill = NA, colour = "black",
                                              linewidth = 0.5),
      axis.line       = ggplot2::element_blank(),
      plot.title      = ggplot2::element_text(size = 12, color = "black",
                                              hjust = 0.5),
      axis.text       = ggplot2::element_text(size = 10, color = "black"),
      axis.title      = ggplot2::element_text(size = 12, color = "black"),
      legend.title    = ggplot2::element_text(size = 10, color = "black"),
      legend.text     = ggplot2::element_text(size = 9,  color = "black"),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )
  )
}


#' Build the likelihood histogram panel
#'
#' @param params data.frame of fitted kinetic parameters.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.diag_likelihood_hist <- function(params) {
  if (is.null(params) || !("fit_likelihood" %in% colnames(params))) {
    return(.diag_placeholder("fit_likelihood\nnot available"))
  }

  ll <- as.numeric(params[["fit_likelihood"]])
  ll <- ll[is.finite(ll)]
  if (length(ll) == 0)
    return(.diag_placeholder("No finite\nlikelihood values"))

  df <- data.frame(ll = ll)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$ll)) +
    ggplot2::geom_histogram(binwidth = 0.02, boundary = 0,
                            fill = "#7f9eba", colour = "white",
                            linewidth = 0.2) +
    ggplot2::labs(
      title = "Likelihood distribution",
      x     = "fit_likelihood",
      y     = "Gene Number"
    ) +
    .theme_diag()
}


#' Build the gene filter summary panel
#'
#' @param params data.frame of fitted kinetic parameters.
#' @param vgenes Named logical vector of selected velocity genes.
#' @param ir_stats data.frame of IR gene statistics.
#' @param n_total Integer. Total number of genes.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.diag_gene_filter <- function(params, vgenes, n_total) {
  n_fitted   <- if (!is.null(params) && "fit_likelihood" %in% colnames(params))
    sum(!is.na(params[["fit_likelihood"]])) else 0L
  n_velocity <- if (!is.null(vgenes)) sum(vgenes, na.rm = TRUE) else 0L

  df <- data.frame(
    category = factor(
      c("Total genes", "Fitted", "Velocity genes"),
      levels = c("Total genes", "Fitted", "Velocity genes")
    ),
    count = c(n_total, n_fitted, n_velocity)
  )

  ggplot2::ggplot(df, ggplot2::aes(
    x    = .data$category,
    y    = .data$count,
    fill = .data$category
  )) +
    ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$count),
                       vjust = -0.4, size = 3, colour = "black") +
    ggplot2::scale_fill_manual(values = c(
      "Total genes"    = "#BFBFBF",
      "Fitted"         = "#7f9eba",
      "Velocity genes" = "#d17c7c"
    )) +
    ggplot2::labs(
      title = "Gene filter summary",
      x     = NULL,
      y     = "Gene Number"
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    .theme_diag() +
    list(ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    ))
}


#' Build the fit_t versus latent_time scatter panel
#'
#' Uses the mean fitted time per cell, averaged across velocity genes, and
#' compares it to the global \code{latent_time}.
#'
#' @param pv A \code{plantvelo} object.
#' @param vgenes Named logical vector of selected velocity genes.
#' @param group_by Character scalar or \code{NULL}. Metadata column used for
#'   colouring.
#' @param input.color Character vector or \code{NULL}. Optional colour palette.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.diag_fit_t_scatter <- function(pv, vgenes, group_by, input.color) {
  fit_t   <- pv@layers[["fit_t"]]
  lt      <- pv@meta.data[["latent_time"]]

  if (is.null(fit_t) || is.null(lt))
    return(.diag_placeholder("fit_t or latent_time\nnot available"))

  # Use only velocity genes when computing per-cell mean fit_t
  if (!is.null(vgenes) && sum(vgenes, na.rm = TRUE) > 0) {
    vg_names <- names(vgenes)[vgenes]
    vg_idx   <- which(rownames(fit_t) %in% vg_names)
    if (length(vg_idx) == 0) vg_idx <- seq_len(nrow(fit_t))
  } else {
    vg_idx <- seq_len(nrow(fit_t))
  }

  mean_ft <- colMeans(make_dense_matrix(fit_t[vg_idx, , drop = FALSE]),
                      na.rm = TRUE)
  n_cells <- length(lt)
  if (length(mean_ft) != n_cells)
    return(.diag_placeholder("Dimension mismatch:\nfit_t vs cells"))

  df <- data.frame(
    mean_fit_t  = as.numeric(mean_ft),
    latent_time = as.numeric(lt)
  )

  # Annotate Pearson correlation
  valid <- is.finite(df$mean_fit_t) & is.finite(df$latent_time)
  r_val <- if (sum(valid) > 5)
    round(stats::cor(df$mean_fit_t[valid], df$latent_time[valid]), 3)
  else NA_real_

  ann_label <- if (is.finite(r_val))
    paste0("r = ", r_val) else ""

  # Colour handling
  if (!is.null(group_by) && group_by %in% colnames(pv@meta.data)) {
    df$color    <- pv@meta.data[[group_by]]
    is_discrete <- !is.numeric(df$color)

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x      = .data$mean_fit_t,
      y      = .data$latent_time,
      colour = .data$color
    )) +
      ggplot2::geom_point(size = 0.5, alpha = 0.6)

    if (is_discrete) {
      if (!is.null(input.color)) {
        p <- p + ggplot2::scale_colour_manual(values = input.color,
                                              name = group_by)
      } else {
        p <- p + ggplot2::scale_colour_discrete(name = group_by)
      }
    } else {
      cols <- input.color %||% c("#440154", "#31688E", "#35B779", "#FDE725")
      p <- p + ggplot2::scale_colour_gradientn(colours = cols, name = group_by)
    }
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data$mean_fit_t,
      y = .data$latent_time
    )) +
      ggplot2::geom_point(size = 0.4, alpha = 0.5, colour = "#7f9eba")
  }

  p <- p +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                         colour = "grey30", linewidth = 0.6,
                         se = FALSE, inherit.aes = FALSE,
                         mapping = ggplot2::aes(
                           x = .data$mean_fit_t,
                           y = .data$latent_time
                         ),
                         data = df[valid, ]) +
    ggplot2::annotate("text",
                      x      = min(df$mean_fit_t, na.rm = TRUE),
                      y      = max(df$latent_time, na.rm = TRUE),
                      label  = ann_label,
                      hjust  = 0, vjust = 1,
                      size   = 3, colour = "grey20") +
    ggplot2::labs(
      title = "Fitted time vs. latent time",
      x     = "mean fit_t (velocity genes)",
      y     = "latent_time"
    ) +
    .theme_diag()

  p
}


#' Build the top-gene likelihood bar chart
#'
#' @param params Data frame of fitted kinetic parameters.
#' @param n_top Integer. Number of top genes shown.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.diag_top_genes <- function(params, n_top) {
  if (is.null(params) || !("fit_likelihood" %in% colnames(params))) {
    return(.diag_placeholder("fit_likelihood\nnot available"))
  }

  df <- data.frame(
    gene             = rownames(params),
    fit_likelihood   = as.numeric(params[["fit_likelihood"]]),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$fit_likelihood), , drop = FALSE]

  .diag_one_bar(
    df,
    n_top,
    fill = "#7f9eba",
    title = "Top fitted genes"
  )
}


#' Build a single likelihood bar chart
#'
#' Internal helper reused by \code{.diag_top_genes()}.
#'
#' @param df data.frame of ranked genes and likelihood values.
#' @param n_top Integer. Number of top genes to display.
#' @param fill Character scalar. Bar fill colour.
#' @param title Character scalar. Panel title.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.diag_one_bar <- function(df, n_top, fill, title) {
  if (nrow(df) == 0)
    return(.diag_placeholder(sub("top.*", "no genes", title)))

  df <- df[order(df$fit_likelihood, decreasing = TRUE), , drop = FALSE]
  df <- head(df, as.integer(n_top))
  df$gene <- factor(df$gene, levels = rev(df$gene))

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data$fit_likelihood,
    y = .data$gene
  )) +
    ggplot2::geom_col(width = 0.7, fill = fill) +
    ggplot2::labs(title = title, x = "fit_likelihood", y = NULL) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    .theme_diag()
}


#' Build a placeholder panel
#'
#' Used when the required data for a diagnostics panel are unavailable.
#'
#' @param msg Character scalar. Message shown in the placeholder panel.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.diag_placeholder <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                      size = 4, colour = "grey50", hjust = 0.5) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey80",
                                           linewidth = 0.4)
    )
}
