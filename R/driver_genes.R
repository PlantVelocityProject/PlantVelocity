# =============================================================================
# PlantVelocity: driver_genes.R
# Driver gene ranking based on spliced and unspliced velocity coherence
# =============================================================================


#' @title Rank Driver Genes by Velocity Coherence
#' @description
#' Filter fitted genes by likelihood, then rank them using velocity coherence
#' computed from the spliced and unspliced modalities.
#'
#' The primary score is controlled by \code{coherence_mode}:
#' \describe{
#'   \item{\code{"s_only"}}{Use spliced coherence only.}
#'   \item{\code{"mean"}}{Average spliced and unspliced coherence.}
#' }
#'
#' @section Output written to \code{@misc}:
#' \code{velocity_gene_stats} contains exactly \code{gene},
#' \code{coherence}, \code{coherence_s}, \code{coherence_u},
#' \code{mean_abs_velocity}, \code{n_active}, and \code{fit_likelihood}.
#'
#' @param pv A \code{plantvelo} object with a velocity graph.
#' @param vkey Character scalar. Spliced velocity layer key.
#' @param modality_s Character scalar. Smoothed spliced layer key.
#' @param modality_u Character scalar. Smoothed unspliced layer key.
#' @param n_top Integer or \code{NULL}. Maximum number of genes to retain.
#' @param min_likelihood Numeric scalar. Absolute likelihood threshold.
#' @param likelihood_quantile Numeric in \eqn{(0, 1]}. Fraction of fitted genes
#'   retained after ranking by likelihood.
#' @param coherence_mode One of \code{"s_only"} or \code{"mean"}.
#' @param self_transitions Logical. Include self transitions when constructing
#'   the transition matrix.
#' @param verbose Logical. Print progress messages.
#'
#' @return A \code{plantvelo} object with
#'   \code{@misc$velocity_gene_stats} populated.
#'
#' @examples
#' \dontrun{
#' pv <- rank_velocity_genes(pv)
#' pv <- rank_velocity_genes(pv, coherence_mode = "mean")
#' }
#' @export
rank_velocity_genes <- function(pv,
                                vkey               = "velocity",
                                modality_s         = "Ms",
                                modality_u         = "Mu",
                                n_top              = NULL,
                                min_likelihood     = 0,
                                likelihood_quantile = 0.3,
                                coherence_mode     = c("s_only", "mean"),
                                self_transitions   = TRUE,
                                verbose            = TRUE) {

  coherence_mode <- match.arg(coherence_mode)

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (!is.numeric(likelihood_quantile) || length(likelihood_quantile) != 1 ||
      likelihood_quantile <= 0 || likelihood_quantile > 1)
    stop("`likelihood_quantile` must be a single numeric value in (0, 1].")

  gkey <- paste0(vkey, "_graph")
  if (is.null(pv@graphs[[gkey]]))
    stop(sprintf(
      "Velocity graph '%s' not found. Run compute_velocity_graph() first.", gkey
    ))

  V_s <- pv@layers[[vkey]]
  if (is.null(V_s))
    stop(sprintf("Velocity layer '%s' not found in @layers.", vkey))

  # ----------------collect modality-specific velocities and smoothed expression----------------

  V_u <- pv@layers[[paste0(vkey, "_u")]]

  Ms <- pv@moments[[modality_s]] %||% t(pv@layers[["spliced"]])
  Mu <- pv@moments[[modality_u]] %||% t(pv@layers[["unspliced"]])

  n_genes    <- nrow(V_s)
  n_cells    <- ncol(V_s)
  gene_names <- rownames(V_s)

  params <- pv@kinetics[["params"]]
  if (!is.null(params)) {
    params <- .validate_dynamics_schema(pv)
  }

  # ----------------extract fit likelihood and valid fitted genes----------------

  fit_ll <- rep(NA_real_, n_genes)
  if (!is.null(params) && "fit_likelihood" %in% colnames(params))
    fit_ll <- params[match(gene_names, rownames(params)), "fit_likelihood"]

  valid_fits <- rep(FALSE, n_genes)
  if (!is.null(params)) {
    valid_fits <- .valid_velocity_fits(params)[
      match(gene_names, rownames(params))
    ]
    valid_fits[is.na(valid_fits)] <- FALSE
  }

  # ----------------Step 1: quality control by fit_likelihood----------------

  # 1a. Global absolute lower bound
  keep_genes <- valid_fits
  if (min_likelihood > 0)
    keep_genes <- keep_genes & (is.finite(fit_ll) & fit_ll >= min_likelihood)

  # 1b. Quantile filter across all fitted genes
  if (likelihood_quantile < 1) {
    eligible <- keep_genes & is.finite(fit_ll)
    if (any(eligible)) {
      threshold <- stats::quantile(
        fit_ll[eligible],
        probs = 1 - likelihood_quantile,
        na.rm = TRUE
      )
      keep_genes <- eligible & fit_ll >= threshold
    } else {
      keep_genes[] <- FALSE
    }

    if (verbose)
      message(sprintf(
        "[QC] Likelihood filter: kept %d / %d fitted genes (top %.0f%%).",
        sum(keep_genes), sum(valid_fits),
        likelihood_quantile * 100
      ))
  }

  # ----------------build transition matrix----------------

  T_mat <- build_transition_matrix(pv, vkey = vkey,
                                   self_transitions = self_transitions)

  # ----------------Step 2: compute modality-specific coherence----------------

  cell_idx <- seq_len(n_cells)

  score_idx <- which(keep_genes)
  .score_modality <- function(V_mat, X_mat) {
    scores <- rep(NA_real_, n_genes)
    if (is.null(V_mat) || is.null(X_mat) || length(score_idx) == 0L) {
      return(scores)
    }

    velocity_idx <- if (!is.null(rownames(V_mat))) {
      match(gene_names[score_idx], rownames(V_mat))
    } else {
      score_idx
    }
    expression_idx <- if (!is.null(colnames(X_mat))) {
      match(gene_names[score_idx], colnames(X_mat))
    } else {
      score_idx
    }
    aligned <- !is.na(velocity_idx) & !is.na(expression_idx) &
      velocity_idx <= nrow(V_mat) & expression_idx <= ncol(X_mat)
    if (any(aligned)) {
      scores[score_idx[aligned]] <- .coherence_score(
        V_mat[velocity_idx[aligned], cell_idx, drop = FALSE],
        X_mat[cell_idx, expression_idx[aligned], drop = FALSE],
        T_mat
      )
    }
    scores
  }

  coh_s <- .score_modality(V_s, Ms)
  coh_u <- .score_modality(V_u, Mu)

  # Primary ranking score determined by coherence_mode
  coherence <- if (coherence_mode == "s_only") {
    coh_s
  } else {
    rowMeans(cbind(coh_s, coh_u), na.rm = TRUE)
  }

  # Auxiliary metrics
  mean_abs_vel <- rowMeans(abs(V_s), na.rm = TRUE)
  n_active     <- rowSums(V_s != 0, na.rm = TRUE)

  # ----------------assemble result table----------------

  result <- data.frame(
    gene              = gene_names,
    coherence         = coherence,
    coherence_s       = coh_s,
    coherence_u       = coh_u,
    mean_abs_velocity = mean_abs_vel,
    n_active          = n_active,
    fit_likelihood    = fit_ll,
    stringsAsFactors  = FALSE
  )

  # Apply likelihood QC mask
  result <- result[keep_genes, , drop = FALSE]

  # Remove invalid rows: coherence missing or all-zero velocity
  result <- result[!is.na(result$coherence) & result$n_active > 0, , drop = FALSE]

  # Step 3: sort by decreasing coherence
  result <- result[order(result$coherence, decreasing = TRUE, na.last = TRUE), ]

  if (!is.null(n_top)) result <- head(result, n_top)

  rownames(result) <- NULL

  # ----------------store results and report----------------

  pv@misc[["velocity_gene_stats"]] <- result

  if (verbose) {
    mode_desc <- if (coherence_mode == "s_only") {
      "coherence_s (spliced only)"
    } else {
      "mean(coherence_s, coherence_u)"
    }
    message(sprintf(
      "Ranked %d genes by %s.",
      nrow(result),
      mode_desc
    ))
    message("Results stored in @misc$velocity_gene_stats.")
  }

  return(pv)
}


#' @title Get Velocity Gene Statistics
#' @description
#' Retrieve the ranked gene statistics table produced by
#' \code{rank_velocity_genes()}.
#'
#' @param pv A \code{plantvelo} object.
#'
#' @return A data frame, or \code{NULL} if
#'   \code{rank_velocity_genes()} has not been run.
#' @examples
#' \dontrun{
#' stats <- get_velocity_gene_stats(pv)
#' }
#' @export
get_velocity_gene_stats <- function(pv) {
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  stats <- pv@misc[["velocity_gene_stats"]]
  if (is.null(stats)) {
    message("No gene stats found. Run rank_velocity_genes() first.")
    return(NULL)
  }

  stats
}


# =============================================================================
# Visualization
# =============================================================================

#' @title Plot Driver Genes
#' @description
#' Plot a positive horizontal lollipop chart of top genes ranked by velocity
#' coherence. Point colour represents \code{fit_likelihood}.
#'
#' @param pv A \code{plantvelo} object.
#' @param n_top Integer. Number of top genes shown.
#' @param title Character scalar or \code{NULL}. Plot title.
#' @param fill_quantile Numeric vector of length two defining likelihood colour
#'   limits when \code{use_quantile = TRUE}.
#' @param use_quantile Logical. Use quantile-based colour limits.
#' @param segment_size Numeric. Lollipop stem width.
#' @param point_size Numeric. Lollipop head size.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' pv <- rank_velocity_genes(pv)
#' plot_driver_genes(pv)
#' }
#' @export
plot_driver_genes <- function(pv,
                              n_top           = 20L,
                              title           = NULL,
                              fill_quantile   = c(0.05, 0.95),
                              use_quantile    = TRUE,
                              segment_size    = 0.5,
                              point_size      = 5) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  if (!requireNamespace("scales", quietly = TRUE))
    stop("Package 'scales' is required. Install with: install.packages('scales')")
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (!is.numeric(fill_quantile) || length(fill_quantile) != 2L ||
      any(!is.finite(fill_quantile)) ||
      fill_quantile[1] < 0 || fill_quantile[2] > 1 ||
      fill_quantile[1] >= fill_quantile[2]) {
    stop("`fill_quantile` must be a numeric vector of length 2 in [0, 1], e.g. c(0.05, 0.95).")
  }

  # ----------------retrieve ranked gene table----------------

  df <- get_velocity_gene_stats(pv)
  if (is.null(df))
    stop("Run rank_velocity_genes() first.")

  df <- df[order(df$coherence, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  df <- head(df, as.integer(n_top))
  df$gene <- factor(df$gene, levels = rev(df$gene))

  # ----------------set fit_likelihood colour limits----------------

  .safe_limits <- function(x, q = c(0.05, 0.95), use_quantile = TRUE) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(c(0, 1))

    lim <- if (use_quantile && length(x) > 2L) {
      as.numeric(stats::quantile(x, probs = q, na.rm = TRUE, names = FALSE))
    } else {
      range(x, na.rm = TRUE)
    }

    if (!all(is.finite(lim))) {
      lim <- range(x, na.rm = TRUE)
    }

    if (lim[1] == lim[2]) {
      eps <- if (lim[1] == 0) 1e-6 else abs(lim[1]) * 0.05
      lim <- c(lim[1] - eps, lim[2] + eps)
    }

    lim
  }

  ll_lim <- .safe_limits(
    df$fit_likelihood,
    q = fill_quantile,
    use_quantile = use_quantile
  )

  # ----------------build plot----------------

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$coherence, y = .data$gene)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$coherence, yend = .data$gene),
      linewidth = segment_size,
      colour = "#9aa7b1"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$fit_likelihood),
      size = point_size
    ) +
    # ggplot2::scale_colour_gradient(
    #   low = "#c1e7ff",
    #   high = "#006db3",
    #   name = "fit_likelihood",
    #   limits = ll_lim,
    #   oob = scales::squish,
    #   na.value = "grey80"
    # ) +
    ggplot2::scale_colour_gradientn(
      colours = c(
        "#D5E3EC",
        "#9DBDD0",
        "#5F8EAA",
        "#234F6D"
      ),
      name = "fit_likelihood",
      limits = ll_lim,
      oob = scales::squish,
      na.value = "grey80"
    ) +
    ggplot2::labs(
      title = title %||% paste0("Top ", n_top, " driver genes"),
      x     = "coherence",
      y     = NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 13, color = "black",
                                              hjust = 0.5),
      axis.text       = ggplot2::element_text(size = 10, color = "black"),
      axis.title      = ggplot2::element_text(size = 12, color = "black"),
      legend.title    = ggplot2::element_text(size = 10, color = "black"),
      legend.text     = ggplot2::element_text(size = 9,  color = "black"),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )

  p
}

# =============================================================================
# Internal helper functions
# =============================================================================

#' Compute velocity coherence score for a single modality
#'
#' Computes the projected neighbourhood velocity for each gene under the
#' transition matrix and correlates it with the observed velocity in the same
#' modality.
#'
#' @param V_mat Velocity matrix of size \eqn{n\_genes \times n\_cells}, or
#'   \code{NULL}.
#' @param X_mat Smoothed expression matrix of size
#'   \eqn{n\_cells \times n\_genes}, or \code{NULL}.
#' @param T_sub Row-normalised transition matrix of size
#'   \eqn{n\_cells \times n\_cells}.
#'
#' @return Numeric vector of coherence scores of length \eqn{n\_genes}. If the
#'   modality is unavailable, returns an all-\code{NA} vector.
#' @keywords internal
.coherence_score <- function(V_mat, X_mat, T_sub) {
  n_genes <- if (!is.null(V_mat)) nrow(V_mat) else
    if (!is.null(X_mat)) ncol(X_mat) else 0L

  na_vec <- rep(NA_real_, n_genes)

  if (is.null(V_mat) || is.null(X_mat)) return(na_vec)
  if (nrow(V_mat) == 0 || nrow(X_mat) == 0) return(na_vec)
  if (ncol(V_mat) != nrow(X_mat) ||
      nrow(T_sub) != nrow(X_mat) || ncol(T_sub) != nrow(X_mat)) {
    stop("Velocity, expression, and transition matrices have incompatible cell dimensions.")
  }

  expression_idx <- seq_len(n_genes)
  if (!is.null(rownames(V_mat)) && !is.null(colnames(X_mat))) {
    expression_idx <- match(rownames(V_mat), colnames(X_mat))
  }

  for (gene_idx in seq_len(n_genes)) {
    x_idx <- expression_idx[gene_idx]
    if (is.na(x_idx) || x_idx > ncol(X_mat)) next

    expression <- as.numeric(X_mat[, x_idx, drop = TRUE])
    velocity <- as.numeric(V_mat[gene_idx, , drop = TRUE])
    projected <- as.numeric(T_sub %*% expression) - expression
    na_vec[gene_idx] <- .cor_colwise(
      matrix(velocity, ncol = 1L),
      matrix(projected, ncol = 1L)
    )[1L]
  }

  na_vec
}


#' Compute Pearson correlation column-wise
#'
#' Efficiently computes Pearson correlations between corresponding columns of
#' two matrices without constructing a full correlation matrix.
#'
#' @param A Numeric matrix of size \eqn{n\_obs \times n\_genes}.
#' @param B Numeric matrix of size \eqn{n\_obs \times n\_genes}.
#'
#' @return Numeric vector of column-wise Pearson correlation coefficients.
#' @keywords internal
.cor_colwise <- function(A, B) {
  A_c   <- sweep(A, 2L, colMeans(A, na.rm = TRUE), `-`)
  B_c   <- sweep(B, 2L, colMeans(B, na.rm = TRUE), `-`)
  ss_A  <- sqrt(colSums(A_c^2, na.rm = TRUE))
  ss_B  <- sqrt(colSums(B_c^2, na.rm = TRUE))
  denom <- ss_A * ss_B
  res   <- colSums(A_c * B_c, na.rm = TRUE) / denom
  res[denom == 0] <- NA_real_
  res
}
