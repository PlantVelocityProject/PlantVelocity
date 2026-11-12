# =============================================================================
# PlantVelocity: driver_genes.R
# Driver gene ranking based on spliced and unspliced velocity coherence
# =============================================================================


#' @title Rank Driver Genes by Velocity Coherence
#' @description
#' Require sufficient unsmoothed U/S evidence and valid fits, inherit
#' \code{@velocity$genes} when present, filter by likelihood, then rank eligible
#' genes by velocity coherence. These are candidate dynamics-associated genes;
#' the score does not establish a causal regulatory role.
#'
#' The primary score is controlled by \code{coherence_mode}:
#' \describe{
#'   \item{\code{"s_only"}}{Use spliced coherence only.}
#'   \item{\code{"mean"}}{Average spliced and unspliced coherence, requiring
#'     finite scores in both modalities.}
#' }
#'
#' @section Output written to \code{@misc}:
#' \code{velocity_gene_stats} contains \code{gene},
#' \code{coherence}, \code{coherence_s}, \code{coherence_u},
#' \code{mean_abs_velocity}, \code{n_active}, \code{fit_likelihood}, and six
#' unsmoothed evidence columns: \code{raw_n_cells}, \code{raw_u_positive},
#' \code{raw_s_positive}, \code{raw_joint_positive}, \code{raw_u_total}, and
#' \code{raw_s_total}. \code{n_active} counts nonzero velocity estimates, not
#' independently detected RNA. This single table records all velocity-layer
#' genes, filter flags, rank and \code{selected}. Ranking thresholds and options
#' are stored in \code{velocity_params$ranking}, alongside the top-level velocity
#' calculation settings. No separate QC table or exclusion-reason column is stored.
#' Numeric missing values are stored as \code{NA}, including all-missing
#' velocity means that otherwise produce \code{NaN}.
#'
#' @param pv A \code{plantvelo} object with a velocity graph.
#' @param vkey Character scalar. Spliced velocity layer key.
#' @param modality_s Character scalar. Smoothed spliced layer key.
#' @param modality_u Character scalar. Smoothed unspliced layer key.
#' @param n_top Integer or \code{NULL}. Maximum number of genes to retain.
#' @param min_likelihood Numeric scalar. Absolute likelihood threshold.
#' @param likelihood_quantile Numeric in \eqn{(0, 1]}. Upper fraction retained
#'   by likelihood among genes passing raw-evidence, velocity and fit filters.
#' @param coherence_mode One of \code{"s_only"} or \code{"mean"}.
#' @param self_transitions Logical. Include self transitions when constructing
#'   the transition matrix.
#' @param verbose Logical. Print progress messages.
#' @param min_u_cells,min_s_cells Non-negative integer minimum numbers of cells
#'   with positive unsmoothed U and S respectively. Defaults are \code{20L}.
#' @param min_joint_cells Non-negative integer minimum number of cells with
#'   both U and S positive. Default \code{10L}.
#' @param min_u_counts,min_s_counts Non-negative finite minimum total unsmoothed
#'   U and S signal. Defaults are \code{20}. Totals use input-layer units, so use
#'   original count layers for count-based interpretation.
#'
#' @details
#' Evidence is read from \code{@layers$unspliced} and \code{@layers$spliced},
#' never from smoothed moments. Missing, negative or non-finite evidence is an
#' error. Named gene/cell axes are matched independently. Defaults are adjustable
#' conservative screening settings, not universal biological cutoffs. Setting
#' individual evidence thresholds to zero disables that requirement.
#' Existing dynamics and velocities are not recomputed or changed.
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
                                verbose            = TRUE,
                                min_u_cells        = 20L,
                                min_s_cells        = 20L,
                                min_joint_cells    = 10L,
                                min_u_counts       = 20,
                                min_s_counts       = 20) {

  coherence_mode <- match.arg(coherence_mode)

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(likelihood_quantile) ||
      likelihood_quantile <= 0 || likelihood_quantile > 1)
    stop("`likelihood_quantile` must be a single numeric value in (0, 1].")
  if (!scalar(min_likelihood) || min_likelihood < 0)
    stop("`min_likelihood` must be finite and non-negative.")
  if (!is.null(n_top) && (!scalar(n_top) || n_top < 0 || n_top != floor(n_top)))
    stop("`n_top` must be NULL or a non-negative integer.")
  for (key in c("vkey", "modality_s", "modality_u")) {
    x <- get(key)
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
      stop(sprintf("`%s` must be one non-empty string.", key))
  }
  raw_thresholds <- .raw_evidence_thresholds(min_u_cells, min_s_cells,
    min_joint_cells, min_u_counts, min_s_counts)
  for (key in c("self_transitions", "verbose")) {
    x <- get(key)
    if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(sprintf("`%s` must be TRUE or FALSE.", key))
  }

  gkey <- paste0(vkey, "_graph")
  if (is.null(pv@graphs[[gkey]]))
    stop(sprintf(
      "Velocity graph '%s' not found. Run compute_velocity_graph() first.", gkey
    ))

  V_s <- pv@layers[[vkey]]
  if (is.null(V_s))
    stop(sprintf("Velocity layer '%s' not found in @layers.", vkey))

  # ----------------collect modality-specific velocities and smoothed expression----------------

  axes <- .velocity_axes(pv)
  gene_names <- rownames(V_s)
  if (is.null(gene_names) && nrow(V_s) == length(axes$gene_ids)) gene_names <- axes$gene_ids
  if (is.null(gene_names) || anyNA(gene_names) || any(!nzchar(gene_names)) ||
      anyDuplicated(gene_names) || !all(gene_names %in% axes$gene_ids))
    stop("Velocity genes must have unique names matching the spliced layer.")
  cells <- axes$cell_ids
  V_s <- .driver_align_matrix(V_s, vkey, gene_names, cells)
  V_u <- pv@layers[[paste0(vkey, "_u")]]
  if (!is.null(V_u)) V_u <- .driver_align_matrix(V_u, paste0(vkey, "_u"), gene_names, cells)
  raw_s <- .driver_align_matrix(pv@layers$spliced, "spliced", gene_names, cells)
  raw_u <- .driver_align_matrix(pv@layers$unspliced, "unspliced", gene_names, cells)
  evidence <- .driver_raw_evidence(raw_u, raw_s)
  Ms <- if (is.null(pv@moments[[modality_s]])) Matrix::t(raw_s) else {
    .driver_align_matrix(pv@moments[[modality_s]], modality_s, gene_names, cells, "cells_by_genes")
  }
  Mu <- if (is.null(pv@moments[[modality_u]])) Matrix::t(raw_u) else {
    .driver_align_matrix(pv@moments[[modality_u]], modality_u, gene_names, cells, "cells_by_genes")
  }
  n_genes <- nrow(V_s)

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

  # Raw independent support is required before relative likelihood selection.
  raw_failures <- .raw_evidence_failures(evidence, raw_thresholds)
  pass_raw <- !Reduce(`|`, raw_failures)
  pass_velocity <- .driver_velocity_mask(pv@velocity$genes, gene_names)
  pass_likelihood <- is.finite(fit_ll) & fit_ll >= min_likelihood
  keep_genes <- valid_fits & pass_raw & pass_velocity & pass_likelihood
  below_quantile <- rep(FALSE, n_genes)
  threshold <- NA_real_
  if (likelihood_quantile < 1) {
    eligible <- keep_genes
    if (any(eligible)) {
      threshold <- stats::quantile(
        fit_ll[eligible],
        probs = 1 - likelihood_quantile,
        na.rm = TRUE
      )
      below_quantile[eligible] <- fit_ll[eligible] < threshold
      keep_genes <- eligible & !below_quantile
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

  T_mat <- if (any(keep_genes)) {
    .driver_align_matrix(build_transition_matrix(pv, vkey = vkey,
      self_transitions = self_transitions), "transition matrix", cells, cells)
  } else NULL

  # ----------------Step 2: compute modality-specific coherence----------------

  score_idx <- which(keep_genes)
  .score_modality <- function(V_mat, X_mat) {
    scores <- rep(NA_real_, n_genes)
    if (is.null(V_mat) || is.null(X_mat) || length(score_idx) == 0L) {
      return(scores)
    }

    scores[score_idx] <- .coherence_score(
      V_mat[score_idx, , drop = FALSE], X_mat[, score_idx, drop = FALSE], T_mat)
    scores
  }

  coh_s <- .score_modality(V_s, Ms)
  coh_u <- .score_modality(V_u, Mu)

  # Primary ranking score determined by coherence_mode
  coherence <- if (coherence_mode == "s_only") {
    coh_s
  } else {
    rowMeans(cbind(coh_s, coh_u))
  }

  # Auxiliary metrics
  mean_abs_vel <- as.numeric(Matrix::rowMeans(abs(V_s), na.rm = TRUE))
  n_active     <- as.numeric(Matrix::rowSums(V_s != 0, na.rm = TRUE))

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

  result <- cbind(result, evidence)
  eligible <- keep_genes & is.finite(coherence) & is.finite(mean_abs_vel) & n_active > 0
  ordered <- which(eligible)
  ordered <- ordered[order(-coherence[ordered], gene_names[ordered])]
  rank <- rep(NA_integer_, n_genes)
  rank[ordered] <- seq_along(ordered)
  selected_idx <- if (is.null(n_top)) ordered else head(ordered, n_top)
  selected <- rep(FALSE, n_genes)
  selected[selected_idx] <- TRUE
  result <- cbind(result, pass_fit = valid_fits, pass_velocity = pass_velocity,
    pass_raw = pass_raw, pass_likelihood = pass_likelihood & !below_quantile,
    rank = rank, selected = selected)
  result <- .normalize_driver_missing(result)
  rownames(result) <- NULL

  # ----------------store results and report----------------

  pv@misc[["velocity_gene_stats"]] <- result
  pv@misc[c("velocity_gene_qc","velocity_gene_qc_params")] <- NULL
  run <- pv@misc$velocity_params
  if (!is.list(run)) run <- list()
  run$qc_version <- NULL
  run$ranking <- list(
    min_u_cells = min_u_cells, min_s_cells = min_s_cells, min_joint_cells = min_joint_cells,
    min_u_counts = min_u_counts, min_s_counts = min_s_counts,
    min_likelihood = min_likelihood, likelihood_quantile = likelihood_quantile,
    likelihood_threshold = threshold, coherence_mode = coherence_mode,
    self_transitions = self_transitions, n_top = n_top, vkey = vkey,
    modality_s = modality_s, modality_u = modality_u)
  pv@misc$velocity_params <- run

  if (verbose) {
    message(sprintf("[QC] Raw U/S evidence: %d / %d genes pass; velocity mask: %d / %d.",
      sum(pass_raw), n_genes, sum(pass_velocity), n_genes))
    mode_desc <- if (coherence_mode == "s_only") {
      "coherence_s (spliced only)"
    } else {
      "mean(coherence_s, coherence_u)"
    }
    message(sprintf(
      "Ranked %d genes by %s.",
      length(selected_idx),
      mode_desc
    ))
    message("Results stored in @misc$velocity_gene_stats.")
  }

  return(pv)
}

.normalize_driver_missing <- function(stats) {
  for (key in names(stats)) {
    x <- stats[[key]]
    if (!is.numeric(x)) next
    missing <- is.nan(x)
    if (any(missing)) {
      x[missing] <- NA_real_
      stats[[key]] <- x
    }
  }
  stats
}


#' @title Get Velocity Gene Statistics
#' @description
#' Retrieve the ranked gene statistics table produced by
#' \code{rank_velocity_genes()}.
#'
#' @param pv A \code{plantvelo} object.
#' @param include_filtered Logical. If \code{TRUE}, return the complete
#'   \code{velocity_gene_stats} table with all genes, evidence and filter flags.
#'   Default \code{FALSE} returns selected genes ordered by rank.
#'
#' @details
#' A stored ranking without the merged selection fields and ranking parameters
#' produces a warning. Re-run \code{rank_velocity_genes()} to update an older
#' saved object; retrieving its table does not retroactively apply filters.
#' Legacy separate QC tables can be read but are not recreated or stored.
#' Numeric \code{NaN} values in saved tables are returned as \code{NA} without
#' modifying the object or recomputing selection.
#'
#' @return A data frame, or \code{NULL} if
#'   \code{rank_velocity_genes()} has not been run.
#' @examples
#' \dontrun{
#' stats <- get_velocity_gene_stats(pv)
#' }
#' @export
get_velocity_gene_stats <- function(pv, include_filtered = FALSE) {
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (!is.logical(include_filtered) || length(include_filtered) != 1L || is.na(include_filtered))
    stop("`include_filtered` must be TRUE or FALSE.")
  stats <- pv@misc$velocity_gene_stats
  merged <- !is.null(stats) && all(c("selected","rank","pass_fit","pass_velocity",
    "pass_raw","pass_likelihood") %in% names(stats)) &&
    is.list(pv@misc$velocity_params) && is.list(pv@misc$velocity_params$ranking)
  if (!merged && include_filtered && !is.null(pv@misc$velocity_gene_qc))
    stats <- pv@misc$velocity_gene_qc
  if (is.null(stats)) {
    message("No gene stats found. Run rank_velocity_genes() first.")
    return(NULL)
  }
  if (!merged)
    warning("Stored driver ranking uses an older schema. Re-run rank_velocity_genes().",
            call. = FALSE)
  stats$exclusion_reason <- NULL
  if (!include_filtered && "selected" %in% names(stats)) {
    stats <- stats[stats$selected %in% TRUE,,drop=FALSE]
    stats <- stats[order(stats$rank,stats$gene,na.last=TRUE),,drop=FALSE]
  }
  rownames(stats) <- NULL
  .normalize_driver_missing(stats)
}


# =============================================================================
# Visualization
# =============================================================================

#' @title Plot Driver Genes
#' @description
#' Plot a positive horizontal lollipop chart of selected driver genes, sorted
#' by \code{coherence} or \code{fit_likelihood} before taking \code{n_top}.
#' The x axis remains coherence; point colour represents fit likelihood.
#'
#' @param pv A \code{plantvelo} object.
#' @param n_top Integer. Number of top genes shown.
#' @param title Character scalar or \code{NULL}. Plot title.
#' @param fill_quantile Numeric vector of length two defining likelihood colour
#'   limits when \code{use_quantile = TRUE}.
#' @param use_quantile Logical. Use quantile-based colour limits.
#' @param segment_size Numeric. Lollipop stem width.
#' @param point_size Numeric. Lollipop head size.
#' @param sort_by Character. Sort descending by \code{"coherence"} (default)
#'   or \code{"fit_likelihood"}. Ties use ascending gene names; missing and
#'   non-finite sorting values come last. Only already selected driver genes
#'   are considered; stored selection flags and ranks are not modified.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' pv <- rank_velocity_genes(pv)
#' plot_driver_genes(pv)
#' plot_driver_genes(pv, sort_by = "fit_likelihood")
#' }
#' @export
plot_driver_genes <- function(pv,
                              n_top           = 20L,
                              title           = NULL,
                              fill_quantile   = c(0.05, 0.95),
                              use_quantile    = TRUE,
                              segment_size    = 0.5,
                              point_size      = 5,
                              sort_by         = c("coherence", "fit_likelihood")) {

  sort_by <- match.arg(sort_by)

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

  sort_values <- df[[sort_by]]
  sort_values[!is.finite(sort_values)] <- NA_real_
  df <- df[order(-sort_values,as.character(df$gene),na.last=TRUE), , drop = FALSE]
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
