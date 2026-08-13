# =============================================================================
# PlantVelocity: latent_time.R
# Shared latent-time inference across multiple genes
# Corresponds to scVelo: tools/_em_model_core.py::latent_time()
# =============================================================================

#' @title Compute Shared Latent Time
#' @description
#' Infer a single shared latent-time axis for each cell by integrating
#' per-gene fitted times from \code{recover_dynamics()}.
#'
#' The algorithm:
#' \enumerate{
#'   \item Filter low-quality genes (\code{fit_likelihood < min_likelihood}).
#'   \item Identify root cell(s) from \code{@meta.data} or \code{@misc}.
#'   \item Re-orient each gene's time axis relative to the root using
#'     \code{root_time()}.
#'   \item Fuse per-gene times into a consensus axis via
#'     \code{compute_shared_time()}.
#'   \item Correct low-confidence cells using KNN smoothing.
#'   \item Normalise the final latent-time vector to \eqn{[0, 1]}.
#' }
#'
#' The result is stored in \code{@meta.data$latent_time}.
#'
#' @param pv A \code{plantvelo} object with \code{@layers$fit_t} and
#'   \code{@kinetics$params} populated by \code{recover_dynamics()}.
#' @param vkey Character scalar. Velocity layer key, used only for reading
#'   \code{@misc[[paste0(vkey, "_params")]]}. Default \code{"velocity"}.
#' @param min_likelihood Numeric. Genes with \code{fit_likelihood} below this
#'   threshold are excluded from latent-time fusion. Default \code{0.1}.
#' @param min_confidence Numeric. Cells whose local consistency score falls
#'   below this threshold are replaced by KNN-smoothed latent time. Default
#'   \code{0.75}.
#' @param root_key Character scalar or \code{NULL}. Column name in
#'   \code{@meta.data} (probability vector) or key in \code{@misc}
#'   (cell index) specifying root cells. \code{NULL} triggers auto-detection.
#'   Default \code{NULL}.
#' @param end_key Character scalar or \code{NULL}. Column name in
#'   \code{@meta.data} containing terminal-state probabilities, used as a weak
#'   time prior. Default \code{NULL}.
#' @param t_max Numeric or \code{NULL}. If provided, rescale the output from
#'   \eqn{[0, 1]} to \eqn{[0, t\_max]}. Default \code{NULL}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{plantvelo} object with \code{@meta.data$latent_time}
#'   added as a numeric vector of length \eqn{n\_cells}.
#'
#' @examples
#' \dontrun{
#' pv <- compute_latent_time(pv)
#' }
#' @export
compute_latent_time <- function(pv,
                                vkey           = "velocity",
                                min_likelihood = 0.1,
                                min_confidence = 0.75,
                                root_key       = NULL,
                                end_key        = NULL,
                                t_max          = NULL,
                                verbose        = TRUE) {

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (is.null(pv@layers[["fit_t"]]))
    stop(paste(
      "latent_time: 'fit_t' layer not found.",
      "Please run recover_dynamics() first."
    ))

  # ----------------read and transpose fitted time matrix----------------
  # @layers$fit_t is stored as genes × cells, whereas root_time() and
  # compute_shared_time() expect cells × genes.
  T_mat   <- t(make_dense_matrix(pv@layers[["fit_t"]]))
  n_cells <- nrow(T_mat)
  n_genes <- ncol(T_mat)

  # Gene names aligned to columns of T_mat
  gene_names <- colnames(T_mat)

  # ----------------filter low-quality genes----------------

  idx_valid <- !is.na(colSums(T_mat))

  params <- pv@kinetics[["params"]]
  if (!is.null(params) &&
      "fit_likelihood" %in% colnames(params) &&
      !is.null(min_likelihood)) {
    ll <- params[gene_names, "fit_likelihood"]
    ll[is.na(ll)] <- 0
    idx_valid <- idx_valid & (ll >= min_likelihood)
  }

  if (sum(idx_valid) < 2) {
    warning("latent_time: fewer than 2 valid genes after filtering. Using all genes.")
    idx_valid <- !is.na(colSums(T_mat))
  }

  T_filt <- T_mat[, idx_valid, drop = FALSE]

  if (verbose)
    message(sprintf(
      "Using %d / %d genes for latent-time fusion (min_likelihood = %.2f).",
      sum(idx_valid), n_genes, min_likelihood
    ))

  # ----------------identify root cell(s)----------------

  root_cell <- .find_root_cell(pv, root_key, T_mat)

  # ----------------re-orient time axes and compute consensus latent time----------------

  if (!is.null(root_cell) && length(root_cell) == 1) {
    result   <- root_time(T_filt, root = root_cell)
    latent_t <- compute_shared_time(result$t_rooted)

  } else if (!is.null(root_cell) && length(root_cell) <= 4) {
    # Multiple roots: compute one latent-time estimate per root and average
    lat_list <- lapply(
      root_cell[seq_len(min(4L, length(root_cell)))],
      function(r) {
        result <- root_time(T_filt, root = r)
        compute_shared_time(result$t_rooted)
      }
    )
    latent_t <- scale_to_01(rowMeans(do.call(cbind, lat_list), na.rm = TRUE))

  } else {
    # No root information available: centre each gene by its minimum fitted time
    # and use the row means as a fallback latent-time estimate.
    if (verbose)
      message("  No root cell found. Using row means of centered fit_t.")
    T_centered <- sweep(T_filt, 2, apply(T_filt, 2, min, na.rm = TRUE), "-")
    latent_t   <- scale_to_01(rowMeans(T_centered, na.rm = TRUE))
  }

  # ----------------apply optional end-point constraint----------------
  # Terminal-state probabilities are used as a weak prior to favour later times
  # for likely fate cells.

  if (!is.null(end_key) && !is.null(pv@meta.data[[end_key]])) {
    fate_probs  <- pv@meta.data[[end_key]]
    fate_probs[is.na(fate_probs)] <- 0
    fate_cells  <- which(fate_probs > 0.5)

    if (length(fate_cells) > 0) {
      fate_list <- lapply(
        fate_cells[seq_len(min(4L, length(fate_cells)))],
        function(f) {
          result <- root_time(T_filt, root = f)
          1 - compute_shared_time(result$t_rooted)
        }
      )
      latent_fate <- scale_to_01(
        rowMeans(do.call(cbind, fate_list), na.rm = TRUE)
      )
      latent_t <- scale_to_01(latent_t + 0.2 * latent_fate)
    }
  }

  # ----------------local consistency correction----------------

  if (!is.null(pv@graphs[["neighbors"]]) && !is.null(min_confidence)) {
    conn <- pv@graphs$neighbors$connectivities
    C    <- Matrix::Matrix(conn, sparse = TRUE)

    # Row-normalised connectivity matrix
    rs <- Matrix::rowSums(C)
    rs[rs == 0] <- 1
    C_norm <- C / rs

    tc <- as.vector(C_norm %*% latent_t)

    # Consistency score: deviation between a cell's latent time and the
    # neighbourhood-smoothed estimate
    max_lt <- max(latent_t, na.rm = TRUE)
    if (max_lt > 0) {
      z      <- sum(latent_t * tc, na.rm = TRUE) /
        max(sum(tc^2, na.rm = TRUE), 1e-10)
      conf   <- (1 - abs(latent_t / max_lt - tc * z / max_lt))^2
      low_conf <- !is.na(conf) & (conf < min_confidence)

      if (any(low_conf)) {
        if (verbose)
          message(sprintf(
            "  Smoothing %d low-confidence cells (min_confidence = %.2f).",
            sum(low_conf), min_confidence
          ))
        # Mask low-confidence source columns and replace by neighbourhood
        # smoothing
        C_masked <- C
        C_masked[, low_conf] <- 0
        rs_m <- Matrix::rowSums(C_masked)
        rs_m[rs_m == 0] <- 1
        C_masked <- C_masked / rs_m
        latent_t <- as.vector(C_masked %*% latent_t)
      }
    }
  }

  # ----------------final normalisation and write-out----------------

  latent_t <- scale_to_01(latent_t)
  if (!is.null(t_max)) latent_t <- latent_t * t_max

  pv@meta.data[["latent_time"]] <- latent_t

  if (verbose)
    message("Added 'latent_time' to @meta.data.")

  return(pv)
}


# =============================================================================
# Internal helper functions
# =============================================================================

#' Identify root-cell indices
#'
#' Search order:
#' \enumerate{
#'   \item User-specified \code{root_key} in \code{@misc} containing integer
#'     cell indices.
#'   \item User-specified \code{root_key} in \code{@meta.data} containing a
#'     root-probability vector.
#'   \item Auto-detection from common root-related fields in
#'     \code{@meta.data}.
#'   \item Final fallback: the cell with the minimum row sum in \code{fit_t}.
#' }
#'
#' @param pv A \code{plantvelo} object.
#' @param root_key Character scalar or \code{NULL}. User-specified root key.
#' @param T_mat Fitted-time matrix of size \eqn{n\_cells \times n\_genes}.
#'
#' @return Integer vector of root-cell indices, or \code{NULL}.
#' @keywords internal
.find_root_cell <- function(pv, root_key, T_mat) {

  # Helper: extract top root-cell indices from a probability vector
  .from_prob_vec <- function(prob_vec) {
    prob_vec[is.na(prob_vec)] <- 0
    if (max(prob_vec, na.rm = TRUE) <= 0) return(NULL)
    top <- order(prob_vec, decreasing = TRUE)[
      seq_len(min(4L, length(prob_vec)))
    ]
    top[prob_vec[top] >= 0.9 * max(prob_vec, na.rm = TRUE)]
  }

  # ----------------user-specified integer index in @misc----------------
  if (!is.null(root_key) && !is.null(pv@misc[[root_key]])) {
    val <- pv@misc[[root_key]]
    if (is.numeric(val)) return(as.integer(val))
  }

  # ----------------user-specified probability vector in @meta.data----------------
  if (!is.null(root_key) && !is.null(pv@meta.data[[root_key]])) {
    result <- .from_prob_vec(pv@meta.data[[root_key]])
    if (!is.null(result)) return(result)
  }

  # ----------------auto-detect common root-related fields----------------
  for (key in c("root_cells", "starting_cells", "root_states_probs")) {
    if (!is.null(pv@meta.data[[key]])) {
      result <- .from_prob_vec(pv@meta.data[[key]])
      if (!is.null(result)) return(result)
    }
  }

  # ----------------minimum row sum in fit_t----------------
  if (!is.null(T_mat)) {
    t_sum <- rowSums(T_mat, na.rm = TRUE)
    return(which.min(t_sum))
  }

  NULL
}
