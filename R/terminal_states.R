# =============================================================================
# PlantVelocity: terminal_states.R
# Terminal-state identification by Markov-chain eigendecomposition
# Corresponds to scVelo: tools/terminal_states.py
# (model-agnostic; reused without conceptual changes)
#
# Mathematical framework:
#   End points: left eigenvectors of T (lambda approximately 1)
#               -> absorbing states where cells tend to remain
#   Root cells: left eigenvectors of backward T
#               -> repelling states from which cells tend to originate
# =============================================================================

#' @title Compute Terminal States (Root Cells and End Points)
#' @description
#' Identify trajectory root cells and end points by eigendecomposition of the
#' velocity transition matrix.
#'
#' \itemize{
#'   \item \strong{End points} (absorbing states): identified from the leading
#'     eigenvectors of the forward transition matrix
#'     (\eqn{\lambda \approx 1}).
#'   \item \strong{Root cells} (repelling states): identified from the leading
#'     eigenvectors of the backward transition matrix.
#' }
#'
#' Results are continuous probability scores in \eqn{[0, 1]}.
#'
#' @section Output written to \code{@meta.data}:
#' \describe{
#'   \item{\code{root_cells}}{Numeric vector. Root-cell probability per cell.}
#'   \item{\code{end_points}}{Numeric vector. End-point probability per cell.}
#' }
#'
#' @param pv A \code{plantvelo} object with velocity graph computed by
#'   \code{compute_velocity_graph()}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#' @param modality Character scalar. Expression layer used in
#'   \code{verify_roots()}. Default \code{"Ms"}.
#' @param groupby Character scalar or \code{NULL}. Column in
#'   \code{@meta.data} used for group-wise analysis. Each group receives its own
#'   within-group normalisation, while eigendecomposition is still performed on
#'   the full graph. \code{NULL} analyses all cells jointly. Default
#'   \code{NULL}.
#' @param groups Character vector or \code{NULL}. Subset of group labels when
#'   \code{groupby} is set. \code{NULL} uses all groups. Default
#'   \code{NULL}.
#' @param self_transitions Logical. Add self-loops before softmax
#'   normalisation. Default \code{FALSE}.
#' @param eps Numeric. Eigenvalue selection threshold
#'   (\eqn{\lambda \geq 1 - \mathrm{eps}}). Default \code{0.001}.
#' @param random_state Integer. Random seed. Default \code{0L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{plantvelo} object with \code{@meta.data$root_cells} and
#'   \code{@meta.data$end_points} added.
#'
#' @examples
#' \dontrun{
#' pv <- compute_terminal_states(pv)
#' }
#' @export
compute_terminal_states <- function(pv,
                                    vkey             = "velocity",
                                    modality         = "Ms",
                                    groupby          = NULL,
                                    groups           = NULL,
                                    self_transitions = FALSE,
                                    eps              = 1e-3,
                                    random_state     = 0L,
                                    verbose          = TRUE) {

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  gkey     <- paste0(vkey, "_graph")
  gkey_neg <- paste0(vkey, "_graph_neg")

  if (is.null(pv@graphs[[gkey]]))
    stop(sprintf(
      "Velocity graph '%s' not found. Run compute_velocity_graph() first.", gkey
    ))

  if (is.null(pv@graphs[["neighbors"]]))
    stop("No neighbor graph found. Run build_neighbor_graph() first.")

  # ----------------read global graph data----------------

  n_cells    <- ncol(pv@layers[["spliced"]])
  dist_mat   <- pv@graphs$neighbors$distances
  T_fwd_full <- Matrix::Matrix(pv@graphs[[gkey]],     sparse = TRUE)
  T_neg_full <- Matrix::Matrix(pv@graphs[[gkey_neg]], sparse = TRUE)

  # Initialise output vectors with zeros
  roots_out <- rep(0, n_cells)
  ends_out  <- rep(0, n_cells)

  # ----------------build full transition matrices and run eigendecomposition----------------
  # Eigendecomposition must be performed on the full graph rather than on
  # per-group subgraphs:
  #   - Subgraph restriction removes cross-group edges and may produce a
  #     substochastic matrix with spectral radius < 1.
  #   - groupby affects only within-group normalisation, not the eigensystem.
  #
  # T_neg stores negative edges. Using scale = -10 amplifies strong reverse
  # edges:
  #   exp(negative * -10) = exp(positive) >> 1
  # which concentrates the stationary distribution on true root cells.
  # If scale = +10 were used instead, exp(negative * 10) would shrink toward
  # zero and the distribution would become nearly uniform.

  T_back <- .softmax_transition(T_neg_full, scale = -10,
                                self_transitions = self_transitions)

  # Repair zero rows in T_back:
  # velocity_graph_neg stores only pi < 0 edges, so cells without reverse
  # neighbours have row sum 0, producing a substochastic matrix and
  # spectral radius < 1. Add self-loops only to zero rows so that such cells
  # become absorbing states while all non-zero rows remain unchanged.
  zero_rows_back <- which(Matrix::rowSums(T_back) == 0)
  if (length(zero_rows_back) > 0) {
    T_back <- T_back + Matrix::sparseMatrix(
      i    = zero_rows_back,
      j    = zero_rows_back,
      x    = rep(1, length(zero_rows_back)),
      dims = c(n_cells, n_cells)
    )
  }

  T_fwd  <- .softmax_transition(T_fwd_full, scale = 10,
                                self_transitions = self_transitions)

  eig_roots <- compute_eigs(T_back, eps = eps, perc = c(2, 98),
                            random_state = random_state)
  eig_ends  <- compute_eigs(T_fwd,  eps = eps, perc = c(2, 98),
                            random_state = random_state)

  n_roots <- ncol(eig_roots$vectors)
  n_ends  <- ncol(eig_ends$vectors)

  if (verbose)
    message(sprintf("  %d root region(s), %d end-point region(s).", n_roots, n_ends))

  # No valid eigenvectors found on either side
  if (n_roots == 0 && n_ends == 0) {
    warning(paste0(
      "compute_terminal_states: eigendecomposition returned no valid eigenvectors ",
      "for either root cells or end points.\n",
      "All @meta.data$root_cells and @meta.data$end_points will be set to 0.\n",
      "Suggestions:\n",
      "  1. Increase eps (for example eps = 0.01 or 0.05) to relax the eigenvalue threshold.\n",
      "  2. Check that compute_velocity_graph() produced a non-trivial velocity graph.\n",
      "  3. Verify that the velocity layer is not all-zero or NA."
    ))
  } else if (n_roots == 0) {
    warning(paste0(
      "compute_terminal_states: no root-cell eigenvectors found. ",
      "@meta.data$root_cells will be all zeros.\n",
      "Try increasing eps or checking the backward transition matrix."
    ))
  } else if (n_ends == 0) {
    warning(paste0(
      "compute_terminal_states: no end-point eigenvectors found. ",
      "@meta.data$end_points will be all zeros.\n",
      "Try increasing eps or checking the forward transition matrix."
    ))
  }

  # ----------------resolve analysis groups----------------

  categories <- if (!is.null(groupby)) {
    all_cats <- levels(as.factor(pv@meta.data[[groupby]]))
    if (!is.null(groups)) intersect(groups, all_cats) else all_cats
  } else {
    list(NULL)
  }

  # ----------------extract group-specific rows, smooth/normalise within group, and write back----------------

  for (cat in categories) {

    if (!is.null(cat)) {
      cell_mask <- pv@meta.data[[groupby]] == cat
      cell_idx  <- which(cell_mask)
    } else {
      cell_mask <- NULL
      cell_idx  <- seq_len(n_cells)
    }

    n_sub <- length(cell_idx)

    # Group-specific connectivity submatrix used to smooth eigenvectors
    conn_sub <- dist_mat[cell_idx, cell_idx, drop = FALSE]

    # ----------------root cells----------------
    if (ncol(eig_roots$vectors) > 0) {
      root_vecs <- eig_roots$vectors[cell_idx, , drop = FALSE]
      roots_raw <- as.vector(Matrix::rowSums(conn_sub %*% root_vecs))
      roots     <- scale_to_01(pmax(roots_raw, 0))
      clip_ub   <- quantile(roots, 0.98, na.rm = TRUE)
      roots     <- scale_to_01(pmin(roots, clip_ub))
      roots     <- verify_roots(pv, roots, cell_idx = cell_idx,
                                modality = modality)
    } else {
      roots <- rep(0, n_sub)
    }

    # ----------------end points----------------
    if (ncol(eig_ends$vectors) > 0) {
      end_vecs <- eig_ends$vectors[cell_idx, , drop = FALSE]
      ends_raw <- as.vector(Matrix::rowSums(conn_sub %*% end_vecs))
      ends     <- scale_to_01(pmax(ends_raw, 0))
      clip_ub  <- quantile(ends, 0.98, na.rm = TRUE)
      ends     <- scale_to_01(pmin(ends, clip_ub))
    } else {
      ends <- rep(0, n_sub)
    }

    # ----------------write back into global output vectors----------------
    roots_out[cell_idx] <- roots
    ends_out[cell_idx]  <- ends

    if (verbose && !is.null(cat))
      message(sprintf("  Processed group '%s' (%d cells).", cat, n_sub))
  }

  # ----------------store results in @meta.data----------------

  pv@meta.data[["root_cells"]] <- roots_out
  pv@meta.data[["end_points"]] <- ends_out

  if (verbose)
    message("Added 'root_cells' and 'end_points' to @meta.data.")

  return(pv)
}


#' @title Eigendecomposition of Transition Matrix
#' @description
#' Compute leading eigenvectors of a row-normalised Markov transition matrix
#' \eqn{T}. Stationary-like states are identified by selecting eigenvectors
#' corresponding to eigenvalues satisfying
#' \eqn{\lambda \geq 1 - \varepsilon}.
#'
#' Uses \pkg{RSpectra} for sparse matrices when available and falls back to
#' base \code{eigen()} otherwise.
#'
#' @param T Sparse or dense transition matrix of size
#'   \eqn{n\_cells \times n\_cells}.
#' @param k Integer. Number of leading eigenvalues to compute. Default
#'   \code{10L}.
#' @param eps Numeric. Eigenvalue selection threshold. Default
#'   \code{1e-3}.
#' @param perc Numeric vector of length 2. Percentile clipping applied to each
#'   eigenvector to reduce noise. Default \code{c(2, 98)}.
#' @param random_state Integer. Random seed. Default \code{0L}.
#'
#' @return Named list with components:
#' \describe{
#'   \item{\code{values}}{Numeric vector of selected eigenvalues.}
#'   \item{\code{vectors}}{Matrix of size \eqn{n\_cells \times n\_selected}.
#'     Eigenvectors are converted to absolute values and each column is scaled
#'     to \eqn{[0, 1]}.}
#' }
#' @export
compute_eigs <- function(T, k = 10L, eps = 1e-3,
                         perc = c(2, 98), random_state = 0L) {
  set.seed(random_state)
  n <- nrow(T)
  k <- min(as.integer(k), n - 2L)

  if (k < 1L)
    return(list(values = numeric(0),
                vectors = matrix(0, n, 0)))

  # Left eigenvectors of T correspond to right eigenvectors of T^T
  T_t <- Matrix::t(T)

  eigvals <- numeric(0)
  eigvecs <- matrix(numeric(0), n, 0)

  tryCatch({
    if (requireNamespace("RSpectra", quietly = TRUE)) {
      res     <- RSpectra::eigs(T_t, k = k, which = "LR",
                                opts = list(tol = 1e-6, maxitr = 1000))
      eigvals <- Re(res$values)
      eigvecs <- Re(res$vectors)
    } else {
      T_dense <- as.matrix(T_t)
      res     <- eigen(T_dense)
      eigvals <- Re(res$values)[seq_len(min(k, n))]
      eigvecs <- Re(res$vectors)[, seq_len(min(k, n)), drop = FALSE]
    }
  }, error = function(e) {
    warning(sprintf(
      paste0(
        "compute_eigs: eigendecomposition failed (%s).\n",
        "Returning empty eigenvectors; terminal states will be set to 0.\n",
        "If using RSpectra, try unloading it to fall back to base::eigen(), ",
        "or reduce k."
      ),
      e$message
    ))
  })

  if (length(eigvals) == 0)
    return(list(values = numeric(0), vectors = matrix(0, n, 0)))

  # Sort in decreasing order
  ord     <- order(eigvals, decreasing = TRUE)
  eigvals <- eigvals[ord]
  eigvecs <- eigvecs[, ord, drop = FALSE]

  # Keep eigenvalues satisfying lambda >= 1 - eps
  idx     <- eigvals >= 1 - eps
  eigvals <- eigvals[idx]
  eigvecs <- abs(eigvecs[, idx, drop = FALSE])

  # Percentile clipping and column scaling to [0, 1]
  if (!is.null(perc) && ncol(eigvecs) > 0) {
    for (j in seq_len(ncol(eigvecs))) {
      lb <- quantile(eigvecs[, j], perc[1] / 100, na.rm = TRUE)
      ub <- quantile(eigvecs[, j], perc[2] / 100, na.rm = TRUE)
      eigvecs[eigvecs[, j] < lb, j] <- 0
      eigvecs[, j] <- pmin(eigvecs[, j], ub)
      mx <- max(eigvecs[, j], na.rm = TRUE)
      if (mx > 0) eigvecs[, j] <- eigvecs[, j] / mx
    }
  }

  list(values = eigvals, vectors = eigvecs)
}


#' @title Verify Root Cell Identification
#' @description
#' Perform a sanity check on the root-cell probability vector.
#'
#' Two strategies are supported:
#' \enumerate{
#'   \item If \code{time_key} is provided, test whether root cells are enriched
#'     at the earliest observed time point. This is useful for time-resolved
#'     plant experiments.
#'   \item Otherwise, use a plasticity score based on gene-count correlation
#'     (scVelo-style). If \code{@kinetics$params} lacks
#'     \code{gene_count_corr}, gene means are used as a fallback proxy.
#' }
#'
#' @param pv A \code{plantvelo} object.
#' @param roots Numeric vector. Root-cell probabilities for the cells indexed
#'   by \code{cell_idx}.
#' @param cell_idx Integer vector. Row indices in \code{@meta.data} that
#'   \code{roots} corresponds to. Default \code{NULL}, meaning all cells.
#' @param modality Character scalar. Smoothed expression key in
#'   \code{@moments} used for plasticity scoring. Default \code{"Ms"}.
#' @param time_key Character scalar or \code{NULL}. Column in
#'   \code{@meta.data} containing sampling or treatment time. Default
#'   \code{NULL}.
#'
#' @return Numeric vector of the same length as \code{roots}, possibly
#'   modified after verification.
#' @export
verify_roots <- function(pv, roots,
                         cell_idx = NULL,
                         modality = "Ms",
                         time_key = NULL) {

  if (is.null(cell_idx)) cell_idx <- seq_len(nrow(pv@meta.data))

  # ----------------strategy A: earliest-time enrichment check----------------
  if (!is.null(time_key) && !is.null(pv@meta.data[[time_key]])) {
    t_obs  <- pv@meta.data[[time_key]][cell_idx]
    t_min  <- min(t_obs, na.rm = TRUE)
    early  <- !is.na(t_obs) & (t_obs == t_min)
    enrich <- mean(roots[early] > 0.5, na.rm = TRUE)

    if (enrich < 0.2)
      message(paste(
        "verify_roots: root cells are not enriched at the earliest time point.",
        "Consider setting root_key manually in compute_latent_time()."
      ))
    return(roots)
  }

  # ----------------strategy B: plasticity score check----------------
  X_full <- if (!is.null(pv@moments[[modality]])) {
    make_dense_matrix(pv@moments[[modality]])
  } else if (!is.null(pv@layers[["spliced"]])) {
    t(make_dense_matrix(pv@layers[["spliced"]]))
  } else {
    NULL
  }

  if (is.null(X_full)) return(roots)

  X <- X_full[cell_idx, , drop = FALSE]

  # Use gene_count_corr if available; otherwise fall back to gene means
  if (!is.null(pv@kinetics[["params"]]) &&
      "gene_count_corr" %in% colnames(pv@kinetics$params)) {
    gcc <- pv@kinetics$params[["gene_count_corr"]]
  } else {
    gcc <- colMeans(X_full, na.rm = TRUE)
  }

  gcc[is.na(gcc)] <- 0
  if (sum(abs(gcc)) == 0) return(roots)

  plasticity <- as.vector(X %*% gcc)
  plasticity <- scale_to_01(plasticity)

  p_ub    <- plasticity > 0.5
  root_ub <- roots > 0.9

  if (sum(p_ub) == 0 || sum(!p_ub) == 0) return(roots)

  n_right  <- sum(root_ub & p_ub)  / sum(p_ub)
  n_false  <- sum(root_ub & !p_ub) / sum(!p_ub)
  n_random <- mean(root_ub, na.rm = TRUE)

  if (n_right > 3 * n_random) {
    roots <- roots * as.numeric(p_ub)
  } else if (n_false > n_random || n_right < n_random) {
    message(paste(
      "verify_roots: uncertain root-cell identification.",
      "Please verify manually with plot_velocity_embedding()."
    ))
  }

  roots
}


#' @title Compute Cell Fate Assignment
#' @description
#' Assign each cell to the most likely terminal-state group using the
#' fundamental matrix \eqn{N = (I - T)^{-1}}.
#'
#' \strong{Warning}: matrix inversion is \eqn{O(n^3)} and is only recommended
#' for datasets with fewer than about 3,000 cells.
#'
#' @param pv A \code{plantvelo} object.
#' @param groupby Character scalar. Column in \code{@meta.data} containing
#'   cell-state or cluster labels. Default \code{"clusters"}.
#' @param disconnected_groups Character vector or \code{NULL}. Groups excluded
#'   from fate computation. Cells in these groups retain their own label.
#'   Default \code{NULL}.
#' @param self_transitions Logical. Add self-loops before transition
#'   normalisation. Default \code{FALSE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{plantvelo} object with \code{@meta.data$cell_fate} and
#'   \code{@meta.data$cell_fate_confidence} added.
#' @export
compute_cell_fate <- function(pv,
                              groupby             = "clusters",
                              disconnected_groups = NULL,
                              self_transitions    = FALSE,
                              verbose             = TRUE) {

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  n_obs <- ncol(pv@layers[["spliced"]])
  if (n_obs > 5000 && verbose)
    message(paste(
      "compute_cell_fate: n_cells =", n_obs,
      "> 5,000. Matrix inversion is O(n^3) and may be slow.",
      "Consider compute_terminal_states() for large datasets."
    ))

  T_mat   <- build_transition_matrix(pv, self_transitions = self_transitions)
  T_dense <- as.matrix(T_mat)

  fate <- tryCatch(
    solve(diag(n_obs) - T_dense),
    error = function(e) {
      if (verbose)
        message("compute_cell_fate: matrix inversion failed; using power iteration.")
      .power_iteration_fundamental(T_dense, n_steps = 20L)
    }
  )

  group_labels <- pv@meta.data[[groupby]]
  if (is.null(group_labels))
    stop(sprintf("groupby key '%s' not found in @meta.data.", groupby))

  cell_fates <- group_labels[apply(fate, 1, which.max)]

  if (!is.null(disconnected_groups)) {
    idx <- group_labels %in% disconnected_groups
    cell_fates[idx] <- group_labels[idx]
  }

  row_sums <- rowSums(fate)
  row_sums[row_sums == 0] <- 1

  pv@meta.data[["cell_fate"]]            <- cell_fates
  pv@meta.data[["cell_fate_confidence"]] <- apply(fate, 1, max) / row_sums

  if (verbose)
    message("Added 'cell_fate' and 'cell_fate_confidence' to @meta.data.")

  return(pv)
}


#' @title Compute Cell Origin Assignment
#' @description
#' Identify the most likely origin state for each cell using the backward
#' transition matrix, analogous to a reverse-direction version of
#' \code{compute_cell_fate()}.
#'
#' @inheritParams compute_cell_fate
#'
#' @return A \code{plantvelo} object with \code{@meta.data$cell_origin} and
#'   \code{@meta.data$cell_origin_confidence} added.
#' @export
compute_cell_origin <- function(pv,
                                groupby             = "clusters",
                                disconnected_groups = NULL,
                                self_transitions    = FALSE,
                                verbose             = TRUE) {

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  n_obs   <- ncol(pv@layers[["spliced"]])
  T_back  <- build_transition_matrix(pv, self_transitions = self_transitions,
                                     backward = TRUE)
  T_dense <- as.matrix(T_back)

  origin <- tryCatch(
    solve(diag(n_obs) - T_dense),
    error = function(e) {
      if (verbose)
        message("compute_cell_origin: matrix inversion failed; using power iteration.")
      .power_iteration_fundamental(T_dense, n_steps = 20L)
    }
  )

  group_labels <- pv@meta.data[[groupby]]
  if (is.null(group_labels))
    stop(sprintf("groupby key '%s' not found in @meta.data.", groupby))

  cell_origins <- group_labels[apply(origin, 1, which.max)]

  if (!is.null(disconnected_groups)) {
    idx <- group_labels %in% disconnected_groups
    cell_origins[idx] <- group_labels[idx]
  }

  row_sums <- rowSums(origin)
  row_sums[row_sums == 0] <- 1

  pv@meta.data[["cell_origin"]]            <- cell_origins
  pv@meta.data[["cell_origin_confidence"]] <- apply(origin, 1, max) / row_sums

  if (verbose)
    message("Added 'cell_origin' and 'cell_origin_confidence' to @meta.data.")

  return(pv)
}


# =============================================================================
# Internal helper functions
# =============================================================================

#' Approximate the fundamental matrix by power iteration
#'
#' Fallback used when direct matrix inversion fails for large matrices:
#' \deqn{N \approx I + T + T^2 + \cdots + T^{n\_steps}}
#'
#' @param T Dense transition matrix.
#' @param n_steps Integer. Number of power-iteration terms. Default
#'   \code{20L}.
#'
#' @return Dense matrix approximating the fundamental matrix.
#' @keywords internal
.power_iteration_fundamental <- function(T, n_steps = 20L) {
  n   <- nrow(T)
  mat <- diag(n)
  Tk  <- diag(n)
  for (i in seq_len(n_steps)) {
    Tk  <- Tk %*% T
    mat <- mat + Tk
  }
  mat
}
