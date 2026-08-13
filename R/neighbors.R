# =============================================================================
# PlantVelocity: neighbor_graph.R
# Build the cell-cell KNN neighbor graph used for local smoothing
# and downstream RNA velocity analysis.
#
# Mathematical framework:
#   1. Compute a K-nearest-neighbor graph from low-dimensional cell
#      representations (typically PCA embeddings)
#   2. Store raw Euclidean neighbor distances as a sparse matrix
#   3. Convert distances into Gaussian-kernel connectivities
#   4. Symmetrise the connectivity matrix and use it as the basis for
#      moment smoothing, velocity graph construction, and latent-time inference
#
# Output is written to @graphs$neighbors, including:
#   - connectivities : symmetric sparse cell-cell weight matrix
#   - distances      : sparse KNN distance matrix
#   - params         : graph construction settings
# =============================================================================

#' @title Build a KNN Neighbor Graph for a plantvelo Object
#' @description
#' Construct a K-nearest-neighbor (KNN) graph from cell embeddings stored in
#' \code{@reductions} (typically PCA). The resulting sparse connectivity and
#' distance matrices are written to \code{@graphs$neighbors} and serve as the
#' required input for \code{compute_moments()}.
#'
#' Two KNN backends are supported:
#' \describe{
#'   \item{\code{"base"}}{Pure-R implementation using \code{stats::dist} for
#'     n \eqn{\le} 5,000 cells, and a batched dot-product strategy for larger
#'     datasets. No extra packages required.}
#'   \item{\code{"bioc"}}{Uses \pkg{BiocNeighbors} (exact KNN via KD-tree /
#'     KMKNN). Recommended for n > 10,000 cells.}
#' }
#'
#' @section Representation selection:
#' \enumerate{
#'   \item If \code{use_pca = TRUE}, the function looks for a reduction named
#'     \code{"pca"} (case-insensitive) in \code{@reductions} and extracts
#'     \code{cell.embeddings}.
#'   \item If PCA is not found, it falls back to \code{@layers$spliced}
#'     (dense conversion).
#' }
#'
#' @section Output written to \code{@graphs$neighbors}:
#' \describe{
#'   \item{\code{connectivities}}{Symmetric sparse matrix (n_cells x n_cells).
#'     Gaussian-kernel weights computed from KNN distances, symmetrised by
#'     averaging.}
#'   \item{\code{distances}}{Sparse matrix (n_cells x n_cells). Raw Euclidean
#'     distances to each KNN.}
#'   \item{\code{params}}{Named list recording: \code{n_neighbors},
#'     \code{method}, \code{n_pcs}, \code{use_pca}.}
#' }
#'
#' @param pv A \code{plantvelo} object created by \code{create_plantvelo()}.
#' @param n_neighbors Integer. Number of nearest neighbors. Default \code{30L}.
#' @param n_pcs Integer. Number of PCA dimensions to use. Default \code{30L}.
#' @param method Character scalar. KNN backend: \code{"base"} (default) or
#'   \code{"bioc"} (requires \pkg{BiocNeighbors}).
#' @param random_state Integer. Random seed for reproducibility. Default \code{0L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{plantvelo} object with \code{@graphs$neighbors} populated.
#'
#' @examples
#' \dontrun{
#' pv <- build_neighbor_graph(pv, n_neighbors = 30)
#' }
#' @export
build_neighbor_graph <- function(pv,
                                  n_neighbors  = 30L,
                                  n_pcs        = 30L,
                                  method       = c("base", "bioc"),
                                  random_state = 0L,
                                  verbose      = TRUE) {

  method <- match.arg(method)

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (!is.numeric(n_neighbors) || length(n_neighbors) != 1 || n_neighbors < 2)
    stop("`n_neighbors` must be a single integer >= 2.")

  if (method == "bioc" && !requireNamespace("BiocNeighbors", quietly = TRUE))
    stop(
      "Package `BiocNeighbors` is required for method = 'bioc'.\n",
      "Install with: BiocManager::install('BiocNeighbors')"
    )

  n_neighbors <- as.integer(n_neighbors)
  n_pcs       <- as.integer(n_pcs)

  # ----------------prepare representation----------------

  if (verbose) message("[1] Extracting cell representation ...")
  set.seed(random_state)

  X_rep <- .get_representation(pv, n_pcs)
  n_obs <- nrow(X_rep)

  # Dimension assertion: nrow(X_rep) must equal n_cells.
  # @layers$spliced is stored as genes × cells, so n_cells = ncol(spliced).
  n_cells_expected <- ncol(pv@layers[["spliced"]])
  if (n_obs != n_cells_expected)
    stop(sprintf(
      "Dimension mismatch: X_rep has %d rows but @layers$spliced has %d cells (ncol). ",
      n_obs, n_cells_expected,
      "Expected nrow(X_rep) == n_cells. ",
      "Check matrix orientation in @reductions or @layers."
    ))

  n_neighbors <- min(n_neighbors, n_obs - 1L)

  if (n_obs > 10000L && method == "base")
    warning(
      "n_cells = ", n_obs, " > 10,000 with method = 'base' may be slow.\n",
      "Consider method = 'bioc': BiocManager::install('BiocNeighbors')"
    )

  if (verbose) message(sprintf(
    "    %d cells | %d dims | %d neighbors | method = '%s'",
    n_obs, ncol(X_rep), n_neighbors, method
  ))

  # ----------------build KNN graph----------------

  if (verbose) message("[2] Computing KNN graph ...")

  result <- if (method == "bioc") {
    .build_bioc_neighbors(X_rep, n_neighbors)
  } else {
    .build_base_neighbors(X_rep, n_neighbors)
  }

  # ----------------store in @graphs slot----------------

  pv@graphs[["neighbors"]] <- list(
    connectivities = result$connectivities,
    distances      = result$distances,
    params = list(
      n_neighbors = n_neighbors,
      method      = method,
      n_pcs       = n_pcs
    )
  )

  if (verbose) message(sprintf(
    "Neighbor graph stored in @graphs$neighbors (%d x %d sparse matrices).",
    n_obs, n_obs
  ))

  return(pv)
}

# =============================================================================
# Internal helper functions
# =============================================================================

#' Extract the cell representation matrix for KNN (PCA only)
#'
#' Looks up the PCA reduction in \code{@reductions} (case-insensitive) and
#' returns the cell embeddings truncated to \code{n_pcs} dimensions.  Stops
#' if PCA is not found.
#' @keywords internal
.get_representation <- function(pv, n_pcs) {
  pca_key <- .find_reduction_key("pca", names(pv@reductions))

  if (is.null(pca_key))
    stop(
      "PCA not found in @reductions. ",
      "Please run RunPCA() on your Seurat object before calling create_plantvelo(), ",
      "or re-run create_plantvelo() after adding PCA results."
    )

  red_obj <- pv@reductions[[pca_key]]
  X <- tryCatch(
    red_obj@cell.embeddings,
    error = function(e) as.matrix(red_obj)
  )

  if (!is.null(n_pcs) && ncol(X) > n_pcs)
    X <- X[, seq_len(n_pcs), drop = FALSE]

  return(X)
}

#' Case-insensitive match for a dimensionality reduction key
#'
#' Returns the first element of \code{available} whose lower-case value
#' matches \code{tolower(name)}, or \code{NULL} if none found.
#' @keywords internal
.find_reduction_key <- function(name, available) {
  hits <- available[tolower(available) == tolower(name)]
  if (length(hits) > 0) return(hits[1])
  NULL
}

#' Base R KNN implementation
#'
#' For n ≤ 5,000 cells uses \code{stats::dist}; for larger datasets uses a
#' batched dot-product expansion of \eqn{\|a - b\|^2 = \|a\|^2 + \|b\|^2 - 2a \cdot b}.
#' Returns a list with \code{distances} (sparse) and \code{connectivities}
#' (Gaussian-kernel, symmetrised).
#' @keywords internal
.build_base_neighbors <- function(X, n_neighbors) {
  n_obs       <- nrow(X)
  knn_dist    <- matrix(0.0, n_obs, n_neighbors)
  knn_indices <- matrix(0L,  n_obs, n_neighbors)

  if (n_obs <= 5000L) {
    # Small dataset: compute the full pairwise distance matrix directly
    D <- as.matrix(stats::dist(X, method = "euclidean"))
    for (i in seq_len(n_obs)) {
      d_i    <- D[i, ]
      d_i[i] <- Inf
      ord    <- order(d_i)[seq_len(n_neighbors)]
      knn_indices[i, ] <- ord
      knn_dist[i, ]    <- d_i[ord]
    }
  } else {
    # Large dataset: use batched dot-product expansion
    # ||a - b||² = ||a||² + ||b||² - 2(a·b)
    X_sq       <- rowSums(X^2)
    batch_size <- 500L
    for (i_start in seq(1L, n_obs, by = batch_size)) {
      i_end   <- min(i_start + batch_size - 1L, n_obs)
      i_batch <- seq(i_start, i_end)
      cross   <- X[i_batch, , drop = FALSE] %*% t(X)
      d_sq    <- X_sq[i_batch] +
                 matrix(X_sq, nrow = length(i_batch), ncol = n_obs, byrow = TRUE) -
                 2 * cross
      d_sq[d_sq < 0] <- 0
      for (a in seq_along(i_batch)) {
        i        <- i_batch[a]
        dv       <- sqrt(d_sq[a, ])
        dv[i]    <- Inf
        ord      <- order(dv)[seq_len(n_neighbors)]
        knn_indices[i, ] <- ord
        knn_dist[i, ]    <- dv[ord]
      }
    }
  }

  list(
    distances      = .knn_to_sparse(knn_dist, knn_indices, n_obs),
    connectivities = .connectivities_gaussian(knn_dist, knn_indices, n_obs)
  )
}

#' BiocNeighbors KNN implementation
#'
#' Delegates to \code{BiocNeighbors::findKNN()} with the KMKNN algorithm.
#' Requires the \pkg{BiocNeighbors} package.
#' @keywords internal
.build_bioc_neighbors <- function(X, n_neighbors) {
  res   <- BiocNeighbors::findKNN(X, k = n_neighbors,
                                    BNPARAM = BiocNeighbors::KmknnParam())
  n_obs <- nrow(X)
  list(
    distances      = .knn_to_sparse(res$distance, res$index, n_obs),
    connectivities = .connectivities_gaussian(res$distance, res$index, n_obs)
  )
}

#' Convert KNN indices and distances to a sparse matrix
#'
#' @param knn_dist  n_obs × k matrix of KNN distances.
#' @param knn_indices n_obs × k matrix of KNN neighbour indices (1-based).
#' @param n_obs Number of observations.
#' @return Sparse n_obs × n_obs distance matrix (dgCMatrix).
#' @keywords internal
.knn_to_sparse <- function(knn_dist, knn_indices, n_obs) {
  k    <- ncol(knn_indices)
  rows <- rep(seq_len(n_obs), each = k)
  cols <- as.vector(t(knn_indices))
  vals <- as.vector(t(knn_dist))
  Matrix::sparseMatrix(i = rows, j = cols, x = vals,
                        dims = c(n_obs, n_obs))
}

#' Build a symmetrised Gaussian-kernel connectivity matrix
#'
#' Computes Gaussian weights \eqn{w_{ij} = \exp(-d_{ij}^2 / \sigma_i^2)}
#' where \eqn{\sigma_i} is the distance to the k-th nearest neighbour of cell
#' i, then symmetrises by averaging: \eqn{(C + C^\top) / 2}.
#' @keywords internal
.connectivities_gaussian <- function(knn_dist, knn_indices, n_obs) {
  k <- ncol(knn_indices)
  # Local bandwidth: distance to the k-th nearest neighbour
  sigma           <- knn_dist[, k]
  sigma[sigma == 0] <- 1
  # Gaussian kernel weights
  weights <- exp(-knn_dist^2 / sigma^2)

  rows <- rep(seq_len(n_obs), each = k)
  cols <- as.vector(t(knn_indices))
  vals <- as.vector(t(weights))

  conn <- Matrix::sparseMatrix(i = rows, j = cols, x = vals,
                                dims = c(n_obs, n_obs))
  # Symmetrise: average both directions
  (conn + Matrix::t(conn)) / 2
}

#' Retain only the top-k edges per row in a connectivity matrix
#'
#' Zeros out all but the k largest values in each row of a sparse connectivity
#' matrix (dense conversion used internally).
#' @keywords internal
.select_top_k <- function(conn, k) {
  conn_dense <- as.matrix(conn)
  for (i in seq_len(nrow(conn_dense))) {
    row_vals <- conn_dense[i, ]
    if (sum(row_vals > 0) > k) {
      thresh                            <- sort(row_vals, decreasing = TRUE)[k + 1]
      conn_dense[i, conn_dense[i, ] <= thresh] <- 0
    }
  }
  Matrix::Matrix(conn_dense, sparse = TRUE)
}

# ====================================================================================
# Internal helper functions used by both build_neighbor_graph() and compute_moments()
# ====================================================================================

#' Check that the neighbor graph has been built
#'
#' Stops with a descriptive error if \code{@graphs$neighbors} is absent or
#' contains neither connectivities nor distances.
#' @keywords internal
.verify_neighbors <- function(pv) {
  if (is.null(pv@graphs[["neighbors"]]))
    stop(paste(
      "No neighbor graph found in @graphs$neighbors.",
      "Please run build_neighbor_graph() first.",
      "Example: pv <- build_neighbor_graph(pv, n_neighbors = 30)"
    ))
  if (is.null(pv@graphs$neighbors$connectivities) &&
      is.null(pv@graphs$neighbors$distances))
    stop("@graphs$neighbors exists but contains no connectivities or distances.")
  invisible(NULL)
}

#' Retrieve the row-normalised connectivity matrix
#'
#' Extracts the connectivity (or distance) matrix from
#' \code{@graphs$neighbors}, optionally truncates to top-k edges, optionally
#' expands to second-order neighbours, and row-normalises so each row sums
#' to 1.  Used internally by \code{compute_moments()}.
#' @keywords internal
.get_connectivities <- function(pv,
                                mode              = "connectivities",
                                n_neighbors       = NULL,
                                recurse_neighbors = FALSE) {
  .verify_neighbors(pv)

  conn <- pv@graphs$neighbors[[mode]]
  if (is.null(conn))
    stop(sprintf("'%s' not found in @graphs$neighbors.", mode))

  conn <- Matrix::Matrix(conn, sparse = TRUE)

  # Optionally truncate to top-k edges
  if (!is.null(n_neighbors))
    conn <- .select_top_k(conn, n_neighbors)

  # Optionally extend to second-order neighbours
  if (recurse_neighbors) {
    conn2  <- conn %*% conn
    conn   <- conn + conn2
    conn@x <- rep(1, length(conn@x))
  }

  # Row-normalise (each row sums to 1)
  row_sums           <- Matrix::rowSums(conn)
  row_sums[row_sums == 0] <- 1
  conn <- conn / row_sums

  conn
}
