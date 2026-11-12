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
#'   \item{\code{"base"}}{Blocked exact Euclidean search without a full
#'     dense cell-by-cell distance matrix. Ties are ordered by cell name.}
#'   \item{\code{"bioc"}}{Uses \pkg{BiocNeighbors} (exact KNN via KD-tree /
#'     KMKNN). Recommended for n > 10,000 cells.}
#' }
#'
#' @section Representation selection:
#' PCA is required. Named embeddings are aligned to spliced cell names;
#' unnamed complete embeddings retain their positional interpretation.
#' Values must be finite. At least two cells are required.
#'
#' @section Result invalidation:
#' Successful reconstruction clears moments, velocity graphs and their dependent
#' results. Fits survive only when existing run records explicitly specify raw
#' expression and no connected-state fitting. Raw counts, base embeddings,
#' custom annotations and independent IR statistics are preserved. RNG is local.
#'
#' @section Output written to \code{@graphs$neighbors}:
#' \describe{
#'   \item{\code{connectivities}}{Symmetric sparse matrix (n_cells x n_cells).
#'     Gaussian-kernel weights computed from KNN distances, symmetrised by
#'     averaging.}
#'   \item{\code{distances}}{Sparse matrix (n_cells x n_cells). Raw Euclidean
#'     distances to each KNN.}
#'   \item{\code{params}}{Named list recording: \code{n_neighbors},
#'     \code{method}, \code{n_pcs}.}
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

  axes <- .pv_axes(pv)
  n_neighbors <- .validate_whole(n_neighbors,"n_neighbors",1L)
  n_pcs <- .validate_whole(n_pcs,"n_pcs",1L)
  if (length(axes$cells)<2L) stop("At least two cells are required.")

  if (method == "bioc" && !requireNamespace("BiocNeighbors", quietly = TRUE))
    stop(
      "Package `BiocNeighbors` is required for method = 'bioc'.\n",
      "Install with: BiocManager::install('BiocNeighbors')"
    )

  n_neighbors <- as.integer(n_neighbors)
  n_pcs       <- as.integer(n_pcs)

  # ----------------prepare representation----------------

  if (verbose) message("[1] Extracting cell representation ...")
  restore_seed <- .local_seed(random_state)
  on.exit(restore_seed(),add=TRUE)

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

  if (n_neighbors>n_obs-1L && verbose) message("Using n_neighbors = n_cells - 1.")
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

  dimnames(result$distances) <- dimnames(result$connectivities) <- list(axes$cells,axes$cells)
  pv <- .invalidate_upstream(pv,neighbors_changed=TRUE)
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

  axes <- .pv_axes(pv)
  X <- .canonical_matrix(X,"PCA",axes$cells,colnames(X) %||% paste0("PC",seq_len(ncol(X))))
  if (!ncol(X) || any(!is.finite(X))) stop("PCA must contain finite values and at least one dimension.")
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
#' Uses a blocked dot-product expansion of \eqn{\|a - b\|^2 = \|a\|^2 + \|b\|^2 - 2a \cdot b}.
#' Returns a list with \code{distances} (sparse) and \code{connectivities}
#' (Gaussian-kernel, symmetrised).
#' @keywords internal
.build_base_neighbors <- function(X, n_neighbors) {
  n <- nrow(X)
  knn_dist <- matrix(0,n,n_neighbors)
  knn_indices <- matrix(0L,n,n_neighbors)
  square <- rowSums(X^2)
  labels <- rownames(X) %||% sprintf("%012d",seq_len(n))
  for (start in seq.int(1L,n,by=256L)) {
    idx <- start:min(n,start+255L)
    d2 <- outer(square[idx],square,"+")-2*tcrossprod(X[idx,,drop=FALSE],X)
    d2 <- pmax(d2,0)
    if (any(!is.finite(d2))) stop("PCA distance calculation overflowed; rescale the representation.")
    for (a in seq_along(idx)) {
      d <- sqrt(d2[a,]); d[idx[a]] <- Inf
      chosen <- order(d,labels)[seq_len(n_neighbors)]
      knn_indices[idx[a],] <- chosen
      knn_dist[idx[a],] <- d[chosen]
    }
  }
  list(distances=.knn_to_sparse(knn_dist,knn_indices,n),
       connectivities=.connectivities_gaussian(knn_dist,knn_indices,n))
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
  labels <- rownames(X)
  for (i in seq_len(n_obs)) {
    ord <- order(res$distance[i,],labels[res$index[i,]])
    res$distance[i,] <- res$distance[i,ord]
    res$index[i,] <- res$index[i,ord]
  }
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
  if (!is.matrix(knn_indices) || !is.matrix(knn_dist) ||
      !identical(dim(knn_indices),dim(knn_dist)) || nrow(knn_indices)!=n_obs ||
      any(!is.finite(knn_indices)) || any(knn_indices!=floor(knn_indices)) ||
      any(knn_indices<1 | knn_indices>n_obs)) stop("Invalid KNN index matrix.")
  .nonnegative_values(knn_dist,"KNN distances")
  if (any(knn_indices==row(knn_indices)) ||
      any(apply(knn_indices,1,anyDuplicated)>0)) stop("KNN indices contain self edges or duplicate indices.")
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
#' matrix without dense conversion. Ties use cell names.
#' @keywords internal
.select_top_k <- function(conn, k) {
  k <- .validate_whole(k,"n_neighbors",1L)
  edges <- Matrix::summary(methods::as(conn,"generalMatrix"))
  edges <- edges[edges$x>0,,drop=FALSE]
  labels <- colnames(conn) %||% sprintf("%012d",seq_len(ncol(conn)))
  rows <- split(seq_len(nrow(edges)),edges$i)
  keep <- unlist(lapply(rows,function(idx) {
    own <- idx[edges$i[idx]==edges$j[idx]]
    other <- idx[edges$i[idx]!=edges$j[idx]]
    c(own,head(other[order(-edges$x[other],labels[edges$j[other]])],k))
  }),use.names=FALSE)
  Matrix::sparseMatrix(i=edges$i[keep],j=edges$j[keep],x=edges$x[keep],
    dims=dim(conn),dimnames=dimnames(conn))
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
#' Extracts the connectivity matrix, or converts distances to Gaussian weights, from
#' \code{@graphs$neighbors}, optionally truncates to top-k edges, optionally
#' expands to second-order neighbours, and row-normalises so each row sums
#' to 1. Empty rows use their own expression. Second-order expansion retains
#' support-uniform weighting and stops above a 512 MiB sparse-output budget.
#' Used internally by \code{compute_moments()}.
#' @keywords internal
.get_connectivities <- function(pv, mode="connectivities", n_neighbors=NULL,
                                recurse_neighbors=FALSE) {
  .verify_neighbors(pv)
  cells <- .pv_axes(pv)$cells
  if (!mode %in% c("connectivities","distances")) stop("Unknown neighbor mode.")
  if (!is.logical(recurse_neighbors) || length(recurse_neighbors)!=1L || is.na(recurse_neighbors))
    stop("recurse_neighbors must be a logical scalar.")
  conn <- .canonical_graph(pv@graphs$neighbors[[mode]],mode,cells)
  if (identical(mode,"distances")) {
    if (!is.null(pv@graphs$neighbors$connectivities)) {
      support <- .canonical_graph(pv@graphs$neighbors$connectivities,"connectivities",cells)
      conn <- .recover_distance_support(conn,support)
    }
    conn <- .distance_connectivity(conn,n_neighbors)
  } else if (!is.null(n_neighbors)) conn <- .select_top_k(conn,n_neighbors)
  if (recurse_neighbors) conn <- .second_order_support(conn)
  out <- .normalize_smoothing(conn)
  dimnames(out) <- list(cells,cells)
  out
}

.recover_distance_support <- function(D,C) {
  e <- Matrix::summary(D)
  support <- Matrix::summary(Matrix::drop0(C))
  support <- support[support$i!=support$j & support$x>0,,drop=FALSE]
  n <- as.double(nrow(D))
  keys <- e$i+(e$j-1)*n
  missing <- !((support$i+(support$j-1)*n) %in% keys)
  if (!any(missing)) return(D)
  reverse <- match(support$j[missing]+(support$i[missing]-1)*n,keys)
  if (anyNA(reverse)) stop("Distance support has ambiguous missing zero edges; rebuild the neighbor graph.")
  # A directed KNN distance is sufficient; symmetrize after applying its
  # original row bandwidth. Adding reverse distances here would change it.
  D
}

.distance_connectivity <- function(D,k=NULL) {
  if (!is.null(k)) k <- .validate_whole(k,"n_neighbors",1L)
  # summary retains explicitly stored zero-distance edges; implicit zeros are absent.
  e <- Matrix::summary(methods::as(D,"generalMatrix"))
  e <- e[e$i!=e$j,,drop=FALSE]
  labels <- colnames(D)
  rows <- split(seq_len(nrow(e)),e$i)
  keep <- integer(); weights <- numeric()
  for (idx in rows) {
    idx <- idx[order(e$x[idx],labels[e$j[idx]])]
    if (!is.null(k)) idx <- head(idx,k)
    sigma <- max(e$x[idx]); if (sigma==0) sigma <- 1
    keep <- c(keep,idx); weights <- c(weights,exp(-(e$x[idx]/sigma)^2))
  }
  W <- Matrix::sparseMatrix(i=e$i[keep],j=e$j[keep],x=weights,
    dims=dim(D),dimnames=dimnames(D))
  (W+Matrix::t(W))/2
}

.second_order_support <- function(C) {
  e <- Matrix::summary(Matrix::drop0(C))
  adjacency <- split(e$j,factor(e$i,levels=seq_len(nrow(C))))
  rows <- vector("list",nrow(C)); count <- 0
  for (i in seq_len(nrow(C))) {
    rows[[i]] <- sort(unique(c(adjacency[[i]],unlist(adjacency[adjacency[[i]]],use.names=FALSE))))
    count <- count+length(rows[[i]])
    if (count*12+(nrow(C)+1)*4 > 512*1024^2)
      stop("Second-order support exceeds 512 MiB; disable recurse_neighbors or reduce n_neighbors.")
  }
  Matrix::sparseMatrix(i=rep(seq_len(nrow(C)),lengths(rows)),j=unlist(rows,use.names=FALSE),
    x=1,dims=dim(C),dimnames=dimnames(C))
}

# =============================================================================
# First-order moments
# =============================================================================

# =============================================================================
# PlantVelocity: moments.R
# Compute first-order KNN-smoothed moments for spliced, unspliced,
# and optionally intron-retained expression layers.
#
# Mathematical framework:
#   Ms = C · S^T
#   Mu = C · U^T
#   Mr = C · R^T
#
# where:
#   - C is the row-normalised cell-cell neighbour connectivity matrix
#   - S, U, and R are stored as genes × cells in @layers
#   - Ms, Mu, and Mr are stored as cells × genes in @moments
#
# This module provides the neighbourhood-smoothing step required before
# kinetics recovery and downstream velocity modelling.
# =============================================================================

#' @title Compute KNN First-Order Moments (Ms, Mu, Mr)
#' @description
#' Smooth the spliced, unspliced, and (optionally) intron-retained count layers
#' by KNN-weighted averaging. The row-normalised connectivity matrix \eqn{C}
#' from \code{@graphs$neighbors} is used as weights:
#'
#' \deqn{M_s = C \cdot S^\top, \quad M_u = C \cdot U^\top, \quad M_r = C \cdot R^\top}
#'
#' where \eqn{S}, \eqn{U}, \eqn{R} are stored as \strong{genes \eqn{\times} cells}
#' in \code{@layers}, and the results are stored as \strong{cells \eqn{\times} genes}
#' in \code{@moments}.
#'
#' The \eqn{M_r} term is the PlantVelocity-specific intron-retained moment and
#' is only computed when \code{include_ir = TRUE} and \code{@layers$ir} exists.
#'
#' @section Prerequisites:
#' \code{build_neighbor_graph()} must be run before this function.
#' The connectivity matrix is read from \code{@graphs$neighbors}.
#'
#' @section Output written to \code{@moments}:
#' \describe{
#'   \item{\code{Ms}}{cells \eqn{\times} genes. Smoothed spliced counts.}
#'   \item{\code{Mu}}{cells \eqn{\times} genes. Smoothed unspliced counts.}
#'   \item{\code{Mr}}{cells \eqn{\times} genes. Smoothed intron-retained counts
#'     retained independently for future IR statistics.}
#' }
#'
#' @param pv A \code{plantvelo} object with \code{@graphs$neighbors} already
#'   populated by \code{build_neighbor_graph()}.
#' @param include_ir Logical. Compute \code{Mr} from \code{@layers$ir} when
#'   available. Old Mr is removed when omitted or absent. Default \code{TRUE}.
#' @param n_neighbors Integer or \code{NULL}. If provided, truncate the
#'   connectivity matrix to the top \code{n_neighbors} edges per cell before
#'   smoothing. \code{NULL} uses all stored edges. Default \code{NULL}.
#' @param mode Character. Which neighbour matrix to use:
#'   \code{"connectivities"} (Gaussian-kernel weights, default) or
#'   \code{"distances"} (nearest edges converted to Gaussian weights).
#' @param recurse_neighbors Logical. Extend to second-order neighbours before
#'   smoothing. Default \code{FALSE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @details
#' All expression and graph axes are aligned by name before multiplication.
#' Successful recomputation clears moment-dependent fits/velocities, all old
#' velocity graphs and downstream results. Explicit raw-source fits survive.
#' No new persistent QC fields are added. Errors leave the input intact.
#' @return A \code{plantvelo} object with \code{@moments} populated.
#'
#' @examples
#' \dontrun{
#' pv <- build_neighbor_graph(pv, n_neighbors = 30)
#' pv <- compute_moments(pv)
#' }
#' @export
compute_moments <- function(pv,
                             include_ir        = TRUE,
                             n_neighbors       = NULL,
                             mode              = c("connectivities", "distances"),
                             recurse_neighbors = FALSE,
                             verbose           = TRUE) {

  mode <- match.arg(mode)
  axes <- .pv_axes(pv)
  if (!is.logical(include_ir) || length(include_ir)!=1L || is.na(include_ir)) stop("include_ir must be a logical scalar.")

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (is.null(pv@layers[["spliced"]]))
    stop("@layers$spliced not found. Please re-run create_plantvelo().")

  if (is.null(pv@layers[["unspliced"]]))
    stop("@layers$unspliced not found. Please re-run create_plantvelo().")

  .verify_neighbors(pv)

  # ----------------retrieve connectivity matrix----------------

  if (verbose) message("[1] Retrieving connectivity matrix ...")

  C <- .get_connectivities(pv,
                            mode              = mode,
                            n_neighbors       = n_neighbors,
                            recurse_neighbors = recurse_neighbors)

  # ----------------compute moments----------------

  if (verbose) message("[2] Computing moments ...")

  # @layers are genes × cells; .smooth_layer() transposes before multiplication
  # and returns cells × genes
  prepare <- function(key) {
    X <- .canonical_matrix(pv@layers[[key]],key,axes$genes,axes$cells)
    .nonnegative_values(X,key)
    X
  }
  Ms <- .smooth_layer(C, prepare("spliced"))
  Mu <- .smooth_layer(C, prepare("unspliced"))
  Mr <- if (include_ir && !is.null(pv@layers$ir)) .smooth_layer(C,prepare("ir")) else NULL
  pv <- .invalidate_upstream(pv,neighbors_changed=FALSE)

  pv@moments[["Ms"]] <- Ms
  pv@moments[["Mu"]] <- Mu
  pv@moments[["Mr"]] <- Mr

  if (verbose) {
    message(sprintf("    Ms: %d cells x %d genes (spliced moments).",
                    nrow(Ms), ncol(Ms)))
    message(sprintf("    Mu: %d cells x %d genes (unspliced moments).",
                    nrow(Mu), ncol(Mu)))
  }

  # ----------------Mr: PlantVelocity intron-retained moments----------------

  if (include_ir) {
    if (!is.null(pv@layers[["ir"]])) {
      # Mr was computed and validated before any object writes.
      if (verbose) message(sprintf(
        "    Mr: %d cells x %d genes (intron-retained moments).",
        nrow(Mr), ncol(Mr)
      ))
    } else {
      if (verbose) message(
        "    Note: @layers$ir not found. Mr skipped. ",
        "Re-run create_plantvelo() with a loom file containing the ",
        "IR count layer for future independent IR statistics."
      )
    }
  }

  if (verbose) message(sprintf(
    "Moments stored in @moments: %s",
    paste(names(pv@moments), collapse = ", ")
  ))

  return(pv)
}

# =============================================================================
# Internal helper functions
# =============================================================================

#' KNN-smooth a single expression layer
#'
#' Multiplies the row-normalised connectivity matrix \eqn{C} (cells × cells)
#' by the transposed layer matrix \eqn{X^\top} (cells × genes) and returns a
#' dense cells × genes smoothed matrix.
#'
#' @param C Row-normalised connectivity matrix (cells × cells, sparse).
#' @param X Expression layer matrix (genes × cells, sparse or dense).
#' @return Dense numeric matrix (cells × genes) with dimnames preserved.
#' @keywords internal
.smooth_layer <- function(C, X) {
  result <- matrix(0,nrow(C),nrow(X),dimnames=list(rownames(C),rownames(X)))
  if (!nrow(X)) return(result)
  for (start in seq.int(1L,nrow(X),by=256L)) {
    idx <- start:min(nrow(X),start+255L)
    result[,idx] <- as.matrix(C %*% Matrix::t(X[idx,,drop=FALSE]))
  }
  if (any(!is.finite(result))) stop("Moment calculation produced non-finite values.")
  result
}

# =============================================================================
# Smoothing normalization
# =============================================================================

.normalize_smoothing <- function(C) {
  total <- Matrix::rowSums(C)
  empty <- which(total==0)
  if (length(empty)) C <- C + Matrix::sparseMatrix(i=empty,j=empty,x=1,
    dims=dim(C),dimnames=dimnames(C))
  Matrix::drop0(Matrix::Diagonal(x=1/Matrix::rowSums(C)) %*% C)
}

# =============================================================================
# Upstream dependency invalidation
# =============================================================================

# Use only existing run records. Unknown provenance is cleared conservatively.
.invalidate_upstream <- function(pv, neighbors_changed=FALSE) {
  fit <- pv@misc$recover_dynamics
  if (!is.list(fit)) fit <- list()
  keep_fit <- isTRUE(fit$use_raw) &&
    (!neighbors_changed || identical(fit$fit_connected_states,FALSE))
  run_keys <- as.character(names(pv@misc))[vapply(pv@misc,function(x) is.list(x) &&
    is.character(x$mode) && length(x$mode)==1L && !is.na(x$mode) &&
    x$mode %in% c("dynamical","stochastic"),logical(1))]
  run_keys <- run_keys[!is.na(run_keys) & nzchar(run_keys) & endsWith(run_keys,"_params")]
  graph_keys <- names(pv@graphs)
  vkeys <- unique(c("velocity",sub("_params$","",run_keys),
    sub("_graph(_neg)?$","",graph_keys[grepl("_graph(_neg)?$",graph_keys)])))
  for (key in vkeys) {
    run <- pv@misc[[paste0(key,"_params")]]
    if (!is.list(run)) run <- list()
    safe <- isTRUE(run$use_raw) && (identical(run$mode,"stochastic") ||
      (identical(run$mode,"dynamical") && keep_fit))
    if (!safe) {
      pv@layers[c(key,paste0(key,c("_u","_r")))] <- NULL
      pv@misc[[paste0(key,"_params")]] <- NULL
      if (identical(key,"velocity")) pv@velocity$genes <- NULL
    }
    pv <- .invalidate_velocity_dependents(pv,key)
  }
  # Known velocity projections are removed by the helper above; also clear
  # unrecorded legacy projections under the reserved default velocity prefix.
  pv@reductions[grep("^velocity_",names(pv@reductions),value=TRUE)] <- NULL
  pv@velocity[c("stream","streamlines","stream_params")] <- NULL
  if (!keep_fit) {
    pv@kinetics$params <- NULL
    pv@layers[c("fit_t","fit_tau")] <- NULL
    pv@misc$recover_dynamics <- NULL
  }
  if (neighbors_changed) pv@moments[c("Ms","Mu","Mr")] <- NULL
  pv
}
