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
#'   available. Silently skipped if the layer is absent. Default \code{TRUE}.
#' @param n_neighbors Integer or \code{NULL}. If provided, truncate the
#'   connectivity matrix to the top \code{n_neighbors} edges per cell before
#'   smoothing. \code{NULL} uses all stored edges. Default \code{NULL}.
#' @param mode Character. Which neighbour matrix to use:
#'   \code{"connectivities"} (Gaussian-kernel weights, default) or
#'   \code{"distances"}.
#' @param recurse_neighbors Logical. Extend to second-order neighbours before
#'   smoothing. Default \code{FALSE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
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
  Ms <- .smooth_layer(C, pv@layers[["spliced"]])
  Mu <- .smooth_layer(C, pv@layers[["unspliced"]])

  pv@moments[["Ms"]] <- Ms
  pv@moments[["Mu"]] <- Mu

  if (verbose) {
    message(sprintf("    Ms: %d cells x %d genes (spliced moments).",
                    nrow(Ms), ncol(Ms)))
    message(sprintf("    Mu: %d cells x %d genes (unspliced moments).",
                    nrow(Mu), ncol(Mu)))
  }

  # ----------------Mr: PlantVelocity intron-retained moments----------------

  if (include_ir) {
    if (!is.null(pv@layers[["ir"]])) {
      Mr <- .smooth_layer(C, pv@layers[["ir"]])
      pv@moments[["Mr"]] <- Mr
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
  # X: genes × cells (dgCMatrix or matrix)
  # C: cells × cells row-normalised connectivity
  # returns: cells × genes dense matrix
  X_sp <- if (inherits(X, "sparseMatrix")) X else Matrix::Matrix(X, sparse = TRUE)
  result <- as.matrix(C %*% Matrix::t(X_sp))
  # Preserve dimnames (rownames = barcodes, colnames = gene names)
  rownames(result) <- colnames(X_sp)   # colnames of X = cell barcodes
  colnames(result) <- rownames(X_sp)   # rownames of X = gene names
  storage.mode(result) <- "double"
  result
}
