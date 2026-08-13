# =============================================================================
# PlantVelocity: velocity_embedding.R
# Velocity embedding projection: project high-dimensional velocity vectors
# into a low-dimensional embedding space (UMAP / PCA / tSNE)
# Corresponds to scVelo: tools/velocity_embedding.py
#
# Core formula:
#   V_emb[i] = sum_j T[i,j] * (E[j] - E[i])
#            = (T %*% E) - E
# where T is the row-normalised Markov transition matrix and E is the
# n_cells x d low-dimensional embedding matrix.
# =============================================================================

#' @title Compute Velocity Embedding
#' @description
#' Project high-dimensional RNA velocity vectors onto a low-dimensional
#' embedding (for example UMAP, PCA, or tSNE) by computing the expected
#' displacement under the velocity-guided Markov transition matrix:
#'
#' \deqn{V_{\mathrm{emb},i} = \sum_j T_{ij}\,(E_j - E_i) = (T E)_i - E_i}
#'
#' where \eqn{T} is the row-normalised transition matrix returned by
#' \code{build_transition_matrix()} and \eqn{E} is the
#' \eqn{n\_cells \times d} embedding matrix.
#'
#' @section Output written to \code{@reductions}:
#' \describe{
#'   \item{\code{<vkey>_<reduction>}}{Numeric matrix
#'     (\eqn{n\_cells \times d}) of projected velocity vectors in embedding
#'     coordinates. Stored in \code{@reductions} alongside PCA/UMAP results,
#'     following the same cell-coordinate convention.}
#' }
#'
#' @param pv A \code{plantvelo} object with a velocity graph already computed
#'   by \code{compute_velocity_graph()}.
#' @param reduction Character scalar. Dimensionality reduction to project into.
#'   Matched case-insensitively against \code{names(pv@reductions)}; a leading
#'   \code{"X_"} prefix is ignored during matching. Default \code{"umap"}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#' @param scale Numeric. Softmax scale parameter used when constructing the
#'   transition matrix. Default \code{10}.
#' @param self_transitions Logical. Add self-loops before softmax
#'   normalisation. Recommended \code{TRUE} to stabilise projected vectors for
#'   cells with weak outgoing transitions. Default \code{TRUE}.
#' @param use_negative_cosines Logical. Subtract the negative-cosine graph to
#'   penalise backward transitions. Default \code{FALSE}.
#' @param weight_diffusion Numeric in \eqn{[0,1]} or \code{NULL}. Weight used to
#'   mix the velocity transition matrix with a diffusion kernel derived from
#'   neighbour connectivities. Default \code{NULL}.
#' @param autoscale Logical. Rescale projected vectors so that the median
#'   vector norm is approximately 1\% of the median inter-cell spacing in the
#'   first embedding dimension. Default \code{TRUE}.
#' @param all_comps Logical. Use all embedding dimensions. If \code{FALSE},
#'   only the first two embedding dimensions are used. Default \code{TRUE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{plantvelo} object with
#'   \code{@reductions[[paste0(vkey, "_", reduction)]]} added as a numeric matrix
#'   of size \eqn{n\_cells \times d}.
#'
#' @examples
#' \dontrun{
#' pv <- compute_velocity_embedding(pv, reduction = "umap")
#' # result stored in pv@reductions$velocity_umap
#' }
#' @export
compute_velocity_embedding <- function(pv,
                                       reduction                = "umap",
                                       vkey                 = "velocity",
                                       scale                = 10,
                                       self_transitions     = TRUE,
                                       use_negative_cosines = FALSE,
                                       weight_diffusion     = NULL,
                                       autoscale            = TRUE,
                                       all_comps            = TRUE,
                                       verbose              = TRUE) {

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  gkey <- paste0(vkey, "_graph")
  if (is.null(pv@graphs[[gkey]]))
    stop(sprintf(
      "Velocity graph '%s' not found. Run compute_velocity_graph() first.", gkey
    ))

  # ----------------retrieve embedding coordinates (n_cells x d)----------------

  E_full <- .get_embedding(pv, reduction)
  n_obs  <- nrow(E_full)
  n_dim  <- if (isTRUE(all_comps)) ncol(E_full) else min(2L, ncol(E_full))
  E      <- E_full[, seq_len(n_dim), drop = FALSE]

  # ----------------build Markov transition matrix----------------

  T_mat <- build_transition_matrix(pv,
                                   vkey                 = vkey,
                                   scale                = scale,
                                   self_transitions     = self_transitions,
                                   use_negative_cosines = use_negative_cosines,
                                   weight_diffusion     = weight_diffusion)

  # ----------------project velocity into embedding space----------------
  #
  # T is row-normalised (row sums = 1), therefore:
  #   V_emb[i] = sum_j T[i,j] * (E[j] - E[i])
  #            = (T %*% E)[i,] - E[i,]
  #
  # Sparse matrix multiplication avoids converting T into a dense matrix
  # (O(nnz x d) instead of O(n^2 x d)).

  V_emb <- as.matrix(T_mat %*% E) - E

  # ----------------autoscale projected vectors----------------
  #
  # Target: median vector norm approximately 1% x median neighbour spacing in the first
  # embedding dimension, matching the logic used in scVelo.

  if (autoscale) {
    norms    <- sqrt(rowSums(V_emb^2, na.rm = TRUE))
    med_norm <- stats::median(norms[norms > 0], na.rm = TRUE)

    if (!is.na(med_norm) && med_norm > 0) {
      emb_scale <- stats::median(abs(diff(sort(E[, 1]))), na.rm = TRUE)
      if (!is.na(emb_scale) && emb_scale > 0) {
        V_emb <- V_emb * (0.01 * emb_scale / med_norm)
      }
    }
  }

  # ----------------restore dimnames----------------

  rownames(V_emb) <- colnames(pv@layers[["spliced"]])
  colnames(V_emb) <- colnames(E)

  # ----------------store in @reductions----------------
  # Velocity embeddings share the same semantics as PCA/UMAP coordinates,
  # so @reductions is a more appropriate destination than @misc.

  out_key <- paste0(vkey, "_", reduction)
  pv@reductions[[out_key]] <- V_emb

  if (verbose)
    message(sprintf(
      "Velocity embedding stored in @reductions$%s (%d x %d).",
      out_key, n_obs, n_dim
    ))

  return(pv)
}


#' @title Get Velocity Embedding Matrix
#' @description
#' Retrieve a previously computed velocity embedding from \code{@reductions}.
#' This is a convenience wrapper around
#' \code{pv@reductions[[paste0(vkey, "_", reduction)]]}.
#'
#' @param pv A \code{plantvelo} object.
#' @param reduction Character scalar. Embedding reduction name. Default
#'   \code{"umap"}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#'
#' @return Numeric matrix of size \eqn{n\_cells \times d}, or \code{NULL} if
#'   the requested velocity embedding has not yet been computed.
#' @examples
#' \dontrun{
#' V_emb <- get_velocity_embedding(pv, reduction = "umap")
#' }
#' @export
get_velocity_embedding <- function(pv, reduction = "umap", vkey = "velocity") {
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  key <- paste0(vkey, "_", reduction)
  mat <- pv@reductions[[key]]

  if (is.null(mat))
    message(sprintf(
      "No velocity embedding for reduction '%s'. Run compute_velocity_embedding() first.",
      reduction
    ))

  mat
}


# =============================================================================
# Internal helper functions
# =============================================================================

#' Retrieve an embedding matrix from \code{@reductions}
#'
#' Supports flexible matching of embedding names:
#' \itemize{
#'   \item Case-insensitive matching, for example \code{"UMAP"},
#'     \code{"umap"}, or \code{"Umap"}
#'   \item With or without an \code{"X_"} prefix, for example
#'     \code{"X_umap"} or \code{"umap"}
#'   \item Seurat-style reduction objects stored in
#'     \code{@reductions[[key]]@cell.embeddings}
#' }
#'
#' If the matched embedding matrix has no column names, default names of the
#' form \code{<BASIS>_1}, \code{<BASIS>_2}, ... are assigned.
#'
#' @param pv A \code{plantvelo} object.
#' @param reduction Character scalar. Requested embedding name.
#'
#' @return Numeric matrix of embedding coordinates.
#' @keywords internal
.get_embedding <- function(pv, reduction) {

  avail <- names(pv@reductions)
  if (length(avail) == 0)
    stop("@reductions is empty. Run dimensionality reduction first.")

  # Compare reduction and available keys after stripping the "X_" prefix
  # and converting to lower case.
  basis_norm <- tolower(sub("^[Xx]_", "", reduction))

  for (key in avail) {
    key_norm <- tolower(sub("^[Xx]_", "", key))
    if (key_norm == basis_norm) {
      object <- pv@reductions[[key]]
      emb <- if (is.matrix(object)) {
        object
      } else if (is.data.frame(object)) {
        as.matrix(object)
      } else if (isS4(object) &&
                 "cell.embeddings" %in% methods::slotNames(object)) {
        methods::slot(object, "cell.embeddings")
      } else if (!is.null(attr(object, "cell.embeddings"))) {
        attr(object, "cell.embeddings")
      } else {
        stop(sprintf(
          "Reduction '%s' does not contain cell embeddings.",
          reduction
        ))
      }
      # Add default column names if the reduction matrix has none
      if (is.null(colnames(emb))) {
        colnames(emb) <- paste0(toupper(basis_norm), "_", seq_len(ncol(emb)))
      }
      return(emb)
    }
  }

  stop(sprintf(
    "Reduction '%s' not found in @reductions. Available: %s",
    reduction, paste(avail, collapse = ", ")
  ))
}
