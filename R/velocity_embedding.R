# =============================================================================
# PlantVelocity: velocity_embedding.R
# Velocity embedding projection: project high-dimensional velocity vectors
# into a low-dimensional embedding space (UMAP / PCA / tSNE)
# Corresponds to scVelo: tools/velocity_embedding.py
#
# Core formula (scVelo's default retain_scale = FALSE):
#   dE_ij = (E[j] - E[i]) / ||E[j] - E[i]||
#   V_emb[i] = sum_j (T[i,j] - mean(T[i, neighbours])) * dE_ij
# The transition diagonal is removed without renormalising the remaining rows.
# =============================================================================

#' @title Compute Velocity Embedding
#' @description
#' Project high-dimensional RNA velocity vectors onto a low-dimensional
#' embedding (for example UMAP, PCA, or tSNE). The
#' \code{projection = "centered_direction"} option uses scVelo's default projection:
#'
#' \deqn{V_{\mathrm{emb},i} = \sum_{j \in N_i}
#'   (T_{ij} - \overline{T}_{i,N_i})\,
#'   \frac{E_j - E_i}{\lVert E_j - E_i \rVert}}
#'
#' Here \eqn{T} is the absolute-row-normalised transition matrix returned by
#' \code{build_transition_matrix()}, with its diagonal removed without further
#' normalisation. \eqn{N_i} contains its non-zero off-diagonal entries.
#' Coincident coordinates contribute zero displacement. Isolated cells retain
#' zero vectors instead of propagating undefined empty-neighbour means.
#'
#' The default \code{projection = "expected_displacement"} restores the matrix projection
#' from PlantVelocity 0.5.1 (commit \code{2a9a4c6}):
#' \deqn{V_{\mathrm{emb}} = T E - E.}
#' This mode uses non-negative row-normalized transition probabilities, retains
#' requested self-transitions, and sets empty transition rows to zero movement. It uses
#' the current transition-kernel construction; it does not restore the entire
#' historical analysis pipeline. Neither mode estimates displacement per unit
#' of biological time.
#'
#' @section Output written to \code{@reductions}:
#' \describe{
#'   \item{\code{<vkey>_<reduction>}}{Numeric matrix
#'     (\eqn{n\_cells \times d}) of projected velocity vectors in embedding
#'     coordinates. Stored in \code{@reductions} alongside PCA/UMAP results,
#'     following the same cell-coordinate convention. Its \code{"projection"}
#'     attribute records the selected method.}
#' }
#'
#' @param pv A \code{plantvelo} object with a velocity graph already computed
#'   by \code{compute_velocity_graph()}.
#' @param reduction Character scalar. Dimensionality reduction to project into.
#'   Matched case-insensitively against \code{names(pv@reductions)}; a leading
#'   \code{"X_"} prefix is ignored during matching. Default \code{"umap"}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#' @param scale Numeric. Exponential kernel scale used when constructing the
#'   transition matrix. Default \code{10}.
#' @param self_transitions Logical. Estimate self-transition confidence before
#'   kernel normalisation. Default \code{TRUE}.
#' @param use_negative_cosines Logical. Include signed reverse weights from
#'   negative cosine similarities, matching scVelo's projection default.
#'   Default \code{NULL}: \code{TRUE} for \code{"centered_direction"} and
#'   \code{FALSE} for \code{"expected_displacement"}. Explicit \code{TRUE}
#'   is incompatible with \code{"expected_displacement"}.
#' @param weight_diffusion Numeric in \eqn{[0,1]} or \code{NULL}. Weight used to
#'   mix an embedding-space Gaussian diffusion kernel with the velocity kernel.
#'   A positive value enables diffusion in \code{reduction}. Default \code{NULL}.
#' @param autoscale Logical. Apply method-specific display scaling. The default
#'   \code{TRUE} uses the Matplotlib-compatible quiver scale for centered
#'   directions and the restored 0.5.1 median-spacing scale for expected
#'   displacements. Use \code{FALSE} to retain unscaled projection values.
#' @param all_comps Logical. Use all embedding dimensions. If \code{FALSE},
#'   only the first two embedding dimensions are used. Default \code{TRUE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#' @param projection Character. \code{"expected_displacement"} (default)
#'   computes the expected coordinate change for one graph step using
#'   non-negative transition probabilities. \code{"centered_direction"}
#'   combines unit neighbour directions with mean-centered transition weights.
#'
#' @section Automatic scaling:
#' For \code{"centered_direction"}, scaling reproduces
#' \code{scvelo.tools.velocity_embedding.quiver_autoscale()}
#' with Matplotlib's default figure geometry and 5 percent axis margins. To match
#' customized Python plotting settings, set \code{options(PlantVelocity.quiver =
#' list(figsize = c(width, height), subplot = c(left, right, bottom, top),
#' margins = c(xmargin, ymargin)))}. Only the first two components determine the
#' scalar; the same scalar is applied to all projected components.
#' For \code{"expected_displacement"}, the restored scale makes the median
#' non-zero norm 1 percent of the median sorted coordinate spacing along the
#' first axis when both statistics are positive. This is a display convention,
#' not a time calibration. Zero fields remain zero. Recomputing a projection
#' clears its cached stream field and paths.
#'
#' @return A \code{plantvelo} object with
#'   \code{@reductions[[paste0(vkey, "_", reduction)]]} added as a numeric matrix
#'   of size \eqn{n\_cells \times d}.
#'
#' @examples
#' \dontrun{
#' pv <- compute_velocity_embedding(pv, reduction = "umap")
#' # result stored in pv@reductions$velocity_umap
#' pv <- compute_velocity_embedding(
#'   pv, projection = "centered_direction", autoscale = FALSE
#' )
#' }
#' @export
compute_velocity_embedding <- function(pv,
                                       reduction            = "umap",
                                       vkey                 = "velocity",
                                       scale                = 10,
                                       self_transitions     = TRUE,
                                       use_negative_cosines = NULL,
                                       weight_diffusion     = NULL,
                                       autoscale            = TRUE,
                                       all_comps            = TRUE,
                                       verbose              = TRUE,
                                       projection = c("expected_displacement", "centered_direction")) {

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")
  projection <- match.arg(projection)
  if (is.null(use_negative_cosines))
    use_negative_cosines <- identical(projection, "centered_direction")
  if (!is.logical(use_negative_cosines) || length(use_negative_cosines) != 1L ||
      is.na(use_negative_cosines))
    stop("`use_negative_cosines` must be NULL, TRUE or FALSE.")
  if (projection == "expected_displacement" && use_negative_cosines)
    stop("`expected_displacement` requires non-negative transition probabilities; set `use_negative_cosines = FALSE`.")

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
  if (!is.numeric(E) || any(!is.finite(E)))
    stop("Embedding coordinates must be a finite numeric matrix.")
  if (n_dim < 2L) stop("Velocity projection requires at least two embedding components.")
  cell_names <- colnames(pv@layers[["spliced"]])
  if (!is.null(cell_names) && !is.null(rownames(E))) {
    if (anyDuplicated(rownames(E)) || !setequal(cell_names, rownames(E)))
      stop("Embedding and expression layers must identify the same cells.")
    E <- E[cell_names, , drop = FALSE]
  }

  # ----------------build Markov transition matrix----------------

  T_mat <- build_transition_matrix(pv,
                                   vkey                 = vkey,
                                   scale                = scale,
                                   self_transitions     = self_transitions,
                                   use_negative_cosines = use_negative_cosines,
                                   weight_diffusion     = weight_diffusion,
                                   basis = if (!is.null(weight_diffusion) &&
                                               weight_diffusion > 0) reduction else NULL)
  if (projection == "expected_displacement" && any(T_mat@x < 0))
    stop("`expected_displacement` requires non-negative transition probabilities.")
  if (nrow(T_mat) != n_obs)
    stop("Velocity graph and embedding must have the same number of cells.")
  if (!is.null(rownames(T_mat)) && !is.null(rownames(E)) &&
      !identical(rownames(T_mat), rownames(E))) {
    if (!setequal(rownames(T_mat), rownames(E)))
      stop("Velocity graph and embedding identify different cells.")
    T_mat <- T_mat[rownames(E), rownames(E), drop = FALSE]
  }

  if (projection == "expected_displacement") {
    # Restored from 0.5.1, commit 2a9a4c6 (projection and autoscaling).
    V_emb <- as.matrix(T_mat %*% E) - E
    # Empty transition rows have no expected movement.
    V_emb[Matrix::rowSums(T_mat) == 0, ] <- 0

    # ----------------autoscale projected vectors----------------
    #
    # Historical scaling: median vector norm is 1% of the median
    # spacing along the first embedding dimension.

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
  } else {
    # ----------------project velocity into embedding space----------------
    #
    # Use sparse row access to avoid an n_cells x n_cells dense allocation.
    diag(T_mat) <- 0
    T_mat <- methods::as(Matrix::drop0(T_mat), "RsparseMatrix")
    V_emb <- matrix(0, nrow = n_obs, ncol = n_dim)
    for (i in seq_len(n_obs)) {
      start <- T_mat@p[i] + 1L
      end <- T_mat@p[i + 1L]
      if (start > end) next
      entries <- start:end
      neighbours <- T_mat@j[entries] + 1L
      displacement <- sweep(E[neighbours, , drop = FALSE], 2L, E[i, ], "-")
      lengths <- sqrt(rowSums(displacement^2))
      lengths[lengths == 0] <- 1
      displacement <- displacement / lengths
      probabilities <- T_mat@x[entries]
      V_emb[i, ] <- colSums(displacement * (probabilities - mean(probabilities)))
    }

    # ----------------autoscale projected vectors----------------
    #
    if (autoscale) {
      quiver_scale <- .velocity_quiver_autoscale(E, V_emb)
      if (is.finite(quiver_scale) && quiver_scale > 0)
        V_emb <- V_emb / (3 * quiver_scale)
    }
  }

  # ----------------restore dimnames----------------

  rownames(V_emb) <- rownames(E) %||% cell_names
  colnames(V_emb) <- colnames(E)
  attr(V_emb, "projection") <- projection

  # ----------------store in @reductions----------------
  # Velocity embeddings share the same semantics as PCA/UMAP coordinates,
  # so @reductions is a more appropriate destination than @misc.

  out_key <- paste0(vkey, "_", reduction)
  pv@reductions[[out_key]] <- V_emb
  cached <- pv@velocity[["stream_params"]]
  basis_norm <- function(x) tolower(sub("^[Xx]_", "", x))
  if (is.null(cached$grid) ||
      (identical(basis_norm(cached$grid$reduction), basis_norm(reduction)) &&
       identical(cached$grid$vkey, vkey))) {
    pv@velocity[c("stream", "streamlines", "stream_params")] <- NULL
  }

  if (verbose)
    message(sprintf(
      "Velocity embedding stored in @reductions$%s (%d x %d).",
      out_key, n_obs, n_dim
    ))

  return(pv)
}


.ensure_velocity_projection <- function(pv, reduction, vkey, projection = NULL) {
  if (!is.null(projection))
    projection <- match.arg(projection, c("expected_displacement", "centered_direction"))
  velocity <- .find_velocity_embedding(pv, reduction, vkey)
  if (is.null(velocity) ||
      (!is.null(projection) && !identical(attr(velocity, "projection"), projection))) {
    pv <- compute_velocity_embedding(pv, reduction = reduction, vkey = vkey,
      projection = projection, verbose = FALSE)
  }
  pv
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

.velocity_quiver_autoscale <- function(embedding, velocity) {
  scale_factor <- max(abs(embedding))
  if (!is.finite(scale_factor) || scale_factor == 0 || nrow(embedding) == 0L)
    return(0)
  geometry <- utils::modifyList(
    list(figsize = c(6.4, 4.8), subplot = c(0.125, 0.9, 0.11, 0.88),
         margins = c(0.05, 0.05)),
    getOption("PlantVelocity.quiver", list())
  )
  if (!is.numeric(geometry$figsize) || length(geometry$figsize) != 2L ||
      any(!is.finite(geometry$figsize)) || any(geometry$figsize <= 0) ||
      !is.numeric(geometry$subplot) || length(geometry$subplot) != 4L ||
      any(!is.finite(geometry$subplot)) ||
      geometry$subplot[2L] <= geometry$subplot[1L] ||
      geometry$subplot[4L] <= geometry$subplot[3L] ||
      !is.numeric(geometry$margins) || length(geometry$margins) != 2L ||
      any(!is.finite(geometry$margins)) || any(geometry$margins <= -0.5))
    stop("Invalid geometry in `options(PlantVelocity.quiver)`.")
  axis_span <- function(values, margin) {
    limits <- range(values)
    magnitude <- max(abs(limits))
    # Matplotlib's linear locator uses nonsingular(..., expander = 0.05).
    if (magnitude < 1e-300) {
      limits <- c(-0.05, 0.05)
    } else if (diff(limits) <= magnitude * 1e-15) {
      limits <- limits + c(-1, 1) * 0.05 * magnitude
    }
    diff(limits) * (1 + 2 * margin)
  }
  normalized <- embedding[, 1:2, drop = FALSE] / scale_factor
  x_span <- axis_span(normalized[, 1L], geometry$margins[1L])
  y_span <- axis_span(normalized[, 2L], geometry$margins[2L])
  aspect <- geometry$figsize[2L] * diff(geometry$subplot[3:4]) /
    (geometry$figsize[1L] * diff(geometry$subplot[1:2]))
  arrow_lengths <- sqrt((velocity[, 1L] / x_span)^2 +
                        (aspect * velocity[, 2L] / y_span)^2)
  1.8 * max(10, sqrt(nrow(embedding))) * mean(arrow_lengths) / scale_factor
}


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
