# =============================================================================
# PlantVelocity: velocity_graph.R
# Velocity transition graph: cosine similarity computation + Markov transition matrix
#
# Core formula:
#   π_{ij} = (X[j] - X[i])^T · V[i] / (‖X[j]-X[i]‖ · ‖V[i]‖)
#
# This file makes no assumption about the source of the velocity
# (the IR-excluded two-state model is supported),
# and only requires that @layers[[vkey]] stores velocity vectors.
# =============================================================================

#' @title Compute Velocity Graph
#' @description
#' Build a directed cell-cell velocity graph based on cosine similarities
#' between each cell's velocity vector and the expression difference to its
#' neighbours.
#'
#' For each cell \eqn{i} and each of its KNN neighbours \eqn{j}:
#' \deqn{\pi_{ij} = \frac{(X_j - X_i)^\top V_i}{\|X_j - X_i\| \cdot \|V_i\|}}
#'
#' The resulting sparse graph stores positive and negative cosine similarities
#' separately and can be further converted into a Markov transition matrix by
#' \code{build_transition_matrix()}.
#'
#' @section Output written to \code{@graphs}:
#' \describe{
#'   \item{\code{<vkey>_graph}}{Sparse \eqn{n\_cells \times n\_cells} matrix of
#'     positive cosine similarities only.}
#'   \item{\code{<vkey>_graph_neg}}{Sparse \eqn{n\_cells \times n\_cells} matrix
#'     of negative cosine similarities only, retaining negative values.}
#' }
#'
#' @section Output written to \code{@meta.data}:
#' \describe{
#'   \item{\code{<vkey>_self_transition}}{Per-cell self-transition probability,
#'     defined as the complement of outward transition confidence.}
#' }
#'
#' @param pv A \code{plantvelo} object with \code{@layers[[vkey]]} populated by
#'   \code{compute_velocity()}.
#' @param vkey Character scalar. Velocity layer key in \code{@layers}. Default
#'   \code{"velocity"}.
#' @param xkey Character scalar. Expression reference layer key in
#'   \code{@moments}, used to compute \eqn{\Delta X}. Default \code{"Ms"}.
#' @param n_recurse_neighbors Integer. KNN recursion depth for expanding the
#'   neighbour set. \code{1L} uses first-order neighbours only; \code{2L}
#'   additionally includes second-order neighbours. Default \code{2L}.
#' @param sqrt_transform Logical or \code{NULL}. Apply a square-root
#'   variance-stabilising transform to both velocity and \eqn{\Delta X}.
#'   \code{NULL} auto-detects from
#'   \code{@misc[[paste0(vkey, "_params")]]$mode} and defaults to
#'   \code{TRUE} for stochastic velocity. Default \code{NULL}.
#' @param gene_subset Character vector, logical vector, or \code{NULL}. Genes
#'   to include in graph construction, intersected with existing eligibility.
#'   Named logical subsets are matched by name; unnamed subsets must be full
#'   length. Unknown names and missing values are errors. \code{NULL} uses the
#'   default layer's \code{@velocity$genes}, or eligibility reconstructed from
#'   dynamical run settings. Unverifiable external/custom layers need an explicit
#'   subset. Default \code{NULL}.
#' @param approx Logical or \code{NULL}. Project \code{X} and \code{V} onto the
#'   top principal components before computing cosines. \code{NULL}
#'   auto-selects \code{TRUE} when the number of genes exceeds 100. Default
#'   \code{NULL}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @details
#' Align gene and cell axes by name and select genes before densifying inputs.
#' Require at least two finite, nonzero velocity genes after selection and
#' missing-column filtering. No missing velocities are imputed. Custom velocity
#' keys do not borrow the default layer mask. Success clears dependent embeddings,
#' stream fields, driver tables, terminal-state results and latent time; failures
#' leave old results intact. Existing cosine/PCA transformations are unchanged.
#'
#' @return A \code{plantvelo} object with updated \code{@graphs} and
#'   \code{@meta.data}.
#'
#' @examples
#' \dontrun{
#' pv <- compute_velocity_graph(pv)
#' }
#' @export
compute_velocity_graph <- function(pv,
                                   vkey                = "velocity",
                                   xkey                = "Ms",
                                   n_recurse_neighbors = 2L,
                                   sqrt_transform      = NULL,
                                   gene_subset         = NULL,
                                   approx              = NULL,
                                   verbose             = TRUE) {

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (is.null(pv@layers[[vkey]]))
    stop(sprintf(
      "Velocity layer '%s' not found. Run compute_velocity() first.", vkey
    ))

  if (is.null(pv@graphs[["neighbors"]]))
    stop("No neighbor graph found. Run build_neighbor_graph() first.")

  axes <- .velocity_axes(pv)
  for (key in c("vkey","xkey")) {
    x <- get(key)
    if (!is.character(x) || length(x)!=1L || is.na(x) || !nzchar(x))
      stop(paste(key,"must be a nonempty string."))
  }
  for (key in c("sqrt_transform","approx","verbose")) {
    x <- get(key)
    if (!is.null(x) && (!is.logical(x) || length(x)!=1L || is.na(x)))
      stop(paste(key,"must be TRUE or FALSE",if(key!="verbose") "or NULL." else "."))
  }
  if (!.finite_scalar(n_recurse_neighbors) || n_recurse_neighbors<1 ||
      n_recurse_neighbors!=floor(n_recurse_neighbors)) stop("n_recurse_neighbors must be a positive integer.")
  all_genes <- axes$gene_ids
  selected <- .resolve_graph_subset(gene_subset,all_genes)
  eligible <- .graph_gene_eligibility(pv,vkey,all_genes,!is.null(gene_subset))
  genes <- all_genes[selected & eligible]
  if (length(genes)<2L) stop("Velocity graph requires at least two eligible genes.")
  V <- t(as.matrix(.subset_canonical_matrix(pv@layers[[vkey]],vkey,genes,
    axes$cell_ids,all_genes,axes$cell_ids)))
  if (!is.null(pv@moments[[xkey]])) {
    X <- as.matrix(.subset_canonical_matrix(pv@moments[[xkey]],xkey,axes$cell_ids,
      genes,axes$cell_ids,all_genes))
  } else {
    X <- t(as.matrix(.subset_canonical_matrix(pv@layers$spliced,"spliced",genes,
      axes$cell_ids,all_genes,axes$cell_ids)))
  }
  usable <- colSums(!is.finite(V))==0 & colSums(!is.finite(X))==0 & colSums(abs(V))>0
  V <- V[,usable,drop=FALSE]
  X <- X[,usable,drop=FALSE]
  if (ncol(V)<2L || nrow(V)<2L) stop("Velocity graph requires at least two usable genes and cells.")
  dist_mat <- .subset_canonical_matrix(pv@graphs$neighbors$distances,"neighbour distances",
    axes$cell_ids,axes$cell_ids)
  values <- if (inherits(dist_mat,"sparseMatrix")) dist_mat@x else as.numeric(dist_mat)
  if (any(!is.finite(values)) || any(values<0)) stop("Neighbour distances must be finite and non-negative.")

  # ----------------build velocity graph----------------

  if (verbose)
    message("Computing velocity graph (cosine similarities)...")

  vkey_mode <- pv@misc[[paste0(vkey, "_params")]]$mode

  vgraph <- VelocityGraph$new(
    X                   = X,
    V                   = V,
    dist_mat            = dist_mat,
    n_recurse_neighbors = n_recurse_neighbors,
    sqrt_transform      = sqrt_transform,
    approx              = approx,
    vkey_mode           = vkey_mode
  )

  if (!any(is.finite(vgraph$V) & vgraph$V!=0)) stop("No cells have a computable velocity direction.")
  vgraph$compute_cosines(verbose = verbose)
  if (!any(vgraph$graph!=0) && !any(vgraph$graph_neg!=0))
    stop("No informative cosine similarities could be computed.")

  # ----------------write results back to plantvelo----------------

  pv <- .invalidate_velocity_dependents(pv,vkey,graph_changed=TRUE)
  pv@meta.data <- .align_cell_metadata(pv@meta.data,axes$cell_ids)
  dimnames(vgraph$graph) <- dimnames(vgraph$graph_neg) <- list(axes$cell_ids,axes$cell_ids)
  pv@graphs[[paste0(vkey, "_graph")]]     <- vgraph$graph
  pv@graphs[[paste0(vkey, "_graph_neg")]] <- vgraph$graph_neg
  pv@meta.data[[paste0(vkey, "_self_transition")]] <- vgraph$self_prob

  if (verbose)
    message(sprintf(
      "Velocity graph stored in @graphs$%s_graph (%d x %d sparse).",
      vkey, nrow(vgraph$graph), ncol(vgraph$graph)
    ))

  return(pv)
}


# =============================================================================
# VelocityGraph R6 class
# =============================================================================

#' VelocityGraph R6 class
#'
#' Internal class that stores velocity graph data and performs cosine
#' similarity computations between velocity vectors and neighbour-state
#' displacement vectors.
#'
#' @keywords internal
VelocityGraph <- R6::R6Class(
  "VelocityGraph",
  public = list(
    #' @field X Expression reference matrix (\eqn{n\_cells \times n\_features}).
    X               = NULL,
    #' @field V Velocity matrix after optional preprocessing and row centering.
    V               = NULL,
    #' @field V_raw Velocity matrix before row centering.
    V_raw           = NULL,
    #' @field indices Neighbour index matrix (\eqn{n\_cells \times k}).
    indices         = NULL,
    #' @field graph Sparse graph of positive cosine similarities.
    graph           = NULL,
    #' @field graph_neg Sparse graph of negative cosine similarities.
    graph_neg       = NULL,
    #' @field self_prob Per-cell self-transition probabilities.
    self_prob       = NULL,
    #' @field n_recurse_neighbors Neighbour recursion depth actually used.
    n_recurse_neighbors = NULL,
    #' @field sqrt_transform Logical. Whether square-root transformation is applied.
    sqrt_transform  = FALSE,

    #' @description Initialise a velocity graph object
    #'
    #' @param X Expression reference matrix (cells × genes), already filtered to
    #'   the selected gene subset.
    #' @param V Velocity matrix (cells × genes), already filtered to the
    #'   selected gene subset.
    #' @param dist_mat Sparse cell-cell distance matrix
    #'   (\eqn{n\_cells \times n\_cells}).
    #' @param n_recurse_neighbors Integer. Neighbour recursion depth. Default
    #'   \code{2L}.
    #' @param sqrt_transform Logical or \code{NULL}. Whether to apply square-root
    #'   transformation.
    #' @param approx Logical or \code{NULL}. Whether to use PCA approximation.
    #' @param vkey_mode Character or \code{NULL}. Velocity mode read from
    #'   \code{@misc[[paste0(vkey, "_params")]]$mode}.
    initialize = function(X, V,
                          dist_mat,
                          n_recurse_neighbors = 2L,
                          sqrt_transform      = NULL,
                          approx              = NULL,
                          vkey_mode           = NULL) {

      # Remove genes that are entirely zero or NA in the velocity matrix
      nans <- apply(V, 2, function(col) all(is.na(col) | col == 0))
      X    <- X[, !nans, drop = FALSE]
      V    <- V[, !nans, drop = FALSE]

      n_genes <- ncol(X)
      n_cells <- nrow(X)

      # ----------------optional PCA approximation----------------
      use_approx <- isTRUE(approx) || (is.null(approx) && n_genes > 100)

      if (use_approx && n_genes > 1) {
        n_pcs   <- max(min(30L, n_genes - 1L, n_cells - 1L), 1L)
        pca_res <- stats::prcomp(X, rank. = n_pcs, scale. = FALSE, center = TRUE)
        PCs     <- pca_res$rotation

        self$X  <- pca_res$x
        V_center <- sweep(V, 2, colMeans(V, na.rm = TRUE), "-")
        self$V   <- V_center %*% PCs
        # Cells with all-zero original velocity remain zero after projection
        self$V[rowSums(abs(V), na.rm = TRUE) == 0, ] <- 0
      } else {
        self$X <- X
        self$V <- V
      }
      self$V_raw <- self$V

      # ----------------optional square-root transformation----------------
      if (is.null(sqrt_transform)) {
        sqrt_transform <- !is.null(vkey_mode) && vkey_mode == "stochastic"
      }
      self$sqrt_transform <- sqrt_transform
      if (sqrt_transform) {
        self$V <- sign(self$V) * sqrt(abs(self$V))
      }

      # ----------------row-center velocity matrix----------------
      # Critical step: remove the global mean direction for each cell.
      row_means <- rowMeans(self$V, na.rm = TRUE)
      self$V    <- sweep(self$V, 1, row_means, "-")

      # ----------------build neighbour index matrix----------------
      self$n_recurse_neighbors <- n_recurse_neighbors %||% 2L
      self$indices <- private$.get_knn_indices(dist_mat)
    },

    #' @description Compute cosine similarities cell by cell and build sparse graphs
    #' @param verbose Logical. Print progress messages.
    compute_cosines = function(verbose = TRUE) {
      n_obs <- nrow(self$X)

      all_vals <- list()
      all_rows <- list()
      all_cols <- list()

      progress_step <- max(floor(n_obs / 10), 1L)

      for (obs_id in seq_len(n_obs)) {
        if (verbose && obs_id %% progress_step == 0) {
          message(sprintf("  %.0f%% done...", obs_id / n_obs * 100))
        }

        # Skip cells with zero velocity
        if (max(abs(self$V[obs_id, ]), na.rm = TRUE) == 0) next

        neighs_idx <- private$.get_iterative_neighbors(
          obs_id, self$n_recurse_neighbors
        )
        if (length(neighs_idx) == 0) next

        # State transition vectors: ΔX = X[j] - X[i]
        dX <- self$X[neighs_idx, , drop = FALSE] -
          matrix(self$X[obs_id, ],
                 nrow = length(neighs_idx), ncol = ncol(self$X),
                 byrow = TRUE)

        if (self$sqrt_transform) {
          dX <- sign(dX) * sqrt(abs(dX))
        }

        val <- .cosine_correlation(dX, self$V[obs_id, ])
        val[is.na(val)] <- 0

        all_vals <- c(all_vals, list(val))
        all_rows <- c(all_rows, list(rep(obs_id, length(neighs_idx))))
        all_cols <- c(all_cols, list(neighs_idx))
      }

      # ----------------return empty graphs if no valid cosines are computed----------------
      if (length(all_vals) == 0) {
        warning("No cosine similarities computed. Check velocity layer.")
        empty <- Matrix::sparseMatrix(i = 1L, j = 1L, x = 0,
                                      dims = c(n_obs, n_obs))
        self$graph     <- empty
        self$graph_neg <- empty
        self$self_prob <- rep(1, n_obs)
        return(invisible(self))
      }

      vals <- unlist(all_vals)
      rows <- unlist(all_rows)
      cols <- unlist(all_cols)
      vals[is.na(vals)] <- 0

      # ----------------build positive and negative sparse graphs----------------
      graph_all <- Matrix::sparseMatrix(
        i = rows, j = cols, x = vals,
        dims = c(n_obs, n_obs)
      )

      # Positive graph (π > 0)
      self$graph      <- graph_all
      self$graph@x[self$graph@x < 0] <- 0
      self$graph      <- Matrix::drop0(self$graph)

      # Negative graph (π < 0, retaining negative values)
      graph_neg       <- graph_all
      graph_neg@x[graph_neg@x > 0] <- 0
      self$graph_neg  <- Matrix::drop0(graph_neg)

      # Self-transition probability = max_conf_98 - cell_conf
      # Weaker outgoing confidence implies stronger self-loop probability.
      confidence     <- apply(self$graph, 1, max)
      p98            <- quantile(confidence, 0.98, na.rm = TRUE)
      self$self_prob <- pmax(p98 - confidence, 0)

      invisible(self)
    }
  ),

  private = list(

    # Extract actual stored KNN edges without densifying the distance matrix
    .get_knn_indices = function(dist_mat) {
      dist_sparse <- Matrix::Matrix(dist_mat, sparse = TRUE)

      if (!methods::is(dist_sparse, "dMatrix")) {
        stop("`dist_mat` must be a numeric cell distance matrix.")
      }

      if (
        methods::is(dist_sparse, "symmetricMatrix") ||
        methods::is(dist_sparse, "triangularMatrix")
      ) {
        dist_sparse <- methods::as(dist_sparse, "generalMatrix")
      }

      dist_row <- methods::as(dist_sparse, "RsparseMatrix")

      n_cells <- nrow(self$X)
      if (nrow(dist_row) != n_cells || ncol(dist_row) != n_cells) {
        stop(sprintf(
          "`dist_mat` must be a %d x %d cell distance matrix; got %d x %d.",
          n_cells, n_cells, nrow(dist_row), ncol(dist_row)
        ))
      }

      neighbor_list <- vector("list", n_cells)

      for (cell_id in seq_len(n_cells)) {
        start <- dist_row@p[cell_id] + 1L
        end   <- dist_row@p[cell_id + 1L]

        if (start > end) {
          neighbor_list[[cell_id]] <- integer()
          next
        }

        positions <- start:end
        neighbors <- dist_row@j[positions] + 1L
        distances <- dist_row@x[positions]

        keep <- (
          neighbors != cell_id &
          is.finite(distances) &
          distances >= 0
        )

        neighbors <- neighbors[keep]
        distances <- distances[keep]
        ordering  <- order(distances, neighbors)

        neighbor_list[[cell_id]] <- neighbors[ordering]
      }

      max_neighbors <- max(lengths(neighbor_list), 0L)

      if (max_neighbors == 0L) {
        stop(
          "Neighbor distance graph has no stored off-diagonal edges."
        )
      }

      indices <- matrix(
        0L,
        nrow = n_cells,
        ncol = max_neighbors
      )

      for (cell_id in seq_len(n_cells)) {
        neighbors <- neighbor_list[[cell_id]]

        if (length(neighbors) > 0L) {
          indices[cell_id, seq_along(neighbors)] <- neighbors
        }
      }

      indices
    },

    # Expand neighbours recursively (first-order or second-order)
    .get_iterative_neighbors = function(obs_id, n_recurse) {
      neighs <- self$indices[obs_id, ]
      neighs <- neighs[neighs > 0]

      if (n_recurse > 1 && length(neighs) > 0) {
        second_order <- unique(as.vector(self$indices[neighs, , drop = FALSE]))
        second_order <- second_order[second_order > 0 & second_order != obs_id]
        neighs <- unique(c(neighs, second_order))
      }

      neighs[neighs != obs_id]
    }
  )
)


# =============================================================================
# Markov transition matrix
# =============================================================================

#' @title Build Velocity Transition Matrix
#' @description
#' Construct the sparse velocity transition kernel used by scVelo. Positive
#' similarities are transformed with \code{expm1(scale * graph)}. Self-transition
#' confidence is estimated from the 98th percentile of outgoing similarities.
#' Negative similarities can either contribute signed reverse weights or
#' non-negative transition probabilities. Rows are normalised by their absolute
#' weight sums; signed rows therefore need not sum algebraically to one.
#'
#' @param pv A \code{plantvelo} object with
#'   \code{@graphs[[paste0(vkey, "_graph")]]} populated by
#'   \code{compute_velocity_graph()}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#' @param scale Numeric. Exponential kernel scale. Larger values yield sharper
#'   transitions. Default \code{10}.
#' @param self_transitions Logical. Estimate self-transition confidence before
#'   applying the exponential kernel. Default \code{TRUE}.
#' @param use_negative_cosines Logical. Give negative similarities signed reverse
#'   weights. If \code{FALSE}, use scVelo's non-negative probability branch.
#'   Default \code{FALSE}.
#' @param weight_diffusion Numeric in \eqn{[0,1]} or \code{NULL}. Weight used to
#'   mix an embedding-space Gaussian diffusion kernel with the velocity kernel.
#'   Only used when \code{basis} is supplied. Default \code{NULL} (zero weight).
#' @param backward Logical. Transpose the complete unnormalised kernel before
#'   row normalisation to obtain reverse transitions. Default \code{FALSE}.
#' @param basis Character scalar or \code{NULL}. Optional embedding used for the
#'   Gaussian distance kernel, following scVelo's \code{basis} argument.
#'   Default \code{NULL}, which leaves diffusion disabled.
#' @param scale_diffusion Positive numeric. Gaussian bandwidth multiplier,
#'   relative to the mean distance of positive transitions. Default \code{1}.
#'
#' @return A sparse transition matrix normalised by absolute row sums
#'   (\code{dgCMatrix}, \eqn{n\_cells \times n\_cells}).
#'   Empty rows remain zero. With \code{use_negative_cosines = TRUE}, the matrix
#'   can contain negative entries and should be used for velocity projection,
#'   rather than sampling transition probabilities.
#'
#' @examples
#' \dontrun{
#' T_mat <- build_transition_matrix(pv)
#' }
#' @export
build_transition_matrix <- function(pv,
                                    vkey                 = "velocity",
                                    scale                = 10,
                                    self_transitions     = TRUE,
                                    use_negative_cosines = FALSE,
                                    weight_diffusion     = NULL,
                                    backward             = FALSE,
                                    basis                = NULL,
                                    scale_diffusion      = 1) {

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  gkey <- paste0(vkey, "_graph")

  T_raw <- pv@graphs[[gkey]]

  if (is.null(T_raw))
    stop(sprintf(
      "Velocity graph '%s' not found. Run compute_velocity_graph() first.", gkey
    ))

  scalar_finite <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x)
  }
  if (!scalar_finite(scale)) stop("`scale` must be a finite numeric scalar.")
  for (flag in list(self_transitions, use_negative_cosines, backward)) {
    if (!is.logical(flag) || length(flag) != 1L || is.na(flag))
      stop("Transition controls must be non-missing logical scalars.")
  }
  if (!is.null(weight_diffusion) &&
      (!scalar_finite(weight_diffusion) || weight_diffusion < 0 || weight_diffusion > 1))
    stop("`weight_diffusion` must be NULL or a number between zero and one.")
  if (!scalar_finite(scale_diffusion) || scale_diffusion <= 0)
    stop("`scale_diffusion` must be a positive finite scalar.")

  cell_order <- colnames(pv@layers[["spliced"]])
  graph <- .as_velocity_sparse(T_raw, cell_order)
  if (self_transitions && nrow(graph) > 0L) {
    confidence <- numeric(nrow(graph))
    if (length(graph@x) > 0L) {
      row_max <- tapply(graph@x, graph@i + 1L, max)
      confidence[as.integer(names(row_max))] <- pmax(0, row_max)
    }
    upper <- stats::quantile(confidence, 0.98, names = FALSE)
    diag(graph) <- pmin(1, pmax(0, upper - confidence))
  }
  graph@x <- expm1(graph@x * scale)
  T_mat <- Matrix::drop0(graph)

  graph_neg <- pv@graphs[[paste0(vkey, "_graph_neg")]]
  if (!is.null(graph_neg)) {
    graph_neg <- .as_velocity_sparse(graph_neg, rownames(graph) %||% cell_order)
    if (!identical(dim(graph_neg), dim(T_mat)))
      stop("Positive and negative velocity graphs must have the same dimensions.")
    if (use_negative_cosines) {
      graph_neg@x <- expm1(-graph_neg@x * scale)
      T_mat <- T_mat - graph_neg
    } else {
      graph_neg@x <- expm1(graph_neg@x * scale)
      T_mat <- Matrix::drop0(T_mat + graph_neg)
      # scVelo adds one only to entries in the combined sparse support.
      T_mat@x <- T_mat@x + 1
    }
  }

  if (backward) T_mat <- Matrix::t(T_mat)
  T_mat <- .normalize_velocity_transition(T_mat)

  if (!is.null(basis)) {
    embedding <- .get_embedding(pv, basis)
    if (nrow(embedding) != nrow(T_mat))
      stop("Diffusion embedding and velocity graph must have the same number of cells.")
    graph_cells <- rownames(T_mat) %||% colnames(pv@layers[["spliced"]])
    if (!is.null(graph_cells) && !is.null(rownames(embedding))) {
      if (anyDuplicated(rownames(embedding)) ||
          !setequal(graph_cells, rownames(embedding)))
        stop("Diffusion embedding and velocity graph identify different cells.")
      embedding <- embedding[graph_cells, , drop = FALSE]
    }
    if (any(!is.finite(embedding))) stop("Diffusion embedding must be finite.")
    edges <- Matrix::summary(Matrix::drop0(T_mat > 0))
    if (nrow(edges) > 0L) {
      distances <- sqrt(rowSums((embedding[edges$i, , drop = FALSE] -
                                embedding[edges$j, , drop = FALSE])^2))
      bandwidth <- scale_diffusion * mean(distances)
      gaussian <- function(multiplier) {
        weights <- if (bandwidth > 0) {
          exp(-0.5 * (distances / (bandwidth * multiplier))^2)
        } else {
          rep(1, length(distances))
        }
        Matrix::sparseMatrix(i = edges$i, j = edges$j, x = weights,
                             dims = dim(T_mat), dimnames = dimnames(T_mat))
      }
      T_mat <- T_mat * gaussian(1)
      if (!is.null(weight_diffusion) && weight_diffusion > 0 && weight_diffusion < 1)
        T_mat <- (1 - weight_diffusion) * T_mat + weight_diffusion * gaussian(0.5)
      T_mat <- .normalize_velocity_transition(T_mat)
    }
  }
  .as_velocity_sparse(T_mat)
}


.as_velocity_sparse <- function(graph, cell_order = NULL) {
  graph <- methods::as(Matrix::Matrix(graph, sparse = TRUE), "generalMatrix")
  graph <- methods::as(methods::as(graph, "dMatrix"), "CsparseMatrix")
  if (nrow(graph) != ncol(graph) || any(!is.finite(graph@x)))
    stop("Velocity graphs must be square numeric matrices with finite entries.")
  if (!is.null(rownames(graph)) && !is.null(colnames(graph))) {
    if (anyDuplicated(rownames(graph)) || anyDuplicated(colnames(graph)) ||
        !setequal(rownames(graph), colnames(graph)))
      stop("Velocity graph rows and columns must identify the same unique cells.")
    target <- cell_order %||% rownames(graph)
    if (anyDuplicated(target) || !setequal(target, rownames(graph)))
      stop("Velocity graphs and expression layers must identify the same cells.")
    if (!identical(rownames(graph), target) || !identical(colnames(graph), target))
      graph <- graph[target, target, drop = FALSE]
  }
  graph
}


.normalize_velocity_transition <- function(kernel) {
  totals <- Matrix::rowSums(abs(kernel))
  totals[totals == 0] <- 1
  Matrix::drop0(kernel / totals)
}


#' @title Simulate Cell Trajectory by Random Walk
#' @description
#' Simulate a directed random walk on the velocity transition matrix, starting
#' from a specified cell.
#'
#' At each step, the next cell is sampled from the row of the transition matrix
#' corresponding to the current cell.
#'
#' @param pv A \code{plantvelo} object with a velocity graph already computed.
#' @param starting_cell Integer. Starting cell index (1-based).
#' @param n_steps Integer. Number of random-walk steps. Default \code{10}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#' @param random_state Integer or \code{NULL}. Random seed. Default
#'   \code{NULL}.
#'
#' @return Integer vector of length \code{n_steps + 1} containing the visited
#'   cell indices along the random-walk path.
#'
#' @examples
#' \dontrun{
#' path <- get_cell_transitions(pv, starting_cell = 1, n_steps = 20)
#' }
#' @export
get_cell_transitions <- function(pv,
                                 starting_cell,
                                 n_steps      = 10,
                                 vkey         = "velocity",
                                 random_state = NULL) {

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (!is.null(random_state)) set.seed(random_state)

  T_mat <- build_transition_matrix(pv, vkey = vkey, self_transitions = FALSE)
  path  <- integer(n_steps + 1)
  path[1] <- starting_cell

  for (step in seq_len(n_steps)) {
    curr  <- path[step]
    probs <- as.vector(T_mat[curr, ])
    probs[is.na(probs) | probs < 0] <- 0
    total <- sum(probs)
    if (total == 0) {
      path[(step + 1):(n_steps + 1)] <- curr
      break
    }
    path[step + 1] <- sample.int(length(probs), 1L,
                                 prob = probs / total)
  }

  path
}


# =============================================================================
# Internal helper functions
# =============================================================================

#' Compute cosine similarity row-wise
#'
#' Computes row-wise cosine similarities between a matrix of neighbour
#' displacement vectors \code{dX} and a single velocity vector \code{v}.
#'
#' @param dX Numeric matrix of neighbour displacement vectors
#'   (\eqn{n\_neighbors \times n\_genes}).
#' @param v Numeric vector of length \eqn{n\_genes}.
#'
#' @return Numeric vector of cosine similarities, one per row of \code{dX}.
#' @keywords internal
.cosine_correlation <- function(dX, v) {
  norm_v <- sqrt(sum(v^2, na.rm = TRUE))
  if (norm_v < 1e-10) return(rep(0, nrow(dX)))

  dot     <- as.vector(dX %*% v)
  norm_dX <- sqrt(rowSums(dX^2, na.rm = TRUE))
  denom   <- norm_dX * norm_v
  denom[denom < 1e-10] <- 1

  dot / denom
}

#' Apply softmax scaling and row normalisation
#'
#' Exponentiates the edge weights of a sparse graph, optionally adds
#' self-transitions on the diagonal, and row-normalises the matrix.
#'
#' @param T_raw Sparse matrix of unnormalised transition scores.
#' @param scale Numeric. Softmax scale parameter. Default \code{10}.
#' @param self_transitions Logical. Add self-loops before normalisation.
#'
#' @return Row-normalised sparse transition matrix.
#' @keywords internal
.softmax_transition <- function(T_raw, scale = 10, self_transitions = FALSE) {
  T_mat   <- T_raw
  T_mat@x <- exp(T_mat@x * scale)
  T_mat   <- Matrix::drop0(T_mat)

  if (self_transitions) {
    T_mat <- T_mat + Matrix::Diagonal(nrow(T_mat), rep(1, nrow(T_mat)))
  }

  rs <- Matrix::rowSums(T_mat)
  rs[rs == 0] <- 1
  T_mat / rs
}

# =============================================================================
# Graph gene selection and eligibility
# =============================================================================

.resolve_graph_subset <- function(subset, genes) {
  if (is.null(subset)) return(rep(TRUE,length(genes)))
  if (anyNA(subset)) stop("gene_subset cannot contain missing/NA values.")
  if (is.character(subset)) {
    if (any(!subset %in% genes)) stop("gene_subset contains unknown genes.")
    return(genes %in% unique(subset))
  }
  if (!is.logical(subset) || !is.null(dim(subset)))
    stop("gene_subset must be a character or logical vector.")
  if (is.null(names(subset))) {
    if (length(subset)!=length(genes)) stop("gene_subset has incompatible length.")
    return(subset)
  }
  if (any(!nzchar(names(subset))) || anyDuplicated(names(subset)) ||
      any(!names(subset) %in% genes)) stop("gene_subset has invalid or unknown names.")
  idx <- match(genes,names(subset))
  out <- rep(FALSE,length(genes))
  out[!is.na(idx)] <- subset[idx[!is.na(idx)]]
  out
}

.graph_gene_eligibility <- function(pv, vkey, genes, explicit) {
  mask <- if (identical(vkey,"velocity")) pv@velocity$genes else NULL
  if (!is.null(mask)) {
    if (!is.null(names(mask)) && any(!names(mask) %in% genes))
      stop("Velocity gene mask has unknown names.")
    return(.driver_velocity_mask(mask,genes))
  }
  run <- pv@misc[[paste0(vkey,"_params")]]
  if (identical(run$mode,"dynamical")) {
    params <- .validate_dynamics_schema(pv)
    cutoff <- if ("min_likelihood" %in% names(run)) run$min_likelihood else .001
    mask <- .vel_genes_filter(params,.valid_velocity_fits(params),cutoff)
    return(.driver_velocity_mask(mask,genes))
  }
  if (!explicit) stop("No verifiable eligibility mask; provide an explicit gene_subset.")
  rep(TRUE,length(genes))
}
