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
#' @param vkey Character scalar. Reserved velocity-layer identifier. Gene
#'   eligibility uses fitted dynamics independently of velocity execution order.
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
#' @param min_gene_coverage Numeric in \eqn{(0,1]}. Minimum fraction of cells
#'   with finite non-negative fitted time per gene. Default \code{0.5}.
#' @param min_cell_coverage Numeric in \eqn{(0,1]}. Minimum fraction of eligible
#'   genes with supported time per cell, with at least two genes required.
#'   Default \code{0.5}. Both coverage defaults are adjustable screening choices.
#'
#' @details
#' Require valid fitted parameters, finite likelihood and at least two genes
#' with sufficient coverage and nonzero time range. No fallback reinstates
#' rejected genes. Missing times remain missing when cell support is insufficient.
#' Explicit roots must be valid indices or probabilities; automatic roots use
#' supported cells and mean gene-normalized times, breaking ties by cell name.
#' Neighbour correction replaces only low-confidence supported targets with
#' weighted times from reliable sources. Cells without reliable neighbours retain
#' their estimate. No QC fields are stored. Errors leave previous results intact.
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
                                vkey = "velocity",
                                min_likelihood = 0.1,
                                min_confidence = 0.75,
                                root_key = NULL,
                                end_key = NULL,
                                t_max = NULL,
                                verbose = TRUE,
                                min_gene_coverage = 0.5,
                                min_cell_coverage = 0.5) {
  if (!inherits(pv, "plantvelo")) stop("`pv` must be a plantvelo object.")
  .validate_likelihood_threshold(min_likelihood)
  for (key in c("min_gene_coverage", "min_cell_coverage")) {
    x <- get(key)
    if (!.finite_scalar(x) || x <= 0 || x > 1) stop(paste(key, "must lie in (0,1]."))
  }
  if (!is.null(min_confidence) && (!.finite_scalar(min_confidence) ||
      min_confidence < 0 || min_confidence > 1)) stop("min_confidence must lie in [0,1] or be NULL.")
  if (!is.null(t_max) && (!.finite_scalar(t_max) || t_max <= 0)) stop("t_max must be finite and positive.")
  for (key in c("vkey", "root_key", "end_key")) {
    x <- get(key)
    if ((key=="vkey" || !is.null(x)) &&
        (!is.character(x) || length(x)!=1L || is.na(x) || !nzchar(x)))
      stop(paste(key, "must be a nonempty string."))
  }
  if (!is.logical(verbose) || length(verbose)!=1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.")
  axes <- .velocity_axes(pv)
  params <- .validate_dynamics_schema(pv)
  ll <- params$fit_likelihood
  valid <- .valid_velocity_fits(params) & is.finite(ll)
  if (!is.null(min_likelihood)) valid <- valid & ll >= min_likelihood
  genes <- axes$gene_ids[valid]
  if (length(genes)<2L) stop("latent_time requires at least two eligible genes; no fallback is used.")
  T <- t(as.matrix(.subset_canonical_matrix(pv@layers$fit_t, "fit_t", genes,
    axes$cell_ids, axes$gene_ids, axes$cell_ids)))
  if (any(is.infinite(T)) || any(T < 0, na.rm=TRUE)) stop("Fitted time must be non-negative and finite or NA.")
  valid <- colSums(is.finite(T)) >= max(2L,ceiling(nrow(T)*min_gene_coverage))
  varying <- vapply(seq_len(ncol(T)), function(j) {
    x <- T[is.finite(T[,j]),j]
    length(x)>1L && diff(range(x))>0
  }, logical(1))
  T <- T[,valid & varying,drop=FALSE]
  if (ncol(T)<2L) stop("latent_time requires at least two genes with sufficient time coverage and range.")
  pv@meta.data <- .align_cell_metadata(pv@meta.data, axes$cell_ids)
  supported <- .time_row_support(T, min_cell_coverage)
  roots <- .find_root_cell(pv,root_key,T,min_cell_coverage)
  merge_roots <- function(indices, explicit=FALSE) {
    times <- list()
    skipped <- 0L
    for (root in indices) {
      if (!supported[root] || sum(is.finite(T[root,]))<2L) { skipped <- skipped+1L; next }
      rooted <- root_time(T,root)$t_rooted
      needed <- max(2L,ceiling(ncol(rooted)*min_cell_coverage))
      time <- compute_shared_time(rooted,min_valid_genes=needed)
      time[!supported] <- NA_real_
      times[[length(times)+1L]] <- time
    }
    if (skipped && explicit && length(times)) warning("Some explicit roots lack supported time and were skipped.")
    if (!length(times)) return(NULL)
    out <- rowMeans(do.call(cbind,times),na.rm=TRUE)
    out[!is.finite(out)] <- NA_real_
    if (length(times)>1L) out <- .scale_finite_time(out)
    out
  }
  latent <- merge_roots(roots,!is.null(root_key))
  if (is.null(latent)) stop("No root has enough supported time observations.")
  if (!is.null(end_key)) {
    if (is.null(pv@meta.data[[end_key]])) stop("end_key not found in cell metadata.")
    probability <- .root_probabilities(pv@meta.data[[end_key]],nrow(T))
    ends <- which(probability>.5 & supported)
    ends <- head(ends[order(-probability[ends],axes$cell_ids[ends])],4L)
    if (length(ends)) {
      end_time <- merge_roots(ends)
      if (!is.null(end_time)) {
        end_time <- .scale_finite_time(1-end_time)
        usable <- is.finite(latent) & is.finite(end_time)
        latent[usable] <- latent[usable] + .2*end_time[usable]
        latent <- .scale_finite_time(latent)
      }
    }
  }
  conn <- pv@graphs$neighbors$connectivities
  if (!is.null(min_confidence) && !is.null(conn)) {
    C <- .subset_canonical_matrix(conn,"neighbour connectivities",axes$cell_ids,axes$cell_ids)
    latent <- .correct_latent_time(latent,C,min_confidence)
  }
  latent <- .scale_finite_time(latent,strict=TRUE)
  if (!is.null(t_max)) latent <- latent*t_max
  pv@meta.data$latent_time <- unname(latent)
  if (verbose) message(sprintf("Added latent_time using %d eligible genes.",ncol(T)))
  pv
}

#' Identify Supported Root-Cell Indices
#' @description Use explicit integer indices or aligned probabilities first,
#'   then common root fields, then supported gene-normalized mean time.
#' @param pv A plantvelo object with cell-aligned metadata.
#' @param root_key Explicit root key or NULL.
#' @param T_mat Cells-by-genes fitted-time matrix, after gene filtering.
#' @param min_cell_coverage Minimum fraction of finite gene times per root cell.
#' @return Integer cell indices; invalid explicit priors raise an error.
#' @keywords internal
.find_root_cell <- function(pv, root_key, T_mat, min_cell_coverage=0.5) {
  cells <- rownames(T_mat)
  if (is.null(cells)) cells <- as.character(seq_len(nrow(T_mat)))
  supported <- .time_row_support(T_mat,min_cell_coverage)
  from_probability <- function(x, restrict=FALSE) {
    x <- .root_probabilities(x,nrow(T_mat))
    idx <- which(x>0 & (!restrict | supported))
    if (!length(idx)) return(integer())
    idx <- idx[x[idx]>=.9*max(x[idx])]
    head(idx[order(-x[idx],cells[idx])],4L)
  }
  if (!is.null(root_key)) {
    x <- pv@misc[[root_key]]
    if (!is.null(x)) {
      if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || any(!is.finite(x)) ||
          any(x!=floor(x) | x<1 | x>nrow(T_mat))) stop("root indices must be valid integer cell indices.")
      return(head(unique(as.integer(x)),4L))
    }
    if (is.null(pv@meta.data[[root_key]])) stop("Explicit root_key is missing.")
    roots <- from_probability(pv@meta.data[[root_key]])
    if (!length(roots)) stop("Explicit root probabilities have no positive support.")
    return(roots)
  }
  for (key in c("root_cells","starting_cells","root_states_probs")) {
    if (!is.null(pv@meta.data[[key]])) {
      if (all(is.na(pv@meta.data[[key]])))
        stop("Automatic root scores are all missing; specify a custom root_key or recompute terminal states.")
      roots <- from_probability(pv@meta.data[[key]],TRUE)
      if (length(roots)) return(roots)
    }
  }
  if (!any(supported)) stop("No supported cells are available for root selection.")
  scaled <- apply(T_mat,2,.scale_finite_time)
  score <- rowMeans(scaled,na.rm=TRUE)
  idx <- which(supported & is.finite(score))
  unname(idx[order(score[idx],cells[idx])][1L])
}

# =============================================================================
# Latent-time coverage, priors and neighbor correction
# =============================================================================

.scale_finite_time <- function(x, strict=FALSE) {
  valid <- is.finite(x)
  out <- rep(NA_real_,length(x))
  if (!any(valid)) {
    if (strict) stop("No supported latent time observations remain.")
    return(out)
  }
  limits <- range(x[valid])
  if (limits[2] <= limits[1]) {
    if (strict) stop("Latent time has no identifiable time range.")
    out[valid] <- 0
  } else out[valid] <- (x[valid]-limits[1])/(limits[2]-limits[1])
  out
}

.time_row_support <- function(T, coverage) {
  rowSums(is.finite(T)) >= max(2L, ceiling(ncol(T)*coverage))
}

.root_probabilities <- function(x, n) {
  if (!is.numeric(x) || !is.null(dim(x)) || length(x)!=n ||
      any(is.infinite(x)) || any(x < 0 | x > 1, na.rm=TRUE))
    stop("Root/end probabilities must match cells and lie in [0,1].")
  x[is.na(x)] <- 0
  x
}

.correct_latent_time <- function(time, C, min_confidence) {
  if (length(dim(C))!=2L || any(dim(C)!=length(time)))
    stop("Neighbour matrix has incompatible cell dimensions.")
  values <- if (inherits(C,"sparseMatrix")) C@x else as.numeric(C)
  if (any(!is.finite(values)) || any(values < 0)) stop("Neighbour weights must be finite and non-negative.")
  valid <- is.finite(time)
  average <- function(source) {
    weights <- as.numeric(C %*% as.numeric(source))
    sums <- as.numeric(C %*% ifelse(source,time,0))
    out <- rep(NA_real_,length(time))
    out[weights>0] <- sums[weights>0]/weights[weights>0]
    out
  }
  tc <- average(valid)
  paired <- valid & is.finite(tc)
  if (!any(paired)) return(time)
  max_t <- max(time[valid])
  if (max_t<=0) return(time)
  z <- sum(time[paired]*tc[paired])/max(sum(tc[paired]^2),1e-10)
  conf <- rep(NA_real_,length(time))
  conf[paired] <- (1-abs(time[paired]/max_t-tc[paired]*z/max_t))^2
  low <- paired & !is.na(conf) & conf < min_confidence
  if (any(low)) {
    replacement <- average(valid & !low & is.finite(conf))
    replace <- low & is.finite(replacement)
    time[replace] <- replacement[replace]
  }
  time
}
