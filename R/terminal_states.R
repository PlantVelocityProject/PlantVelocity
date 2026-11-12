# =============================================================================
# PlantVelocity: terminal_states.R
# Terminal-state identification by Markov-chain eigendecomposition
# Relative near-unit spectral support and explicit-target absorption.
# Scores and kernels differ from previous releases and are not scVelo parity.
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
#' Results are relative spectral scores in \eqn{[0, 1]}, not calibrated
#' probabilities. Isolated, unselected or constant-score cells receive NA.
#' The forward and backward kernels are shared with build_transition_matrix().
#' Strongly connected blocks are solved separately and their left eigenvectors
#' are lifted through downstream blocks. Repeated stationary directions and
#' near-unit slow modes are scored using their complete subspace. Sparse
#' failures or incomplete spectra error; no large dense fallback is used.
#'
#' @section Output written to \code{@meta.data}:
#' \describe{
#'   \item{\code{root_cells}}{Numeric vector. Relative root-cell score per cell.}
#'   \item{\code{end_points}}{Numeric vector. Relative end-point score per cell.}
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
#'   within-component, within-group normalisation. Cross-group transition edges
#'   are retained. \code{NULL} analyses all cells jointly. Default
#'   \code{NULL}.
#' @param groups Character vector or \code{NULL}. Subset of group labels when
#'   \code{groupby} is set. \code{NULL} uses all groups. Default
#'   \code{NULL}.
#' @param self_transitions Logical. Add self-loops using the shared transition
#'   kernel. Default \code{FALSE}.
#' @param eps Numeric. Eigenvalue selection threshold
#'   (\eqn{|\lambda-1| \leq \mathrm{eps}}), in (0,1). Default \code{0.001}.
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
compute_terminal_states <- function(pv, vkey="velocity", modality="Ms",
                                    groupby=NULL, groups=NULL, self_transitions=FALSE,
                                    eps=1e-3, random_state=0L, verbose=TRUE) {
  axes <- .pv_axes(pv); .validate_eps(eps)
  .validate_whole(random_state,"random_state",0L)
  pv@meta.data <- .align_cell_metadata(pv@meta.data,axes$cells)
  if (!is.null(groups) && is.null(groupby)) stop("groups requires groupby.")
  labels <- rep("all",length(axes$cells))
  if (!is.null(groupby)) {
    if (!is.character(groupby) || length(groupby)!=1L || is.na(groupby) ||
        is.null(pv@meta.data[[groupby]])) stop("groupby key is missing or invalid.")
    labels <- as.character(pv@meta.data[[groupby]])
    if (!is.null(groups) && (!is.character(groups) || anyNA(groups) || any(!groups %in% labels)))
      stop("groups contains unknown group labels.")
  }
  C <- .get_connectivities(pv)
  forward <- .canonical_graph(build_transition_matrix(pv,vkey=vkey,
    self_transitions=self_transitions),"forward transition",axes$cells)
  backward <- .canonical_graph(build_transition_matrix(pv,vkey=vkey,
    self_transitions=self_transitions,backward=TRUE),"backward transition",axes$cells)
  roots <- .terminal_scores(backward,C,eps,random_state)
  ends <- .terminal_scores(forward,C,eps,random_state)
  component <- .graph_components(forward,strong=FALSE)
  selected <- !is.na(labels) & (is.null(groups) | labels %in% groups)
  selected <- which(selected)
  normalize <- function(values) {
    result <- rep(NA_real_,length(values))
    keys <- interaction(component,labels,drop=TRUE)
    for (idx in split(selected,keys[selected])) {
      idx <- idx[is.finite(values[idx])]
      result[idx] <- .normalize_terminal_score(values[idx])
    }
    result
  }
  roots <- normalize(roots); ends <- normalize(ends)
  if (any(is.finite(roots))) roots <- verify_roots(pv,roots,modality=modality)
  if (anyNA(roots) || anyNA(ends))
    message("Terminal scores are NA for isolated, unselected or uninformative cells/groups.")
  pv <- .clear_terminal_dependents(pv)
  pv@meta.data$root_cells <- roots
  pv@meta.data$end_points <- ends
  if (verbose) message("Updated root_cells and end_points (relative spectral scores).")
  pv
}


#' @title Eigendecomposition of Transition Matrix
#' @description
#' Compute leading eigenvectors of a row-normalised Markov transition matrix
#' \eqn{T}. Stationary-like states are identified by selecting eigenvectors
#' corresponding to eigenvalues satisfying
#' \eqn{|\lambda-1| \leq \varepsilon} with negligible imaginary part.
#'
#' Uses base \code{eigen()} for at most 128 cells and \pkg{RSpectra} for
#' larger kernels. Larger failures do not fall back to dense decomposition.
#' Input must be finite, non-negative and row stochastic. Raw eigenvectors
#' are checked for residual error before display processing. RNG is preserved.
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
#'     to \eqn{[0, 1]} when perc is not NULL. These processed columns are no
#'     longer raw eigenvectors.}
#' }
#' @export
compute_eigs <- function(T, k=10L, eps=1e-3, perc=c(2,98), random_state=0L) {
  T <- .validate_markov(T)
  k <- .validate_whole(k,"k",1L); .validate_eps(eps)
  if (!is.null(perc) && (!is.numeric(perc) || length(perc)!=2L ||
      any(!is.finite(perc)) || perc[1]<0 || perc[2]>100 || perc[1]>perc[2]))
    stop("perc must be NULL or two ordered percentiles in [0,100].")
  restore_seed <- .local_seed(random_state); on.exit(restore_seed(),add=TRUE)
  request <- if (nrow(T)<=128L) min(k,nrow(T)) else min(k,nrow(T)-2L)
  raw <- .raw_eigs(T,request)
  raw$values <- head(raw$values,request)
  raw$vectors <- raw$vectors[,seq_len(length(raw$values)),drop=FALSE]
  res <- .select_unit_eigs(T,raw,eps)
  res$vectors <- abs(res$vectors)
  if (!is.null(perc) && ncol(res$vectors)) for (j in seq_len(ncol(res$vectors))) {
    x <- res$vectors[,j]
    bound <- stats::quantile(x,perc/100,names=FALSE)
    x[x<bound[1]] <- 0; x <- pmin(x,bound[2])
    if (max(x)>0) x <- x/max(x)
    res$vectors[,j] <- x
  }
  res
}


#' @title Verify Root Cell Identification
#' @description
#' Apply an expression/sampling-time heuristic to relative root-cell scores.
#' Axes are aligned by cell and gene names. Missing scores are preserved;
#' non-finite expression or constant plasticity skips the optional correction.
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
#' @param roots Numeric vector. Relative root-cell scores for the cells indexed
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
verify_roots <- function(pv, roots, cell_idx=NULL, modality="Ms", time_key=NULL) {
  axes <- .pv_axes(pv)
  meta <- .align_cell_metadata(pv@meta.data,axes$cells)
  if (is.null(cell_idx)) cell_idx <- seq_along(axes$cells)
  if (!is.numeric(cell_idx) || any(!is.finite(cell_idx)) || any(cell_idx!=floor(cell_idx)) ||
      any(cell_idx<1 | cell_idx>length(axes$cells)) || anyDuplicated(cell_idx) ||
      !is.numeric(roots) || length(roots)!=length(cell_idx) || any(is.infinite(roots)) ||
      any(roots<0 | roots>1,na.rm=TRUE)) stop("Invalid root scores or cell_idx.")
  valid <- is.finite(roots)
  if (!any(valid)) return(roots)
  if (!is.null(time_key)) {
    if (!is.character(time_key) || length(time_key)!=1L || is.na(time_key) || is.null(meta[[time_key]]))
      stop("time_key must identify an existing metadata column.")
    time <- meta[[time_key]][cell_idx]
    if (!is.numeric(time) || !any(is.finite(time))) {
      message("verify_roots: no finite sampling times; correction skipped.")
      return(roots)
    }
    early <- valid & is.finite(time) & time==min(time[is.finite(time)])
    if (any(early) && mean(roots[early]>.5)<.2)
      message("verify_roots: roots are not enriched at the earliest sampling time.")
    return(roots)
  }
  X <- if (!is.null(pv@moments[[modality]]))
    .canonical_matrix(pv@moments[[modality]],modality,axes$cells,axes$genes) else
    Matrix::t(.canonical_matrix(pv@layers$spliced,"spliced",axes$genes,axes$cells))
  values <- if (inherits(X,"sparseMatrix")) X@x else as.numeric(X)
  if (any(!is.finite(values))) {
    message("verify_roots: non-finite expression; correction skipped."); return(roots)
  }
  gcc <- Matrix::colMeans(X)
  params <- pv@kinetics$params
  if (!is.null(params) && "gene_count_corr" %in% names(params)) {
    if (is.null(rownames(params)) || anyDuplicated(rownames(params)) ||
        any(!rownames(params) %in% axes$genes)) stop("gene_count_corr has invalid gene names.")
    gcc <- numeric(length(axes$genes)); gcc[match(rownames(params),axes$genes)] <- params$gene_count_corr
  }
  if (any(!is.finite(gcc)) || !any(gcc!=0)) return(roots)
  plasticity <- as.numeric(X[cell_idx,,drop=FALSE]%*%gcc)
  if (diff(range(plasticity[valid]))<=1e-12) {
    message("verify_roots: no expression variation; correction skipped."); return(roots)
  }
  plasticity <- (plasticity-min(plasticity[valid]))/diff(range(plasticity[valid]))
  high <- valid & plasticity>.5; low <- valid & !high
  if (!any(high) || !any(low)) return(roots)
  root_high <- roots>.9
  right <- mean(root_high[high]); wrong <- mean(root_high[low]); baseline <- mean(root_high[valid])
  if (right>3*baseline) roots[valid] <- roots[valid]*as.numeric(high[valid])
  else if (wrong>baseline || right<baseline) message("verify_roots: root identification needs biological verification.")
  roots
}


#' @title Compute Target-Conditioned Cell Fate
#' @description
#' Assign cells to specified terminal groups using absorbing Markov-chain
#' probabilities. Solve (I-Q)B=R with sparse systems; never invert the full
#' transition matrix. Targets are absorbing for this calculation only.
#'
#' @param pv A plantvelo object with a velocity graph.
#' @param groupby Metadata key used to group selected target cells. Default
#'   "clusters". Target labels must be non-missing and non-empty.
#' @param disconnected_groups Character vector of excluded group labels, or
#'   NULL. These cells receive NA. Edges into them consume probability mass;
#'   remaining edges are not renormalized. Cannot contain target cells.
#' @param self_transitions Logical. Add self-transitions before normalization.
#' @param verbose Logical. Print progress messages.
#' @param vkey Velocity graph key. Default "velocity".
#' @param target_cells Required target cell names or a complete logical mask.
#'   Named masks must cover every cell exactly once. NULL errors; terminal
#'   scores are not automatically thresholded to choose targets.
#' @param tol Numeric tolerance in (0,1) for probability and solver checks.
#'   Default 1e-8. This is a numerical control, not statistical confidence.
#' @param max_iter Positive integer. Iteration limit if sparse factorization
#'   fails. Default 10000. Iteration also checks remaining transient mass;
#'   failure to converge errors without updating results.
#'
#' @details
#' Probabilities are summed within each target group before label selection.
#' Confidence is the largest unconditioned target-group absorption probability.
#' Non-target closed classes and excluded cells can consume probability mass,
#' so target probabilities may sum to less than one. They are not rescaled.
#' Unreachable/excluded cells have NA label and confidence. Tied groups (within
#' tol) have an NA label and retain the largest probability as confidence.
#' Target cells have their own group label and confidence one. No probability
#' matrix, target mask or additional QC fields are persisted.
#'
#' @return A plantvelo object updating only cell_fate and cell_fate_confidence
#'   in metadata (with metadata aligned to the spliced cell axis).
#' @examples
#' \dontrun{
#' pv <- compute_cell_fate(pv, target_cells = terminal_ids, groupby = "clusters")
#' }
#' @export
compute_cell_fate <- function(pv, groupby="clusters", disconnected_groups=NULL,
                              self_transitions=FALSE, verbose=TRUE,
                              vkey="velocity", target_cells=NULL, tol=1e-8,
                              max_iter=10000L) {
  .compute_absorption(pv,groupby,disconnected_groups,self_transitions,verbose,
    vkey,target_cells,tol,max_iter,backward=FALSE)
}


#' @title Compute Target-Conditioned Cell Origin
#' @description
#' Assign cells to specified root groups using absorption on the backward
#' graph. The kernel is transposed and row-normalized by build_transition_matrix;
#' this is a reverse-graph origin estimate, not a calibrated biological ancestry
#' probability. All target, exclusion, tie and convergence rules match fate.
#' @inheritParams compute_cell_fate
#' @return A plantvelo object updating cell_origin and cell_origin_confidence
#'   in metadata. No full probability matrix is stored.
#' @examples
#' \dontrun{
#' pv <- compute_cell_origin(pv, target_cells = root_ids, groupby = "clusters")
#' }
#' @export
compute_cell_origin <- function(pv, groupby="clusters", disconnected_groups=NULL,
                                self_transitions=FALSE, verbose=TRUE,
                                vkey="velocity", target_cells=NULL, tol=1e-8,
                                max_iter=10000L) {
  .compute_absorption(pv,groupby,disconnected_groups,self_transitions,verbose,
    vkey,target_cells,tol,max_iter,backward=TRUE)
}

# =============================================================================
# Sparse graph components
# =============================================================================

# Iterative Kosaraju traversal; no recursion or dense adjacency matrices.
.graph_components <- function(T,strong=TRUE) {
  n <- nrow(T)
  if (!n) return(integer())
  e <- Matrix::summary(Matrix::drop0(methods::as(T,"generalMatrix")))
  e <- e[e$x>0 & e$i!=e$j,,drop=FALSE]
  forward <- split(e$j,factor(e$i,levels=seq_len(n)))
  reverse <- split(e$i,factor(e$j,levels=seq_len(n)))
  seen <- logical(n); finish <- integer(n); done <- 0L
  stack <- integer(n); cursor <- integer(n)
  if (!strong) {
    forward <- lapply(seq_len(n),function(i) c(forward[[i]],reverse[[i]]))
  }
  component <- integer(n); number <- 0L
  for (root in seq_len(n)) {
    if (seen[root]) next
    number <- number+1L; depth <- 1L; stack[depth] <- root; seen[root] <- TRUE
    component[root] <- number
    while (depth>0L) {
      node <- stack[depth]
      if (cursor[node]<length(forward[[node]])) {
        cursor[node] <- cursor[node]+1L
        child <- forward[[node]][cursor[node]]
        if (!seen[child]) {
          depth <- depth+1L; stack[depth] <- child; seen[child] <- TRUE
          component[child] <- number
        }
      } else {
        done <- done+1L; finish[done] <- node; depth <- depth-1L
      }
    }
  }
  if (!strong) return(component)
  component[] <- 0L; number <- 0L
  for (root in rev(finish)) {
    if (component[root]!=0L) next
    number <- number+1L; depth <- 1L; stack[depth] <- root; component[root] <- number
    while (depth>0L) {
      node <- stack[depth]; depth <- depth-1L
      for (child in reverse[[node]]) if (component[child]==0L) {
        component[child] <- number; depth <- depth+1L; stack[depth] <- child
      }
    }
  }
  component
}

.closed_class_cells <- function(T,external=rep(FALSE,nrow(T))) {
  component <- .graph_components(T)
  if (!length(component)) return(logical())
  closed <- rep(TRUE,max(component))
  closed[component[external]] <- FALSE
  e <- Matrix::summary(Matrix::drop0(T))
  cross <- e$x>0 & component[e$i]!=component[e$j]
  closed[component[e$i[cross]]] <- FALSE
  closed[component]
}

# =============================================================================
# Checked terminal spectrum and subspace scores
# =============================================================================

.validate_markov <- function(T,allow_zero=FALSE,tol=1e-8) {
  if (length(dim(T))!=2L || nrow(T)!=ncol(T) || !nrow(T))
    stop("Transition kernel must be a non-empty square matrix.")
  if (!((is.matrix(T) && is.numeric(T)) || methods::is(T,"dMatrix") || methods::is(T,"iMatrix")))
    stop("Transition kernel must be numeric.")
  .nonnegative_values(T,"Transition kernel")
  totals <- Matrix::rowSums(T)
  if (any(abs(totals-1)>tol & !(allow_zero & totals==0)))
    stop("Transition kernel must be row stochastic (or have explicit zero rows).")
  methods::as(Matrix::Matrix(T,sparse=TRUE),"generalMatrix")
}

.validate_eps <- function(eps) {
  if (!.finite_scalar(eps) || eps<=0 || eps>=1) stop("eps must lie in (0,1).")
}

.raw_eigs <- function(T,k) {
  n <- nrow(T)
  if (n<=128L) {
    res <- eigen(as.matrix(Matrix::t(T)))
  } else {
    if (k>n-2L) stop("Cannot resolve complete near-unit spectrum within the sparse eigenvalue budget.")
    res <- tryCatch(RSpectra::eigs(Matrix::t(T),k=k,which="LR",
      opts=list(tol=1e-10,maxitr=10000,ncv=min(n,max(2L*k+1L,64L)))),error=function(e)
        stop("Sparse eigendecomposition failed: ",conditionMessage(e)))
    if (length(res$values)<k || (!is.null(res$nconv) && res$nconv<k))
      stop("Sparse eigendecomposition did not converge for all requested vectors.")
  }
  ord <- order(Re(res$values),decreasing=TRUE)
  list(values=res$values[ord],vectors=res$vectors[,ord,drop=FALSE])
}

.select_unit_eigs <- function(T,res,eps) {
  valid <- abs(res$values-1)<=eps & abs(Im(res$values))<=1e-8
  values <- Re(res$values[valid]); vectors <- Re(res$vectors[,valid,drop=FALSE])
  if (any(!is.finite(values)) || any(!is.finite(vectors))) stop("Non-finite eigensystem.")
  if (ncol(vectors)) for (j in seq_len(ncol(vectors))) {
    v <- vectors[,j]
    residual <- sqrt(sum((as.numeric(Matrix::t(T)%*%v)-values[j]*v)^2)) /
      max(sqrt(sum(v^2)),.Machine$double.eps)
    if (!is.finite(residual) || residual>1e-6) stop("Eigenvector residual exceeds 1e-6.")
  }
  rownames(vectors) <- rownames(T)
  list(values=values,vectors=vectors)
}

.block_unit_eigs <- function(T,eps,count) {
  n <- nrow(T)
  if (n<=128L) return(.select_unit_eigs(T,.raw_eigs(T,n),eps))
  budget <- max(10L,count)
  repeat {
    requested <- min(budget+1L,n-2L)
    res <- .raw_eigs(T,requested)
    space <- .select_unit_eigs(T,res,eps)
    if (ncol(space$vectors)<count) stop("Sparse spectrum missed repeated stationary eigenvalues.")
    # LR orders by real part. Rejected complex candidates can precede valid
    # real ones; only a real-part boundary certifies the unrequested spectrum.
    if (min(Re(res$values))<1-eps) return(space)
    if (budget>=128L || requested==n-2L)
      stop("Cannot confirm complete near-unit eigenspace within the 128-column budget.")
    budget <- min(128L,budget*2L)
  }
}

.valid_terminal_lift <- function(shifted,solution,rhs) {
  if (inherits(solution,"error") || any(!is.finite(solution))) return(FALSE)
  residual <- as.matrix(shifted%*%solution)-rhs
  all(is.finite(residual)) && all(sqrt(colSums(residual^2))<=
    1e-8*pmax(1,sqrt(colSums(rhs^2))))
}

.complete_unit_space <- function(T,eps) {
  n <- nrow(T)
  component <- .graph_components(T)
  closed <- .closed_class_cells(T)
  count <- length(unique(component[closed]))
  if (count>128L) stop("More than 128 closed classes; terminal spectrum exceeds the per-component budget.")
  if (n<=128L) return(.block_unit_eigs(T,eps,count)$vectors)
  if (max(component)==1L) return(.block_unit_eigs(T,eps,count)$vectors)

  # Kosaraju labels condensation components in source-to-sink order. The
  # spectrum of this block triangular matrix is the union of block spectra.
  blocks <- split(seq_len(n),component)
  spectra <- vector("list",length(blocks)); values <- numeric()
  columns <- vector("list",length(blocks))
  for (b in seq_along(blocks)) {
    idx <- blocks[[b]]; Q <- T[idx,idx,drop=FALSE]
    # Infinity-norm bound: this entire block lies outside the selected disk.
    if (max(Matrix::rowSums(Q))<1-eps) next
    spectra[[b]] <- .block_unit_eigs(Q,eps,as.integer(all(closed[idx])))
    k <- length(spectra[[b]]$values)
    if (length(values)+k>128L)
      stop("Complete near-unit eigenspace exceeds the 128-column budget.")
    columns[[b]] <- length(values)+seq_len(k)
    values <- c(values,spectra[[b]]$values)
  }
  vectors <- matrix(0,n,length(values),dimnames=list(rownames(T),NULL))
  edges <- Matrix::summary(Matrix::drop0(T))
  cross <- component[edges$i]!=component[edges$j] & edges$x>0
  edges <- edges[cross,,drop=FALSE]
  if (any(component[edges$i]>component[edges$j]))
    stop("Internal terminal component ordering is invalid.")
  predecessors <- split(edges$i,factor(component[edges$j],levels=seq_along(blocks)))
  for (b in seq_along(blocks)) {
    idx <- blocks[[b]]; own <- columns[[b]]
    if (length(own)) vectors[idx,own] <- spectra[[b]]$vectors
    previous <- unique(predecessors[[b]])
    if (!length(previous) || !length(values)) next
    rhs <- as.matrix(Matrix::t(T[previous,idx,drop=FALSE]) %*%
      vectors[previous,,drop=FALSE])
    needed <- which(colSums(abs(rhs))>0)
    if (!length(needed)) next
    Q <- Matrix::t(T[idx,idx,drop=FALSE])
    # One sparse factorization per distinct eigenvalue in this block; this
    # retains transient near-unit modes instead of using stationary modes only.
    for (cols in split(needed,match(values[needed],unique(values[needed])))) {
      lambda <- values[cols[1L]]
      shifted <- lambda*Matrix::Diagonal(length(idx))-Q
      target <- rhs[,cols,drop=FALSE]
      solution <- tryCatch(Matrix::solve(shifted,target),error=identity)
      same <- which(abs(spectra[[b]]$values-lambda)<=1e-10)
      if (!.valid_terminal_lift(shifted,solution,target) && length(same)) {
        # Repeated slow eigenvalues can give a compatible singular equation.
        # A bordered solve fixes the free local nullspace coordinates, already
        # represented by this block's independent eigenvectors. Incompatible
        # forcing (a defective repeated mode) must still fail explicitly.
        local <- qr(spectra[[b]]$vectors[,same,drop=FALSE],tol=1e-8)
        if (local$rank!=length(same)) stop("Terminal eigenspace is numerically rank deficient.")
        basis <- Matrix::Matrix(qr.Q(local),sparse=TRUE)
        system <- rbind(cbind(shifted,basis),cbind(Matrix::t(basis),
          Matrix::sparseMatrix(i=integer(),j=integer(),dims=rep(ncol(basis),2))))
        target <- rbind(target,matrix(0,ncol(basis),length(cols)))
        solution <- tryCatch(Matrix::solve(system,target),error=identity)
        if (!inherits(solution,"error"))
          solution <- as.matrix(solution)[seq_along(idx),,drop=FALSE]
      }
      if (!.valid_terminal_lift(shifted,solution,rhs[,cols,drop=FALSE]))
        stop("Cannot lift terminal eigenspace: incompatible repeated eigenvalue or excessive residual.")
      vectors[idx,cols] <- as.matrix(solution)
    }
  }
  # Lifting must satisfy the original full operator, not just diagonal blocks.
  .select_unit_eigs(T,list(values=values,vectors=vectors),eps)$vectors
}

.subspace_score <- function(vectors) {
  if (!ncol(vectors)) return(rep(NA_real_,nrow(vectors)))
  qr_space <- qr(vectors,tol=1e-8,LAPACK=FALSE)
  if (qr_space$rank!=ncol(vectors)) stop("Terminal eigenspace is numerically rank deficient.")
  U <- qr.Q(qr_space)
  sqrt(rowSums(U^2))
}

.terminal_scores <- function(T,C,eps,seed) {
  restore_seed <- .local_seed(seed); on.exit(restore_seed(),add=TRUE)
  T <- .validate_markov(T,allow_zero=TRUE)
  off <- T; diag(off) <- 0; off <- Matrix::drop0(off)
  active <- Matrix::rowSums(off)+Matrix::colSums(off)>0
  out <- rep(NA_real_,nrow(T))
  if (!any(active)) stop("Transition graph has no informative non-self edges.")
  empty <- which(Matrix::rowSums(T)==0)
  if (length(empty)) T <- T+Matrix::sparseMatrix(i=empty,j=empty,x=1,dims=dim(T))
  component <- .graph_components(T,strong=FALSE)
  for (idx in split(which(active),component[active])) {
    sub <- T[idx,idx,drop=FALSE]
    score <- .subspace_score(.complete_unit_space(sub,eps))
    if (any(!is.finite(score))) stop("No valid near-unit eigenspace for terminal analysis.")
    W <- .normalize_smoothing(C[idx,idx,drop=FALSE])
    out[idx] <- as.numeric(W%*%score)
  }
  out
}

.normalize_terminal_score <- function(x) {
  if (!length(x)) return(x)
  if (diff(range(x))<=1e-12*max(1,max(abs(x)))) return(rep(NA_real_,length(x)))
  x <- (x-min(x))/diff(range(x))
  x <- pmin(x,stats::quantile(x,.98,names=FALSE))
  if (diff(range(x))<=1e-12) return(rep(NA_real_,length(x)))
  (x-min(x))/diff(range(x))
}

# =============================================================================
# Explicit-target absorption helpers
# =============================================================================

# Target-conditioned absorption. All matrices and component labels are transient.
.resolve_target_cells <- function(target_cells,cells) {
  x <- target_cells
  if (is.null(x)) stop("target_cells is required; specify terminal/root cell names or a complete logical mask.")
  if (!is.null(dim(x)) || anyNA(x)) stop("target_cells must be cell names or a non-missing logical mask.")
  if (is.character(x)) {
    if (!length(x) || anyDuplicated(x) || any(!x %in% cells))
      stop("target_cells contains empty, duplicate or unknown cells.")
    return(cells %in% x)
  }
  if (!is.logical(x) || length(x)!=length(cells))
    stop("target_cells must be cell names or a complete logical mask.")
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x)) || !setequal(names(x),cells))
      stop("target_cells mask must name all cells exactly once.")
    x <- x[match(cells,names(x))]
  }
  if (!any(x)) stop("target_cells selects no cells.")
  unname(x)
}

.target_reachable <- function(T,targets,excluded) {
  e <- Matrix::summary(Matrix::drop0(T))
  e <- e[e$x>0 & !excluded[e$i] & !excluded[e$j],,drop=FALSE]
  incoming <- split(e$i,factor(e$j,levels=seq_len(nrow(T))))
  reachable <- targets
  queue <- integer(nrow(T)); initial <- which(targets)
  queue[seq_along(initial)] <- initial; last <- length(initial); first <- 1L
  while (first<=last) {
    node <- queue[first]; first <- first+1L
    for (parent in incoming[[node]]) if (!reachable[parent]) {
      reachable[parent] <- TRUE; last <- last+1L; queue[last] <- parent
    }
  }
  reachable
}

.valid_absorption_solution <- function(B,Q,R,tol) {
  if (any(!is.finite(B)) || any(B < -tol | B > 1+tol)) return(FALSE)
  residual <- max(abs(as.matrix(Q%*%B)+R-B)) / max(1,max(abs(R)))
  is.finite(residual) && residual<=tol
}

.iterate_absorption <- function(Q,R,tol,max_iter) {
  B <- matrix(0,nrow(R),ncol(R)); survival <- rep(1,nrow(Q))
  for (step in seq_len(max_iter)) {
    next_B <- as.matrix(Q%*%B)+R
    survival <- as.numeric(Q%*%survival)
    change <- max(abs(next_B-B)); B <- next_B
    if (max(survival)<=tol && change<=tol && .valid_absorption_solution(B,Q,R,tol))
      return(B)
  }
  stop("Absorption iteration did not converge within max_iter; no results were written.")
}

.absorption_assignment <- function(T,labels,targets,excluded,tol,max_iter) {
  T <- .validate_markov(T,allow_zero=TRUE,tol=tol)
  n <- nrow(T)
  answer <- rep(NA_character_,n); confidence <- rep(NA_real_,n)
  answer[targets] <- labels[targets]; confidence[targets] <- 1
  group_names <- sort(unique(labels[targets]))
  remaining <- which(!targets & !excluded)
  trapped <- logical(n)
  if (length(remaining)) {
    # Detect boundary escape with strictly positive edges, not a probability cutoff.
    edges <- Matrix::summary(Matrix::drop0(T))
    escape <- logical(n)
    escape[edges$i[edges$x>0 & (targets[edges$j] | excluded[edges$j])]] <- TRUE
    trapped[remaining] <- .closed_class_cells(T[remaining,remaining,drop=FALSE],escape[remaining])
  }
  reachable <- .target_reachable(T,targets,excluded)
  transient <- which(!targets & !excluded & !trapped & reachable)
  if (!length(transient)) return(list(label=answer,confidence=confidence))
  component <- .graph_components(T,strong=FALSE)
  for (idx in split(transient,component[transient])) {
    terminal <- which(targets & component==component[idx[1L]])
    Q <- T[idx,idx,drop=FALSE]
    A <- Matrix::Diagonal(nrow(Q))-Q
    factorization <- tryCatch(withCallingHandlers(Matrix::lu(A),
      warning=function(w) stop(conditionMessage(w))),error=function(e) NULL)
    boundary <- T[idx,terminal,drop=FALSE]
    candidate <- sort(unique(labels[terminal]))
    best <- numeric(length(idx)); second <- rep(-Inf,length(idx))
    winning <- rep(NA_character_,length(idx)); total <- numeric(length(idx))
    for (start in seq.int(1L,length(candidate),by=16L)) {
      groups <- candidate[start:min(length(candidate),start+15L)]
      membership <- Matrix::sparseMatrix(i=which(labels[terminal] %in% groups),
        j=match(labels[terminal][labels[terminal] %in% groups],groups),x=1,
        dims=c(length(terminal),length(groups)))
      R <- as.matrix(boundary%*%membership)
      B <- if (is.null(factorization)) NULL else tryCatch(
        withCallingHandlers(as.matrix(Matrix::solve(factorization,R)),
          warning=function(w) stop(conditionMessage(w))),error=function(e) NULL)
      if (is.null(B) || !.valid_absorption_solution(B,Q,R,tol))
        B <- .iterate_absorption(Q,R,tol,max_iter)
      B <- pmin(pmax(B,0),1)
      total <- total+rowSums(B)
      for (j in seq_along(groups)) {
        p <- B[,j]
        better <- is.na(winning) | p>best
        second <- ifelse(better,best,pmax(second,p))
        # With the first candidate there is no competing target group yet.
        second[is.na(winning)] <- -Inf
        best[better] <- p[better]; winning[better] <- groups[j]
      }
    }
    if (any(!is.finite(total)) || any(total>1+tol))
      stop("Absorption probabilities exceed total mass one; no results were written.")
    supported <- best>0
    confidence[idx[supported]] <- best[supported]
    unique_best <- supported & best-second>tol
    answer[idx[unique_best]] <- winning[unique_best]
  }
  list(label=answer,confidence=confidence)
}

.compute_absorption <- function(pv,groupby,disconnected_groups,self_transitions,
                                verbose,vkey,target_cells,tol,max_iter,backward) {
  axes <- .pv_axes(pv)
  if (!is.character(groupby) || length(groupby)!=1L || is.na(groupby) || !nzchar(groupby))
    stop("groupby must be a non-empty metadata key.")
  if (!is.character(vkey) || length(vkey)!=1L || is.na(vkey) || !nzchar(vkey)) stop("vkey must be a non-empty string.")
  if (!.finite_scalar(tol) || tol<=0 || tol>=1) stop("tol must lie in (0,1).")
  max_iter <- .validate_whole(max_iter,"max_iter",1L)
  if (!is.logical(verbose) || length(verbose)!=1L || is.na(verbose)) stop("verbose must be a logical scalar.")
  meta <- .align_cell_metadata(pv@meta.data,axes$cells)
  labels <- meta[[groupby]]
  if (is.null(labels) || !is.atomic(labels) || !is.null(dim(labels))) stop("groupby metadata labels are missing or invalid.")
  labels <- as.character(labels)
  targets <- .resolve_target_cells(target_cells,axes$cells)
  if (anyNA(labels[targets]) || any(!nzchar(labels[targets]))) stop("Target cells must have non-missing, non-empty group labels.")
  if (!is.null(disconnected_groups) && (!is.character(disconnected_groups) ||
      anyNA(disconnected_groups) || any(!disconnected_groups %in% labels)))
    stop("disconnected_groups contains unknown group labels.")
  excluded <- labels %in% disconnected_groups
  if (any(targets & excluded)) stop("Target cells overlap excluded disconnected_groups.")
  T <- .canonical_graph(build_transition_matrix(pv,vkey=vkey,
    self_transitions=self_transitions,backward=backward),"transition",axes$cells)
  result <- .absorption_assignment(T,labels,targets,excluded,tol,max_iter)
  key <- if (backward) "cell_origin" else "cell_fate"
  meta[[key]] <- result$label
  meta[[paste0(key,"_confidence")]] <- result$confidence
  pv@meta.data <- meta
  if (verbose) message("Updated ",key," with target-group absorption assignments.")
  pv
}

# =============================================================================
# Terminal result invalidation
# =============================================================================

.clear_terminal_dependents <- function(pv) {
  fields <- c("latent_time","cell_fate","cell_fate_confidence","cell_origin","cell_origin_confidence")
  pv@meta.data[intersect(fields,names(pv@meta.data))] <- NULL
  .clear_velocity_gene_results(pv)
}
