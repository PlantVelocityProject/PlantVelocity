# Name alignment and shared unsmoothed evidence for fitting and driver ranking.
# Keep count matrices sparse; smoothing is never used to count supporting cells.
.driver_align_matrix <- function(x, label, genes, cells,
                                 orientation = c("genes_by_cells", "cells_by_genes")) {
  orientation <- match.arg(orientation)
  if (is.null(x)) stop(sprintf("Required layer '%s' is missing.", label), call. = FALSE)
  if (length(dim(x)) != 2L ||
      !((is.matrix(x) && is.numeric(x)) || methods::is(x, "dMatrix") || methods::is(x, "iMatrix"))) {
    stop(sprintf("Layer '%s' must be a numeric matrix.", label), call. = FALSE)
  }
  index <- function(observed, expected, size, axis) {
    if (is.null(observed)) {
      if (size != length(expected))
        stop(sprintf("Layer '%s' has incompatible %s dimensions.", label, axis), call. = FALSE)
      return(seq_along(expected))
    }
    if (anyNA(observed) || any(!nzchar(observed)) || anyDuplicated(observed))
      stop(sprintf("Layer '%s' has invalid %s names.", label, axis), call. = FALSE)
    idx <- match(expected, observed)
    if (anyNA(idx)) stop(sprintf("Layer '%s' has missing %s names.", label, axis), call. = FALSE)
    idx
  }
  by_gene <- identical(orientation, "genes_by_cells")
  rows <- index(rownames(x), if (by_gene) genes else cells, nrow(x),
                if (by_gene) "gene" else "cell")
  cols <- index(colnames(x), if (by_gene) cells else genes, ncol(x),
                if (by_gene) "cell" else "gene")
  x <- x[rows, cols, drop = FALSE]
  dimnames(x) <- if (by_gene) list(genes, cells) else list(cells, genes)
  x
}

.raw_evidence_thresholds <- function(min_u_cells = 20L, min_s_cells = 20L,
                                     min_joint_cells = 10L,
                                     min_u_counts = 20, min_s_counts = 20) {
  thresholds <- as.list(environment())
  for (key in names(thresholds)) {
    x <- thresholds[[key]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 ||
        (grepl("cells$", key) && x != floor(x))) {
      stop(sprintf("`%s` must be a finite non-negative %s.", key,
        if (grepl("cells$", key)) "integer" else "number"), call. = FALSE)
    }
  }
  thresholds
}

# Work on one block at a time, including for dense input layers.
.raw_evidence_counts <- function(u, s) {
  for (key in c("u", "s")) {
    x <- get(key)
    values <- if (methods::is(x, "Matrix")) methods::slot(x, "x") else x
    if (any(!is.finite(values)))
      stop("Raw spliced/unspliced evidence must be finite.", call. = FALSE)
    if (any(values < 0))
      stop("Raw spliced/unspliced evidence must be non-negative.", call. = FALSE)
  }
  detected_u <- u > 0
  detected_s <- s > 0
  u_total <- as.numeric(Matrix::rowSums(u))
  s_total <- as.numeric(Matrix::rowSums(s))
  if (any(!is.finite(u_total)) || any(!is.finite(s_total)))
    stop("Raw spliced/unspliced totals must be finite.", call. = FALSE)
  list(
    raw_u_positive = as.numeric(Matrix::rowSums(detected_u)),
    raw_s_positive = as.numeric(Matrix::rowSums(detected_s)),
    raw_joint_positive = as.numeric(Matrix::rowSums(detected_u & detected_s)),
    raw_u_total = u_total,
    raw_s_total = s_total
  )
}

.raw_evidence_failures <- function(evidence, thresholds) {
  list(
    raw_u_cells = evidence$raw_u_positive < thresholds$min_u_cells,
    raw_s_cells = evidence$raw_s_positive < thresholds$min_s_cells,
    raw_joint_cells = evidence$raw_joint_positive < thresholds$min_joint_cells,
    raw_u_counts = evidence$raw_u_total < thresholds$min_u_counts,
    raw_s_counts = evidence$raw_s_total < thresholds$min_s_counts
  )
}

.driver_raw_evidence <- function(u, s) {
  keys <- c("raw_u_positive", "raw_s_positive", "raw_joint_positive",
            "raw_u_total", "raw_s_total")
  evidence <- stats::setNames(lapply(keys, function(key) numeric(nrow(u))), keys)
  if (nrow(u)) for (start in seq.int(1L, nrow(u), by = 256L)) {
    idx <- seq.int(start, min(start + 255L, nrow(u)))
    block <- .raw_evidence_counts(u[idx, , drop = FALSE], s[idx, , drop = FALSE])
    for (key in keys) evidence[[key]][idx] <- block[[key]]
  }
  as.data.frame(c(list(raw_n_cells = rep(ncol(u), nrow(u))), evidence), row.names = NULL)
}

# Return only a temporary named mask; write no QC fields to the pv object.
.raw_evidence_filter <- function(pv, genes, min_u_cells = 20L, min_s_cells = 20L,
                                 min_joint_cells = 10L,
                                 min_u_counts = 20, min_s_counts = 20) {
  thresholds <- .raw_evidence_thresholds(min_u_cells, min_s_cells,
    min_joint_cells, min_u_counts, min_s_counts)
  axes <- .velocity_axes(pv)
  aligned_block <- function(layer, label, selected) {
    if (!is.null(layer) && is.null(rownames(layer))) {
      if (length(dim(layer)) != 2L || nrow(layer) != length(axes$gene_ids))
        stop(sprintf("Layer '%s' has incompatible gene dimensions.", label), call. = FALSE)
      # An unnamed full axis follows canonical order, not the block's length.
      layer <- layer[match(selected, axes$gene_ids), , drop = FALSE]
      rownames(layer) <- selected
    }
    .driver_align_matrix(layer, label, selected, axes$cell_ids)
  }
  keep <- stats::setNames(logical(length(genes)), genes)
  if (length(genes)) for (start in seq.int(1L, length(genes), by = 256L)) {
    idx <- seq.int(start, min(start + 255L, length(genes)))
    u <- aligned_block(pv@layers$unspliced, "unspliced", genes[idx])
    s <- aligned_block(pv@layers$spliced, "spliced", genes[idx])
    evidence <- .raw_evidence_counts(u, s)
    keep[idx] <- !Reduce(`|`, .raw_evidence_failures(evidence, thresholds))
  }
  keep
}

.driver_velocity_mask <- function(mask, genes) {
  if (is.null(mask)) return(rep(TRUE, length(genes)))
  if (!is.logical(mask) || !is.null(dim(mask)))
    stop("@velocity$genes must be a logical vector.", call. = FALSE)
  if (!is.null(names(mask))) {
    if (anyNA(names(mask)) || any(!nzchar(names(mask))) || anyDuplicated(names(mask)))
      stop("@velocity$genes must have valid unique gene names.", call. = FALSE)
    mask <- mask[match(genes, names(mask))]
  } else if (length(mask) != length(genes)) {
    stop("Unnamed @velocity$genes must match the velocity gene dimension.", call. = FALSE)
  }
  !is.na(mask) & mask
}
