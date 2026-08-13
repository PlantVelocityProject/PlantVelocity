# =============================================================================
# PlantVelocity: ir_statistics.R
# Standalone IR-associated signal analysis: validate aligned IR, unspliced,
# and spliced layers; compute cell- and gene-level statistics; and provide
# post-hoc QC, cell-landscape, and gene-association plots.
# IR remains excluded from two-state dynamics and velocity inference.
#
# Core formulas:
#   library_total[i] = sum_g (IR[g,i] + U[g,i] + S[g,i])
#   IR_fraction[i]   = sum_g IR[g,i] / library_total[i]
#   IR_norm[g,i]     = log1p(scale_factor * IR[g,i] / library_total[i])
# Genes are retained by minimum detected-cell and total-count thresholds.
# =============================================================================

.validate_ir_layer <- function(layer, name) {
  if (!(is.matrix(layer) || inherits(layer, "Matrix"))) {
    stop(sprintf("`@layers$%s` must be a matrix-like object.", name))
  }

  values <- if (inherits(layer, "sparseMatrix")) {
    layer@x
  } else {
    as.numeric(layer)
  }
  if (any(!is.finite(values)) || any(values < 0)) {
    stop(sprintf(
      "`@layers$%s` must contain finite non-negative values.",
      name
    ))
  }

  invisible(layer)
}

.validate_ir_layers <- function(pv) {
  if (!inherits(pv, "plantvelo")) {
    stop("`pv` must be a plantvelo object.")
  }

  keys <- c("ir", "unspliced", "spliced")
  missing <- keys[vapply(
    keys,
    function(key) is.null(pv@layers[[key]]),
    logical(1)
  )]
  if (length(missing)) {
    stop(sprintf(
      "Missing required layer(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  layers <- lapply(keys, function(key) {
    layer <- pv@layers[[key]]
    .validate_ir_layer(layer, key)
    layer
  })
  names(layers) <- keys

  reference_dimnames <- dimnames(layers$ir)
  aligned <- vapply(
    layers[-1L],
    function(layer) identical(dimnames(layer), reference_dimnames),
    logical(1)
  )
  if (length(reference_dimnames) != 2L ||
      is.null(reference_dimnames[[1L]]) ||
      is.null(reference_dimnames[[2L]]) ||
      any(!aligned)) {
    stop(paste(
      "IR, unspliced, and spliced layers must have identical",
      "gene and cell names."
    ))
  }

  layers
}

.normalize_ir_sparse <- function(ir, library_total, scale_factor) {
  sparse_ir <- methods::as(Matrix::Matrix(ir, sparse = TRUE), "dgCMatrix")
  weights <- ifelse(library_total > 0, scale_factor / library_total, 0)
  normalized <- sparse_ir %*% Matrix::Diagonal(x = weights)
  normalized <- methods::as(normalized, "dgCMatrix")
  normalized@x <- log1p(normalized@x)
  dimnames(normalized) <- dimnames(ir)
  normalized
}

.validate_ir_stat_arguments <- function(min_cells,
                                        min_total_count,
                                        scale_factor,
                                        n_cells) {
  if (is.null(min_cells)) {
    min_cells <- max(10L, ceiling(0.01 * n_cells))
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1L ||
      !is.finite(min_cells) || min_cells < 1 ||
      min_cells != as.integer(min_cells)) {
    stop("`min_cells` must be a positive whole-number scalar or NULL.")
  }
  if (!is.numeric(min_total_count) || length(min_total_count) != 1L ||
      !is.finite(min_total_count) || min_total_count < 0) {
    stop("`min_total_count` must be a non-negative numeric scalar.")
  }
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be a positive numeric scalar.")
  }

  list(
    min_cells = as.integer(min_cells),
    min_total_count = as.numeric(min_total_count),
    scale_factor = as.numeric(scale_factor)
  )
}

#' Compute Standalone IR-Excluded Statistics
#'
#' Compute descriptive cell- and gene-level statistics from a gene-level IR
#' count layer. The calculation reads intron-excluded unspliced and spliced
#' counts only to construct library totals and fractions. It does not modify
#' the supplied object or any two-state dynamics result.
#'
#' @param pv A `plantvelo` object with aligned `ir`, `unspliced`, and `spliced`
#'   layers.
#' @param min_cells Minimum number of cells with a positive IR count required
#'   to retain a gene. `NULL` uses `max(10, ceiling(0.01 * n_cells))`.
#' @param min_total_count Minimum total IR count required to retain a gene.
#' @param scale_factor Positive library-size scale used for normalized IR
#'   abundance.
#'
#' @return A standalone `plantvelo_ir_stats` object containing cell metrics,
#'   gene metrics, a sparse normalized IR matrix, and aligned identifiers.
#'
#' @examples
#' \dontrun{
#' ir_stats <- compute_ir_stats(pv)
#' }
#'
#' @export
compute_ir_stats <- function(pv,
                             min_cells = NULL,
                             min_total_count = 20,
                             scale_factor = 1e4) {
  layers <- .validate_ir_layers(pv)
  ir <- layers$ir
  unspliced <- layers$unspliced
  spliced <- layers$spliced
  n_cells <- ncol(ir)
  if (n_cells < 1L) {
    stop("IR statistics require at least one cell.")
  }

  arguments <- .validate_ir_stat_arguments(
    min_cells,
    min_total_count,
    scale_factor,
    n_cells
  )

  ir_depth <- as.numeric(Matrix::colSums(ir))
  library_total <- ir_depth +
    as.numeric(Matrix::colSums(unspliced)) +
    as.numeric(Matrix::colSums(spliced))
  ir_fraction <- rep(NA_real_, n_cells)
  valid_library <- library_total > 0
  ir_fraction[valid_library] <-
    ir_depth[valid_library] / library_total[valid_library]
  ir_detected <- as.integer(Matrix::colSums(ir > 0))

  gene_total <- as.numeric(Matrix::rowSums(ir))
  gene_detected <- as.numeric(Matrix::rowSums(ir > 0))
  keep <- gene_detected >= arguments$min_cells &
    gene_total >= arguments$min_total_count

  retained_ir <- ir[keep, , drop = FALSE]
  normalized_ir <- .normalize_ir_sparse(
    retained_ir,
    library_total,
    arguments$scale_factor
  )
  burden <- if (nrow(normalized_ir)) {
    as.numeric(Matrix::colMeans(normalized_ir))
  } else {
    rep(NA_real_, n_cells)
  }
  burden[!valid_library] <- NA_real_

  retained_total <- gene_total[keep] +
    as.numeric(Matrix::rowSums(unspliced[keep, , drop = FALSE])) +
    as.numeric(Matrix::rowSums(spliced[keep, , drop = FALSE]))
  gene_fraction <- rep(NA_real_, sum(keep))
  valid_genes <- retained_total > 0
  gene_fraction[valid_genes] <-
    gene_total[keep][valid_genes] / retained_total[valid_genes]

  cell_metrics <- data.frame(
    cell = colnames(ir),
    ir_depth_cell = ir_depth,
    ir_fraction_cell = ir_fraction,
    ir_detected_genes = ir_detected,
    ir_burden_norm = burden,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  gene_metrics <- data.frame(
    gene = rownames(ir)[keep],
    ir_detection_rate_gene = gene_detected[keep] / n_cells,
    ir_fraction_gene = gene_fraction,
    ir_total_count = gene_total[keep],
    ir_mean_abundance_norm = if (nrow(normalized_ir)) {
      as.numeric(Matrix::rowMeans(normalized_ir))
    } else {
      numeric()
    },
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  structure(
    list(
      cell_metrics = cell_metrics,
      gene_metrics = gene_metrics,
      normalized_ir = normalized_ir,
      cells = colnames(ir),
      genes = rownames(ir)[keep]
    ),
    class = "plantvelo_ir_stats"
  )
}

# =============================================================================
# Shared IR plotting utilities
# =============================================================================

.require_ir_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for IR plots.")
  }
}

.validate_ir_stats <- function(ir_stats) {
  if (!inherits(ir_stats, "plantvelo_ir_stats")) {
    stop("`ir_stats` must be returned by compute_ir_stats().")
  }
  required <- c(
    "cell_metrics", "gene_metrics", "normalized_ir", "cells", "genes"
  )
  if (!identical(names(ir_stats), required)) {
    stop("`ir_stats` has an invalid schema.")
  }
  invisible(ir_stats)
}

.validate_ir_stats_for_pv <- function(pv, ir_stats) {
  .validate_ir_stats(ir_stats)
  layers <- .validate_ir_layers(pv)
  normalized_genes_compatible <- if (length(ir_stats$genes)) {
    identical(ir_stats$genes, rownames(ir_stats$normalized_ir))
  } else {
    nrow(ir_stats$normalized_ir) == 0L
  }
  compatible <-
    identical(ir_stats$cells, colnames(layers$ir)) &&
    identical(ir_stats$cells, ir_stats$cell_metrics$cell) &&
    identical(ir_stats$cells, colnames(ir_stats$normalized_ir)) &&
    identical(ir_stats$genes, ir_stats$gene_metrics$gene) &&
    normalized_genes_compatible &&
    all(ir_stats$genes %in% rownames(layers$ir))
  if (!compatible) {
    stop("`ir_stats` is not compatible with `pv`.")
  }
  invisible(layers)
}

.ir_theme <- function() {
  ggplot2::theme_classic() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 13,
      color = "black",
      hjust = 0.5
    ),
    axis.text = ggplot2::element_text(size = 10, color = "black"),
    axis.title = ggplot2::element_text(size = 12, color = "black"),
    legend.title = ggplot2::element_text(size = 10, color = "black"),
    legend.text = ggplot2::element_text(size = 9, color = "black"),
    legend.key.size = grid::unit(0.4, "cm")
    )
}

.combine_ir_panels <- function(panels, ncol) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(panels)
  }
  patchwork::wrap_plots(panels, ncol = ncol)
}

.add_ir_gene_labels <- function(plot, data, x, y, rank, n) {
  if (is.null(n) || n <= 0L || !nrow(data)) {
    return(plot)
  }
  ordered <- order(data[[rank]], decreasing = TRUE, na.last = NA)
  label_data <- data[utils::head(ordered, n), , drop = FALSE]
  mapping <- ggplot2::aes(
    x = .data[[x]],
    y = .data[[y]],
    label = .data$gene
  )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot + ggrepel::geom_text_repel(
      data = label_data,
      mapping = mapping,
      inherit.aes = FALSE,
      size = 3,
      max.overlaps = Inf
    )
  } else {
    plot + ggplot2::geom_text(
      data = label_data,
      mapping = mapping,
      inherit.aes = FALSE,
      size = 3,
      vjust = -0.5
    )
  }
}

# =============================================================================
# IR signal QC and overview plots
# =============================================================================

.build_ir_qc_panels <- function(pv, ir_stats, label_top = 5L) {
  .require_ir_ggplot2()

  layers <- .validate_ir_stats_for_pv(
    pv,
    ir_stats
  )

  total_rna <- as.numeric(
    Matrix::colSums(layers$ir)
  ) +
    as.numeric(
      Matrix::colSums(layers$unspliced)
    ) +
    as.numeric(
      Matrix::colSums(layers$spliced)
    )

  cell_data <- ir_stats$cell_metrics

  cell_data$total_rna <- total_rna[
    match(
      cell_data$cell,
      colnames(layers$ir)
    )
  ]

  depth_rho <- suppressWarnings(
    stats::cor(
      cell_data$total_rna,
      cell_data$ir_depth_cell,
      method = "spearman",
      use = "complete.obs"
    )
  )

  # ---------------------------------------------------------------------------
  # Panel A: sequencing-depth dependence
  # ---------------------------------------------------------------------------

  depth <- ggplot2::ggplot(
    cell_data,
    ggplot2::aes(
      x = log10(.data$total_rna + 1),
      y = log10(.data$ir_depth_cell + 1)
    )
  ) +
    ggplot2::geom_bin_2d(
      bins = 35
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c(
        "#D5E3EC",
        "#9DBDD0",
        "#5F8EAA",
        "#234F6D"
      ),
      name = "Cell density"
    ) +
    ggplot2::labs(
      title = "IR signal depth dependence",
      subtitle = sprintf(
        "Spearman rho = %.2f",
        depth_rho
      ),
      x = "log10 total RNA counts + 1",
      y = "log10 total IR signal + 1"
    ) +
    .ir_theme()

  # ---------------------------------------------------------------------------
  # Panel B: gene-level prevalence
  # ---------------------------------------------------------------------------

  gene_data <- ir_stats$gene_metrics

  prevalence <- ggplot2::ggplot(
    gene_data,
    ggplot2::aes(
      x = .data$ir_detection_rate_gene,
      y = .data$ir_mean_abundance_norm,
      size = .data$ir_total_count,
      colour = .data$ir_mean_abundance_norm
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.8
    ) +
    ggplot2::scale_colour_gradientn(
      colours = c(
        "#D5E3EC",
        "#9DBDD0",
        "#5F8EAA",
        "#234F6D"
      ),
      name = "Mean normalized\nIR signal"
    ) +
    ggplot2::labs(
      title = "Gene-level IR signal prevalence",
      x = "IR signal detection rate",
      y = "Mean normalized IR signal",
      size = "Total IR signal"
    ) +
    .ir_theme()

  prevalence <- .add_ir_gene_labels(
    prevalence,
    gene_data,
    "ir_detection_rate_gene",
    "ir_mean_abundance_norm",
    "ir_mean_abundance_norm",
    label_top
  )

  # ---------------------------------------------------------------------------
  # Panel C: cell-level sparsity
  # ---------------------------------------------------------------------------

  sparsity <- ggplot2::ggplot(
    cell_data,
    ggplot2::aes(
      x = .data$ir_detected_genes
    )
  ) +
    ggplot2::geom_histogram(
      bins = 30,
      fill = "#456A8A",
      colour = "white"
    ) +
    ggplot2::labs(
      title = "Cell-level IR signal sparsity",
      x = "Detected IR genes",
      y = "Cells"
    ) +
    .ir_theme()

  # ---------------------------------------------------------------------------
  # Panel D: cell-level IR signal fraction
  # ---------------------------------------------------------------------------

  fraction <- ggplot2::ggplot(
    cell_data,
    ggplot2::aes(
      x = .data$ir_fraction_cell
    )
  ) +
    ggplot2::geom_density(
      fill = "#B7C9D6",
      colour = "#456A8A",
      linewidth = 0.6,
      alpha = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = "Cell-level IR signal fraction",
      x = "IR signal fraction",
      y = "Density"
    ) +
    .ir_theme()

  list(
    depth = depth,
    sparsity = sparsity,
    prevalence = prevalence,
    fraction = fraction
  )
}


#' Plot IR-associated Signal QC and Overview
#'
#' Build the four Page 1 panels for sequencing-depth dependence,
#' gene-level prevalence, cell-level sparsity, and cell-level
#' IR-associated signal fraction.
#'
#' @param pv A `plantvelo` object used to retrieve aligned signal layers.
#' @param ir_stats A standalone result returned by
#'   \code{\link{compute_ir_stats}()}.
#' @param label_top Number of genes with high normalized IR-associated signal
#'   labelled in the prevalence panel. Use `NULL` or zero to disable labels.
#'
#' @return A `patchwork` object when patchwork is installed, otherwise a named
#'   list of four ggplot objects.
#'
#' @examples
#' \dontrun{
#' ir_stats <- compute_ir_stats(pv)
#' plot_ir_qc(pv, ir_stats)
#' }
#'
#' @export
plot_ir_qc <- function(pv, ir_stats, label_top = 5L) {
  panels <- .build_ir_qc_panels(
    pv,
    ir_stats,
    label_top = label_top
  )

  .combine_ir_panels(
    panels,
    ncol = 2L
  )
}

# =============================================================================
# Cell-level IR landscape plots
# =============================================================================

.build_ir_cell_panels <- function(pv,
                                  ir_stats,
                                  reduction = "umap",
                                  point_size = 0.25,
                                  alpha = 0.90,
                                  raster_points = FALSE,
                                  raster_dpi = 1200,
                                  input_color = c(
                                    "#EDF3F7",
                                    "#C9DCE8",
                                    "#91B8CF",
                                    "#578DAE",
                                    "#234F6D"
                                  )) {
  .require_ir_ggplot2()
  if (!is.logical(raster_points) || length(raster_points) != 1L ||
      is.na(raster_points)) {
    stop("`raster_points` must be TRUE or FALSE.")
  }
  if (!is.numeric(raster_dpi) || length(raster_dpi) != 1L ||
      !is.finite(raster_dpi) || raster_dpi <= 0) {
    stop("`raster_dpi` must be one positive number.")
  }
  if (raster_points && !requireNamespace("ggrastr", quietly = TRUE)) {
    stop(
      paste0(
        "Rasterization requires package 'ggrastr'. ",
        "Install it with: install.packages('ggrastr')"
      )
    )
  }
  .validate_ir_stats_for_pv(pv, ir_stats)
  embedding <- .get_embedding(pv, reduction)
  if (ncol(embedding) < 2L || is.null(rownames(embedding))) {
    stop(sprintf(
      "Reduction '%s' must have cell names and at least two dimensions.",
      reduction
    ))
  }
  cell_index <- match(rownames(embedding), ir_stats$cells)
  if (anyNA(cell_index)) {
    stop("Reduction cells do not align with `ir_stats$cells`.")
  }

  data <- ir_stats$cell_metrics

  .embedding_plot <- function(metric, legend) {
    plot <- plot_velocity_embedding(
      pv = pv,
      reduction = reduction,
      group_by = data[[metric]][cell_index],
      show_arrows = FALSE,
      point_size = point_size,
      alpha = alpha,
      raster_points = raster_points,
      raster_dpi = raster_dpi,
      input_color = input_color,
      title = ""
    )
    colour_scale <- plot$scales$get_scales("colour")
    if (!is.null(colour_scale)) {
      colour_scale$name <- legend
    }
    plot
  }

  list(
    fraction = .embedding_plot(
      "ir_fraction_cell",
      "IR signal fraction"
    ),
    burden = .embedding_plot(
      "ir_burden_norm",
      "Normalized IR burden"
    ),
    detected = .embedding_plot(
      "ir_detected_genes",
      "Detected IR genes"
    )
  )
}

#' Plot the Cell-Level IR Landscape
#'
#' Color a cell embedding by IR fraction, normalized burden, and detected IR
#' genes, then summarize their pairwise correlations. Velocity arrows are not
#' included because IR is excluded from the inferred vector field.
#'
#' @param pv A `plantvelo` object containing the requested reduction.
#' @param ir_stats A standalone result returned by
#'   \code{\link{compute_ir_stats}()}.
#' @param reduction Name of a reduction in `pv@reductions`.
#' @param point_size Non-negative cell-point size forwarded to
#'   \code{\link{plot_velocity_embedding}()}.
#' @param alpha Cell-point opacity forwarded to
#'   \code{\link{plot_velocity_embedding}()}.
#' @param raster_points Logical; rasterize embedding point layers with
#'   `ggrastr`. Axes, legends, titles, and the correlation panel remain vector
#'   graphics.
#' @param raster_dpi Positive internal resolution for rasterized point layers.
#' @param input_color Continuous colour palette forwarded to
#'   \code{\link{plot_velocity_embedding}()}.
#'
#' @return A `patchwork` object when patchwork is installed, otherwise a named
#'   list of four ggplot objects.
#'
#' @examples
#' \dontrun{
#' ir_stats <- compute_ir_stats(pv)
#' plot_ir_cell_landscape(pv, ir_stats, reduction = "umap")
#' }
#'
#' @export
plot_ir_cell_landscape <- function(pv,
                                   ir_stats,
                                   reduction = "umap",
                                   point_size = 0.25,
                                   alpha = 0.90,
                                   raster_points = FALSE,
                                   raster_dpi = 1200,
                                   input_color = c(
                                     "#EDF3F7",
                                     "#C9DCE8",
                                     "#91B8CF",
                                     "#578DAE",
                                     "#234F6D"
                                   )) {
  panels <- .build_ir_cell_panels(
    pv,
    ir_stats,
    reduction = reduction,
    point_size = point_size,
    alpha = alpha,
    raster_points = raster_points,
    raster_dpi = raster_dpi,
    input_color = input_color
  )
  .combine_ir_panels(panels, ncol = 3)
}

# =============================================================================
# Gene-level IR association plots
# =============================================================================

.validate_ir_gene_n_top <- function(n_top) {
  if (!is.numeric(n_top) || is.complex(n_top) || length(n_top) != 1L ||
      !is.finite(n_top) || n_top < 1L ||
      n_top > .Machine$integer.max || n_top != floor(n_top)) {
    stop("`n_top` must be a positive whole-number scalar.")
  }
  as.integer(n_top)
}

.validate_ir_latent_bins <- function(latent_bins) {
  if (!is.numeric(latent_bins) || is.complex(latent_bins) ||
      length(latent_bins) != 1L || !is.finite(latent_bins) ||
      latent_bins < 2L || latent_bins > .Machine$integer.max ||
      latent_bins != floor(latent_bins)) {
    stop("`latent_bins` must be a whole-number scalar of at least 2.")
  }
  as.integer(latent_bins)
}

.select_ir_genes <- function(ir_stats, genes, n_top) {
  available <- ir_stats$genes
  if (is.null(genes)) {
    order_index <- order(
      ir_stats$gene_metrics$ir_fraction_gene,
      decreasing = TRUE,
      na.last = NA
    )
    return(utils::head(
      ir_stats$gene_metrics$gene[order_index],
      n_top
    ))
  }

  if (!is.character(genes) || !length(genes)) {
    stop("`genes` must be a non-empty character vector or NULL.")
  }
  genes <- unique(genes)
  missing <- setdiff(genes, available)
  if (length(missing) == length(genes)) {
    stop("None of the requested genes are available in `ir_stats`.")
  }
  if (length(missing)) {
    warning(sprintf(
      "Ignoring unavailable genes: %s.",
      paste(missing, collapse = ", ")
    ))
  }
  intersect(genes, available)
}

.get_ir_latent_time <- function(pv, ir_stats) {
  latent <- pv@meta.data[["latent_time"]]
  if (is.null(latent) || !is.numeric(latent) ||
      is.null(rownames(pv@meta.data))) {
    return(NULL)
  }
  latent <- latent[match(ir_stats$cells, rownames(pv@meta.data))]
  valid <- is.finite(latent)
  if (sum(valid) < 2L || diff(range(latent[valid])) <= 0) {
    return(NULL)
  }

  latent
}

.latent_ir_summary <- function(pv, ir_stats, genes, latent_bins) {
  latent <- .get_ir_latent_time(pv, ir_stats)
  if (is.null(latent)) {
    return(NULL)
  }
  valid <- is.finite(latent)

  breaks <- seq(
    min(latent[valid]),
    max(latent[valid]),
    length.out = latent_bins + 1L
  )
  bins <- cut(
    latent,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
  matrix_data <- as.matrix(
    ir_stats$normalized_ir[genes, , drop = FALSE]
  )
  observed_bins <- sort(unique(bins[!is.na(bins)]))
  summaries <- lapply(observed_bins, function(bin) {
    rowMeans(matrix_data[, bins == bin, drop = FALSE])
  })
  mean_matrix <- do.call(cbind, summaries)
  colnames(mean_matrix) <- paste0("bin_", observed_bins)
  z_matrix <- t(scale(t(mean_matrix)))
  z_matrix[!is.finite(z_matrix)] <- 0
  peak <- max.col(z_matrix, ties.method = "first")
  z_matrix[order(peak), , drop = FALSE]
}

.build_ir_gene_plot <- function(pv,
                                ir_stats,
                                plot_type = c(
                                  "ranking",
                                  "heatmap",
                                  "trend",
                                  "driver"
                                ),
                                genes = NULL,
                                n_top = 20L,
                                latent_bins = 20L) {
  .require_ir_ggplot2()
  .validate_ir_stats_for_pv(pv, ir_stats)
  plot_type <- match.arg(plot_type)
  if (is.null(genes) || plot_type == "driver") {
    n_top <- .validate_ir_gene_n_top(n_top)
  }
  if (!length(ir_stats$genes)) {
    stop("No usable IR genes are available for the gene association plot.")
  }

  selected <- .select_ir_genes(ir_stats, genes, n_top)

  if (plot_type == "ranking") {
    gene_data <- ir_stats$gene_metrics[
      match(selected, ir_stats$gene_metrics$gene),
      ,
      drop = FALSE
    ]
    gene_data <- gene_data[
      order(gene_data$ir_fraction_gene, decreasing = TRUE),
      ,
      drop = FALSE
    ]
    gene_data$gene <- factor(
      gene_data$gene,
      levels = rev(gene_data$gene)
    )

    return(
      ggplot2::ggplot(
        gene_data,
        ggplot2::aes(x = .data$ir_fraction_gene, y = .data$gene)
      ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0,
            xend = .data$ir_fraction_gene,
            yend = .data$gene
          ),
          colour = "grey70"
        ) +
        ggplot2::geom_point(
          ggplot2::aes(
            colour = .data$ir_detection_rate_gene,
            size = .data$ir_total_count
          )
        ) +
        ggplot2::scale_colour_gradientn(
          colours = c(
            "#D5E3EC",
            "#9DBDD0",
            "#5F8EAA",
            "#234F6D"
          ),
          name = "IR signal\ndetection rate"
        ) +
        ggplot2::scale_size_continuous(
          range = c(1.5, 5),
          name = "Total IR signal"
        ) +
        ggplot2::labs(
          title = "Top IR-associated genes",
          x = "Gene-level IR signal fraction",
          y = NULL,
          size = "Total IR signal"
        ) +
        .ir_theme()
    )
  }

  if (plot_type == "heatmap") {
    latent_bins <- .validate_ir_latent_bins(latent_bins)
    latent_matrix <- .latent_ir_summary(
      pv,
      ir_stats,
      selected,
      latent_bins
    )
    if (is.null(latent_matrix)) {
      stop("`latent_time` is unavailable or unusable for this plot.")
    }
  }

  if (plot_type == "trend") {
    latent_bins <- .validate_ir_latent_bins(latent_bins)
    latent <- .get_ir_latent_time(pv, ir_stats)
    if (is.null(latent)) {
      stop("`latent_time` is unavailable or unusable for this plot.")
    }
  }

  if (plot_type == "heatmap") {
    heatmap_data <- as.data.frame(
      as.table(latent_matrix),
      stringsAsFactors = FALSE
    )
    names(heatmap_data) <- c("gene", "latent_bin", "z_score")
    heatmap_data$latent_bin <- factor(
      heatmap_data$latent_bin,
      levels = colnames(latent_matrix),
      ordered = TRUE
    )
    return(
      ggplot2::ggplot(
        heatmap_data,
        ggplot2::aes(
          x = .data$latent_bin,
          y = .data$gene,
          fill = .data$z_score
        )
      ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
          low = "#355FA3",
          mid = "white",
          high = "#B33A3A",
          midpoint = 0,
          name = "Mean IR signal\n(within-gene z-score)"
        ) +
        ggplot2::labs(
          title = "IR signal across latent time",
          x = "Latent-time bin",
          y = NULL
        ) +
        .ir_theme() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    )
  }

  if (plot_type == "trend") {
    trend_colours <- c(
      "#3B4CC0",  # blue
      "#648FFF",  # light blue
      "#2F7F9D",  # steel blue
      "#00A6A6",  # teal
      "#009E73",  # green
      "#59A14F",  # leaf green
      "#A6C84C",  # yellow green
      "#E6AB02",  # ochre
      "#F28E2B",  # orange
      "#E76F51",  # coral
      "#D62728",  # red
      "#B23A48",  # dark red
      "#D45087",  # magenta
      "#CC79A7",  # pink purple
      "#9467BD",  # purple
      "#6F4C9B",  # dark purple
      "#8C564B",  # brown
      "#C49A6C",  # tan
      "#4D4D4D",  # dark grey
      "#76B7B2"   # muted cyan
    )

    gene_colours <- stats::setNames(
      rep(trend_colours, length.out = length(selected)),
      selected
    )

    normalized <- as.matrix(
      ir_stats$normalized_ir[selected, , drop = FALSE]
    )
    trend_data <- do.call(rbind, lapply(seq_along(selected), function(index) {
      data.frame(
        gene = selected[index],
        latent_time = latent,
        value = normalized[index, ],
        stringsAsFactors = FALSE
      )
    }))
    return(
      ggplot2::ggplot(
        trend_data,
        ggplot2::aes(
          x = .data$latent_time,
          y = .data$value,
          colour = .data$gene
        )
      ) +
        ggplot2::geom_point(alpha = 0.15, size = 0.5) +
        ggplot2::stat_summary_bin(
          fun = mean,
          geom = "line",
          bins = latent_bins,
          linewidth = 0.8
        ) +
        ggplot2::facet_wrap(~gene, scales = "free_y") +
        ggplot2::labs(
          title = "Post-hoc relationship",
          x = "Latent time",
          y = "Normalized IR signal",
          colour = "Gene"
        ) +
        ggplot2::scale_x_continuous(
          breaks = c(0, 0.5, 1),
          labels = c("0", "0.5", "1.0")
        ) +
        ggplot2::scale_colour_manual(
          values = gene_colours
        ) +
        .ir_theme() +
        ggplot2::theme(legend.position = "none",
                       strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_text(size = 11))
    )
  }

  drivers <- pv@misc[["velocity_gene_stats"]]
  required_driver_columns <- c("gene", "coherence", "fit_likelihood")
  if (is.null(drivers) ||
      any(!required_driver_columns %in% colnames(drivers))) {
    stop(
      "`velocity_gene_stats` with gene, coherence, and fit_likelihood is required."
    )
  }
  driver_data <- merge(
    ir_stats$gene_metrics,
    drivers[, required_driver_columns, drop = FALSE],
    by = "gene"
  )
  if (!nrow(driver_data)) {
    stop("No overlapping genes are available for the driver association.")
  }
  driver_plot <- ggplot2::ggplot(
    driver_data,
    ggplot2::aes(
      x = .data$coherence,
      y = .data$ir_fraction_gene,
      colour = .data$fit_likelihood
    )
  ) +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::scale_colour_gradientn(
      colours = c(
        "#D5E3EC",
        "#9DBDD0",
        "#5F8EAA",
        "#234F6D"
      ),
      name = "Fit likelihood"
    ) +
    ggplot2::labs(
      title = "Post-hoc IR-driver association",
      x = "Velocity coherence",
      y = "Gene-level IR signal fraction"
    ) +
    .ir_theme()

  driver_order <- order(
    driver_data$coherence,
    decreasing = TRUE,
    na.last = NA
  )
  top_driver_genes <- utils::head(
    driver_data$gene[driver_order],
    n_top
  )
  overlap <- intersect(selected, top_driver_genes)
  if (length(overlap)) {
    label_data <- driver_data[
      driver_data$gene %in% overlap,
      ,
      drop = FALSE
    ]
    driver_plot <- .add_ir_gene_labels(
      driver_plot,
      label_data,
      "coherence",
      "ir_fraction_gene",
      "coherence",
      length(overlap)
    )
  }

  driver_plot
}

#' Plot Gene-Level IR and Post-Hoc Dynamics Associations
#'
#' Draw one selected gene-level IR plot using existing two-state latent time or
#' velocity-driver statistics when requested. These plots are descriptive
#' post-hoc analyses and do not alter dynamics results.
#'
#' @param pv A `plantvelo` object.
#' @param ir_stats A standalone result returned by
#'   \code{\link{compute_ir_stats}()}.
#' @param plot_type One of `"ranking"`, `"heatmap"`, `"trend"`, or `"driver"`.
#' @param genes Optional character vector of genes to highlight or display.
#' @param n_top Positive number of top IR genes selected when `genes = NULL`.
#' @param latent_bins Number of deterministic latent-time bins used by
#'   `"heatmap"` and `"trend"`.
#'
#' @return A single `ggplot` object.
#'
#' @examples
#' \dontrun{
#' ir_stats <- compute_ir_stats(pv)
#' plot_ir_gene_association(pv, ir_stats)
#' }
#'
#' @export
plot_ir_gene_association <- function(pv,
                                     ir_stats,
                                     plot_type = c(
                                       "ranking",
                                       "heatmap",
                                       "trend",
                                       "driver"
                                     ),
                                     genes = NULL,
                                     n_top = 20L,
                                     latent_bins = 20L) {
  .build_ir_gene_plot(
    pv,
    ir_stats,
    plot_type = plot_type,
    genes = genes,
    n_top = n_top,
    latent_bins = latent_bins
  )
}
