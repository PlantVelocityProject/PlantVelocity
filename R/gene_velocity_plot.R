# =============================================================================
# PlantVelocity: gene_velocity_plot.R
# Gene-wise RNA velocity colouring on existing cell embeddings.
# Display signed velocity estimates without recomputing transition graphs or
# projecting vectors. Missing estimates remain distinct from genuine zeros.
# =============================================================================

#' Plot Gene-wise Velocity on a Cell Embedding
#'
#' Colour cells in an existing embedding by the signed RNA velocity of each
#' requested gene. This displays scalar velocity estimates, not projected
#' arrows, and requires neither a velocity graph nor a velocity embedding.
#'
#' @param pv A \code{plantvelo} object containing a velocity layer and cell
#'   embedding.
#' @param genes Non-empty character vector of gene names. Repeated names are
#'   displayed once, in first-occurrence order. Unknown genes cause an error.
#' @param reduction Embedding name, matched case-insensitively; a leading
#'   \code{X_} is ignored. The first two dimensions are displayed.
#' @param vkey Name of the genes-by-cells velocity layer, such as
#'   \code{"velocity"} or \code{"velocity_u"}.
#' @param ncol Positive integer specifying the number of panel columns.
#'   Use \code{1L} to stack genes vertically.
#' @param scale_mode Either \code{"per_gene"} (each gene has its own symmetric
#'   colour limits) or \code{"shared"} (one symmetric range across all requested
#'   genes and displayed cells). Independent ranges emphasize within-gene
#'   patterns; identical colours then need not represent identical velocities.
#' @param limits Optional finite numeric vector \code{c(lower, upper)}, with
#'   \code{lower < 0 < upper}. Overrides automatic limits for all panels.
#'   Out-of-range values are colour-clipped without changing the plotted data.
#' @param point_size Non-negative point size.
#' @param alpha Point opacity between zero and one.
#' @param raster_points Whether to rasterize points using \pkg{ggrastr}.
#' @param raster_dpi Positive raster resolution in dots per inch.
#' @param input_color Three colours for negative velocity, zero, and positive
#'   velocity, respectively. Missing values are always shown in \code{grey80}.
#'
#' @details
#' Cells are aligned by name rather than matrix position. Embeddings may
#' contain a subset of velocity-layer cells, but all displayed cells must be
#' present in that layer. Only requested velocity rows are converted to a
#' dense matrix; the input object is not modified.
#'
#' Non-finite estimates are displayed as missing, never as zero. When
#' \code{pv@misc[[paste0(vkey, "_params")]]$mode} is \code{"dynamical"}
#' and \code{pv@kinetics$params$fit_model} is available, genes absent from the
#' fit table or not marked \code{"fitted"} are also displayed as missing.
#' A layer ending in \code{_u} falls back to its base layer's mode metadata.
#' In \code{"stochastic"} mode, fit status is ignored. Without mode metadata,
#' only genes explicitly marked \code{"unfitted"} are masked. Quality flags
#' in \code{pv@velocity$genes} are not used to discard estimates.
#'
#' A gene with no finite displayed estimates produces a warning. Automatic
#' limits fall back to \code{c(-1, 1)} for all-zero or entirely missing data.
#' Zero remains the colour midpoint, including with asymmetric fixed limits.
#'
#' @return A \pkg{ggplot2} plot for one gene, or a \pkg{patchwork} composition
#'   for multiple genes. Each panel's data contains \code{cell}, \code{x},
#'   \code{y}, \code{gene}, and \code{velocity}; missing or unfitted estimates
#'   have \code{NA} velocity.
#' @examples
#' \dontrun{
#' plot_gene_velocity(pv, genes = "AT1G01010")
#' p <- plot_gene_velocity(pv, genes = c("AT1G01010", "AT1G01020"),
#'                         ncol = 1L, scale_mode = "shared")
#' ggplot2::ggsave("gene_velocity.png", p, width = 4, height = 7, dpi = 300)
#' }
#' @export
plot_gene_velocity <- function(pv,
                               genes,
                               reduction = "umap",
                               vkey = "velocity",
                               ncol = 3L,
                               scale_mode = c("per_gene", "shared"),
                               limits = NULL,
                               point_size = 0.25,
                               alpha = 0.9,
                               raster_points = FALSE,
                               raster_dpi = 600,
                               input_color = c("#2166AC", "#F7F7F7", "#B2182B")) {
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.", call. = FALSE)
  valid_strings <- function(x) {
    is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x))
  }
  if (!valid_strings(genes))
    stop("`genes` must be a non-empty character vector without missing names.",
         call. = FALSE)
  genes <- unique(genes)
  for (key in c("reduction", "vkey")) {
    value <- get(key)
    if (!valid_strings(value) || length(value) != 1L)
      stop(sprintf("`%s` must be one non-empty string.", key), call. = FALSE)
  }
  scale_mode <- match.arg(scale_mode)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(ncol) || ncol < 1 || ncol != floor(ncol))
    stop("`ncol` must be a positive integer.", call. = FALSE)
  if (!scalar(point_size) || point_size < 0)
    stop("`point_size` must be finite and non-negative.", call. = FALSE)
  if (!scalar(alpha) || alpha < 0 || alpha > 1)
    stop("`alpha` must be between zero and one.", call. = FALSE)
  if (!is.logical(raster_points) || length(raster_points) != 1L ||
      is.na(raster_points))
    stop("`raster_points` must be TRUE or FALSE.", call. = FALSE)
  if (!scalar(raster_dpi) || raster_dpi <= 0)
    stop("`raster_dpi` must be finite and positive.", call. = FALSE)
  if (!valid_strings(input_color) || length(input_color) != 3L)
    stop("`input_color` must contain three colours (negative, zero, positive).",
         call. = FALSE)
  if (!is.null(limits) &&
      (!is.numeric(limits) || length(limits) != 2L ||
       any(!is.finite(limits)) || limits[1] >= 0 || limits[2] <= 0))
    stop("`limits` must be finite c(lower, upper) with lower < 0 < upper.",
         call. = FALSE)

  prepared <- .gene_velocity_plot_data(pv, genes, reduction, vkey)
  axis_labels <- sub("umap", "UMAP", colnames(prepared$embedding), ignore.case = TRUE)
  values <- prepared$values
  symmetric_limits <- function(x) {
    finite <- x[is.finite(x)]
    extent <- if (length(finite)) max(abs(finite)) else 0
    if (extent == 0) extent <- 1
    c(-extent, extent)
  }
  shared_limits <- if (!is.null(limits)) limits else if (scale_mode == "shared") {
    symmetric_limits(values)
  } else NULL
  missing_genes <- genes[rowSums(is.finite(values)) == 0L]
  if (length(missing_genes))
    warning("No finite velocity estimates for: ", paste(missing_genes, collapse = ", "),
            ". These cells are shown in grey.", call. = FALSE)

  panels <- lapply(seq_along(genes), function(i) {
    dat <- data.frame(cell = rownames(prepared$embedding),
                      x = prepared$embedding[, 1], y = prepared$embedding[, 2],
                      gene = genes[i], velocity = as.numeric(values[i, ]),
                      row.names = NULL)
    panel_limits <- if (is.null(shared_limits)) symmetric_limits(dat$velocity) else shared_limits
    points <- ggplot2::geom_point(size = point_size, alpha = alpha)
    if (raster_points) points <- ggrastr::rasterise(points, dpi = raster_dpi)
    ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data$y,
                                   colour = .data$velocity)) +
      points +
      ggplot2::scale_colour_gradient2(
        low = input_color[1], mid = input_color[2], high = input_color[3],
        midpoint = 0, limits = panel_limits, oob = scales::squish,
        na.value = "grey80", name = "Velocity"
      ) +
      ggplot2::labs(title = genes[i], x = axis_labels[1],
                    y = axis_labels[2]) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, color = "black", hjust = 0.5),
        axis.title = ggplot2::element_text(size = 12, color = "black"),
        axis.text = ggplot2::element_text(size = 10, color = "black"),
        legend.title = ggplot2::element_text(size = 10, color = "black"),
        legend.text = ggplot2::element_text(size = 9, color = "black"),
        legend.key.size = grid::unit(0.4, "cm")
      )
  })
  if (length(panels) == 1L) return(panels[[1]])
  patchwork::wrap_plots(panels, ncol = ncol,
                        guides = if (is.null(shared_limits)) "keep" else "collect")
}

.gene_velocity_plot_data <- function(pv, genes, reduction, vkey) {
  velocity <- pv@layers[[vkey]]
  if (is.null(velocity))
    stop(sprintf("Velocity layer '%s' not found. Run compute_velocity() first.", vkey),
         call. = FALSE)
  numeric_matrix <- function(x) {
    (is.matrix(x) && is.numeric(x)) || methods::is(x, "dMatrix") ||
      methods::is(x, "iMatrix")
  }
  if (!numeric_matrix(velocity))
    stop("The velocity layer must be a numeric genes-by-cells matrix.", call. = FALSE)
  check_names <- function(x, label) {
    if (is.null(x) || !length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
      stop(label, " must have non-empty, unique names without missing values.",
           call. = FALSE)
  }
  check_names(rownames(velocity), "Velocity genes")
  check_names(colnames(velocity), "Velocity cells")
  unknown <- setdiff(genes, rownames(velocity))
  if (length(unknown))
    stop("Unknown genes in velocity layer: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  embedding <- .get_embedding(pv, reduction)
  if (!numeric_matrix(embedding) || ncol(embedding) < 2L || nrow(embedding) == 0L)
    stop("The embedding must be numeric with cells and at least two dimensions.",
         call. = FALSE)
  check_names(rownames(embedding), "Embedding cells")
  embedding <- as.matrix(embedding[, 1:2, drop = FALSE])
  if (any(!is.finite(embedding)))
    stop("Displayed embedding coordinates must be finite.", call. = FALSE)
  unknown <- setdiff(rownames(embedding), colnames(velocity))
  if (length(unknown))
    stop("Embedding cells absent from velocity layer: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  values <- as.matrix(velocity[genes, rownames(embedding), drop = FALSE])
  values[!is.finite(values)] <- NA_real_

  metadata <- pv@misc[[paste0(vkey, "_params")]]
  if (is.null(metadata) && grepl("_u$", vkey))
    metadata <- pv@misc[[paste0(sub("_u$", "", vkey), "_params")]]
  mode <- if (is.list(metadata)) metadata$mode else NULL
  params <- pv@kinetics$params
  if (!identical(mode, "stochastic") && is.data.frame(params) &&
      "fit_model" %in% names(params)) {
    check_names(rownames(params), "Kinetic fit genes")
    status <- as.character(params$fit_model[match(genes, rownames(params))])
    unfitted <- if (identical(mode, "dynamical")) {
      is.na(status) | status != "fitted"
    } else {
      !is.na(status) & status == "unfitted"
    }
    values[unfitted, ] <- NA_real_
  }
  list(values = values, embedding = embedding)
}
