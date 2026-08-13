# =============================================================================
# PlantVelocity: velocity_grid.R
#
# Plot locally averaged RNA velocity arrows on a regular embedding grid.
# =============================================================================

.scale_velocity_grid_arrows <- function(grid,
                                        coords,
                                        arrow_scale,
                                        arrow_min_length,
                                        arrow_max_length) {
  required <- c("x", "y", "dx", "dy", "valid")
  if (!is.data.frame(grid) || !all(required %in% names(grid))) {
    stop("`grid` must contain x, y, dx, dy, and valid columns.",
         call. = FALSE)
  }
  if (!is.numeric(grid$x) || !is.numeric(grid$y) ||
      !is.numeric(grid$dx) || !is.numeric(grid$dy) ||
      !is.logical(grid$valid) || length(grid$valid) != nrow(grid)) {
    stop("Grid coordinates and vectors must be numeric, and `valid` must be logical.",
         call. = FALSE)
  }
  if ((!is.matrix(coords) && !is.data.frame(coords)) ||
      ncol(coords) < 2L) {
    stop("`coords` must be a numeric matrix or data.frame with at least two columns.",
         call. = FALSE)
  }
  coords <- as.matrix(coords)
  if (!is.numeric(coords)) {
    stop("`coords` must be a numeric matrix or data.frame with at least two columns.",
         call. = FALSE)
  }

  magnitude <- sqrt(grid$dx^2 + grid$dy^2)
  keep <- !is.na(grid$valid) & grid$valid &
    is.finite(grid$x) & is.finite(grid$y) &
    is.finite(grid$dx) & is.finite(grid$dy) &
    is.finite(magnitude) & magnitude > 1e-12
  grid <- grid[keep, , drop = FALSE]
  magnitude <- magnitude[keep]

  if (nrow(grid) == 0L) {
    return(data.frame(
      x = numeric(),
      y = numeric(),
      xend = numeric(),
      yend = numeric(),
      magnitude = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  embedding_span <- max(
    diff(range(coords[, 1L], finite = TRUE)),
    diff(range(coords[, 2L], finite = TRUE))
  )
  if (!is.finite(embedding_span) || embedding_span <= 0) {
    embedding_span <- 1
  }

  magnitude_cap <- stats::quantile(
    magnitude,
    probs = 0.95,
    names = FALSE,
    na.rm = TRUE
  )
  if (!is.finite(magnitude_cap) || magnitude_cap <= 1e-12) {
    magnitude_cap <- max(magnitude)
  }

  normalized <- pmin(magnitude, magnitude_cap) / magnitude_cap
  visible_length <- embedding_span * (
    arrow_min_length +
      (arrow_max_length - arrow_min_length) * sqrt(normalized)
  ) * arrow_scale

  data.frame(
    x = grid$x,
    y = grid$y,
    xend = grid$x + grid$dx / magnitude * visible_length,
    yend = grid$y + grid$dy / magnitude * visible_length,
    magnitude = magnitude,
    stringsAsFactors = FALSE
  )
}


.align_velocity_grid_cells <- function(embedding, velocity_embedding) {
  rowname_state <- function(x) {
    identifiers <- rownames(x)
    if (is.null(identifiers)) {
      return(list(present = FALSE, valid = FALSE, values = NULL))
    }
    list(
      present = TRUE,
      valid = length(identifiers) == nrow(x) &&
        all(!is.na(identifiers)) && all(nzchar(identifiers)) &&
        !anyDuplicated(identifiers),
      values = identifiers
    )
  }

  embedding_rows <- rowname_state(embedding)
  velocity_rows <- rowname_state(velocity_embedding)
  if (embedding_rows$present != velocity_rows$present) {
    stop(
      "Embedding and velocity embedding must either both have row names or both have none.",
      call. = FALSE
    )
  }
  if (embedding_rows$present &&
      (!embedding_rows$valid || !velocity_rows$valid)) {
    stop(
      "Embedding and velocity embedding row names must be complete and unique.",
      call. = FALSE
    )
  }
  if (embedding_rows$present) {
    if (!setequal(embedding_rows$values, velocity_rows$values)) {
      stop(
        "Embedding and velocity embedding row names identify different cells.",
        call. = FALSE
      )
    }
    velocity_embedding <- velocity_embedding[
      embedding_rows$values,
      ,
      drop = FALSE
    ]
  } else if (nrow(embedding) != nrow(velocity_embedding)) {
    stop(
      "Embedding and velocity embedding must have the same number of rows.",
      call. = FALSE
    )
  }

  velocity_embedding
}


#' @title Plot Velocity Embedding on a Regular Grid
#'
#' @description
#' Plot cells in a two-dimensional embedding overlaid with locally averaged
#' RNA velocity arrows on a regular grid. Cell points and grid arrows can be
#' rasterized independently while axes, legends, and text remain vector
#' graphics.
#'
#' @param pv A \code{plantvelo} object.
#' @param reduction Character scalar. Embedding name, such as \code{"umap"},
#'   \code{"pca"}, or \code{"tsne"}. Default \code{"umap"}.
#' @param vkey Character scalar. Velocity layer key. Default \code{"velocity"}.
#' @param group_by Character scalar naming a metadata column, a numeric vector,
#'   or \code{NULL}, using the same colour semantics as
#'   \code{plot_velocity_embedding()}. Default \code{NULL}.
#' @param n_grid Integer. Number of regular grid positions along each axis.
#'   Default \code{30L}.
#' @param n_neighbors Integer or \code{NULL}. Number of nearest cells used to
#'   estimate velocity at each grid position. If \code{NULL}, an automatic
#'   value is used.
#' @param smooth Positive numeric. Gaussian bandwidth as a multiple of the mean
#'   grid spacing. Default \code{0.5}.
#' @param min_mass Numeric. Lower-bound control for the interpolated
#'   velocity-strength mask. Default \code{1}.
#' @param cutoff_perc Numeric in \eqn{[0,100]}. Percentile threshold applied to
#'   local velocity length. Default \code{5}.
#' @param x_pad,y_pad Non-negative numeric. Fractional padding added to the
#'   embedding ranges before constructing the grid. Default \code{0.01}.
#' @param arrow_scale Non-negative numeric. Global arrow-length multiplier.
#'   Default \code{0.8}.
#' @param arrow_min_length,arrow_max_length Non-negative numeric. Minimum and
#'   maximum visible arrow lengths as fractions of the total embedding span.
#'   The maximum must be at least the minimum. Defaults are \code{0.008} and
#'   \code{0.035}.
#' @param arrow_size Non-negative numeric. Arrow-head length in cm. Default
#'   \code{0.06}.
#' @param arrow_linewidth Non-negative numeric. Arrow shaft linewidth. Default
#'   \code{0.15}.
#' @param arrow_alpha Numeric in \eqn{[0,1]}. Arrow opacity. Default \code{0.90}.
#' @param arrow_angle Numeric strictly between 0 and 90. Arrow-head half-angle
#'   in degrees. Default \code{20}.
#' @param arrow_type Character. Arrow-head type, either \code{"closed"} or
#'   \code{"open"}. Default \code{"closed"}.
#' @param point_size Non-negative numeric. Cell-point size. Default \code{0.25}.
#' @param alpha Numeric in \eqn{[0,1]}. Cell-point opacity. Default \code{0.90}.
#' @param point_shape Finite numeric. Cell-point shape. Default \code{16}.
#' @param raster_points Logical. Rasterize the cell-point layer with
#'   \pkg{ggrastr}. Default \code{TRUE}.
#' @param raster_arrows Logical. Rasterize the grid-arrow layer with
#'   \pkg{ggrastr}. Default \code{FALSE}.
#' @param raster_dpi Positive numeric. Resolution in dots per inch for
#'   rasterized layers. Default \code{1200}.
#' @param input_color Character vector or \code{NULL}. Custom discrete or
#'   continuous colour palette for cells. Default \code{NULL}.
#' @param title Character scalar or \code{NULL}. Plot title. Default
#'   \code{NULL}.
#' @param legend Logical. Whether to show the legend. Default \code{TRUE}.
#' @param fixed_aspect Logical. Whether to use an equal x/y coordinate ratio.
#'   Default \code{FALSE}.
#'
#' @details
#' Grid velocities are estimated with the same Gaussian nearest-neighbour
#' interpolation and masking used by \code{plot_velocity_stream()}, but this
#' function draws one arrow per valid grid position and does not calculate or
#' cache streamlines. Arrow magnitudes are clipped at their 95th percentile,
#' square-root normalized, and mapped to the requested visible length range.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plot_velocity_embedding_grid <- function(
    pv,
    reduction = "umap",
    vkey = "velocity",
    group_by = NULL,
    n_grid = 30L,
    n_neighbors = NULL,
    smooth = 0.5,
    min_mass = 1,
    cutoff_perc = 5,
    x_pad = 0.01,
    y_pad = 0.01,
    arrow_scale = 0.8,
    arrow_min_length = 0.008,
    arrow_max_length = 0.035,
    arrow_size = 0.06,
    arrow_linewidth = 0.15,
    arrow_alpha = 0.90,
    arrow_angle = 20,
    arrow_type = c("closed", "open"),
    point_size = 0.25,
    alpha = 0.90,
    point_shape = 16,
    raster_points = TRUE,
    raster_arrows = FALSE,
    raster_dpi = 1200,
    input_color = NULL,
    title = NULL,
    legend = TRUE,
    fixed_aspect = FALSE) {

  arrow_type <- match.arg(arrow_type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required. Install it with: install.packages('ggplot2')",
      call. = FALSE
    )
  }
  if (!inherits(pv, "plantvelo")) {
    stop("`pv` must be a plantvelo object.", call. = FALSE)
  }

  controls <- .validate_stream_controls(
    n_grid = n_grid,
    n_neighbors = n_neighbors,
    smooth = smooth,
    min_mass = min_mass,
    cutoff_perc = cutoff_perc,
    x_pad = x_pad,
    y_pad = y_pad,
    density = 1,
    max_length = 1,
    integration_direction = "both"
  )

  scalar_finite <- function(x) {
    is.numeric(x) && !is.logical(x) && length(x) == 1L && is.finite(x)
  }
  scalar_logical <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }

  if (!is.character(reduction) || length(reduction) != 1L ||
      is.na(reduction) || !nzchar(reduction)) {
    stop("`reduction` must be one non-empty character value.", call. = FALSE)
  }
  if (!is.character(vkey) || length(vkey) != 1L ||
      is.na(vkey) || !nzchar(vkey)) {
    stop("`vkey` must be one non-empty character value.", call. = FALSE)
  }
  if (!scalar_finite(arrow_scale) || arrow_scale < 0) {
    stop("`arrow_scale` must be a non-negative finite scalar.", call. = FALSE)
  }
  if (!scalar_finite(arrow_min_length) || arrow_min_length < 0 ||
      !scalar_finite(arrow_max_length) ||
      arrow_max_length < arrow_min_length) {
    stop(
      "`arrow_min_length` and `arrow_max_length` must be non-negative finite scalars, and `arrow_max_length` must be greater than or equal to `arrow_min_length`.",
      call. = FALSE
    )
  }
  if (!scalar_finite(arrow_size) || arrow_size < 0) {
    stop("`arrow_size` must be a non-negative finite scalar.", call. = FALSE)
  }
  if (!scalar_finite(arrow_linewidth) || arrow_linewidth < 0) {
    stop("`arrow_linewidth` must be a non-negative finite scalar.", call. = FALSE)
  }
  if (!scalar_finite(arrow_alpha) || arrow_alpha < 0 || arrow_alpha > 1) {
    stop("`arrow_alpha` must be a finite scalar between 0 and 1.", call. = FALSE)
  }
  if (!scalar_finite(arrow_angle) || arrow_angle <= 0 || arrow_angle >= 90) {
    stop("`arrow_angle` must be a finite scalar strictly between 0 and 90.",
         call. = FALSE)
  }
  if (!scalar_finite(point_size) || point_size < 0) {
    stop("`point_size` must be a non-negative finite scalar.", call. = FALSE)
  }
  if (!scalar_finite(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a finite scalar between 0 and 1.", call. = FALSE)
  }
  if (!scalar_finite(point_shape)) {
    stop("`point_shape` must be a finite numeric scalar.", call. = FALSE)
  }
  if (!scalar_logical(raster_points)) {
    stop("`raster_points` must be a single non-missing logical value.",
         call. = FALSE)
  }
  if (!scalar_logical(raster_arrows)) {
    stop("`raster_arrows` must be a single non-missing logical value.",
         call. = FALSE)
  }
  if (!scalar_finite(raster_dpi) || raster_dpi <= 0) {
    stop("`raster_dpi` must be a positive finite scalar.", call. = FALSE)
  }
  if (!scalar_logical(legend)) {
    stop("`legend` must be a single non-missing logical value.", call. = FALSE)
  }
  if (!scalar_logical(fixed_aspect)) {
    stop("`fixed_aspect` must be a single non-missing logical value.",
         call. = FALSE)
  }
  if ((raster_points || raster_arrows) &&
      !requireNamespace("ggrastr", quietly = TRUE)) {
    stop(
      "Rasterization requires package 'ggrastr'. Install it with: install.packages('ggrastr')",
      call. = FALSE
    )
  }

  embedding <- .get_embedding(pv, reduction)
  if (!is.matrix(embedding) && !is.data.frame(embedding)) {
    stop(
      "Embedding returned by `.get_embedding()` must be a matrix or data.frame.",
      call. = FALSE
    )
  }
  embedding <- as.matrix(embedding)
  if (!is.numeric(embedding)) {
    stop("Embedding must be numeric.", call. = FALSE)
  }
  if (nrow(embedding) == 0L) {
    stop("The selected embedding contains no cells.", call. = FALSE)
  }
  if (ncol(embedding) < 2L) {
    stop("Embedding must have at least two columns.", call. = FALSE)
  }

  velocity <- .find_velocity_embedding(pv, reduction, vkey)

  if (is.null(velocity)) {
    pv <- compute_velocity_embedding(
      pv,
      reduction = reduction,
      vkey = vkey,
      verbose = FALSE
    )
    velocity <- .find_velocity_embedding(pv, reduction, vkey)
  }
  if (is.null(velocity)) {
    stop("Velocity embedding could not be retrieved after computation.",
         call. = FALSE)
  }
  if (!is.matrix(velocity) && !is.data.frame(velocity)) {
    stop("Velocity embedding must be a matrix or data.frame.", call. = FALSE)
  }
  velocity <- as.matrix(velocity)
  if (!is.numeric(velocity)) {
    stop("Velocity embedding must be numeric.", call. = FALSE)
  }
  velocity <- .align_velocity_grid_cells(embedding, velocity)
  if (nrow(velocity) != nrow(embedding) || ncol(velocity) < 2L) {
    stop(
      "Velocity embedding has incompatible dimensions with the selected embedding.",
      call. = FALSE
    )
  }

  grid_data <- .estimate_stream_grid(
    coords = embedding[, 1:2, drop = FALSE],
    velocity = velocity[, 1:2, drop = FALSE],
    n_grid = controls$n_grid,
    n_neighbors = controls$n_neighbors,
    smooth = controls$smooth,
    min_mass = controls$min_mass,
    cutoff_perc = controls$cutoff_perc,
    x_pad = controls$x_pad,
    y_pad = controls$y_pad
  )
  arrow_data <- .scale_velocity_grid_arrows(
    grid = grid_data,
    coords = embedding,
    arrow_scale = arrow_scale,
    arrow_min_length = arrow_min_length,
    arrow_max_length = arrow_max_length
  )

  cell_df <- data.frame(
    x = as.numeric(embedding[, 1L]),
    y = as.numeric(embedding[, 2L]),
    stringsAsFactors = FALSE
  )
  color_vals <- .resolve_color_by(pv, group_by)
  is_discrete <- !is.null(color_vals) && !is.numeric(color_vals)
  discrete_palette <- NULL

  if (!is.null(color_vals)) {
    if (length(color_vals) != nrow(cell_df)) {
      stop(
        "Resolved `group_by` values must have length equal to the number of cells.",
        call. = FALSE
      )
    }
    if (is_discrete) {
      normalized_colors <- .normalize_discrete_colors(
        color_vals,
        palette = input_color
      )
      cell_df$color <- normalized_colors$values
      discrete_palette <- normalized_colors$palette
    } else {
      cell_df$color <- as.numeric(color_vals)
    }
  }

  legend_name <- if (is.character(group_by) && length(group_by) == 1L) {
    group_by
  } else {
    NULL
  }

  plot <- ggplot2::ggplot()
  point_layer <- function(mapping, fixed_colour = NULL) {
    arguments <- list(
      data = cell_df,
      mapping = mapping,
      size = point_size,
      alpha = alpha,
      shape = point_shape,
      stroke = 0,
      inherit.aes = FALSE
    )
    if (!is.null(fixed_colour)) {
      arguments$colour <- fixed_colour
    }
    if (raster_points) {
      arguments$raster.dpi <- raster_dpi
      do.call(ggrastr::geom_point_rast, arguments)
    } else {
      do.call(ggplot2::geom_point, arguments)
    }
  }

  if (is.null(color_vals)) {
    plot <- plot + point_layer(
      mapping = ggplot2::aes(x = .data$x, y = .data$y),
      fixed_colour = "#AAAAAA"
    )
  } else {
    plot <- plot + point_layer(
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$y,
        colour = .data$color
      )
    )

    if (is_discrete) {
      legend_guide <- ggplot2::guide_legend(
        override.aes = list(
          size = 3,
          alpha = 1,
          shape = point_shape,
          stroke = 0
        )
      )
      if (!is.null(discrete_palette)) {
        plot <- plot + ggplot2::scale_colour_manual(
          values = discrete_palette,
          name = legend_name,
          drop = FALSE,
          guide = legend_guide
        )
      } else {
        plot <- plot + ggplot2::scale_colour_discrete(
          name = legend_name,
          drop = FALSE,
          guide = legend_guide
        )
      }
    } else {
      gradient_colors <- input_color %||%
        c("#440154", "#31688E", "#35B779", "#FDE725")
      plot <- plot + ggplot2::scale_colour_gradientn(
        colours = gradient_colors,
        name = legend_name
      )
    }
  }

  if (nrow(arrow_data) == 0L) {
    warning(
      "No grid arrows passed the velocity field filters.",
      call. = FALSE
    )
  } else {
    arrow_layer <- ggplot2::geom_segment(
      data = arrow_data,
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$y,
        xend = .data$xend,
        yend = .data$yend
      ),
      colour = "grey0",
      alpha = arrow_alpha,
      linewidth = arrow_linewidth,
      lineend = "round",
      linejoin = "mitre",
      arrow = grid::arrow(
        length = grid::unit(arrow_size, "cm"),
        angle = arrow_angle,
        type = arrow_type,
        ends = "last"
      ),
      show.legend = FALSE,
      inherit.aes = FALSE
    )
    if (raster_arrows) {
      arrow_layer <- ggrastr::rasterise(arrow_layer, dpi = raster_dpi)
    }
    plot <- plot + arrow_layer
  }

  dims <- paste0(toupper(sub("^[Xx]_", "", reduction)), "_", 1:2)
  plot <- plot +
    ggplot2::labs(title = title, x = dims[1L], y = dims[2L]) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = 0.025)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = 0.025)
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.45
      ),
      axis.line = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = 13,
        colour = "black",
        hjust = 0.5
      ),
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title = ggplot2::element_text(size = 12, colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(1.5, "mm"),
      legend.title = ggplot2::element_text(size = 10, colour = "black"),
      legend.text = ggplot2::element_text(size = 9, colour = "black"),
      legend.key.size = grid::unit(0.4, "cm"),
      legend.spacing.y = grid::unit(0.02, "cm"),
      legend.box.spacing = grid::unit(2, "mm"),
      plot.margin = ggplot2::margin(
        t = 4,
        r = 5,
        b = 4,
        l = 4,
        unit = "mm"
      )
    )

  if (fixed_aspect) {
    plot <- plot + ggplot2::coord_fixed(ratio = 1, clip = "off")
  }
  if (!legend) {
    plot <- plot + ggplot2::theme(legend.position = "none")
  }

  plot
}
