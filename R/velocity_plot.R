# =============================================================================
# PlantVelocity: velocity_plot.R
#
# Plotting functions for PlantVelocity:
#   1. Plot cell embeddings with RNA velocity arrows
#   2. Plot 2-D phase portraits for spliced/unspliced dynamics
#   3. Plot latent time, root cells, and end-point probabilities
# =============================================================================

#' @title Plot Velocity Embedding
#'
#' @description
#' Plot cells in a two-dimensional embedding overlaid with RNA velocity arrows.
#' Cell points can be rasterized at high resolution while arrows, axes, legends,
#' and text remain vector graphics.
#'
#' @param pv A plantvelo object.
#' @param reduction Embedding name, such as "umap", "pca", or "tsne".
#' @param vkey Velocity layer key.
#' @param group_by Metadata column, numeric vector, or NULL.
#' @param show_arrows Logical. Whether to draw velocity arrows.
#'
#' @param density Numeric scalar in \code{[0, 1]}. Fraction of valid, non-zero velocity
#'   vectors displayed as arrows. For example, 0.005 displays approximately
#'   0.5 percent of valid velocity vectors. A value of 1 displays all valid
#'   arrows, while 0 displays no arrows.
#'
#' @param arrow_scale Global arrow-length multiplier.
#' @param arrow_min_length Minimum visible arrow length as a fraction of the
#'   total embedding span.
#' @param arrow_max_length Maximum visible arrow length as a fraction of the
#'   total embedding span.
#' @param arrow_size Arrow-head length in cm.
#' @param arrow_linewidth Arrow shaft linewidth.
#' @param arrow_alpha Arrow opacity.
#' @param arrow_angle Arrow-head half-angle in degrees. Smaller values produce
#'   narrower and more directional arrowheads.
#' @param arrow_type Either "closed" or "open".
#' @param arrow_colour Fixed arrow colour. NULL uses group colours for a
#'   discrete group_by variable and black otherwise.
#'
#' @param point_size Cell-point size.
#' @param alpha Cell-point opacity.
#' @param point_shape Cell-point shape.
#'
#' @param raster_points Whether to rasterize the cell-point layer.
#' @param raster_arrows Whether to rasterize the arrow layer.
#' @param raster_dpi Internal rasterization resolution.
#'
#' @param input_color Custom discrete or continuous colour palette.
#' @param title Plot title.
#' @param legend Whether to show the legend.
#' @param fixed_aspect Whether to use an equal x/y coordinate ratio.
#' @param seed Non-negative integer random seed used for reproducible arrow
#'   subsampling.
#'
#' @return A ggplot object.
#'
#' @export
plot_velocity_embedding <- function(
    pv,
    reduction   = "umap",
    vkey        = "velocity",
    group_by    = NULL,
    show_arrows = TRUE,

    # ==========================================================
    # Arrow density
    # ==========================================================
    density = 0.10,

    # ==========================================================
    # Arrow parameters
    # ==========================================================
    arrow_scale      = 0.6,
    arrow_min_length = 0.008,
    arrow_max_length = 0.035,
    arrow_size       = 0.08,
    arrow_linewidth  = 0.15,
    arrow_alpha      = 0.90,
    arrow_angle      = 20,
    arrow_type       = c("closed", "open"),
    arrow_colour     = NULL,

    # ==========================================================
    # Point parameters
    # ==========================================================
    point_size  = 0.25,
    alpha       = 0.90,
    point_shape = 16,

    # ==========================================================
    # Rasterization parameters
    # ==========================================================
    raster_points = TRUE,
    raster_arrows = FALSE,
    raster_dpi    = 1200,

    # ==========================================================
    # Other plot parameters
    # ==========================================================
    input_color  = NULL,
    title        = NULL,
    legend       = TRUE,
    fixed_aspect = FALSE,

    # Reproducible arrow sampling
    seed = 42L
) {

  # ============================================================
  # 1. Package and object checks
  # ============================================================

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required. ",
      "Install it with: install.packages('ggplot2')",
      call. = FALSE
    )
  }

  if ((isTRUE(raster_points) || isTRUE(raster_arrows)) &&
      !requireNamespace("ggrastr", quietly = TRUE)) {
    stop(
      "Rasterization requires package 'ggrastr'. ",
      "Install it with: install.packages('ggrastr')",
      call. = FALSE
    )
  }

  if (!inherits(pv, "plantvelo")) {
    stop(
      "`pv` must be a plantvelo object.",
      call. = FALSE
    )
  }

  arrow_type <- match.arg(arrow_type)

  # ============================================================
  # 2. Parameter validation
  # ============================================================

  if (!is.logical(show_arrows) ||
      length(show_arrows) != 1L ||
      is.na(show_arrows)) {
    stop(
      "`show_arrows` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.numeric(density) ||
      length(density) != 1L ||
      !is.finite(density) ||
      density < 0 ||
      density > 1) {
    stop(
      "`density` must be one finite number between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      seed < 0 ||
      seed > .Machine$integer.max ||
      seed != floor(seed)) {
    stop(
      "`seed` must be one non-negative integer no greater than ",
      "`.Machine$integer.max`.",
      call. = FALSE
    )
  }

  seed <- as.integer(seed)

  if (!is.numeric(point_size) ||
      length(point_size) != 1L ||
      !is.finite(point_size) ||
      point_size < 0) {
    stop(
      "`point_size` must be a non-negative finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(alpha) ||
      length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha < 0 ||
      alpha > 1) {
    stop(
      "`alpha` must be between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(arrow_scale) ||
      length(arrow_scale) != 1L ||
      !is.finite(arrow_scale) ||
      arrow_scale < 0) {
    stop(
      "`arrow_scale` must be a non-negative finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(arrow_size) ||
      length(arrow_size) != 1L ||
      !is.finite(arrow_size) ||
      arrow_size < 0) {
    stop(
      "`arrow_size` must be a non-negative finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(arrow_linewidth) ||
      length(arrow_linewidth) != 1L ||
      !is.finite(arrow_linewidth) ||
      arrow_linewidth < 0) {
    stop(
      "`arrow_linewidth` must be a non-negative finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(arrow_alpha) ||
      length(arrow_alpha) != 1L ||
      !is.finite(arrow_alpha) ||
      arrow_alpha < 0 ||
      arrow_alpha > 1) {
    stop(
      "`arrow_alpha` must be between 0 and 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(arrow_angle) ||
      length(arrow_angle) != 1L ||
      !is.finite(arrow_angle) ||
      arrow_angle <= 0 ||
      arrow_angle >= 90) {
    stop(
      "`arrow_angle` must be one finite number between 0 and 90.",
      call. = FALSE
    )
  }

  if (!is.numeric(arrow_min_length) ||
      !is.numeric(arrow_max_length) ||
      length(arrow_min_length) != 1L ||
      length(arrow_max_length) != 1L ||
      !is.finite(arrow_min_length) ||
      !is.finite(arrow_max_length) ||
      arrow_min_length < 0 ||
      arrow_max_length < arrow_min_length) {
    stop(
      "`arrow_min_length` and `arrow_max_length` must be valid, ",
      "and `arrow_max_length` must be greater than or equal to ",
      "`arrow_min_length`.",
      call. = FALSE
    )
  }

  if (!is.numeric(raster_dpi) ||
      length(raster_dpi) != 1L ||
      !is.finite(raster_dpi) ||
      raster_dpi <= 0) {
    stop(
      "`raster_dpi` must be one positive number.",
      call. = FALSE
    )
  }

  # ============================================================
  # 3. Helper for reproducible sampling
  #
  # Restore the original global RNG state after sampling so that
  # plotting does not affect downstream random-number operations.
  # ============================================================

  sample_with_seed <- function(x, size, seed) {

    seed_existed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )

    if (seed_existed) {
      old_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }

    on.exit(
      {
        if (seed_existed) {

          assign(
            ".Random.seed",
            old_seed,
            envir = .GlobalEnv
          )

        } else if (exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )) {

          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )

    set.seed(seed)

    sample(
      x,
      size = size,
      replace = FALSE
    )
  }

  # ============================================================
  # 4. Retrieve embedding coordinates
  # ============================================================

  E <- .get_embedding(
    pv,
    reduction
  )

  if (!is.matrix(E) && !is.data.frame(E)) {
    stop(
      "Embedding returned by `.get_embedding()` must be ",
      "a matrix or data.frame.",
      call. = FALSE
    )
  }

  E <- as.matrix(E)

  if (ncol(E) < 2L) {
    stop(
      "Embedding must have at least two columns.",
      call. = FALSE
    )
  }

  if (nrow(E) == 0L) {
    stop(
      "The selected embedding contains no cells.",
      call. = FALSE
    )
  }

  dims <- paste0(
    toupper(
      sub(
        "^[Xx]_",
        "",
        reduction
      )
    ),
    "_",
    1:2
  )

  df <- data.frame(
    x = as.numeric(E[, 1]),
    y = as.numeric(E[, 2]),
    stringsAsFactors = FALSE
  )

  # ============================================================
  # 5. Retrieve or calculate velocity embedding
  # ============================================================

  reduction_clean <- tolower(
    sub(
      "^[Xx]_",
      "",
      reduction
    )
  )

  emb_key <- paste0(
    vkey,
    "_",
    reduction_clean
  )

  alternative_emb_key <- paste0(
    vkey,
    "_",
    reduction
  )

  V_emb <- pv@reductions[[emb_key]] %||%
    pv@reductions[[alternative_emb_key]]

  if (is.null(V_emb) && isTRUE(show_arrows) && density > 0) {

    pv <- compute_velocity_embedding(
      pv,
      reduction = reduction,
      vkey    = vkey,
      verbose = FALSE
    )

    V_emb <- pv@reductions[[emb_key]] %||%
      pv@reductions[[alternative_emb_key]]
  }

  if (!is.null(V_emb)) {

    V_emb <- as.matrix(V_emb)

    if (nrow(V_emb) != nrow(E) ||
        ncol(V_emb) < 2L) {
      stop(
        "Velocity embedding has incompatible dimensions ",
        "with the selected embedding.",
        call. = FALSE
      )
    }
  }

  # ============================================================
  # 6. Resolve colour mapping
  # ============================================================

  color_vals <- .resolve_color_by(
    pv,
    group_by
  )

  is_discrete <- !is.null(color_vals) &&
    !is.numeric(color_vals)

  discrete_palette <- NULL

  if (!is.null(color_vals)) {

    if (length(color_vals) != nrow(df)) {
      stop(
        "Resolved `group_by` values must have length equal ",
        "to the number of cells.",
        call. = FALSE
      )
    }

    if (is_discrete) {

      normalized_colors <- .normalize_discrete_colors(
        color_vals,
        palette = input_color
      )

      df$color <- normalized_colors$values
      discrete_palette <- normalized_colors$palette

    } else {

      df$color <- as.numeric(color_vals)
    }
  }

  legend_name <- if (
    is.character(group_by) &&
    length(group_by) == 1L
  ) {
    group_by
  } else {
    NULL
  }

  # ============================================================
  # 7. Initialize plot
  # ============================================================

  p <- ggplot2::ggplot()

  # ============================================================
  # 8. Internal helper for adding point layers
  # ============================================================

  add_point_layer <- function(
    mapping,
    fixed_colour = NULL
  ) {

    point_arguments <- list(
      data        = df,
      mapping     = mapping,
      size        = point_size,
      alpha       = alpha,
      shape       = point_shape,
      stroke      = 0,
      inherit.aes = FALSE
    )

    if (!is.null(fixed_colour)) {
      point_arguments$colour <- fixed_colour
    }

    if (isTRUE(raster_points)) {

      point_arguments$raster.dpi <- raster_dpi

      do.call(
        ggrastr::geom_point_rast,
        point_arguments
      )

    } else {

      do.call(
        ggplot2::geom_point,
        point_arguments
      )
    }
  }

  # ============================================================
  # 9. Add cell-point layer
  # ============================================================

  if (is.null(color_vals)) {

    p <- p +
      add_point_layer(
        mapping = ggplot2::aes(
          x = .data$x,
          y = .data$y
        ),
        fixed_colour = "#AAAAAA"
      )

  } else {

    p <- p +
      add_point_layer(
        mapping = ggplot2::aes(
          x      = .data$x,
          y      = .data$y,
          colour = .data$color
        )
      )

    if (is_discrete) {

      legend_guide <- ggplot2::guide_legend(
        override.aes = list(
          size   = 3,
          alpha  = 1,
          shape  = point_shape,
          stroke = 0
        )
      )

      if (!is.null(discrete_palette)) {

        p <- p +
          ggplot2::scale_colour_manual(
            values = discrete_palette,
            name   = legend_name,
            drop   = FALSE,
            guide  = legend_guide
          )

      } else {

        p <- p +
          ggplot2::scale_colour_discrete(
            name  = legend_name,
            drop  = FALSE,
            guide = legend_guide
          )
      }

    } else {

      gradient_colors <- input_color %||%
        c(
          "#440154",
          "#31688E",
          "#35B779",
          "#FDE725"
        )

      p <- p +
        ggplot2::scale_colour_gradientn(
          colours = gradient_colors,
          name    = legend_name
        )
    }
  }

  # ============================================================
  # 10. Calculate and add velocity arrows
  # ============================================================

  arrow_indices <- integer(0)

  if (isTRUE(show_arrows) &&
      !is.null(V_emb) &&
      density > 0) {

    epsilon <- 1e-12

    # ----------------------------------------------------------
    # Retrieve all velocity vectors
    # ----------------------------------------------------------

    vx_all <- as.numeric(
      V_emb[, 1]
    )

    vy_all <- as.numeric(
      V_emb[, 2]
    )

    magnitude_all <- sqrt(
      vx_all^2 + vy_all^2
    )

    # ----------------------------------------------------------
    # Identify finite, non-zero velocity vectors
    # ----------------------------------------------------------

    valid_velocity <- is.finite(vx_all) &
      is.finite(vy_all) &
      is.finite(magnitude_all) &
      magnitude_all > epsilon

    valid_idx <- which(
      valid_velocity
    )

    n_valid <- length(
      valid_idx
    )

    if (n_valid > 0L) {

      # --------------------------------------------------------
      # Convert density into arrow number
      # --------------------------------------------------------

      n_draw <- as.integer(
        ceiling(n_valid * density)
      )

      n_draw <- max(
        1L,
        min(n_draw, n_valid)
      )

      # --------------------------------------------------------
      # Reproducible random sampling
      # --------------------------------------------------------

      arrow_indices <- if (n_draw < n_valid) {

        sample_with_seed(
          x    = valid_idx,
          size = n_draw,
          seed = seed
        )

      } else {

        valid_idx
      }

      arr_idx <- arrow_indices

      vx <- vx_all[arr_idx]
      vy <- vy_all[arr_idx]

      magnitude <- magnitude_all[arr_idx]

      # --------------------------------------------------------
      # Determine embedding span
      # --------------------------------------------------------

      x_range <- range(
        E[, 1],
        finite = TRUE
      )

      y_range <- range(
        E[, 2],
        finite = TRUE
      )

      embedding_span <- max(
        diff(x_range),
        diff(y_range)
      )

      if (!is.finite(embedding_span) ||
          embedding_span <= 0) {
        embedding_span <- 1
      }

      # --------------------------------------------------------
      # Calculate the magnitude reference from all valid vectors
      #
      # This keeps arrow-length scaling stable when density or
      # seed changes.
      # --------------------------------------------------------

      magnitude_valid <- magnitude_all[
        valid_idx
      ]

      magnitude_cap <- stats::quantile(
        magnitude_valid,
        probs = 0.95,
        names = FALSE,
        na.rm = TRUE
      )

      if (!is.finite(magnitude_cap) ||
          magnitude_cap <= epsilon) {

        magnitude_cap <- max(
          magnitude_valid,
          na.rm = TRUE
        )
      }

      magnitude_clipped <- pmin(
        magnitude,
        magnitude_cap
      )

      magnitude_normalized <-
        magnitude_clipped / magnitude_cap

      magnitude_normalized[
        !is.finite(magnitude_normalized)
      ] <- 0

      # --------------------------------------------------------
      # Convert velocity magnitude into visible arrow length
      # --------------------------------------------------------

      minimum_length <- embedding_span *
        arrow_min_length

      maximum_length <- embedding_span *
        arrow_max_length

      visible_length <- minimum_length +
        (maximum_length - minimum_length) *
        sqrt(magnitude_normalized)

      visible_length <- visible_length *
        arrow_scale

      # --------------------------------------------------------
      # Unit direction vectors
      # --------------------------------------------------------

      ux <- vx / magnitude
      uy <- vy / magnitude

      dx <- ux * visible_length
      dy <- uy * visible_length

      df_arr <- data.frame(
        x    = E[arr_idx, 1],
        y    = E[arr_idx, 2],
        xend = E[arr_idx, 1] + dx,
        yend = E[arr_idx, 2] + dy,
        stringsAsFactors = FALSE
      )

      # --------------------------------------------------------
      # Arrow-head specification
      # --------------------------------------------------------

      arrow_specification <- grid::arrow(
        length = grid::unit(
          arrow_size,
          "cm"
        ),
        angle = arrow_angle,
        type  = arrow_type,
        ends  = "last"
      )

      # --------------------------------------------------------
      # Determine arrow colours
      # --------------------------------------------------------

      use_group_arrow_colors <-
        is.null(arrow_colour) &&
        !is.null(color_vals) &&
        is_discrete

      if (use_group_arrow_colors) {

        df_arr$color <- df$color[
          arr_idx
        ]

        arrow_layer <- ggplot2::geom_segment(
          data = df_arr,
          mapping = ggplot2::aes(
            x      = .data$x,
            y      = .data$y,
            xend   = .data$xend,
            yend   = .data$yend,
            colour = .data$color
          ),
          arrow       = arrow_specification,
          linewidth   = arrow_linewidth,
          alpha       = arrow_alpha,
          lineend     = "round",
          linejoin    = "mitre",
          show.legend = FALSE,
          inherit.aes = FALSE
        )

      } else {

        fixed_arrow_colour <- arrow_colour %||%
          "black"

        arrow_layer <- ggplot2::geom_segment(
          data = df_arr,
          mapping = ggplot2::aes(
            x    = .data$x,
            y    = .data$y,
            xend = .data$xend,
            yend = .data$yend
          ),
          colour      = fixed_arrow_colour,
          arrow       = arrow_specification,
          linewidth   = arrow_linewidth,
          alpha       = arrow_alpha,
          lineend     = "round",
          linejoin    = "mitre",
          show.legend = FALSE,
          inherit.aes = FALSE
        )
      }

      # --------------------------------------------------------
      # Optionally rasterize arrow layer
      # --------------------------------------------------------

      if (isTRUE(raster_arrows)) {

        arrow_layer <- ggrastr::rasterise(
          arrow_layer,
          dpi = raster_dpi
        )
      }

      p <- p + arrow_layer
    }
  }

  # ============================================================
  # 11. Labels and publication-style theme
  # ============================================================

  p <- p +
    ggplot2::labs(
      title = title,
      x     = dims[1],
      y     = dims[2]
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = 0.025
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = 0.025
      )
    ) +
    ggplot2::theme_classic(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill      = NA,
        colour    = "black",
        linewidth = 0.45
      ),

      axis.line = ggplot2::element_blank(),

      plot.title = ggplot2::element_text(
        size   = 13,
        colour = "black",
        hjust  = 0.5
      ),

      axis.text = ggplot2::element_text(
        size   = 10,
        colour = "black"
      ),

      axis.title = ggplot2::element_text(
        size   = 12,
        colour = "black"
      ),

      axis.ticks = ggplot2::element_line(
        colour    = "black",
        linewidth = 0.35
      ),

      axis.ticks.length = grid::unit(
        1.5,
        "mm"
      ),

      legend.title = ggplot2::element_text(
        size   = 10,
        colour = "black"
      ),

      legend.text = ggplot2::element_text(
        size   = 9,
        colour = "black"
      ),

      legend.key.size = grid::unit(
        0.4,
        "cm"
      ),

      legend.spacing.y = grid::unit(
        0.02,
        "cm"
      ),

      legend.box.spacing = grid::unit(
        2,
        "mm"
      ),

      plot.margin = ggplot2::margin(
        t    = 4,
        r    = 5,
        b    = 4,
        l    = 4,
        unit = "mm"
      )
    )

  if (isTRUE(fixed_aspect)) {

    p <- p +
      ggplot2::coord_fixed(
        ratio = 1,
        clip  = "off"
      )
  }

  if (!isTRUE(legend)) {

    p <- p +
      ggplot2::theme(
        legend.position = "none"
      )
  }

  # Store sampling information in the plot object
  attr(p, "velocity_arrow_indices") <- arrow_indices
  attr(p, "velocity_arrow_number")  <- length(arrow_indices)
  attr(p, "velocity_arrow_density") <- density
  attr(p, "velocity_arrow_seed")    <- seed

  return(p)
}

# =============================================================================
# plot_phase_portrait()
# =============================================================================

#' @title Plot Phase Portrait
#' @description
#' Plot the spliced-versus-unspliced phase portrait for a single gene.
#' Optionally overlays the fitted kinetic trajectory recovered by
#' \code{recover_dynamics()}.
#'
#' Expression values are preferentially read from \code{@moments}; if the
#' requested smoothed modalities are absent, the function falls back to the raw
#' layers in \code{@layers}.
#'
#' @param pv A \code{plantvelo} object.
#' @param gene Character scalar. Gene name.
#' @param modality_s Character scalar. Smoothed spliced layer key in
#'   \code{@moments}. Default \code{"Ms"}.
#' @param modality_u Character scalar. Smoothed unspliced layer key in
#'   \code{@moments}. Default \code{"Mu"}.
#' @param group_by Character scalar. Colouring scheme. Supported values include
#'   \code{"state"} (requires \code{@layers$fit_t}), \code{"latent_time"}, or
#'   any column in \code{@meta.data}. Default \code{"state"}.
#' @param show_fit Logical. Overlay the fitted kinetic curve. Requires
#'   \code{@kinetics$params}. Default \code{TRUE}.
#' @param n_curve Integer. Number of points used to draw the fitted kinetic
#'   curve. Default \code{300L}.
#' @param input_color Character vector or \code{NULL}. Custom colour palette.
#'   \itemize{
#'     \item For \code{group_by = "state"}, named or unnamed discrete colours
#'       override the default induction/repression palette.
#'     \item For continuous \code{group_by} values, the palette is passed to
#'       \code{ggplot2::scale_colour_gradientn()}.
#'     \item For discrete metadata columns, the palette is passed to
#'       \code{ggplot2::scale_colour_manual()}.
#'     \item \code{NULL} uses built-in defaults.
#'   }
#' @param fit_color Character scalar. Colour of the fitted kinetic curve. Default \code{"grey30"}.
#' @param point_size Numeric. Point size. Default \code{2}.
#' @param alpha Numeric. Point opacity. Default \code{0.7}.
#' @param title Character scalar or \code{NULL}. Plot title. Default
#'   \code{NULL}, which uses the gene name.
#' @param legend Logical. Show the legend. Default \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_phase_portrait(pv, gene = "AT1G01010")
#' plot_phase_portrait(
#'   pv,
#'   gene = "AT1G01010",
#'   input_color = c(induction = "navy", repression = "firebrick")
#' )
#' }
#' @export
plot_phase_portrait <- function(pv,
                                gene,
                                modality_s  = "Ms",
                                modality_u  = "Mu",
                                group_by    = "state",
                                show_fit    = TRUE,
                                n_curve     = 300L,
                                input_color = NULL,
                                fit_color = "grey30",
                                point_size  = 1,
                                alpha       = 0.7,
                                title       = NULL,
                                legend      = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  all_genes <- rownames(pv@layers[["spliced"]])
  g_idx     <- which(all_genes == gene)
  if (length(g_idx) == 0)
    stop(sprintf("Gene '%s' not found in @layers$spliced.", gene))

  # ----------------extract spliced and unspliced expression----------------
  S <- if (!is.null(pv@moments[[modality_s]])) {
    pv@moments[[modality_s]][, g_idx]
  } else {
    pv@layers[["spliced"]][g_idx, ]
  }

  U <- if (!is.null(pv@moments[[modality_u]])) {
    pv@moments[[modality_u]][, g_idx]
  } else {
    pv@layers[["unspliced"]][g_idx, ]
  }

  df <- data.frame(s = as.numeric(S), u = as.numeric(U))

  params <- pv@kinetics[["params"]]
  fit_t <- if (identical(group_by, "state")) {
    .get_fit_time_matrix(pv)
  } else {
    NULL
  }
  uses_dynamics <-
    !is.null(fit_t) ||
    (isTRUE(show_fit) && !is.null(params))
  if (uses_dynamics && !is.null(params)) {
    params <- .validate_dynamics_schema(pv)
  }

  # ----------------resolve colour mapping----------------
  if (identical(group_by, "state") && !is.null(fit_t)) {
    t_g <- fit_t[, gene]
    valid_switch <- !is.null(params) && gene %in% rownames(params) &&
      identical(as.character(params[gene, "fit_model"]), "fitted") &&
      is.finite(params[gene, "fit_t_"]) && params[gene, "fit_t_"] >= 0
    finite_time <- is.finite(t_g)
    state <- rep("unknown", length(t_g))
    if (valid_switch && any(finite_time)) {
      state[finite_time] <- ifelse(
        t_g[finite_time] < params[gene, "fit_t_"],
        "induction", "repression"
      )
    }

    phase_levels <- if (valid_switch && any(finite_time)) {
      c("induction", "repression", if (any(!finite_time)) "unknown")
    } else {
      "unknown"
    }
    df$color <- factor(state, levels = phase_levels)

    default_state <- c(
      induction    = "#2166AC",
      repression   = "#D6604D",
      unknown      = "#888888"
    )
    state_colors <- input_color %||% default_state
    if (identical(phase_levels, "unknown")) {
      state_colors <- c(unknown = "#888888")
    } else if ("unknown" %in% phase_levels && !is.null(input_color)) {
      if (is.null(names(state_colors))) {
        if (length(state_colors) >= 2L &&
            length(state_colors) < length(phase_levels)) {
          state_colors <- c(state_colors, "#888888")
        }
      } else if (!("unknown" %in% names(state_colors))) {
        state_colors <- c(state_colors, unknown = "#888888")
      }
    }

    norm <- .normalize_discrete_colors(
      df$color,
      palette = state_colors
    )
    df$color <- norm$values

    color_scale <- ggplot2::scale_colour_manual(
      values = norm$palette,
      name   = "State",
      drop   = TRUE
    )

  } else if (group_by == "latent_time" && !is.null(pv@meta.data$latent_time)) {
    df$color <- pv@meta.data$latent_time
    pal_grad <- input_color %||% c("#440154", "#31688E", "#35B779", "#FDE725")

    color_scale <- ggplot2::scale_colour_gradientn(
      colours = pal_grad,
      name    = "Latent time"
    )

  } else if (group_by %in% colnames(pv@meta.data)) {
    df$color <- pv@meta.data[[group_by]]

    if (is.numeric(df$color)) {
      pal_grad <- input_color %||% c("#440154", "#31688E", "#35B779", "#FDE725")
      color_scale <- ggplot2::scale_colour_gradientn(
        colours = pal_grad,
        name    = group_by
      )
    } else {
      norm <- .normalize_discrete_colors(df$color, palette = input_color)
      df$color <- norm$values

      color_scale <- if (!is.null(norm$palette)) {
        ggplot2::scale_colour_manual(
          values = norm$palette,
          name   = group_by,
          drop   = FALSE
        )
      } else {
        ggplot2::scale_colour_discrete(
          limits = levels(df$color),
          name   = group_by,
          drop   = FALSE
        )
      }
    }

  } else {
    df$color    <- "#888888"
    color_scale <- NULL
  }

  # ----------------build scatter plot----------------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$s, y = .data$u, colour = .data$color)
  ) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::labs(
      title = title %||% gene,
      x     = "spliced",
      y     = "unspliced"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.border    = ggplot2::element_rect(fill = NA, colour = "black",
                                              linewidth = 0.5),
      axis.line       = ggplot2::element_blank(),
      plot.title      = ggplot2::element_text(size = 13, color = "black",
                                              hjust = 0.5),
      axis.text       = ggplot2::element_text(size = 10, color = "black"),
      axis.title      = ggplot2::element_text(size = 12, color = "black"),
      legend.title    = ggplot2::element_text(size = 10, color = "black"),
      legend.text     = ggplot2::element_text(size = 9,  color = "black"),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )

  if (!is.null(color_scale)) {
    p <- p + color_scale
  }

  # ----------------overlay fitted kinetic curve----------------
  if (show_fit && !is.null(params) && gene %in% rownames(params)) {

    curve_df <- .phase_curve(gene, params, n_curve)
    if (!is.null(curve_df)) {
      p <- p + ggplot2::geom_path(
        data        = curve_df,
        mapping     = ggplot2::aes(x = .data$s, y = .data$u),
        colour      = fit_color,
        linewidth   = 0.8,
        inherit.aes = FALSE
      )
    }
  }

  if (!legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  p
}


# =============================================================================
# plot_latent_time() / plot_root_end_points()
# =============================================================================

#' @title Plot Latent Time
#' @description
#' Convenience wrapper around \code{plot_velocity_embedding()} that colours
#' cells by \code{latent_time} and disables velocity arrows.
#'
#' @param pv A \code{plantvelo} object.
#' @param reduction Character scalar. Embedding to plot. Default
#'   \code{"umap"}.
#' @param ... Further arguments passed to \code{plot_velocity_embedding()}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' plot_latent_time(pv)
#' }
#' @export
plot_latent_time <- function(pv, reduction = "umap", ...) {
  plot_velocity_embedding(pv, reduction = reduction,
                          group_by = "latent_time",
                          show_arrows = FALSE, ...)
}


#' @title Plot Root Cells or End Points
#' @description
#' Convenience wrapper around \code{plot_velocity_embedding()} that colours
#' cells by \code{root_cells} or \code{end_points} values stored in
#' \code{@meta.data}.
#'
#' @param pv A \code{plantvelo} object.
#' @param type Character scalar. \code{"root"} (default) or \code{"end"}.
#' @param reduction Character scalar. Embedding to plot. Default
#'   \code{"umap"}.
#' @param ... Further arguments passed to \code{plot_velocity_embedding()}.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' plot_root_end_points(pv, type = "root")
#' plot_root_end_points(pv, type = "end")
#' }
#' @export
plot_root_end_points <- function(pv, type = "root", reduction = "umap", ...) {
  col <- if (type == "root") "root_cells" else "end_points"
  if (is.null(pv@meta.data[[col]]))
    stop(sprintf(
      "'%s' not found in @meta.data. Run compute_terminal_states() first.", col
    ))
  plot_velocity_embedding(pv, reduction = reduction,
                          group_by = col,
                          show_arrows = FALSE, ...)
}

# =============================================================================
# Internal helper functions
# =============================================================================

#' Normalise discrete colour inputs
#'
#' Standardises a discrete colour vector and an optional palette so that colour
#' assignments follow a consistent ordering:
#' \enumerate{
#'   \item For factors, preserve the original factor levels.
#'   \item For non-factors, sort unique values alphabetically.
#' }
#'
#' If \code{palette} is unnamed, colours are assigned by level order. If it is
#' named, it is reordered to match the resolved levels.
#'
#' @param x Discrete vector to be converted to a factor.
#' @param palette Character vector or \code{NULL}. Optional colour palette.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{values}}{Factor with standardised levels.}
#'   \item{\code{palette}}{Resolved named palette or \code{NULL}.}
#'   \item{\code{levels}}{Character vector of resolved levels.}
#' }
#' @keywords internal
.normalize_discrete_colors <- function(x, palette = NULL) {
  # Seurat-like ordering:
  # 1) factor: keep factor levels
  # 2) non-factor: sort alphabetically
  if (is.factor(x)) {
    lev <- levels(x)
    x_fac <- factor(as.character(x), levels = lev)
  } else {
    x_chr <- as.character(x)
    lev <- sort(unique(x_chr))
    x_fac <- factor(x_chr, levels = lev)
  }

  # No custom palette supplied
  if (is.null(palette)) {
    return(list(values = x_fac, palette = NULL, levels = lev))
  }

  # Unnamed palette: assign by level order
  if (is.null(names(palette))) {
    if (length(palette) < length(lev)) {
      stop("`palette` has fewer colours than the number of discrete groups.")
    }
    pal <- palette[seq_along(lev)]
    names(pal) <- lev
  } else {
    # Named palette: reorder by unified levels
    miss <- setdiff(lev, names(palette))
    if (length(miss) > 0) {
      stop(sprintf(
        "Named `palette` is missing colours for: %s",
        paste(miss, collapse = ", ")
      ))
    }
    pal <- palette[lev]
  }

  list(values = x_fac, palette = pal, levels = lev)
}

#' Resolve \code{group_by} values from a plantvelo object
#'
#' Supports three cases:
#' \enumerate{
#'   \item \code{NULL}: return \code{NULL}.
#'   \item Numeric vector of length equal to the number of cells: return it
#'     directly.
#'   \item Character scalar matching a column in \code{@meta.data}: return that
#'     metadata column.
#' }
#'
#' If a character scalar does not match any metadata column, a warning is
#' issued and \code{NULL} is returned.
#'
#' @param pv A \code{plantvelo} object.
#' @param group_by Grouping specification passed by the user.
#'
#' @return A vector of grouping values, or \code{NULL}.
#' @keywords internal
.resolve_color_by <- function(pv, group_by) {
  if (is.null(group_by)) return(NULL)

  # Direct numeric vector
  if (is.numeric(group_by) && length(group_by) == ncol(pv@layers[["spliced"]]))
    return(group_by)

  # @meta.data column name
  if (is.character(group_by) && length(group_by) == 1) {
    if (group_by %in% colnames(pv@meta.data))
      return(pv@meta.data[[group_by]])
    warning(sprintf("'%s' not found in @meta.data. Ignoring group_by.", group_by))
    return(NULL)
  }

  NULL
}


#' Compute a fitted kinetic curve for the phase portrait
#'
#' Uniformly samples time points and evaluates the two-state ODE solutions
#' using fitted kinetic parameters.
#'
#' @param gene Character scalar. Gene name.
#' @param params data.frame of fitted kinetic parameters.
#' @param n_curve Integer. Number of time points used for the fitted curve.
#'
#' @return A \code{data.frame} with columns \code{u}, \code{s}, and
#'   \code{phase}, or
#'   \code{NULL} if the fitted curve cannot be constructed.
#' @keywords internal
.phase_curve <- function(gene, params, n_curve) {
  if (is.null(params) || is.null(rownames(params)) ||
      !(gene %in% rownames(params))) {
    return(NULL)
  }
  .reject_velocity_legacy_schema(params)

  required <- c(
    "fit_model", "fit_alpha", "fit_beta", "fit_gamma", "fit_t_",
    "fit_scaling"
  )
  if (any(!required %in% colnames(params))) return(NULL)

  numeric_required <- required[-1L]
  if (any(vapply(
    numeric_required,
    function(name) !is.null(dim(params[[name]])),
    logical(1)
  ))) {
    return(NULL)
  }

  parameters <- params[gene, , drop = FALSE]
  if (!identical(as.character(parameters$fit_model), "fitted")) {
    return(NULL)
  }

  values <- lapply(parameters[numeric_required], function(value) {
    if (is.list(value) && is.null(dim(value)) && length(value) == 1L) {
      value[[1L]]
    } else {
      value
    }
  })
  if (any(lengths(values) != 1L) ||
      any(!vapply(values, is.numeric, logical(1))) ||
      any(!vapply(values, function(value) is.null(dim(value)), logical(1)))) {
    return(NULL)
  }
  if (any(!is.finite(unlist(values))) ||
      values$fit_alpha < 0 || values$fit_beta <= 0 ||
      values$fit_gamma <= 0 || values$fit_t_ < 0 ||
      values$fit_scaling <= 0 ||
      !is.numeric(n_curve) || length(n_curve) != 1L ||
      !is.finite(n_curve) || n_curve < 1L) {
    return(NULL)
  }

  time <- seq(0, values$fit_t_ * 2, length.out = as.integer(n_curve))
  vectorized <- vectorize_2state(
    time,
    values$fit_t_,
    values$fit_alpha,
    values$fit_beta,
    values$fit_gamma
  )
  u_curve <- u_solution(
    vectorized$tau,
    vectorized$u0,
    vectorized$alpha,
    values$fit_beta
  )
  s_curve <- s_solution(
    vectorized$tau,
    vectorized$s0,
    vectorized$u0,
    vectorized$alpha,
    values$fit_beta,
    values$fit_gamma
  )
  if (any(!is.finite(u_curve)) || any(!is.finite(s_curve))) return(NULL)
  data.frame(
    u = pmax(as.numeric(u_curve) * values$fit_scaling, 0),
    s = pmax(as.numeric(s_curve), 0),
    phase = ifelse(time < values$fit_t_, "induction", "repression")
  )
}
