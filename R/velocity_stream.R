# =============================================================================
# PlantVelocity: velocity_stream.R
#
# Algorithm:
#   1. Read two-dimensional cell embeddings and velocity embeddings.
#   2. Interpolate velocities onto a regular grid with Gaussian weights.
#   3. Mask grid points using velocity strength and local velocity length.
#   4. Trace deterministic paths in both directions through the valid field.
#   5. Render the resulting paths with ggplot.
# =============================================================================


.validate_stream_controls <- function(n_grid,
                                      n_neighbors,
                                      smooth,
                                      min_mass,
                                      cutoff_perc,
                                      x_pad,
                                      y_pad,
                                      density,
                                      max_length,
                                      integration_direction) {
  scalar_finite <- function(x) {
    is.numeric(x) && !is.logical(x) && length(x) == 1L && is.finite(x)
  }

  if (!scalar_finite(n_grid) || n_grid < 2 ||
      n_grid > .Machine$integer.max || n_grid != floor(n_grid))
    stop("`n_grid` must be an integer greater than or equal to 2.")
  if (!is.null(n_neighbors) &&
      (!scalar_finite(n_neighbors) ||
       n_neighbors < 1 || n_neighbors > .Machine$integer.max ||
       n_neighbors != floor(n_neighbors)))
    stop("`n_neighbors` must be NULL or a positive integer.")
  if (!scalar_finite(smooth) || smooth <= 0)
    stop("`smooth` must be a positive finite scalar.")
  if (!scalar_finite(density) || density <= 0)
    stop("`density` must be a positive finite scalar.")
  if (!scalar_finite(max_length) || max_length <= 0)
    stop("`max_length` must be a positive finite scalar.")
  if (!scalar_finite(min_mass))
    stop("`min_mass` must be a finite scalar.")
  if (!scalar_finite(cutoff_perc) ||
      cutoff_perc < 0 || cutoff_perc > 100)
    stop("`cutoff_perc` must be between 0 and 100.")
  if (!scalar_finite(x_pad) || x_pad < 0)
    stop("`x_pad` must be a non-negative finite scalar.")
  if (!scalar_finite(y_pad) || y_pad < 0)
    stop("`y_pad` must be a non-negative finite scalar.")

  integration_direction <- match.arg(
    integration_direction,
    choices = c("both", "forward", "backward")
  )

  list(
    n_grid = as.integer(n_grid),
    n_neighbors = if (is.null(n_neighbors)) NULL else as.integer(n_neighbors),
    smooth = as.numeric(smooth),
    min_mass = as.numeric(min_mass),
    cutoff_perc = as.numeric(cutoff_perc),
    x_pad = as.numeric(x_pad),
    y_pad = as.numeric(y_pad),
    density = as.numeric(density),
    max_length = as.numeric(max_length),
    integration_direction = integration_direction
  )
}


.estimate_stream_grid <- function(coords,
                                  velocity,
                                  n_grid,
                                  n_neighbors,
                                  smooth,
                                  min_mass,
                                  cutoff_perc,
                                  x_pad,
                                  y_pad) {
  is_numeric_matrix <- function(x) {
    is.matrix(x) && typeof(x) %in% c("double", "integer") &&
      !is.factor(x) && is.null(attr(x, "levels"))
  }

  if (!is_numeric_matrix(coords) || !is_numeric_matrix(velocity) ||
      ncol(coords) < 2L ||
      ncol(velocity) < 2L || nrow(coords) != nrow(velocity))
    stop("`coords` and `velocity` must be aligned numeric matrices with at least two columns.")

  coords <- coords[, 1:2, drop = FALSE]
  velocity <- velocity[, 1:2, drop = FALSE]
  keep <- apply(is.finite(coords), 1L, all) &
    apply(is.finite(velocity), 1L, all)
  coords <- coords[keep, , drop = FALSE]
  velocity <- velocity[keep, , drop = FALSE]
  n_cells <- nrow(coords)

  if (n_cells < 5L)
    stop("At least five cells with finite coordinates and velocities are required.")

  xr <- range(coords[, 1L])
  yr <- range(coords[, 2L])
  if (diff(xr) == 0 || diff(yr) == 0)
    stop("Both coordinate axes must have non-zero ranges.")

  x_range <- xr + c(-1, 1) * diff(xr) * x_pad
  y_range <- yr + c(-1, 1) * diff(yr) * y_pad
  gx <- seq(x_range[1L], x_range[2L], length.out = as.integer(n_grid))
  gy <- seq(y_range[1L], y_range[2L], length.out = as.integer(n_grid))
  grid <- expand.grid(x = gx, y = gy, KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)

  k <- if (is.null(n_neighbors)) {
    max(1L, floor(n_cells / 50L))
  } else {
    min(as.integer(n_neighbors), n_cells)
  }
  bandwidth <- mean(c(diff(gx)[1L], diff(gy)[1L])) * smooth
  n_pts <- nrow(grid)
  bytes_per_row <- n_cells * 8 * 4
  chunk_size <- as.integer(max(
    1L,
    min(500L, floor(64 * 1024^2 / bytes_per_row))
  ))
  dx <- rep(NA_real_, n_pts)
  dy <- rep(NA_real_, n_pts)
  mass <- rep(NA_real_, n_pts)
  local_length <- rep(NA_real_, n_pts)

  for (chunk_start in seq(1L, n_pts, by = chunk_size)) {
    chunk_end <- min(chunk_start + chunk_size - 1L, n_pts)
    idx <- chunk_start:chunk_end
    d2 <- outer(grid$x[idx], coords[, 1L], `-`)^2 +
      outer(grid$y[idx], coords[, 2L], `-`)^2

    for (i in seq_along(idx)) {
      d2_row <- d2[i, ]
      kth <- sort.int(d2_row, partial = k, na.last = NA)[k]
      candidates <- which(d2_row <= kth)
      nn <- candidates[order(d2_row[candidates], candidates)][seq_len(k)]
      distance <- sqrt(d2[i, nn])
      weights <- stats::dnorm(distance, sd = bandwidth)
      mass[idx[i]] <- sum(weights)
      dx[idx[i]] <- sum(weights * velocity[nn, 1L]) / max(1, mass[idx[i]])
      dy[idx[i]] <- sum(weights * velocity[nn, 2L]) / max(1, mass[idx[i]])
      local_length[idx[i]] <- sum(colMeans(abs(velocity[nn, , drop = FALSE])))
    }
  }

  speed <- sqrt(dx^2 + dy^2)
  finite_speed <- speed[is.finite(speed)]
  speed_cutoff <- min(10^(min_mass - 6), max(finite_speed) * 0.9)
  length_cutoff <- stats::quantile(
    local_length[is.finite(local_length)],
    probs = cutoff_perc / 100,
    names = FALSE
  )
  valid <- is.finite(speed) & speed >= speed_cutoff &
    is.finite(local_length) & local_length >= length_cutoff
  dx[!valid] <- NA_real_
  dy[!valid] <- NA_real_

  out <- data.frame(
    x = grid$x,
    y = grid$y,
    dx = dx,
    dy = dy,
    mass = mass,
    local_length = local_length,
    valid = valid,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$y, out$x), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n_neighbors") <- k
  out
}


.make_stream_field <- function(grid) {
  required <- c("x", "y", "dx", "dy", "valid")
  if (!is.data.frame(grid) || !all(required %in% names(grid)))
    stop("`grid` must contain x, y, dx, dy, and valid columns.")

  if (!is.numeric(grid$x) || !is.numeric(grid$y) ||
      any(!is.finite(grid$x)) || any(!is.finite(grid$y)))
    stop("`grid$x` and `grid$y` must contain only finite numeric coordinates.")

  gx <- sort(unique(grid$x))
  gy <- sort(unique(grid$y))
  nx <- length(gx)
  ny <- length(gy)

  if (nx < 2L || ny < 2L)
    stop("`grid` must contain at least two unique x and y coordinates.")
  if (nrow(grid) != nx * ny)
    stop("`grid` must contain every x/y coordinate combination exactly once.")
  if (anyDuplicated(grid[c("x", "y")]))
    stop("`grid` cannot contain duplicated x/y coordinate pairs.")

  is_regular_spacing <- function(values) {
    expected <- seq(min(values), max(values), length.out = length(values))
    spacing <- diff(range(values)) / (length(values) - 1L)
    coord_scale <- max(abs(c(values, expected)))
    tolerance <- max(
      8 * .Machine$double.eps * max(coord_scale, .Machine$double.xmin),
      1e-8 * abs(spacing)
    )
    all(abs(values - expected) <= tolerance)
  }
  if (!is_regular_spacing(gx) || !is_regular_spacing(gy))
    stop("`grid` x and y coordinates must each be evenly spaced.")

  grid <- grid[order(grid$y, grid$x), , drop = FALSE]
  x_span <- diff(range(gx))
  y_span <- diff(range(gy))

  if (!is.finite(x_span) || !is.finite(y_span) ||
      x_span <= 0 || y_span <= 0)
    stop("`grid` must describe a complete two-dimensional grid with non-zero spans.")

  list(
    gx = gx,
    gy = gy,
    dx = matrix(as.numeric(grid$dx), nrow = nx, ncol = ny),
    dy = matrix(as.numeric(grid$dy), nrow = nx, ncol = ny),
    valid = matrix(as.logical(grid$valid), nrow = nx, ncol = ny),
    x_min = gx[1L],
    y_min = gy[1L],
    x_span = x_span,
    y_span = y_span,
    nx = nx,
    ny = ny
  )
}


.sample_stream_field <- function(field, point) {
  valid_field <- is.list(field) &&
    all(c("dx", "dy", "valid", "x_span", "y_span", "nx", "ny") %in%
        names(field)) &&
    is.numeric(field$nx) && length(field$nx) == 1L &&
    is.numeric(field$ny) && length(field$ny) == 1L &&
    is.finite(field$nx) && is.finite(field$ny) &&
    field$nx >= 2L && field$ny >= 2L &&
    field$nx == floor(field$nx) && field$ny == floor(field$ny) &&
    is.numeric(field$x_span) && length(field$x_span) == 1L &&
    is.numeric(field$y_span) && length(field$y_span) == 1L &&
    is.finite(field$x_span) && is.finite(field$y_span) &&
    field$x_span > 0 && field$y_span > 0
  if (!valid_field)
    return(NULL)

  nx <- as.integer(field$nx)
  ny <- as.integer(field$ny)
  if (!is.matrix(field$dx) || !is.matrix(field$dy) || !is.matrix(field$valid) ||
      !identical(dim(field$dx), c(nx, ny)) ||
      !identical(dim(field$dy), c(nx, ny)) ||
      !identical(dim(field$valid), c(nx, ny)))
    return(NULL)

  if (!is.numeric(point) || length(point) != 2L || any(!is.finite(point)) ||
      any(point < 0 | point > 1))
    return(NULL)
  point <- unname(as.numeric(point))

  x_index <- point[1L] * (nx - 1L) + 1
  y_index <- point[2L] * (ny - 1L) + 1
  ix <- min(as.integer(floor(x_index)), nx - 1L)
  iy <- min(as.integer(floor(y_index)), ny - 1L)
  tx <- x_index - ix
  ty <- y_index - iy

  corner_valid <- c(
    field$valid[ix, iy], field$valid[ix + 1L, iy],
    field$valid[ix, iy + 1L], field$valid[ix + 1L, iy + 1L]
  )
  if (!isTRUE(all(corner_valid)))
    return(NULL)

  bilinear <- function(values) {
    values[ix, iy] * (1 - tx) * (1 - ty) +
      values[ix + 1L, iy] * tx * (1 - ty) +
      values[ix, iy + 1L] * (1 - tx) * ty +
      values[ix + 1L, iy + 1L] * tx * ty
  }
  vx <- bilinear(field$dx)
  vy <- bilinear(field$dy)
  speed <- sqrt(vx^2 + vy^2)
  direction <- unname(c(vx / field$x_span, vy / field$y_span))
  direction_norm <- sqrt(sum(direction^2))

  if (!is.finite(vx) || !is.finite(vy) || !is.finite(speed) ||
      !is.finite(direction_norm) || direction_norm <= .Machine$double.eps)
    return(NULL)

  c(
    direction_x = direction[1L] / direction_norm,
    direction_y = direction[2L] / direction_norm,
    speed = unname(speed)
  )
}


.rk4_stream_step <- function(field, point, step_size, direction_sign) {
  if (!is.numeric(step_size) || length(step_size) != 1L ||
      !is.finite(step_size) || step_size <= 0 ||
      !is.numeric(direction_sign) || length(direction_sign) != 1L ||
      !is.finite(direction_sign) || !(direction_sign %in% c(-1, 1)))
    return(NULL)

  direction_at <- function(location) {
    sampled <- .sample_stream_field(field, location)
    if (is.null(sampled))
      return(NULL)
    direction_sign * unname(sampled[c("direction_x", "direction_y")])
  }

  k1 <- direction_at(point)
  if (is.null(k1)) return(NULL)
  k2 <- direction_at(point + 0.5 * step_size * k1)
  if (is.null(k2)) return(NULL)
  k3 <- direction_at(point + 0.5 * step_size * k2)
  if (is.null(k3)) return(NULL)
  k4 <- direction_at(point + step_size * k3)
  if (is.null(k4)) return(NULL)

  next_point <- point + step_size * (k1 + 2 * k2 + 2 * k3 + k4) / 6
  endpoint <- .sample_stream_field(field, next_point)
  if (is.null(endpoint))
    return(NULL)

  list(point = as.numeric(next_point), speed = unname(endpoint["speed"]))
}


.adaptive_stream_step <- function(field, point, base_step, direction_sign) {
  if (!is.numeric(base_step) || length(base_step) != 1L ||
      !is.finite(base_step) || base_step <= 0)
    return(NULL)

  start <- .sample_stream_field(field, point)
  if (is.null(start))
    return(NULL)

  minimum_step <- base_step / 8
  step_size <- base_step
  repeat {
    candidate <- .rk4_stream_step(field, point, step_size, direction_sign)
    if (is.null(candidate))
      return(NULL)

    endpoint <- .sample_stream_field(field, candidate$point)
    if (is.null(endpoint))
      return(NULL)

    cosine <- sum(
      start[c("direction_x", "direction_y")] *
        endpoint[c("direction_x", "direction_y")]
    )
    angle <- acos(min(1, max(-1, cosine)))
    if (!is.finite(angle) || angle <= pi / 6 || step_size <= minimum_step)
      return(candidate)

    step_size <- max(minimum_step, step_size / 2)
  }
}


.trace_stream_direction <- function(field,
                                    seed,
                                    direction_sign,
                                    length_limit,
                                    occupied) {
  if (!is.list(field) || !is.numeric(field$nx) || !is.numeric(field$ny) ||
      length(field$nx) != 1L || length(field$ny) != 1L ||
      !is.finite(field$nx) || !is.finite(field$ny) ||
      field$nx < 2L || field$ny < 2L ||
      !is.numeric(length_limit) || length(length_limit) != 1L ||
      !is.finite(length_limit) || length_limit <= 0)
    return(NULL)

  nx <- as.integer(field$nx)
  ny <- as.integer(field$ny)
  base_step <- 0.5 * min(1 / (nx - 1L), 1 / (ny - 1L))
  start <- .sample_stream_field(field, seed)
  if (is.null(start))
    return(NULL)

  cell_key <- function(point) {
    ix <- min(as.integer(floor(point[1L] * (nx - 1L))) + 1L, nx - 1L)
    iy <- min(as.integer(floor(point[2L] * (ny - 1L))) + 1L, ny - 1L)
    paste(ix, iy, sep = ":")
  }
  is_occupied <- function(point) {
    if (!is.matrix(occupied) || nrow(occupied) == 0L || ncol(occupied) == 0L)
      return(FALSE)
    ox <- min(as.integer(floor(point[1L] * nrow(occupied))) + 1L,
              nrow(occupied))
    oy <- min(as.integer(floor(point[2L] * ncol(occupied))) + 1L,
              ncol(occupied))
    isTRUE(occupied[ox, oy])
  }

  path <- matrix(c(seed[1L], seed[2L], unname(start["speed"])),
                 ncol = 3L,
                 dimnames = list(NULL, c("x", "y", "speed")))
  visited <- cell_key(seed)
  travelled <- 0
  raw_steps <- ceiling(length_limit / (base_step / 8)) + 1
  max_steps <- if (is.finite(raw_steps)) {
    as.integer(max(1, min(raw_steps, 1000000)))
  } else {
    1000000L
  }

  for (step_number in seq_len(max_steps)) {
    next_step <- .adaptive_stream_step(
      field = field,
      point = path[nrow(path), c("x", "y")],
      base_step = base_step,
      direction_sign = direction_sign
    )
    if (is.null(next_step) || any(!is.finite(next_step$point)) ||
        any(next_step$point < 0 | next_step$point > 1))
      break

    segment_length <- sqrt(sum((next_step$point - path[nrow(path), c("x", "y")])^2))
    if (!is.finite(segment_length) || segment_length <= .Machine$double.eps ||
        travelled + segment_length > length_limit || is_occupied(next_step$point))
      break

    next_key <- cell_key(next_step$point)
    if (!identical(next_key, visited[length(visited)])) {
      if (length(visited) >= 4L &&
          next_key %in% visited[seq_len(length(visited) - 3L)])
        break
      visited <- c(visited, next_key)
    }

    path <- rbind(path, c(next_step$point, next_step$speed))
    travelled <- travelled + segment_length
  }

  storage.mode(path) <- "double"
  path
}


.empty_streamlines <- function() {
  data.frame(
    line_id = integer(),
    point_id = integer(),
    x = numeric(),
    y = numeric(),
    speed = numeric(),
    arrow = logical(),
    stringsAsFactors = FALSE
  )
}


.is_valid_stream_grid_cache <- function(grid) {
  required_columns <- c(
    "x", "y", "dx", "dy", "mass", "local_length", "valid"
  )
  numeric_columns <- setdiff(required_columns, "valid")
  if (!is.data.frame(grid) || nrow(grid) == 0L ||
      anyDuplicated(names(grid)) ||
      !all(required_columns %in% names(grid)))
    return(FALSE)
  if (!all(vapply(
    grid[numeric_columns],
    function(column) is.numeric(column) && !is.logical(column),
    logical(1L)
  )))
    return(FALSE)
  if (!is.logical(grid$valid) || anyNA(grid$valid) ||
      any(!is.finite(grid$x)) || any(!is.finite(grid$y)))
    return(FALSE)
  if (!all(is.finite(grid$dx[grid$valid])) ||
      !all(is.finite(grid$dy[grid$valid])))
    return(FALSE)

  field <- tryCatch(
    .make_stream_field(grid),
    error = function(error) NULL
  )
  if (is.null(field) || !is.list(field) ||
      !all(c("dx", "dy", "valid", "nx", "ny") %in% names(field)))
    return(FALSE)
  expected_dimensions <- c(field$nx, field$ny)
  is.matrix(field$dx) && is.matrix(field$dy) && is.matrix(field$valid) &&
    identical(dim(field$dx), expected_dimensions) &&
    identical(dim(field$dy), expected_dimensions) &&
    identical(dim(field$valid), expected_dimensions)
}


.is_valid_streamline_cache <- function(streamlines) {
  required_columns <- c("line_id", "point_id", "x", "y", "speed", "arrow")
  numeric_columns <- setdiff(required_columns, "arrow")
  if (!is.data.frame(streamlines) || anyDuplicated(names(streamlines)) ||
      !all(required_columns %in% names(streamlines)))
    return(FALSE)
  if (!all(vapply(
    streamlines[numeric_columns],
    function(column) is.numeric(column) && !is.logical(column),
    logical(1L)
  )))
    return(FALSE)
  if (!is.logical(streamlines$arrow) || anyNA(streamlines$arrow))
    return(FALSE)
  if (nrow(streamlines) == 0L)
    return(TRUE)
  if (any(!is.finite(streamlines$line_id)) ||
      any(!is.finite(streamlines$point_id)) ||
      any(!is.finite(streamlines$x)) ||
      any(!is.finite(streamlines$y)) ||
      any(!is.finite(streamlines$speed)) ||
      any(streamlines$line_id <= 0) ||
      any(streamlines$line_id != floor(streamlines$line_id)) ||
      any(streamlines$point_id != floor(streamlines$point_id)) ||
      any(streamlines$speed < 0))
    return(FALSE)

  line_runs <- rle(streamlines$line_id)
  if (any(line_runs$lengths < 4L) || anyDuplicated(line_runs$values) ||
      !all(line_runs$values == seq_along(line_runs$values)))
    return(FALSE)
  run_ends <- cumsum(line_runs$lengths)
  run_starts <- c(1L, head(run_ends, -1L) + 1L)
  all(vapply(
    seq_along(run_starts),
    function(run_id) {
      indices <- run_starts[run_id]:run_ends[run_id]
      point_ids <- streamlines$point_id[indices]
      arrows <- streamlines$arrow[indices]
      x_steps <- diff(streamlines$x[indices])
      y_steps <- diff(streamlines$y[indices])
      all(point_ids == seq_along(indices)) &&
        sum(arrows) == 1L && !arrows[length(arrows)] &&
        all(x_steps != 0 | y_steps != 0)
    },
    logical(1L)
  ))
}


.stream_grid_params <- function(reduction, vkey, controls) {
  list(
    reduction = reduction,
    vkey = vkey,
    n_grid = controls$n_grid,
    n_neighbors = controls$n_neighbors,
    smooth = controls$smooth,
    min_mass = controls$min_mass,
    cutoff_perc = controls$cutoff_perc,
    x_pad = controls$x_pad,
    y_pad = controls$y_pad
  )
}


.stream_path_params <- function(controls) {
  list(
    density = controls$density,
    max_length = controls$max_length,
    integration_direction = controls$integration_direction
  )
}


.find_velocity_embedding <- function(pv, reduction, vkey) {
  reduction_suffix <- tolower(sub("^[Xx]_", "", reduction))
  key_prefix <- paste0(vkey, "_")
  reduction_names <- names(pv@reductions)
  prefix_match <- startsWith(tolower(reduction_names), tolower(key_prefix))
  candidate_names <- reduction_names[prefix_match]
  if (length(candidate_names) == 0L)
    return(NULL)

  candidate_suffixes <- tolower(sub("^[Xx]_", "", substr(
    candidate_names,
    nchar(key_prefix) + 1L,
    nchar(candidate_names)
  )))
  candidate_names <- candidate_names[candidate_suffixes == reduction_suffix]
  if (length(candidate_names) == 0L)
    return(NULL)

  requested_name <- paste0(vkey, "_", reduction)
  if (requested_name %in% candidate_names)
    return(pv@reductions[[requested_name]])
  if (length(candidate_names) == 1L)
    return(pv@reductions[[candidate_names]])

  stop(sprintf(
    "Velocity embedding for reduction '%s' is ambiguous: %s.",
    reduction,
    paste(candidate_names, collapse = ", ")
  ))
}


.trace_velocity_streams <- function(grid,
                                    density,
                                    max_length,
                                    integration_direction) {
  scalar_finite <- function(x) {
    is.numeric(x) && !is.logical(x) && length(x) == 1L && is.finite(x)
  }
  if (!scalar_finite(density) || density <= 0)
    stop("`density` must be a positive finite scalar.")
  if (!scalar_finite(max_length) || max_length <= 0)
    stop("`max_length` must be a positive finite scalar.")
  integration_direction <- match.arg(
    integration_direction,
    choices = c("both", "forward", "backward")
  )

  if (is.data.frame(grid) && "valid" %in% names(grid) &&
      !any(grid$valid %in% TRUE))
    return(.empty_streamlines())

  field <- .make_stream_field(grid)
  if (!any(field$valid %in% TRUE))
    return(.empty_streamlines())

  raw_occupancy_n <- ceiling(30 * density)
  if (!is.finite(raw_occupancy_n) || raw_occupancy_n < 1 ||
      raw_occupancy_n > .Machine$integer.max)
    stop("`density` produces an occupancy grid that is too large to represent.")
  max_occupancy_elements <- 25000000
  max_occupancy_n <- floor(sqrt(max_occupancy_elements))
  if (raw_occupancy_n > max_occupancy_n)
    stop("`density` exceeds the occupancy grid element budget; lower `density`.")
  occupancy_n <- as.integer(raw_occupancy_n)
  occupied <- matrix(FALSE, nrow = occupancy_n, ncol = occupancy_n)

  valid_grid <- grid[grid$valid %in% TRUE, , drop = FALSE]
  finite_mass <- is.finite(valid_grid$mass)
  candidates <- valid_grid[order(
    !finite_mass,
    -ifelse(finite_mass, valid_grid$mass, 0),
    valid_grid$y,
    valid_grid$x
  ), , drop = FALSE]

  occupancy_index <- function(point) {
    bounded <- pmin(pmax(as.numeric(point), 0), 1)
    pmin(
      pmax(as.integer(floor(bounded * occupancy_n)) + 1L, 1L),
      occupancy_n
    )
  }
  seed_path <- function(seed, sampled) {
    path <- matrix(
      c(seed[1L], seed[2L], unname(sampled["speed"])),
      ncol = 3L,
      dimnames = list(NULL, c("x", "y", "speed"))
    )
    storage.mode(path) <- "double"
    path
  }
  normalise_path <- function(path, fallback) {
    if (is.null(path))
      return(fallback)
    path <- as.matrix(path[, c("x", "y", "speed"), drop = FALSE])
    storage.mode(path) <- "double"
    colnames(path) <- c("x", "y", "speed")
    path
  }

  lines <- list()
  line_id <- 0L
  for (candidate_id in seq_len(nrow(candidates))) {
    candidate <- candidates[candidate_id, , drop = FALSE]
    seed <- c(
      (candidate$x - field$x_min) / field$x_span,
      (candidate$y - field$y_min) / field$y_span
    )
    seed <- pmin(pmax(as.numeric(seed), 0), 1)
    seed_cell <- occupancy_index(seed)
    if (isTRUE(occupied[seed_cell[1L], seed_cell[2L]]))
      next

    sampled <- .sample_stream_field(field, seed)
    if (is.null(sampled))
      next
    initial_path <- seed_path(seed, sampled)

    if (identical(integration_direction, "both")) {
      backward <- normalise_path(
        .trace_stream_direction(field, seed, -1, max_length / 2, occupied),
        initial_path
      )
      forward <- normalise_path(
        .trace_stream_direction(field, seed, 1, max_length / 2, occupied),
        initial_path
      )
      backward <- backward[seq.int(nrow(backward), 1L), , drop = FALSE]
      path <- rbind(
        backward[seq_len(max(0L, nrow(backward) - 1L)), , drop = FALSE],
        forward
      )
    } else if (identical(integration_direction, "backward")) {
      path <- normalise_path(
        .trace_stream_direction(field, seed, -1, max_length, occupied),
        initial_path
      )
      path <- path[seq.int(nrow(path), 1L), , drop = FALSE]
    } else {
      path <- normalise_path(
        .trace_stream_direction(field, seed, 1, max_length, occupied),
        initial_path
      )
    }
    storage.mode(path) <- "double"

    if (nrow(path) < 4L)
      next
    segment_lengths <- sqrt(rowSums((path[-1L, c("x", "y"), drop = FALSE] -
                                       path[-nrow(path), c("x", "y"), drop = FALSE])^2))
    path_length <- sum(segment_lengths)
    if (!is.finite(path_length) || path_length < sqrt(2) / occupancy_n)
      next

    arrow_point <- which(cumsum(segment_lengths) >= path_length / 2)[1L]
    arrow_point <- min(max(as.integer(arrow_point), 1L), nrow(path) - 1L)
    arrow <- rep(FALSE, nrow(path))
    arrow[arrow_point] <- TRUE

    line_id <- line_id + 1L
    lines[[line_id]] <- data.frame(
      line_id = rep.int(line_id, nrow(path)),
      point_id = seq_len(nrow(path)),
      x = field$x_min + path[, "x"] * field$x_span,
      y = field$y_min + path[, "y"] * field$y_span,
      speed = path[, "speed"],
      arrow = arrow,
      stringsAsFactors = FALSE
    )

    occupied_cells <- t(vapply(
      seq_len(nrow(path)),
      function(point_id) occupancy_index(path[point_id, c("x", "y")]),
      integer(2L)
    ))
    occupied[cbind(occupied_cells[, 1L], occupied_cells[, 2L])] <- TRUE
  }

  if (length(lines) == 0L)
    return(.empty_streamlines())
  streamlines <- do.call(rbind, lines)
  rownames(streamlines) <- NULL
  streamlines
}


#' @title Compute Velocity Stream Field
#' @description
#' Interpolate cell-level velocity embeddings onto a regular two-dimensional
#' grid using Gaussian neighbourhood weights, then trace deterministic stream
#' paths through the interpolated field in normalized axis coordinates.
#'
#' Grid validity is determined by two masks: interpolated velocity strength and
#' local velocity length. Invalid regions receive \code{NA} velocity components,
#' which prevents paths from crossing unsupported parts of the embedding.
#'
#' @param pv A \code{plantvelo} object.
#' @param reduction Character scalar. Embedding key in \code{@reductions}.
#'   Default \code{"umap"}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#' @param n_grid Integer. Side length of the interpolation grid. The total
#'   number of grid points is \eqn{n\_grid^2}. Default \code{50L}.
#' @param n_neighbors Integer or \code{NULL}. Number of nearest cells used for
#'   each grid point. If \code{NULL}, an automatic value is used. Default
#'   \code{NULL}.
#' @param smooth Numeric. Gaussian bandwidth as a multiple of the mean grid
#'   spacing. Default \code{0.5}.
#' @param min_mass Numeric. Controls the lower bound used by the interpolated
#'   velocity-strength mask. Default \code{1}.
#' @param cutoff_perc Numeric in \eqn{[0,100]}. Percentile threshold applied to
#'   local velocity length. Default \code{5}.
#' @param x_pad,y_pad Numeric. Fractional padding added to the embedding ranges
#'   before constructing the interpolation grid. Default \code{0.01}.
#' @param density Positive numeric. Controls occupancy-grid resolution for
#'   seed spacing; larger values permit denser paths. Default \code{2}.
#' @param max_length Positive numeric. Maximum path length measured in
#'   normalized axis coordinates. Default \code{4}.
#' @param integration_direction Character. Integrate in both directions,
#'   forward, or backward. One of \code{"both"}, \code{"forward"}, or
#'   \code{"backward"}. Default \code{"both"}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @section Output written to \code{@velocity}:
#' \describe{
#'   \item{\code{stream}}{A \code{data.frame} with columns \code{x},
#'     \code{y}, \code{dx}, \code{dy}, \code{mass}, \code{local_length}, and
#'     \code{valid}. Each row represents one regular-grid point; \code{dx} and
#'     \code{dy} are \code{NA} where \code{valid} is false.}
#'   \item{\code{streamlines}}{A \code{data.frame} with columns \code{line_id},
#'     \code{point_id}, \code{x}, \code{y}, \code{speed}, and \code{arrow} for
#'     the traced paths and their midpoint arrow positions.}
#'   \item{\code{stream_params}}{A list containing separate \code{grid} and
#'     \code{path} parameter lists used to validate cached results.}
#' }
#'
#' @return The input \code{plantvelo} object with \code{@velocity$stream},
#'   \code{@velocity$streamlines}, and \code{@velocity$stream_params} updated.
#'
#' @examples
#' \dontrun{
#' pv <- compute_velocity_stream(pv)
#' field <- get_velocity_stream_field(pv)
#' }
#' @export
compute_velocity_stream <- function(pv,
                                    reduction = "umap",
                                    vkey = "velocity",
                                    n_grid = 50L,
                                    n_neighbors = NULL,
                                    smooth = 0.5,
                                    min_mass = 1,
                                    cutoff_perc = 5,
                                    x_pad = 0.01,
                                    y_pad = 0.01,
                                    density = 2,
                                    max_length = 4,
                                    integration_direction = c("both", "forward", "backward"),
                                    verbose = TRUE) {

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  controls <- .validate_stream_controls(
    n_grid = n_grid,
    n_neighbors = n_neighbors,
    smooth = smooth,
    min_mass = min_mass,
    cutoff_perc = cutoff_perc,
    x_pad = x_pad,
    y_pad = y_pad,
    density = density,
    max_length = max_length,
    integration_direction = integration_direction
  )

  embedding <- .get_embedding(pv, reduction)
  velocity_embedding <- .find_velocity_embedding(pv, reduction, vkey)
  if (is.null(velocity_embedding)) {
    pv <- compute_velocity_embedding(
      pv,
      reduction = reduction,
      vkey = vkey,
      verbose = FALSE
    )
    velocity_embedding <- .find_velocity_embedding(pv, reduction, vkey)
  }
  if (is.null(velocity_embedding))
    stop("Could not obtain velocity embedding. Run compute_velocity_embedding() first.")

  rowname_state <- function(x) {
    identifiers <- rownames(x)
    if (is.null(identifiers))
      return(list(present = FALSE, valid = FALSE, values = NULL))
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
    stop("Embedding and velocity embedding must either both have row names or both have none.")
  } else if (embedding_rows$present &&
             (!embedding_rows$valid || !velocity_rows$valid)) {
    stop("Embedding and velocity embedding row names must be complete and unique.")
  } else if (embedding_rows$present) {
    if (!setequal(embedding_rows$values, velocity_rows$values))
      stop("Embedding and velocity embedding row names identify different cells.")
    velocity_embedding <- velocity_embedding[embedding_rows$values, , drop = FALSE]
  } else if (nrow(embedding) != nrow(velocity_embedding)) {
    stop("Embedding and velocity embedding must have the same number of rows.")
  }

  coords <- as.matrix(embedding[, 1:2, drop = FALSE])
  velocity <- as.matrix(velocity_embedding[, 1:2, drop = FALSE])
  grid <- .estimate_stream_grid(
    coords = coords,
    velocity = velocity,
    n_grid = controls$n_grid,
    n_neighbors = controls$n_neighbors,
    smooth = controls$smooth,
    min_mass = controls$min_mass,
    cutoff_perc = controls$cutoff_perc,
    x_pad = controls$x_pad,
    y_pad = controls$y_pad
  )
  streamlines <- .trace_velocity_streams(
    grid = grid,
    density = controls$density,
    max_length = controls$max_length,
    integration_direction = controls$integration_direction
  )

  pv@velocity[["stream"]] <- grid
  pv@velocity[["streamlines"]] <- streamlines
  pv@velocity[["stream_params"]] <- list(
    grid = .stream_grid_params(reduction, vkey, controls),
    path = .stream_path_params(controls)
  )

  if (isTRUE(verbose)) {
    message(sprintf(
      "Stream field: %d x %d grid, %d valid points, %d streamlines.",
      controls$n_grid,
      controls$n_grid,
      sum(grid$valid %in% TRUE),
      length(unique(streamlines$line_id))
    ))
  }

  return(pv)

}


#' @title Get Velocity Stream Field
#' @description
#' Retrieve only the precomputed regular-grid velocity field stored in
#' \code{@velocity$stream}. This function does not return the cached paths or
#' parameter lists.
#'
#' @param pv A \code{plantvelo} object.
#' @param reduction Character scalar. Embedding key used when computing the
#'   stream field. Currently retained for interface consistency. Default
#'   \code{"umap"}.
#' @param vkey Character scalar. Velocity layer key. Currently retained for
#'   interface consistency. Default \code{"velocity"}.
#'
#' @return A regular-grid \code{data.frame} with columns \code{x} and \code{y}
#'   for coordinates; \code{dx} and \code{dy} for interpolated velocity;
#'   \code{mass} for cumulative Gaussian weight; \code{local_length} for the
#'   local velocity-length statistic; and \code{valid} for the combined mask.
#'   Returns \code{NULL} with a message if no cached field is available.
#' @examples
#' \dontrun{
#' field <- get_velocity_stream_field(pv)
#' }
#' @export
get_velocity_stream_field <- function(pv,
                                      reduction = "umap",
                                      vkey      = "velocity") {
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  key   <- "stream"
  field <- pv@velocity[[key]]

  if (is.null(field)) {
    message("No stream field found. Run compute_velocity_stream() first.")
    return(NULL)
  }

  field
}


.scale_stream_linewidths <- function(speed, min_linewidth, max_linewidth) {
  midpoint <- (min_linewidth + max_linewidth) / 2
  linewidth <- rep(midpoint, length(speed))
  finite_speed <- speed[is.finite(speed)]

  if (length(finite_speed) == 0L) {
    return(linewidth)
  }

  limits <- stats::quantile(
    finite_speed,
    probs = c(0.05, 0.95),
    names = FALSE
  )
  if (!all(is.finite(limits)) || limits[2L] <= limits[1L]) {
    return(linewidth)
  }

  finite <- is.finite(speed)
  scaled <- (speed[finite] - limits[1L]) / (limits[2L] - limits[1L])
  scaled <- pmin(1, pmax(0, scaled))
  linewidth[finite] <- min_linewidth +
    sqrt(scaled) * (max_linewidth - min_linewidth)
  linewidth
}


.stream_arrow_segments <- function(streamlines, arrow_fraction) {
  starts <- which(streamlines$arrow)
  starts <- starts[
    starts < nrow(streamlines) &
      streamlines$line_id[starts + 1L] == streamlines$line_id[starts]
  ]
  if (length(starts) == 0L || arrow_fraction <= 0) {
    return(data.frame(
      x = numeric(),
      y = numeric(),
      xend = numeric(),
      yend = numeric()
    ))
  }

  n_starts <- length(starts)
  target <- min(
    n_starts,
    max(1L, floor(n_starts * arrow_fraction + 0.5))
  )
  selected <- unique(round(seq.int(1L, n_starts, length.out = target)))
  starts <- starts[selected]

  data.frame(
    x = streamlines$x[starts],
    y = streamlines$y[starts],
    xend = streamlines$x[starts + 1L],
    yend = streamlines$y[starts + 1L],
    stringsAsFactors = FALSE
  )
}


#' @title Plot Velocity Stream
#' @description
#' Plot RNA velocity as continuous stream lines in a two-dimensional embedding.
#'
#' The line width is scaled by local velocity speed. Arrowheads are assigned to
#' a fraction of paths by selecting evenly spaced indices in deterministic path
#' order, with 70 percent of paths selected by default. Invalid field regions stop
#' integration so that disconnected embedding regions remain visually separated.
#'
#' Cell points are plotted underneath the stream lines and can be coloured by a
#' metadata column in \code{@meta.data}.
#'
#' @param pv A \code{plantvelo} object.
#' @param reduction Character scalar. Embedding key in \code{@reductions}.
#'   Default \code{"umap"}.
#' @param vkey Character scalar. Velocity layer key. Default
#'   \code{"velocity"}.
#' @param group_by Character scalar or \code{NULL}. Column name in
#'   \code{@meta.data} used for colouring cells. Default \code{NULL}.
#' @param n_grid Integer. Side length of the interpolation grid. Default
#'   \code{50L}.
#' @param n_neighbors Integer or \code{NULL}. Number of nearest cells used for
#'   each grid point. If \code{NULL}, an automatic value is used. Default
#'   \code{NULL}.
#' @param smooth Numeric. Gaussian bandwidth as a multiple of the mean grid
#'   spacing. Default \code{0.5}.
#' @param min_mass Numeric. Controls the lower bound used by the interpolated
#'   velocity-strength mask. Default \code{1}.
#' @param cutoff_perc Numeric in \eqn{[0,100]}. Percentile threshold applied to
#'   local velocity length. Default \code{5}.
#' @param x_pad,y_pad Numeric. Fractional padding added to the embedding ranges
#'   before constructing the interpolation grid. Default \code{0.01}.
#' @param density Positive numeric. Controls occupancy-grid resolution for
#'   seed spacing; larger values permit denser paths. Default \code{2}.
#' @param max_length Positive numeric. Maximum path length measured in
#'   normalized axis coordinates. Default \code{4}.
#' @param integration_direction Character. Integrate in both directions,
#'   forward, or backward. One of \code{"both"}, \code{"forward"}, or
#'   \code{"backward"}. Default \code{"both"}.
#' @param stream_color Character scalar. Stream line colour. Default
#'   \code{"grey0"}.
#' @param stream_alpha Numeric. Stream line opacity. Default \code{0.8}.
#' @param stream_min_linewidth Non-negative numeric. Minimum line width used by
#'   the local-speed mapping. Together with \code{stream_linewidth}, must satisfy
#'   \code{stream_min_linewidth <= stream_linewidth}. Default \code{0.10}.
#' @param stream_linewidth Non-negative numeric. Maximum line width used by the
#'   local-speed mapping. Default \code{0.35}.
#' @param arrow_fraction Numeric in \eqn{[0,1]}. Fraction of paths that receive
#'   an arrow. Selection uses evenly spaced indices in deterministic path order;
#'   \code{0} selects none and \code{1} selects all paths. Default \code{0.7}.
#' @param arrow_size Non-negative numeric. Arrow-head size in cm. Default
#'   \code{0.10}.
#' @param arrow_linewidth Non-negative numeric. Fixed arrow line width,
#'   independent of local velocity speed. Default \code{0.15}.
#' @param arrow_angle Numeric. Arrow-head angle in degrees. Must satisfy
#'   \code{0 < arrow_angle < 90}. Default \code{20}.
#' @param arrow_type Character. Arrow-head type, either \code{"open"} or
#'   \code{"closed"}. Default \code{"closed"}.
#' @param point_size Numeric. Cell point size. Default \code{0.30}.
#' @param point_shape Numeric. Cell point shape passed to the point layer.
#'   Default \code{16}.
#' @param alpha Numeric. Cell point opacity. Default \code{0.55}.
#' @param input_color Character vector or \code{NULL}. Custom colour palette
#'   for cell colouring. Uses the same semantics as
#'   \code{plot_velocity_embedding()}. Default \code{NULL}.
#' @param raster_points Logical. Rasterize cell points with \pkg{ggrastr} while
#'   keeping stream paths as vectors. If \code{TRUE}, \pkg{ggrastr} is required.
#'   Default \code{TRUE}.
#' @param raster_arrows Logical. Rasterize arrow segments with \pkg{ggrastr}
#'   while keeping stream paths as vectors. If \code{TRUE}, \pkg{ggrastr} is
#'   required. Default \code{FALSE}.
#' @param raster_dpi Positive numeric. Resolution in dots per inch used for
#'   rasterized points and arrows. Default \code{1200}.
#' @param title Character scalar or \code{NULL}. Plot title. Default
#'   \code{NULL}.
#' @param legend Logical. Show the legend. Default \code{TRUE}.
#' @param recompute Logical. If \code{FALSE}, a grid with matching grid
#'   parameters is reused. Matching path parameters also reuse cached paths;
#'   otherwise only paths are retraced. A grid-parameter mismatch recomputes
#'   both levels. If \code{TRUE}, both levels are recomputed. Default
#'   \code{TRUE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @details
#' Finite local speeds are clipped to their 5th and 95th percentiles, then
#' square-root transformed into the bounds set by
#' \code{stream_min_linewidth} and \code{stream_linewidth}. If the finite speed
#' distribution is absent or degenerate, all widths use the midpoint of those
#' bounds; non-finite speeds also use the midpoint. Arrow paths are selected by
#' evenly spaced indices in deterministic path order.
#'
#' If either \code{raster_points} or \code{raster_arrows} is \code{TRUE},
#' \pkg{ggrastr} is required. Because \code{raster_points} defaults to
#' \code{TRUE}, default calls require \pkg{ggrastr}. Setting both options to
#' \code{FALSE} removes this requirement and uses vector layers for points and
#' arrows; stream paths remain vector layers in all cases.
#'
#' @section Cache reuse:
#' With \code{recompute = FALSE}, matching grid parameters reuse the cached
#' regular grid. If path parameters also match, cached paths are reused;
#' otherwise only the paths are retraced. Different grid parameters cause both
#' the grid and paths to be recomputed. With \code{recompute = TRUE}, both cache
#' levels are recomputed unconditionally.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plot_velocity_stream(pv, group_by = "cell_type")
#' plot_velocity_stream(pv, group_by = "cell_type", density = 2.5, max_length = 3)
#'
#' p <- plot_velocity_stream(
#'   pv,
#'   group_by = "cell_type",
#'   arrow_fraction = 0.5,
#'   raster_points = TRUE
#' )
#' ggplot2::ggsave(
#'   "cell_type.velocity_stream.png",
#'   plot = p,
#'   width = 110,
#'   height = 85,
#'   units = "mm",
#'   dpi = 600
#' )
#' }
#' @export
plot_velocity_stream <- function(pv,
                                 reduction = "umap",
                                 vkey = "velocity",
                                 group_by = NULL,
                                 n_grid = 50L,
                                 n_neighbors = NULL,
                                 smooth = 0.5,
                                 min_mass = 1,
                                 cutoff_perc = 5,
                                 x_pad = 0.01,
                                 y_pad = 0.01,
                                 density = 2,
                                 max_length = 4,
                                 integration_direction = c("both", "forward", "backward"),
                                 stream_color = "grey0",
                                 stream_alpha = 0.8,
                                 stream_min_linewidth = 0.10,
                                 stream_linewidth = 0.35,
                                 arrow_fraction = 0.7,
                                 arrow_size = 0.10,
                                 arrow_linewidth = 0.15,
                                 arrow_angle = 20,
                                 arrow_type = c("closed", "open"),
                                 point_size = 0.30,
                                 point_shape = 16,
                                 alpha = 0.55,
                                 input_color = NULL,
                                 raster_points = TRUE,
                                 raster_arrows = FALSE,
                                 raster_dpi = 1200,
                                 title = NULL,
                                 legend = TRUE,
                                 recompute = TRUE,
                                 verbose = TRUE) {

  arrow_type <- match.arg(arrow_type)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  controls <- .validate_stream_controls(
    n_grid = n_grid,
    n_neighbors = n_neighbors,
    smooth = smooth,
    min_mass = min_mass,
    cutoff_perc = cutoff_perc,
    x_pad = x_pad,
    y_pad = y_pad,
    density = density,
    max_length = max_length,
    integration_direction = integration_direction
  )
  scalar_finite <- function(x) {
    is.numeric(x) && !is.logical(x) && length(x) == 1L && is.finite(x)
  }
  scalar_logical <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
  }
  if (!scalar_finite(stream_alpha) || stream_alpha < 0 || stream_alpha > 1)
    stop("`stream_alpha` must be a finite scalar between 0 and 1.")
  if (!scalar_finite(alpha) || alpha < 0 || alpha > 1)
    stop("`alpha` must be a finite scalar between 0 and 1.")
  if (!scalar_finite(stream_min_linewidth) || stream_min_linewidth < 0)
    stop("`stream_min_linewidth` must be a non-negative finite scalar.")
  if (!scalar_finite(stream_linewidth) || stream_linewidth < 0)
    stop("`stream_linewidth` must be a non-negative finite scalar.")
  if (stream_min_linewidth > stream_linewidth)
    stop("`stream_min_linewidth` must not exceed `stream_linewidth`.")
  if (!scalar_finite(arrow_fraction) || arrow_fraction < 0 || arrow_fraction > 1)
    stop("`arrow_fraction` must be a finite scalar between 0 and 1.")
  if (!scalar_finite(arrow_size) || arrow_size < 0)
    stop("`arrow_size` must be a non-negative finite scalar.")
  if (!scalar_finite(arrow_linewidth) || arrow_linewidth < 0)
    stop("`arrow_linewidth` must be a non-negative finite scalar.")
  if (!scalar_finite(arrow_angle) || arrow_angle <= 0 || arrow_angle >= 90)
    stop("`arrow_angle` must be a finite scalar strictly between 0 and 90.")
  if (!scalar_finite(point_size) || point_size < 0)
    stop("`point_size` must be a non-negative finite scalar.")
  if (!scalar_finite(point_shape))
    stop("`point_shape` must be a finite numeric scalar.")
  if (!scalar_logical(raster_points))
    stop("`raster_points` must be a single non-missing logical value.")
  if (!scalar_logical(raster_arrows))
    stop("`raster_arrows` must be a single non-missing logical value.")
  if (!scalar_finite(raster_dpi) || raster_dpi <= 0)
    stop("`raster_dpi` must be a positive finite scalar.")
  if ((raster_points || raster_arrows) &&
      !requireNamespace("ggrastr", quietly = TRUE)) {
    stop(
      "Rasterization requires package 'ggrastr'. Install it with: install.packages('ggrastr')",
      call. = FALSE
    )
  }
  if (!scalar_logical(recompute))
    stop("`recompute` must be a single non-missing logical value.")
  if (!scalar_logical(legend))
    stop("`legend` must be a single non-missing logical value.")
  if (!scalar_logical(verbose))
    stop("`verbose` must be a single non-missing logical value.")

  E <- .get_embedding(pv, reduction)
  dims <- paste0(toupper(sub("^[Xx]_", "", reduction)), "_", 1:2)

  cell_df <- data.frame(
    x = as.numeric(E[, 1]),
    y = as.numeric(E[, 2]),
    stringsAsFactors = FALSE
  )

  color_vals  <- .resolve_color_by(pv, group_by)
  if (!is.null(color_vals) && length(color_vals) != nrow(cell_df)) {
    stop(
      "Resolved cell colours must have one value per embedding row.",
      call. = FALSE
    )
  }
  is_discrete <- !is.null(color_vals) && !is.numeric(color_vals)
  if (!is.null(color_vals)) cell_df$color <- color_vals

  requested_grid_params <- .stream_grid_params(reduction, vkey, controls)
  requested_path_params <- .stream_path_params(controls)
  cached_grid <- pv@velocity[["stream"]]
  cached_paths <- pv@velocity[["streamlines"]]
  cached_params <- pv@velocity[["stream_params"]]
  reuse_grid <- !recompute && .is_valid_stream_grid_cache(cached_grid) &&
    is.list(cached_params) &&
    identical(cached_params$grid, requested_grid_params)

  if (!reuse_grid) {
    pv <- compute_velocity_stream(
      pv = pv,
      reduction = reduction,
      vkey = vkey,
      n_grid = controls$n_grid,
      n_neighbors = controls$n_neighbors,
      smooth = controls$smooth,
      min_mass = controls$min_mass,
      cutoff_perc = controls$cutoff_perc,
      x_pad = controls$x_pad,
      y_pad = controls$y_pad,
      density = controls$density,
      max_length = controls$max_length,
      integration_direction = controls$integration_direction,
      verbose = verbose
    )
    stream_df <- pv@velocity[["streamlines"]]
  } else {
    grid_df <- cached_grid
    if (.is_valid_streamline_cache(cached_paths) &&
        identical(cached_params$path, requested_path_params)) {
      stream_df <- cached_paths
      if (verbose)
        message("Using cached stream field and paths.")
    } else {
      stream_df <- .trace_velocity_streams(
        grid = grid_df,
        density = controls$density,
        max_length = controls$max_length,
        integration_direction = controls$integration_direction
      )
      if (verbose)
        message("Using cached stream field and recomputed paths.")
    }
  }

  p <- ggplot2::ggplot()
  point_layer <- function(mapping, colour = NULL) {
    arguments <- list(
      data = cell_df,
      mapping = mapping,
      size = point_size,
      alpha = alpha,
      shape = point_shape,
      stroke = 0,
      inherit.aes = FALSE
    )
    if (!is.null(colour))
      arguments$colour <- colour

    if (raster_points) {
      arguments$raster.dpi <- raster_dpi
      do.call(ggrastr::geom_point_rast, arguments)
    } else {
      do.call(ggplot2::geom_point, arguments)
    }
  }

  if (is.null(color_vals)) {
    p <- p + point_layer(
      mapping = ggplot2::aes(x = .data$x, y = .data$y),
      colour = "#BFBFBF"
    )
  } else if (is_discrete) {
    p <- p + point_layer(
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                             colour = .data$color)
    )
    if (!is.null(input_color)) {
      p <- p + ggplot2::scale_colour_manual(
        values = input_color,
        name = group_by,
        guide  = ggplot2::guide_legend(
          override.aes = list(
            size  = 3.2,
            alpha = 1)
        )
      )
    } else {
      p <- p + ggplot2::scale_colour_discrete(
        name = group_by,
        guide = ggplot2::guide_legend(
          override.aes = list(
            size  = 3.2,
            alpha = 1
          )
          )
        )
    }
  } else {
    cols <- input_color %||% c("#440154", "#31688E", "#35B779", "#FDE725")
    p <- p + point_layer(
      mapping = ggplot2::aes(x = .data$x, y = .data$y,
                             colour = .data$color)
    ) +
      ggplot2::scale_colour_gradientn(colours = cols, name = group_by)
  }

  if (nrow(stream_df) == 0L) {
    warning(
      "No streamlines passed the field and length filters.",
      call. = FALSE
    )
  } else {
    stream_df$linewidth <- .scale_stream_linewidths(
      stream_df$speed,
      stream_min_linewidth,
      stream_linewidth
    )
    arrow_df <- .stream_arrow_segments(stream_df, arrow_fraction)

    p <- p +
      ggplot2::geom_path(
        data = stream_df,
        mapping = ggplot2::aes(
          x = .data$x,
          y = .data$y,
          group = .data$line_id,
          linewidth = .data$linewidth
        ),
        colour = stream_color,
        alpha = stream_alpha,
        lineend = "round",
        inherit.aes = FALSE
      ) +
      ggplot2::scale_linewidth_identity()

    if (nrow(arrow_df) > 0L) {
      arrow_layer <- ggplot2::geom_segment(
        data = arrow_df,
        mapping = ggplot2::aes(
          x = .data$x,
          y = .data$y,
          xend = .data$xend,
          yend = .data$yend
        ),
        colour = stream_color,
        alpha = stream_alpha,
        linewidth = arrow_linewidth,
        arrow = grid::arrow(
          length = grid::unit(arrow_size, "cm"),
          angle = arrow_angle,
          type = arrow_type
        ),
        inherit.aes = FALSE
      )
      if (raster_arrows)
        arrow_layer <- ggrastr::rasterise(arrow_layer, dpi = raster_dpi)
      p <- p + arrow_layer
    }
  }

  p <- p +
    ggplot2::labs(title = title, x = dims[1], y = dims[2]) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.border    = ggplot2::element_rect(fill = NA, colour = "black",
                                              linewidth = 0.5),
      axis.line       = ggplot2::element_blank(),
      plot.title      = ggplot2::element_text(size = 13, colour = "black",
                                              hjust = 0.5),
      axis.text       = ggplot2::element_text(size = 10, colour = "black"),
      axis.title      = ggplot2::element_text(size = 12, colour = "black"),
      legend.title    = ggplot2::element_text(size = 10, colour = "black"),
      legend.text     = ggplot2::element_text(size = 9,  colour = "black"),
      legend.key.size = ggplot2::unit(0.4, "cm")
    )

  if (!legend)
    p <- p + ggplot2::theme(legend.position = "none")

  p
}
