# =============================================================================
# PlantVelocity: compute_velocity.R
# RNA velocity computation (S4 plantvelo API)
#   vS = beta * U(t) - gamma * S(t)
#   vU = (alpha * o - beta * U(t)) * scaling
# =============================================================================

#' @title Compute RNA Velocity
#' @description
#' Compute per-cell, per-gene RNA velocity vectors from kinetic parameters
#' recovered by \code{recover_dynamics()}.
#'
#' Two modes are supported:
#' \describe{
#'   \item{\code{"dynamical"}}{Two-state U/S velocity. Requires
#'     \code{recover_dynamics()}.}
#'   \item{\code{"stochastic"}}{Fallback residual velocity estimated by
#'     per-gene regression. No prior kinetic fitting is required.}
#' }
#'
#' If \code{mode = NULL}, the velocity mode is automatically inferred from
#' \code{@kinetics$params}. A valid fitted-dynamics schema containing only
#' \code{"unfitted"} genes produces correctly shaped all-\code{NA} dynamical
#' velocity layers rather than falling back to stochastic velocity.
#' Dynamical mode accepts only the \code{"ir_excluded_2state_v1"} schema
#' written by \code{recover_dynamics()}; the per-gene \code{fit_model} field
#' identifies fitted and unfitted genes. Legacy fits must be refitted before
#' velocity computation.
#'
#' @section Output written to \code{@layers}:
#' \describe{
#'   \item{\code{velocity}}{genes \eqn{\times} cells. Primary velocity
#'     \eqn{dS/dt}.}
#'   \item{\code{velocity_u}}{genes \eqn{\times} cells. Unspliced velocity
#'     \eqn{dU/dt}.}
#' }
#'
#' @section Output written to \code{@velocity}:
#' \describe{
#'   \item{\code{genes}}{Named logical vector. \code{TRUE} indicates
#'     high-quality velocity genes passing the \code{min_likelihood} threshold
#'     and internal scaling filters.}
#' }
#'
#' @param pv A \code{plantvelo} object. For dynamical modes,
#'   \code{@kinetics$params} must be populated by \code{recover_dynamics()}.
#' @param mode Character scalar or \code{NULL}. Supported values are
#'   \code{"dynamical"} and \code{"stochastic"}. If \code{NULL}, dynamical
#'   mode is used when \code{@kinetics$params} exists; otherwise stochastic
#'   mode is used. Default \code{NULL}.
#' @param min_likelihood Numeric. Minimum fitted likelihood required for a gene
#'   to be retained as a velocity gene in dynamical modes. Default
#'   \code{0.001}.
#' @param min_r2 Numeric. Reserved for future use; currently not applied in
#'   dynamical modes. Default \code{0.01}.
#' @param use_raw Logical. Use raw count layers instead of smoothed moments.
#'   Default \code{FALSE}.
#' @param vkey Character scalar. Name of the primary velocity layer written to
#'   \code{@layers}. Default \code{"velocity"}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{plantvelo} object with updated \code{@layers},
#'   \code{@velocity}, and \code{@misc} slots.
#'
#' @examples
#' \dontrun{
#' pv <- compute_velocity(pv)
#' }
#' @export
compute_velocity <- function(pv,
                             mode           = NULL,
                             min_likelihood = 0.001,
                             min_r2         = 0.01,
                             use_raw        = FALSE,
                             vkey           = "velocity",
                             verbose        = TRUE) {

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  params <- pv@kinetics[["params"]]

  # ----------------auto-detect velocity mode----------------

  if (is.null(mode)) {
    if (!is.null(params)) {
      params <- .validate_dynamics_schema(pv)
      mode <- "dynamical"
    } else {
      mode <- "stochastic"
      if (verbose)
        message(paste(
          "No fitted parameters found in @kinetics$params.",
          "Using stochastic mode.",
          "Run recover_dynamics() for the full dynamical model."
        ))
    }
  } else if (identical(mode, "dynamical")) {
    params <- .validate_dynamics_schema(pv)
  }

  if (!identical(mode, "dynamical") && !identical(mode, "stochastic")) {
    stop(sprintf(
      "Unknown mode '%s'. Use 'dynamical' or 'stochastic'.",
      mode
    ))
  }

  empty_dynamical_fit <- identical(mode, "dynamical") &&
    all(params[["fit_model"]] == "unfitted")
  if (!empty_dynamical_fit && identical(mode, "dynamical") && !use_raw &&
      (is.null(pv@moments[["Ms"]]) || is.null(pv@moments[["Mu"]]))) {
    stop("@moments$Ms/Mu not found. Please run compute_moments() first.")
  }

  if (verbose)
    message(sprintf("Computing velocities (mode = '%s')...", mode))

  # ----------------compute velocities by mode----------------

  if (identical(mode, "dynamical")) {
    pv <- .vel_2state(pv, min_likelihood, use_raw, vkey, verbose)
  } else if (identical(mode, "stochastic")) {
    pv <- .vel_stochastic(pv, use_raw, vkey, verbose)
  }
  pv@layers[[paste0(vkey, "_r")]] <- NULL

  # ----------------store run parameters----------------

  pv@misc[[paste0(vkey, "_params")]] <- list(
    mode           = mode,
    min_likelihood = min_likelihood,
    min_r2         = min_r2,
    use_raw        = use_raw
  )

  if (verbose)
    message(sprintf(
      "Velocity stored in @layers$%s. Velocity genes in @velocity$genes.",
      vkey
    ))

  return(pv)
}

# =============================================================================
# Internal helper functions
# =============================================================================

#' Reject legacy velocity parameter schemas
#'
#' @param params Kinetics parameter table, or \code{NULL}.
#' @return Invisibly returns \code{NULL}; errors when legacy columns are found.
#' @keywords internal
.reject_velocity_legacy_schema <- function(params) {
  if (is.null(params)) return(invisible(NULL))

  legacy_columns <- c(
    "fit_alpha_u", "fit_alpha_r", "fit_t_ur", "fit_t_rs",
    "fit_delta", "fit_eta", "fit_gamma_r", "fit_std_r", "fit_r0"
  )
  detected <- intersect(legacy_columns, colnames(params))
  if (length(detected) > 0L) {
    stop(sprintf(
      paste0(
        "Detected legacy dynamics columns (%s). Re-run recover_dynamics() ",
        "to create the '%s' schema; legacy parameters are not ",
        "converted automatically."
      ),
      paste(detected, collapse = ", "), .DYNAMICS_MODEL_VERSION
    ))
  }

  invisible(NULL)
}

#' Retrieve canonical velocity matrix axes
#'
#' @param pv A \code{plantvelo} object.
#' @return A list containing canonical \code{cell_ids} and \code{gene_ids}.
#' @keywords internal
.velocity_axes <- function(pv) {
  spliced <- pv@layers[["spliced"]]
  if (is.null(spliced) || length(dim(spliced)) != 2L) {
    stop("@layers$spliced must be a two-dimensional matrix.")
  }

  gene_ids <- rownames(spliced)
  cell_ids <- colnames(spliced)
  valid_names <- function(values) {
    !is.null(values) && length(values) > 0L && !anyNA(values) &&
      all(nzchar(values)) && !anyDuplicated(values)
  }
  if (!valid_names(gene_ids)) {
    stop("@layers$spliced must have unique, non-missing gene row names.")
  }
  if (!valid_names(cell_ids)) {
    stop("@layers$spliced must have unique, non-missing cell column names.")
  }

  list(cell_ids = cell_ids, gene_ids = gene_ids)
}

#' Validate fitted dynamics for velocity computation
#'
#' @param pv A \code{plantvelo} object.
#' @return The validated \code{@kinetics$params} data frame.
#' @keywords internal
.validate_dynamics_schema <- function(pv) {
  axes <- .velocity_axes(pv)
  params <- pv@kinetics[["params"]]
  if (is.null(params)) {
    stop(sprintf(
      paste0(
        "@kinetics$params is missing for schema '%s'. ",
        "Please run recover_dynamics() again."
      ),
      .DYNAMICS_MODEL_VERSION
    ))
  }
  if (!is.data.frame(params)) {
    stop(sprintf(
      paste0(
        "@kinetics$params must be a data frame using schema '%s'. ",
        "Please run recover_dynamics() again."
      ),
      .DYNAMICS_MODEL_VERSION
    ))
  }

  .reject_velocity_legacy_schema(params)

  recovery_metadata <- pv@misc[["recover_dynamics"]]
  version <- if (is.list(recovery_metadata)) {
    recovery_metadata[["model_version"]]
  } else {
    NULL
  }
  if (!identical(version, .DYNAMICS_MODEL_VERSION)) {
    stop(sprintf(
      paste0(
        "Incompatible dynamics schema: model_version must be exactly '%s'. ",
        "Please run recover_dynamics() again."
      ),
      .DYNAMICS_MODEL_VERSION
    ))
  }

  expected_columns <- c(.DEFAULT_PARS_NAMES, "fit_model")
  if (!identical(colnames(params), expected_columns)) {
    stop(sprintf(
      paste0(
        "Incompatible dynamics schema '%s'; expected columns: %s. ",
        "Please run recover_dynamics() again."
      ),
      .DYNAMICS_MODEL_VERSION, paste(expected_columns, collapse = ", ")
    ))
  }

  fit_model <- as.character(params[["fit_model"]])
  valid_fit_status <- c("fitted", "unfitted")
  if (anyNA(fit_model) || any(!fit_model %in% valid_fit_status)) {
    invalid_models <- unique(fit_model[is.na(fit_model) |
                                         !fit_model %in% valid_fit_status])
    invalid_models[is.na(invalid_models)] <- "NA"
    stop(sprintf(
      paste0(
        "Incompatible dynamics schema '%s'; invalid fit_model values: %s. ",
        "Please run recover_dynamics() again."
      ),
      .DYNAMICS_MODEL_VERSION, paste(invalid_models, collapse = ", ")
    ))
  }
  params[["fit_model"]] <- fit_model

  all_genes <- axes$gene_ids
  param_genes <- rownames(params)
  invalid_gene_names <- is.null(param_genes) || anyNA(param_genes) ||
    any(!nzchar(param_genes)) || anyDuplicated(param_genes)
  if (invalid_gene_names || length(param_genes) != length(all_genes) ||
      !setequal(param_genes, all_genes)) {
    stop(sprintf(
      paste0(
        "@kinetics$params using schema '%s' is not aligned to the genes ",
        "in @layers$spliced. Please run recover_dynamics() again."
      ),
      .DYNAMICS_MODEL_VERSION
    ))
  }
  params <- params[match(all_genes, param_genes), , drop = FALSE]

  params
}

#' Identify valid fitted genes
#'
#' @param params Validated kinetics parameter table.
#' @return Logical vector with one value per parameter-table row.
#' @keywords internal
.valid_velocity_fits <- function(params) {
  n <- if (is.data.frame(params)) nrow(params) else 0L
  required <- c(
    "fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_scaling"
  )
  if (!is.data.frame(params) ||
      any(!c("fit_model", required) %in% colnames(params))) {
    return(rep(FALSE, n))
  }

  valid <- !is.na(params[["fit_model"]]) &
    params[["fit_model"]] == "fitted"
  for (column in required) {
    values <- params[[column]]
    if (!is.numeric(values)) return(rep(FALSE, n))
    valid <- valid & is.finite(values)
  }
  valid & params[["fit_scaling"]] != 0
}

#' Retrieve and align an expression matrix (cells × genes)
#'
#' Prefer the smoothed moment in \code{@moments} specified by
#' \code{moment_key}. If unavailable, fall back to the raw layer in
#' \code{@layers} specified by \code{layer_key}, transposing it from
#' genes \eqn{\times} cells to cells \eqn{\times} genes. Named axes are
#' strictly matched and reordered to canonical spliced-layer names; unnamed
#' axes are accepted only at the exact canonical dimensions.
#'
#' Returns \code{NULL} if neither source is available.
#'
#' @param pv A \code{plantvelo} object.
#' @param moment_key Character scalar. Name of the matrix in \code{@moments}.
#' @param layer_key Character scalar. Name of the matrix in \code{@layers}.
#' @param use_raw Logical. Whether to force use of raw layers.
#' @param required Logical. Whether absence of both sources is an error.
#'
#' @return A dense numeric matrix of cells \eqn{\times} genes, or \code{NULL}.
#' @keywords internal
.get_expr_matrix <- function(pv, moment_key, layer_key, use_raw,
                             required = FALSE) {
  axes <- .velocity_axes(pv)
  if (!use_raw && !is.null(pv@moments[[moment_key]])) {
    aligned <- .align_internal_layer(
      pv@moments[[moment_key]], paste0("@moments$", moment_key),
      axes$cell_ids, axes$gene_ids,
      source_orientation = "cells_by_genes", required = TRUE
    )
    return(make_dense_matrix(aligned))
  }
  if (!is.null(pv@layers[[layer_key]])) {
    aligned <- .align_internal_layer(
      pv@layers[[layer_key]], paste0("@layers$", layer_key),
      axes$cell_ids, axes$gene_ids,
      source_orientation = "genes_by_cells", required = TRUE
    )
    return(make_dense_matrix(aligned))
  }
  if (required) {
    if (use_raw) {
      stop(sprintf("Required raw layer @layers$%s is missing.", layer_key))
    }
    stop(sprintf(
      "Neither @moments$%s nor @layers$%s is available.", moment_key, layer_key
    ))
  }
  NULL
}

#' Retrieve and align fitted latent time
#'
#' @param pv A \code{plantvelo} object.
#' @return A dense cells-by-genes fitted-time matrix, or \code{NULL}.
#' @keywords internal
.get_fit_time_matrix <- function(pv) {
  fit_t <- pv@layers[["fit_t"]]
  if (is.null(fit_t)) return(NULL)

  axes <- .velocity_axes(pv)
  aligned <- .align_internal_layer(
    fit_t, "@layers$fit_t", axes$cell_ids, axes$gene_ids,
    source_orientation = "genes_by_cells", required = TRUE
  )
  make_dense_matrix(aligned)
}

#' Filter velocity genes by likelihood and scaling
#'
#' Apply a lower bound on \code{fit_likelihood} and remove genes with extreme
#' \code{fit_scaling} values based on the 5th and 95th percentiles.
#'
#' Returns a named logical vector with length equal to the number of genes in
#' \code{params}.
#'
#' @param params data.frame of fitted kinetic parameters.
#' @param gene_subset Logical vector indicating genes with valid fitted
#'   parameters.
#' @param min_likelihood Numeric. Minimum allowed fitted likelihood.
#'
#' @return Named logical vector indicating retained velocity genes.
#' @keywords internal
.vel_genes_filter <- function(params, gene_subset, min_likelihood) {
  vgenes <- gene_subset

  # Likelihood filter
  if ("fit_likelihood" %in% colnames(params) && !is.null(min_likelihood)) {
    ll <- params[["fit_likelihood"]]
    ll[is.na(ll)] <- 0
    vgenes <- vgenes & (ll > min_likelihood)
  }

  # Extreme scaling filter
  if ("fit_scaling" %in% colnames(params)) {
    sc <- params[["fit_scaling"]]
    fitted_sc <- sc[gene_subset & is.finite(sc)]
    if (length(fitted_sc) > 0L) {
      lb <- quantile(fitted_sc, 0.05, na.rm = TRUE)
      ub <- quantile(fitted_sc, 0.95, na.rm = TRUE)
      lb <- min(lb, 0.03)
      ub <- max(ub, 3.0)
      vgenes <- vgenes & is.finite(sc) & sc >= lb & sc <= ub
    }
  }

  names(vgenes) <- rownames(params)
  vgenes
}

#' Compute two-state dynamical velocity
#'
#' Implements the scVelo-compatible U/S velocity equations, consistent with
#' \code{BaseDynamics$get_vt()} and \code{get_wt()} in the
#' \code{use_ir = FALSE} branch:
#' \itemize{
#'   \item \code{vt = β·U(t) - γ·S_obs}
#'   \item \code{wt = (α·o - β·U(t)) · scaling}
#' }
#'
#' If \code{fit_t} is available, \code{U(t)} is computed from the ODE solution.
#' Otherwise, the scaled observed unspliced expression is used directly.
#'
#' @param pv A \code{plantvelo} object with two-state kinetic parameters in
#'   \code{@kinetics$params}.
#' @param min_likelihood Numeric. Minimum likelihood threshold used for
#'   velocity gene filtering.
#' @param use_raw Logical. Use raw layers instead of smoothed moments.
#' @param vkey Character scalar. Base name for the output velocity layers.
#' @param verbose Logical. Print progress messages.
#'
#' @return A \code{plantvelo} object with \code{@layers$velocity},
#'   \code{@layers$velocity_u}, and \code{@velocity$genes} updated.
#' @keywords internal
.vel_2state <- function(pv, min_likelihood, use_raw, vkey, verbose) {
  params <- .validate_dynamics_schema(pv)

  axes       <- .velocity_axes(pv)
  all_genes  <- axes$gene_ids
  cell_ids   <- axes$cell_ids
  n_genes   <- length(all_genes)
  n_cells   <- length(cell_ids)

  gene_subset <- .valid_velocity_fits(params)
  vt_mat <- matrix(NA_real_, n_genes, n_cells,
                   dimnames = list(all_genes, cell_ids))
  wt_mat <- matrix(NA_real_, n_genes, n_cells,
                   dimnames = list(all_genes, cell_ids))

  if (!any(gene_subset)) {
    if (!all(params[["fit_model"]] == "unfitted")) {
      stop("No fitted genes found. Please run recover_dynamics() first.")
    }
    pv@layers[[vkey]]               <- vt_mat
    pv@layers[[paste0(vkey, "_u")]] <- wt_mat
    pv@layers[[paste0(vkey, "_r")]] <- NULL
    pv@velocity[["genes"]] <- .vel_genes_filter(
      params, gene_subset, min_likelihood
    )
    if (verbose) {
      message(sprintf(
        "  0 / 0 velocity genes selected (min_likelihood = %.4f).",
        min_likelihood
      ))
    }
    return(pv)
  }

  fitted_genes <- rownames(params)[gene_subset]

  U <- .get_expr_matrix(pv, "Mu", "unspliced", use_raw, required = TRUE)
  S <- .get_expr_matrix(pv, "Ms", "spliced",   use_raw, required = TRUE)
  fit_t <- .get_fit_time_matrix(pv)

  for (gene in fitted_genes) {
    g <- match(gene, all_genes)

    alpha_g   <- params[gene, "fit_alpha"]
    beta_g    <- params[gene, "fit_beta"]
    gamma_g   <- params[gene, "fit_gamma"]
    scaling_g <- params[gene, "fit_scaling"]
    t__g      <- params[gene, "fit_t_"]

    u_obs <- U[, g]
    s_obs <- S[, g]

    t_g      <- if (!is.null(fit_t)) fit_t[, g] else NULL
    has_time <- !is.null(t_g) && any(is.finite(t_g))

    if (has_time) {
      # Curve-based velocity
      vz  <- vectorize_2state(t_g, t__g, alpha_g, beta_g, gamma_g)
      ut  <- u_solution(vz$tau, vz$u0, vz$alpha, beta_g)
      o_g <- as.integer(t_g < t__g)

      vt_mat[g, ] <- beta_g * ut - gamma_g * s_obs
      wt_mat[g, ] <- (alpha_g * o_g - beta_g * ut) * scaling_g

    } else {
      # Residual velocity
      u_sc <- u_obs / scaling_g
      o_g  <- rep(1L, n_cells)

      vt_mat[g, ] <- beta_g * u_sc - gamma_g * s_obs
      wt_mat[g, ] <- (alpha_g * o_g - beta_g * u_sc) * scaling_g
    }
  }

  pv@layers[[vkey]]               <- vt_mat
  pv@layers[[paste0(vkey, "_u")]] <- wt_mat
  pv@layers[[paste0(vkey, "_r")]] <- NULL

  vgenes <- .vel_genes_filter(params, gene_subset, min_likelihood)
  pv@velocity[["genes"]] <- vgenes

  if (verbose) {
    n_vg <- sum(vgenes, na.rm = TRUE)
    message(sprintf("  %d / %d velocity genes selected (min_likelihood = %.4f).",
                    n_vg, sum(gene_subset), min_likelihood))
  }

  pv
}

#' Compute stochastic (residual) velocity
#'
#' Estimate per-gene \eqn{\gamma} from a no-intercept ordinary least squares
#' regression of \code{Mu} on \code{Ms}, and return the residual velocity
#' \eqn{Mu - \gamma \cdot Ms}. This mode is used as a fallback when no kinetic
#' parameters are available.
#'
#' @param pv A \code{plantvelo} object.
#' @param use_raw Logical. Use raw layers instead of smoothed moments.
#' @param vkey Character scalar. Name of the output velocity layer.
#' @param verbose Logical. Print progress messages.
#'
#' @return A \code{plantvelo} object with stochastic velocity stored in
#'   \code{@layers[[vkey]]} and velocity gene flags stored in
#'   \code{@velocity$genes}.
#' @keywords internal
.vel_stochastic <- function(pv, use_raw, vkey, verbose) {
  axes      <- .velocity_axes(pv)
  all_genes <- axes$gene_ids
  cell_ids  <- axes$cell_ids
  n_genes   <- length(all_genes)
  n_cells   <- length(cell_ids)

  Ms <- .get_expr_matrix(pv, "Ms", "spliced", use_raw, required = TRUE)
  Mu <- .get_expr_matrix(pv, "Mu", "unspliced", use_raw, required = TRUE)

  # No-intercept regression: gamma = sum(Ms * Mu) / sum(Ms^2) for each gene
  ss_xx <- colSums(Ms^2,   na.rm = TRUE)
  ss_xy <- colSums(Ms * Mu, na.rm = TRUE)
  gamma  <- ifelse(ss_xx > 0, ss_xy / ss_xx, 0)

  # Residual velocity (genes × cells)
  vt_mat <- t(Mu - sweep(Ms, 2, gamma, "*"))
  rownames(vt_mat) <- all_genes
  colnames(vt_mat) <- cell_ids

  pv@layers[[vkey]] <- vt_mat

  # Velocity genes: gamma must be positive and finite
  vgenes <- gamma > 0 & is.finite(gamma)
  names(vgenes) <- all_genes
  pv@velocity[["genes"]] <- vgenes

  if (verbose)
    message(sprintf("  %d stochastic velocity genes (gamma > 0).",
                    sum(vgenes, na.rm = TRUE)))

  pv
}
