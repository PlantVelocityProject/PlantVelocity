# =============================================================================
# PlantVelocity: recover_dynamics.R
# Full EM-based dynamics recovery module
# (base_dynamics.R and dynamics_recovery.R are integrated here)
#
# File structure:
#   Part 1 ── BaseDynamics R6 base class
#              (single-gene data container + fitting interface)
#              Corresponds to scVelo: tools/_em_model_utils.py::BaseDynamics
#              PlantVelocity extension: support for the R layer
#              (r, std_r, weights_r)
#   Part 2 ── DynamicsRecovery R6 class
#              (single-gene EM core, inherits from BaseDynamics)
#              Corresponds to scVelo: tools/_em_model_core.py::DynamicsRecovery
#              PlantVelocity implementation of two-state U/S kinetics
#   Part 3 ── recover_dynamics() S4 public API
#              (multi-gene dispatcher)
#              Corresponds to scVelo:
#              tools/_em_model_core.py::recover_dynamics()
#   Part 4 ── Internal helper functions
#              (.make_view / .resolve_genes / .fit_single_gene)
# =============================================================================


# =============================================================================
# Part 1: BaseDynamics R6 base class
# =============================================================================

#' BaseDynamics R6 class
#'
#' Container for single-gene U/S data and the current parameter estimates.
#' Provides the full computational interface required by the EM fitting
#' procedure.
#'
#' The two-state model is
#' \deqn{dU/dt = \alpha o - \beta U,}
#' \deqn{dS/dt = \beta U - \gamma S.}
#'
#' The subclass \code{DynamicsRecovery} inherits from this class and implements
#' the complete EM workflow.
#'
#' @export
BaseDynamics <- R6::R6Class(
  "BaseDynamics",
  public = list(
    # ----------------data fields----------------
    #' @field gene Gene name.
    gene    = NULL,
    #' @field u Smoothed unspliced expression vector (or raw values when
    #'   \code{use_raw = TRUE}).
    u       = NULL,
    #' @field s Smoothed spliced expression vector.
    s       = NULL,
    #' @field use_raw Logical. Whether raw unsmoothed counts are used.
    use_raw = FALSE,

    # ----------------parameter fields updated during EM----------------
    #' @field alpha Shared transcription rate \eqn{\alpha} into U.
    alpha   = NULL,
    #' @field beta Splicing rate \eqn{\beta}.
    beta    = NULL,
    #' @field gamma Degradation rate \eqn{\gamma}.
    gamma   = NULL,
    #' @field scaling U/S scaling factor.
    scaling = NULL,
    #' @field t_ Induction-to-repression transcription switch time.
    t_      = NULL,
    # ----------------weights and standard deviations----------------
    #' @field weights Logical mask of valid cells after filtering zero-valued
    #'   and extreme-expression observations.
    weights        = NULL,
    #' @field weights_upper Logical mask of high-expression cells used for
    #'   log-likelihood evaluation.
    weights_upper  = NULL,
    #' @field weights_outer Logical mask of cells outside the main trajectory,
    #'   used in differential dynamics tests.
    weights_outer  = NULL,
    #' @field std_u Standard deviation of the U layer, used for loss
    #'   normalisation.
    std_u   = 1,
    #' @field std_s Standard deviation of the S layer.
    std_s   = 1,

    # ----------------EM state----------------
    #' @field t Global latent time vector (\eqn{n\_cells}).
    t       = NULL,
    #' @field tau Within-phase relative time (induction or repression).
    tau     = NULL,
    #' @field o Cell state indicator (\code{1 = induction}, \code{0 = repression}).
    o       = NULL,
    #' @field loss Numeric vector of loss values across iterations.
    loss    = NULL,
    #' @field pars Matrix of parameter history
    #'   (\eqn{n\_parameters \times n\_iterations}).
    pars    = NULL,
    #' @field likelihood Final fitted likelihood value.
    likelihood  = NULL,
    #' @field varx Residual variance.
    varx    = NULL,
    #' @field recoverable Logical. Whether enough informative data points are
    #'   available for fitting.
    recoverable = TRUE,

    #' @field steady_state_ratio Steady-state U/S rate ratio used for
    #'   regularisation.
    steady_state_ratio = NULL,
    #' @field u0_ Value of \eqn{U_0} at the switch time (initial value of the
    #'   repression phase).
    u0_     = NULL,
    #' @field s0_ Value of \eqn{S_0} at the switch time.
    s0_     = NULL,

    # ----------------Nelder-Mead optimisation settings----------------
    #' @field max_iter Integer. Maximum number of EM iterations.
    max_iter    = 10L,
    #' @field perc Truncation percentile. Default \code{99}.
    perc        = 99,
    #' @field optim_method Optimisation method. Default \code{"Nelder-Mead"}.
    optim_method = "Nelder-Mead",
    #' @field optim_control Control list passed to \code{optim()}.
    optim_control = NULL,

    # ----------------neighbour connectivity for state smoothing----------------
    #' @field connectivities Neighbour connectivity matrix (sparse), used for
    #'   soft state assignment.
    connectivities = NULL,
    #' @field fit_connected_states Logical. Whether neighbour smoothing is used.
    fit_connected_states = TRUE,

    # ----------------------------------------------------------------
    # Constructor
    # ----------------------------------------------------------------

    #' @description Initialise a single-gene dynamics fitting object
    #'
    #' @param adata PlantVelocity internal \code{adata}-like object.
    #' @param gene Gene name (character scalar) or gene index (integer).
    #' @param u Optional U vector supplied directly, bypassing \code{adata}
    #'   extraction.
    #' @param s Optional S vector supplied directly.
    #' @param use_raw Logical. Use raw counts (\code{spliced}/\code{unspliced})
    #'   instead of smoothed layers (\code{Ms}/\code{Mu}).
    #' @param perc Truncation percentile. Default \code{99}.
    #' @param max_iter Integer. Maximum number of EM iterations. Default
    #'   \code{10L}.
    #' @param fit_scaling Logical. Whether to fit the U/S scaling ratio.
    #' @param fit_connected_states Logical. Whether to use neighbour connectivity
    #'   smoothing.
    #' @param fit_basal_transcription Logical. Whether to subtract basal
    #'   transcription background.
    #' @param init_vals Optional multipliers for EM initial values. \code{NULL}
    #'   uses defaults.
    initialize = function(adata, gene,
                          u = NULL, s = NULL,
                          use_raw = FALSE,
                          perc = 99,
                          max_iter = 10L,
                          fit_scaling = TRUE,
                          fit_connected_states = TRUE,
                          fit_basal_transcription = FALSE,
                          init_vals = NULL) {
      self$gene    <- gene
      self$perc    <- perc
      self$max_iter <- as.integer(max_iter)
      self$use_raw  <- use_raw

      # Run max_iter / 5 iterations per parameter block
      # (analogous to scVelo simplex_kwargs)
      self$optim_control <- list(
        maxit = max(as.integer(max_iter / 5), 2L)
      )

      # ----------------resolve gene index----------------
      if (is.character(gene)) {
        gene_idx <- which(adata$var_names == gene)
        if (length(gene_idx) == 0) stop(sprintf("Gene '%s' not found in adata.", gene))
        gene_idx <- gene_idx[1]
      } else {
        gene_idx <- as.integer(gene)
      }

      # ----------------read U/S data----------------
      if (is.null(u) || is.null(s)) {
        u_layer <- if (use_raw || is.null(adata$layers$Mu)) "unspliced" else "Mu"
        s_layer <- if (use_raw || is.null(adata$layers$Ms)) "spliced"   else "Ms"
        u <- as.vector(make_dense_matrix(
          adata$layers[[u_layer]][, gene_idx, drop = FALSE]
        ))
        s <- as.vector(make_dense_matrix(
          adata$layers[[s_layer]][, gene_idx, drop = FALSE]
        ))
      }
      self$u <- u
      self$s <- s

      # ----------------basal transcription correction----------------
      if (fit_basal_transcription) {
        u0_basal <- min(self$u, na.rm = TRUE)
        s0_basal <- min(self$s, na.rm = TRUE)
        self$u <- self$u - u0_basal
        self$s <- self$s - s0_basal
      }

      # Initialise parameters as NULL
      self$alpha <- self$beta <- self$gamma <- NULL
      self$scaling <- self$t_ <- NULL

      # Neighbour connectivity matrix (row-normalised)
      self$fit_connected_states <- fit_connected_states
      if (fit_connected_states && !is.null(adata$uns$neighbors)) {
        conn <- Matrix::Matrix(adata$uns$neighbors$connectivities, sparse = TRUE)
        rs   <- Matrix::rowSums(conn)
        rs[rs == 0] <- 1
        self$connectivities <- conn / rs
      }

      # ----------------initialise weights----------------
      tryCatch(
        self$initialize_weights(),
        warning = function(w) {
          self$recoverable <- FALSE
          message(sprintf("BaseDynamics: gene '%s' could not be initialized: %s",
                          self$gene, conditionMessage(w)))
        }
      )
      private$init_vals_mult <- init_vals
    },

    #' @description Initialise cell-weight masks
    #'
    #' Filters cells with \code{U = 0} or \code{S = 0} and truncates
    #' high-expression outliers. Also computes per-layer standard deviations
    #' for loss normalisation.
    #' @param weighted Logical. Whether to exclude high-expression outliers
    #'   when constructing the cell-weight mask.
    initialize_weights = function(weighted = TRUE) {
      nonzero_s <- self$s > 0
      nonzero_u <- self$u > 0

      weights <- nonzero_s & nonzero_u & is.finite(self$s) & is.finite(self$u)
      self$recoverable <- sum(weights) > 2
      if (!self$recoverable) return(invisible(NULL))

      if (weighted) {
        ub_s <- quantile(self$s[weights], self$perc / 100, na.rm = TRUE)
        ub_u <- quantile(self$u[weights], self$perc / 100, na.rm = TRUE)
        if (ub_s > 0) weights <- weights & (self$s <= ub_s)
        if (ub_u > 0) weights <- weights & (self$u <= ub_u)
      }

      self$weights <- weights
      u_w <- self$u[weights]
      s_w <- self$s[weights]

      self$std_u <- sd(u_w)
      self$std_s <- sd(s_w)
      if (self$std_u == 0 || is.na(self$std_u)) self$std_u <- 1
      if (self$std_s == 0 || is.na(self$std_s)) self$std_s <- 1

      # High-expression mask used for log-likelihood computation
      self$weights_upper <- weights
      if (any(weights)) {
        u_max <- max(u_w)
        s_max <- max(s_w)
        w_upper <- (self$u > u_max / 3) & (self$s > s_max / 3)
        self$weights_upper <- weights & w_upper
      }

      invisible(self)
    },

    # ----------------------------------------------------------------
    # Data accessors
    # ----------------------------------------------------------------

    #' @description Get the effective cell-weight mask
    #' @param weighted One of \code{"outer"}, \code{"upper"},
    #'   \code{TRUE}, or \code{FALSE}.
    #' @param weights_cluster Optional logical mask for a cluster subset.
    get_weights = function(weighted = TRUE, weights_cluster = NULL) {
      w <- switch(
        as.character(weighted),
        "outer" = self$weights_outer %||% self$weights,
        "upper" = self$weights_upper,
        "TRUE"  = self$weights,
        "FALSE" = rep(TRUE, length(self$u))
      )
      if (is.null(w)) w <- rep(TRUE, length(self$u))
      if (!is.null(weights_cluster) && length(w) == length(weights_cluster)) {
        w <- w & weights_cluster
      }
      w
    },

    #' @description Get scaled and weighted \code{(u, s)} data
    #' @param scaling Scaling factor. \code{NULL} uses \code{self$scaling}.
    #' @param weighted Weight mode.
    #' @param weights_cluster Optional cluster mask.
    #' @return Named list with components \code{u} and \code{s}.
    get_reads = function(scaling = NULL, weighted = NULL,
                         weights_cluster = NULL) {
      sc <- if (is.null(scaling)) self$scaling %||% 1 else scaling
      u  <- self$u / sc
      s  <- self$s

      if (!is.null(weighted) || !is.null(weights_cluster)) {
        w <- self$get_weights(weighted, weights_cluster)
        u <- u[w]; s <- s[w]
      }
      list(u = u, s = s)
    },

    #' @description Get the current parameter set
    #' @param alpha,beta,gamma,scaling,t_
    #'   Optional parameter overrides. \code{NULL} uses the corresponding
    #'   value stored in the object.
    #' @return Named list containing \code{alpha}, \code{beta}, \code{gamma},
    #'   \code{scaling}, and \code{t_}. Missing inputs fall back to stored
    #'   values.
    get_vars = function(alpha = NULL, beta = NULL, gamma = NULL,
                        scaling = NULL, t_ = NULL) {
      list(
        alpha   = alpha   %||% self$alpha   %||% 1,
        beta    = beta    %||% self$beta    %||% 1,
        gamma   = gamma   %||% self$gamma   %||% 1,
        scaling = scaling %||% self$scaling %||% 1,
        t_      = t_      %||% self$t_      %||% 1
      )
    },

    # ----------------------------------------------------------------
    # Residuals and loss functions
    # ----------------------------------------------------------------

    #' @description Compute normalised residuals between predictions and observations
    #'
    #' Uses the two-state U/S formulas.
    #'
    #' @param ... Optional parameter overrides
    #'   (\code{alpha}, \code{beta}, \code{gamma}, \code{scaling},
    #'   \code{t_}).
    #' @param weighted Weight mode.
    #' @param weights_cluster Optional logical mask for a cluster subset.
    #' @return Named list with components \code{udiff}, \code{rdiff},
    #'   \code{sdiff}, and \code{reg}.
    get_dists = function(..., weighted = TRUE, weights_cluster = NULL) {
      v <- self$get_vars(...)
      reads <- self$get_reads(v$scaling, weighted, weights_cluster)
      u <- reads$u; s <- reads$s
      n <- length(u)

      # Keep time slicing consistent with get_reads(): if weighting is applied,
      # self$t must be sliced accordingly. Otherwise vectorize_*state() returns
      # length = n_cells while u/s have length = n_weighted, which causes
      # recycling warnings in R.
      t_full <- self$t
      if (!is.null(t_full) && length(t_full) == length(self$u)) {
        if (!is.null(weighted) || !is.null(weights_cluster)) {
          w <- self$get_weights(weighted, weights_cluster)
          t_used <- t_full[w]
        } else {
          t_used <- t_full
        }
      } else {
        t_used <- NULL
      }
      # Fall back to a constant time if lengths still mismatch
      if (is.null(t_used) || length(t_used) != n) {
        fallback <- v$t_ / 2
        t_used   <- rep(fallback, n)
      }

      vz <- vectorize_2state(
        t_used,
        v$t_, v$alpha, v$beta, v$gamma
      )
      ut <- u_solution(vz$tau, vz$u0, vz$alpha, v$beta)
      st <- s_solution(vz$tau, vz$s0, vz$u0, vz$alpha, v$beta, v$gamma)

      udiff <- (ut - u) / self$std_u * v$scaling
      rdiff <- rep(0, n)
      sdiff <- (st - s) / self$std_s

      # Regularisation term based on steady-state ratio
      reg <- 0
      if (!is.null(self$steady_state_ratio)) {
        reg <- (v$gamma / v$beta - self$steady_state_ratio) *
          s / self$std_s
      }

      list(udiff = udiff, rdiff = rdiff, sdiff = sdiff, reg = reg)
    },

    #' @description Compute the total squared distance
    #' @details
    #' \code{distx = udiff^2 + rdiff^2 + sdiff^2} (plus \code{reg^2} when
    #' regularisation is enabled).
    #' @param ... Arguments passed to \code{get_dists()}.
    #' @param regularize Logical. Whether to include the regularisation term.
    get_distx = function(..., regularize = TRUE) {
      d <- self$get_dists(...)
      distx <- d$udiff^2 + d$rdiff^2 + d$sdiff^2
      if (regularize) distx <- distx + d$reg^2
      distx
    },

    #' @description Get total sum of squared errors (SSE)
    #' @param ... Arguments passed to \code{get_distx()}.
    get_se = function(...) sum(self$get_distx(...), na.rm = TRUE),

    #' @description Get mean squared error (MSE)
    #' @param ... Arguments passed to \code{get_distx()}.
    get_mse = function(...) mean(self$get_distx(...), na.rm = TRUE),

    #' @description Loss function used during EM optimisation
    #' @details Equal to the total sum of squared errors.
    #' @param ... Arguments passed to \code{get_se()}.
    get_loss = function(...) self$get_se(...),

    #' @description Compute the log-likelihood under a noise model
    #' @param ... Optional parameter overrides passed to \code{get_dists()}.
    #' @param varx Variance. \code{NULL} estimates it from the residuals.
    #' @param noise_model One of \code{"normal"} or \code{"laplace"}.
    #' @param weighted Weight mode passed to \code{get_dists()}.
    get_loglikelihood = function(..., varx = NULL, noise_model = "normal",
                                 weighted = "upper") {
      d <- self$get_dists(..., weighted = weighted)
      distx <- d$udiff^2 + d$rdiff^2 + d$sdiff^2 + d$reg^2

      eucl_distx <- sqrt(distx)
      n <- max(length(distx) - length(self$u) * 0.01, 2)

      if (is.null(varx)) {
        varx <- mean(distx, na.rm = TRUE) -
          mean(sign(d$sdiff) * eucl_distx, na.rm = TRUE)^2
      }
      if (varx == 0) varx <- 1e-10

      if (noise_model == "normal") {
        loglik <- -1 / 2 / n * sum(distx, na.rm = TRUE) / varx
        loglik <- loglik - 1 / 2 * log(2 * pi * varx)
      } else if (noise_model == "laplace") {
        loglik <- -1 / sqrt(2) / n * sum(eucl_distx, na.rm = TRUE) / sqrt(varx)
        loglik <- loglik - 1 / 2 * log(2 * varx)
      } else {
        stop(sprintf("Unsupported noise model: '%s'", noise_model))
      }
      loglik
    },

    #' @description Get the likelihood value
    #' @details Computed as \code{exp(log-likelihood)}.
    #' @param ... Arguments passed to \code{get_loglikelihood()}.
    get_likelihood = function(...) exp(self$get_loglikelihood(...)),

    #' @description Compute the residual variance
    #' @param ... Optional parameter overrides passed to \code{get_dists()}.
    #' @param weighted Weight mode. Default \code{"upper"}.
    get_variance = function(..., weighted = "upper") {
      d     <- self$get_dists(..., weighted = weighted)
      distx <- d$udiff^2 + d$rdiff^2 + d$sdiff^2
      mean(distx, na.rm = TRUE) -
        mean(sign(d$sdiff) * sqrt(distx), na.rm = TRUE)^2
    },

    # ----------------------------------------------------------------
    # Velocity vector calculations
    # ----------------------------------------------------------------

    #' @description Predict the fitted \code{U(t)} trajectory
    #' @param ... Optional parameter overrides passed to \code{get_vars()}.
    get_ut = function(...) {
      v <- self$get_vars(...)
      vz <- vectorize_2state(self$t, v$t_, v$alpha, v$beta, v$gamma)
      u_solution(vz$tau, vz$u0, vz$alpha, v$beta)
    },

    #' @description Predict the fitted \code{S(t)} trajectory
    #' @param ... Optional parameter overrides passed to \code{get_vars()}.
    get_st = function(...) {
      v <- self$get_vars(...)
      vz <- vectorize_2state(self$t, v$t_, v$alpha, v$beta, v$gamma)
      s_solution(vz$tau, vz$s0, vz$u0, vz$alpha, v$beta, v$gamma)
    },

    #' @description Compute the \eqn{dS/dt} velocity vector
    #' @details
    #' \deqn{v_t = \beta \cdot U(t) - \gamma \cdot S(t)}
    #' @param ... Optional parameter overrides passed to the trajectory methods.
    get_vt = function(...) {
      v <- self$get_vars(...)
      ut <- self$get_ut(...)
      st <- self$get_st(...)
      v$beta * ut - v$gamma * st
    },

    #' @description Compute the \eqn{dU/dt} velocity vector
    #' @details
    #' \deqn{w_t = [\alpha \cdot o - \beta \cdot U(t)] \cdot scaling}
    #' @param ... Optional parameter overrides passed to the trajectory methods.
    get_wt = function(...) {
      v  <- self$get_vars(...)
      o  <- if (!is.null(self$o)) self$o else rep(1, length(self$u))
      ut <- self$get_ut(...)
      (v$alpha * o - v$beta * ut) * v$scaling
    }
  ),

  private = list(
    init_vals_mult = NULL
  )
)


# =============================================================================
# Part 2: DynamicsRecovery R6 class
# =============================================================================

#' DynamicsRecovery R6 class
#'
#' Inherits from \code{BaseDynamics} and implements the complete EM-based
#' dynamics recovery workflow for a single gene.
#'
#' Handles parameter initialisation, multi-step coordinate-descent updates,
#' and parameter refinement.
#' When \code{refit_time = TRUE}, each proposal projects cells onto a bounded
#' theoretical trajectory grid using standardised Euclidean distance in U/S
#' space.
#'
#' EM workflow inside \code{fit()}:
#' \preformatted{
#' Pretraining:
#'   fit_t_and_alpha() -> fit_scaling_() -> fit_rates() -> fit_t_()
#'
#' Main loop:
#'   fit_all()   # alpha, beta, gamma, switch time, and optional scaling
#'
#' Finalisation:
#'   compute tau / o / likelihood / variance directly
#' }
#'
#' @export
DynamicsRecovery <- R6::R6Class(
  "DynamicsRecovery",
  inherit = BaseDynamics,
  public = list(

    #' @field fit_scaling Logical or positive numeric scalar controlling
    #'   whether the U/S scaling ratio is fitted or fixed.
    fit_scaling        = TRUE,
    #' @field fit_steady_states Logical. Whether to explicitly fit steady states.
    fit_steady_states  = TRUE,
    #' @field refit_time Logical. Whether to reassign time during EM.
    #'   \code{TRUE} corresponds to EM mode; \code{FALSE} keeps time fixed.
    refit_time         = FALSE,
    #' @field assignment_mode Time assignment mode. \code{NULL} defaults to
    #'   \code{"projection"}; no other modes are supported.
    assignment_mode    = NULL,
    #' @field steady_state_prior Optional logical prior mask for steady-state cells.
    steady_state_prior = NULL,
    #' @field pval_steady P-value from the steady-state bimodality test.
    pval_steady        = 1,
    #' @field steady_u Steady-state mean of the U layer.
    steady_u           = NULL,
    #' @field steady_s Steady-state mean of the S layer.
    steady_s           = NULL,
    #' @field high_pars_resolution Logical. Whether to store a high-resolution
    #'   parameter trajectory across iterations.
    high_pars_resolution = FALSE,

    # ----------------------------------------------------------------
    # Constructor
    # ----------------------------------------------------------------

    #' @description Initialise a \code{DynamicsRecovery} object
    #'
    #' @param adata PlantVelocity internal \code{adata}-like object.
    #' @param gene Gene name or gene index.
    #' @param load_pars Logical. Load existing fitted parameters from
    #'   \code{adata$var} when available.
    #' @param fit_scaling Logical or positive numeric scalar. \code{TRUE} fits
    #'   scaling, \code{FALSE} fixes it at one, and a numeric value fixes it
    #'   at that value. Default \code{TRUE}.
    #' @param fit_steady_states Logical. Whether to fit steady states.
    #'   Default \code{TRUE}.
    #' @param refit_time Logical. Whether parameter proposals trigger projection
    #'   of cells onto the fitted trajectory. Default \code{TRUE}.
    #' @param steady_state_prior Optional logical prior mask for steady-state cells.
    #' @param high_pars_resolution Logical. Whether to record a detailed
    #'   parameter trajectory. Default \code{FALSE}.
    #' @param ... Additional arguments passed to
    #'   \code{BaseDynamics$initialize()}.
    initialize = function(adata, gene,
                          load_pars          = FALSE,
                          fit_scaling        = TRUE,
                          fit_steady_states  = TRUE,
                          refit_time         = FALSE,
                          steady_state_prior = NULL,
                          high_pars_resolution = FALSE,
                          ...) {
      # Call parent constructor to load data and initialise weights
      super$initialize(adata, gene, ...)

      valid_fit_scaling <- isTRUE(fit_scaling) || identical(fit_scaling, FALSE) ||
        (is.numeric(fit_scaling) && length(fit_scaling) == 1L &&
           is.finite(fit_scaling) && fit_scaling > 0)
      if (!valid_fit_scaling) {
        stop("fit_scaling must be TRUE, FALSE, or a positive numeric scalar.")
      }
      self$fit_scaling        <- fit_scaling
      self$fit_steady_states  <- fit_steady_states
      self$refit_time         <- isTRUE(refit_time)
      self$steady_state_prior <- steady_state_prior
      self$high_pars_resolution <- high_pars_resolution

      if (!self$recoverable) return(invisible(self))

      if (isTRUE(load_pars)) {
        self$load_pars(adata, gene)
      } else {
        self$initialize_params()
      }
    },

    # ----------------------------------------------------------------
    # Parameter initialisation
    # ----------------------------------------------------------------

    #' @description Estimate EM initial parameters from data statistics
    #'
    #' Corresponds to \code{DynamicsRecovery.initialize()} in scVelo.
    #' PlantVelocity extends this step with initial estimates for
    #' \eqn{\delta}, \eqn{\eta}, and \eqn{\gamma_R}.
    initialize_params = function() {
      u   <- self$u;  s <- self$s
      w   <- self$weights
      u_w <- u[w];    s_w <- s[w]

      # ----------------estimate scaling----------------
      self$std_u <- sd(u_w);  self$std_s <- sd(s_w)
      if (self$std_u == 0 || is.na(self$std_u)) self$std_u <- 1
      if (self$std_s == 0 || is.na(self$std_s)) self$std_s <- 1

      if (isTRUE(self$fit_scaling)) {
        scaling <- self$std_u / self$std_s
      } else if (is.numeric(self$fit_scaling)) {
        scaling <- self$fit_scaling
      } else {
        scaling <- 1
      }
      scaling <- max(as.numeric(scaling)[1], 1e-6)
      u_sc <- u / scaling;  u_w_sc <- u_w / scaling

      # ----------------initialise beta = 1 and estimate gamma from high-quantile regression----------------
      perc_high <- 0.98
      w_s  <- s_w >= quantile(s_w, perc_high, na.rm = TRUE)
      w_u  <- u_w_sc >= quantile(u_w_sc, perc_high, na.rm = TRUE)
      w_g  <- if (!is.null(self$steady_state_prior)) w_s | self$steady_state_prior[w] else w_s

      beta  <- 1
      gamma <- linreg_slope(
        knn_convolve(u_w_sc, w_g),
        knn_convolve(s_w,    w_g)
      ) + 1e-6

      # Adaptive adjustment to avoid extremely small or large gamma
      if (gamma < 0.05 / scaling) gamma <- gamma * 1.2
      if (gamma > 1.5 / scaling)  gamma <- gamma / 1.2

      gamma <- max(gamma, 1e-6)

      # ----------------initialise shared alpha----------------
      u_inf <- mean(u_w_sc[w_u | w_s], na.rm = TRUE)
      s_inf <- mean(s_w[w_s],           na.rm = TRUE)
      alpha <- u_inf * beta

      # ----------------bimodality test for steady-state detection----------------
      bm   <- private$.test_bimodality(u_w_sc)
      bm_s <- private$.test_bimodality(s_w)

      self$pval_steady <- max(bm$pval, bm_s$pval)
      self$steady_u    <- bm$high_mean   * scaling
      self$steady_s    <- bm_s$high_mean

      if (self$pval_steady < 1e-3) {
        # Significant steady state: use steady-state constraints to refine alpha
        u_inf <- mean(c(u_inf, bm$high_mean), na.rm = TRUE)
        alpha <- gamma * s_inf
        beta  <- alpha / u_inf
      }

      # ----------------store initial parameters and switch state----------------
      t_ <- max(tau_inv_u(u_inf, 0, alpha, beta), 0.1)
      self$alpha   <- alpha
      self$beta    <- beta
      self$gamma   <- gamma
      self$scaling <- max(scaling, 1e-6)
      self$t_      <- t_

      self$u0_ <- u_solution(t_, 0, alpha, beta)
      self$s0_ <- s_solution(t_, 0, 0, alpha, beta, gamma)

      # ----------------initialise time assignment and loss----------------
      self$pars <- matrix(
        c(self$alpha, self$beta, self$gamma, self$t_, self$scaling),
        ncol = 1,
        dimnames = list(
          c("alpha","beta","gamma","t_","scaling"),
          NULL
        )
      )

      # Simple initial time assignment: uniform over [0, t_]
      n_cells <- length(self$u)
      self$t   <- seq(0, t_, length.out = n_cells)
      self$tau <- self$t
      self$o   <- as.integer(self$t < self$t_)
      self$loss <- numeric(0)
      self$loss <- c(self$loss, self$get_loss())

      # ----------------optimise scaling----------------
      if (isTRUE(self$fit_scaling)) {
        self$initialize_scaling(sight = 0.5)
        self$initialize_scaling(sight = 0.1)
      }

      self$steady_state_ratio <- self$gamma / self$beta

      invisible(self)
    },

    #' @description Perform a local line search around the current scaling value
    #' @param sight Numeric. Relative search range around the current value.
    initialize_scaling = function(sight = 0.5) {
      z_vals <- self$scaling + seq(-1, 1, length.out = 4) * self$scaling * sight
      for (z in z_vals) {
        self$update(scaling = z, beta = self$beta / self$scaling * z)
      }
      invisible(self)
    },

    #' @description Load existing fitted parameters from \code{adata$var}
    #' @details Only the current two-state schema is accepted.
    #' @param adata PlantVelocity internal \code{adata}-like object.
    #' @param gene Gene name or index.
    load_pars = function(adata, gene) {
      if (is.character(gene)) {
        idx <- which(adata$var_names == gene)[1]
      } else {
        idx <- as.integer(gene)
      }

      required_columns <- c(.DEFAULT_PARS_NAMES, "fit_model")
      missing_columns <- setdiff(required_columns, names(adata$var))
      if (length(missing_columns) > 0L) {
        stop(sprintf(
          paste0(
            "Cannot load fitted dynamics: incompatible or legacy schema; ",
            "missing columns: %s. Expected model_version '%s'."
          ),
          paste(missing_columns, collapse = ", "), .DYNAMICS_MODEL_VERSION
        ))
      }

      get_value <- function(column) {
        value <- adata$var[[column]][idx]
        if (length(value) == 1L) value else NULL
      }
      version_values <- c(
        adata$model_version %||% NULL,
        adata$uns[["recover_dynamics"]][["model_version"]] %||% NULL,
        if ("model_version" %in% names(adata$var)) {
          get_value("model_version")
        } else {
          NULL
        }
      )
      version_values <- as.character(version_values)
      if (length(version_values) == 0L || anyNA(version_values) ||
          any(version_values != .DYNAMICS_MODEL_VERSION)) {
        stop(sprintf(
          "Cannot load fitted dynamics: model_version must be '%s'.",
          .DYNAMICS_MODEL_VERSION
        ))
      }
      fitted_model <- as.character(get_value("fit_model"))
      if (length(fitted_model) != 1L || is.na(fitted_model) ||
          !identical(fitted_model, "fitted")) {
        stop("Cannot load fitted dynamics: fit_model must be 'fitted'.")
      }

      get_par <- function(key, required = TRUE) {
        v <- adata$var[[paste0("fit_", key)]]
        value <- v[idx]
        valid <- length(value) == 1L && is.numeric(value) && is.finite(value)
        if (!valid && required) {
          stop(sprintf(
            "Cannot load fitted dynamics: fit_%s must be finite for gene '%s'.",
            key, self$gene
          ))
        }
        if (valid) value else NULL
      }

      loaded_scaling <- max(get_par("scaling"), 1e-6)
      self$scaling <- if (identical(self$fit_scaling, FALSE)) {
        1
      } else if (is.numeric(self$fit_scaling)) {
        if (length(self$fit_scaling) != 1L ||
            !is.finite(self$fit_scaling) || self$fit_scaling <= 0) {
          stop("fit_scaling must be TRUE, FALSE, or a positive numeric scalar.")
        }
        self$fit_scaling
      } else {
        loaded_scaling
      }
      self$alpha   <- max(get_par("alpha"), 0)
      self$beta    <- max(get_par("beta"), 1e-6)
      self$gamma   <- max(get_par("gamma"), 1e-6)
      self$t_      <- max(get_par("t_"), 0.01)
      self$u0_ <- u_solution(self$t_,  0, self$alpha, self$beta)
      self$s0_ <- s_solution(self$t_,  0, 0, self$alpha, self$beta, self$gamma)

      self$pars <- matrix(
        c(self$alpha, self$beta, self$gamma, self$t_, self$scaling),
        ncol = 1,
        dimnames = list(
          c("alpha","beta","gamma","t_","scaling"),
          NULL
        )
      )

      n_cells <- length(self$u)
      self$t   <- seq(0, self$t_, length.out = n_cells)
      self$tau <- self$t
      self$o   <- as.integer(self$t < self$t_)
      self$steady_state_ratio <- self$gamma / self$beta
      self$loss <- c(self$get_loss())

      invisible(self)
    },

    # ----------------------------------------------------------------
    # Main EM workflow
    # ----------------------------------------------------------------

    #' @description Run the full EM-based dynamics fitting workflow
    #'
    #' Workflow:
    #' pretraining -> main loop -> finalisation
    #'
    #' @param assignment_mode Time assignment mode. \code{NULL} selects
    #'   \code{"projection"}; other modes are not supported.
    fit = function(assignment_mode = NULL) {
      assignment_mode <- assignment_mode %||% "projection"
      if (!identical(assignment_mode, "projection")) {
        stop(sprintf(
          "Unsupported assignment_mode '%s'; only 'projection' is available.",
          assignment_mode
        ))
      }
      self$assignment_mode <- assignment_mode

      if (self$max_iter <= 0) {
        self$update()
        return(invisible(self))
      }

      private$.pretrain()
      private$.run_em(assignment_mode)
      private$.finalize()

      invisible(self)
    },

    # ----------------------------------------------------------------
    # Single-parameter and multi-parameter optimisation
    # ----------------------------------------------------------------

    #' @description Optimise \eqn{\alpha} while keeping other parameters fixed
    #' @param sight Relative search range. Default \code{0.5}.
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_alpha = function(sight = 0.5, ...) {
      val <- self$alpha
      private$.fit_single_param(
        param_name = "alpha",
        init_val   = val,
        sight      = sight,
        fn         = function(x) self$get_mse(alpha = max(x[1], 0), ...),
        ...
      )
    },

    #' @description Optimise \eqn{\beta} while keeping other parameters fixed
    #' @param sight Relative search range. Default \code{0.5}.
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_beta = function(sight = 0.5, ...) {
      private$.fit_single_param("beta", self$beta, sight,
                                function(x) self$get_mse(beta = max(x[1], 1e-6), ...), ...)
    },

    #' @description Optimise \eqn{\gamma} while keeping other parameters fixed
    #' @param sight Relative search range. Default \code{0.5}.
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_gamma = function(sight = 0.5, ...) {
      private$.fit_single_param("gamma", self$gamma, sight,
                                function(x) self$get_mse(gamma = max(x[1], 1e-6), ...), ...)
    },

    #' @description Jointly optimise \code{(t_, alpha)}
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_t_and_alpha = function(...) {
      alpha_vals <- self$alpha + seq(-1, 1, length.out = 5) * self$alpha / 10
      for (a in alpha_vals) self$update(alpha = a)

      x0 <- c(self$t_, self$alpha)
      fn <- function(x) self$get_mse(t_ = max(x[1], 0.01),
                                     alpha = max(x[2], 0), ...)
      res <- private$.run_optim(fn, x0)
      self$update(t_ = res[1], alpha = res[2])
    },

    #' @description Jointly optimise shared transcription and processing rates
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_rates = function(...) {
      x0 <- c(self$alpha, self$gamma)
      fn <- function(x) self$get_mse(alpha = max(x[1], 0),
                                     gamma = max(x[2], 1e-6), ...)
      res <- private$.run_optim(fn, x0, tol = 1e-2)
      self$update(alpha = res[1], gamma = res[2])
    },

    #' @description Optimise the transcription switch time \code{t_}
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_t_ = function(...) {
      fn <- function(x) self$get_mse(t_ = max(x[1], 0.01), ...)
      res <- private$.run_optim(fn, self$t_)
      self$update(t_ = res[1])
    },

    #' @description Jointly optimise all rates \code{(alpha, beta, gamma)}
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_rates_all = function(...) {
      x0 <- c(self$alpha, self$beta, self$gamma)
      fn <- function(x) self$get_mse(alpha = max(x[1], 0),
                                     beta = max(x[2], 1e-6),
                                     gamma = max(x[3], 1e-6), ...)
      res <- private$.run_optim(fn, x0, tol = 1e-2)
      self$update(alpha = res[1], beta = res[2], gamma = res[3])
    },

    #' @description Jointly optimise \code{(t_, alpha, beta, gamma)}
    #' @details Main optimisation loop for the two-state model.
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_t_and_rates = function(...) {
      x0 <- c(self$t_, self$alpha, self$beta, self$gamma)
      fn <- function(x) self$get_mse(
        t_ = max(x[1], 0.01), alpha = max(x[2], 0),
        beta = max(x[3], 1e-6), gamma = max(x[4], 1e-6), ...
      )
      res <- private$.run_optim(fn, x0, tol = 1e-2)
      self$update(t_ = res[1], alpha = res[2], beta = res[3], gamma = res[4])
    },

    #' @description Jointly optimise \code{(t_, beta, scaling)}
    #' @details Scaling optimisation step.
    #' @param ... Additional parameter overrides used in the loss calculation.
    fit_scaling_ = function(...) {
      x0 <- c(self$t_, self$beta, self$scaling)
      fn <- function(x) self$get_mse(
        t_ = max(x[1], 0.01), beta = max(x[2], 1e-6),
        scaling = max(x[3], 1e-6), ...
      )
      res <- private$.run_optim(fn, x0)
      self$update(t_ = res[1], beta = res[2], scaling = res[3])
    },

    # ----------------------------------------------------------------
    # Core updater: accept/reject strategy
    # ----------------------------------------------------------------

    #' @description Parameter update controller
    #'
    #' A proposed parameter update is accepted only when its finite loss does
    #' not exceed the current loss. Also includes automatic adjustment of
    #' \code{t_} to avoid unrealistically distant switch times.
    #'
    #' Corresponds to \code{DynamicsRecovery.update()} in scVelo.
    #'
    #' @param t Optional finite, non-negative global time vector with one value
    #'   per cell. \code{NULL} uses trajectory projection when
    #'   \code{refit_time = TRUE}.
    #' @param t_ Transcription switch time.
    #' @param alpha,beta,gamma Optional kinetic
    #'   parameter updates.
    #'   \code{NULL} keeps the current value.
    #' @param scaling Optional scaling update.
    #' @param adjust_t_ Logical. Whether to automatically adjust \code{t_}.
    #'   Default \code{TRUE}.
    #' @return Logical. \code{TRUE} if the proposed update was accepted.
    update = function(t = NULL, t_ = NULL, alpha = NULL,
                      beta = NULL, gamma = NULL, scaling = NULL,
                      adjust_t_ = TRUE) {
      loss_prev <- if (length(self$loss) > 0) tail(self$loss, 1) else 1e9

      # Resolve proposed parameters, using current values for NULL inputs
      v <- self$get_vars(
        alpha = alpha, beta = beta, gamma = gamma, scaling = scaling, t_ = t_
      )
      v$t_ <- max(v$t_, 0.01)
      v$alpha <- max(v$alpha, 0)
      v$beta <- max(v$beta, 1e-6)
      v$gamma <- max(v$gamma, 1e-6)
      v$scaling <- max(v$scaling, 1e-6)

      if (!is.null(t) &&
          (length(t) != length(self$u) || any(!is.finite(t)) || any(t < 0))) {
        stop("t must contain one finite, non-negative value per cell.")
      }

      # Reassign time under the proposed parameters
      t_new <- if (is.null(t)) {
        private$.assign_time(v, force = FALSE)$t
      } else {
        t
      }

      # Compute loss under the proposed parameters
      loss_new <- private$.compute_loss_with(v, t_new)
      perform_update <- isTRUE(is.finite(loss_new) && loss_new <= loss_prev)

      # Automatic t_ adjustment to avoid overly distant switch times
      if (adjust_t_ && !is.null(self$o) && any(self$o == 1)) {
        on_cells <- self$o == 1
        if (any(on_cells)) {
          alt_t_ <- max(t_new[on_cells])
          alt_t_ <- alt_t_ + max(t_new) / length(t_new) * sum(t_new == v$t_)

          if (alt_t_ > 0 && alt_t_ < v$t_) {
            v_alt <- v;  v_alt$t_ <- alt_t_
            t_alt <- private$.assign_time(v_alt, force = FALSE)$t
            loss_alt <- private$.compute_loss_with(v_alt, t_alt)
            comparison_loss <- if (is.finite(loss_new)) {
              min(loss_new, loss_prev)
            } else {
              loss_prev
            }

            if (is.finite(loss_alt) && loss_alt <= comparison_loss) {
              v <- v_alt;  t_new <- t_alt;  loss_new <- loss_alt
              perform_update <- TRUE
            }
          }
        }
      }

      # Accept update
      if (perform_update) {
        self$alpha   <- v$alpha
        self$beta    <- v$beta
        self$gamma   <- v$gamma
        self$scaling <- v$scaling
        self$t_      <- v$t_

        self$t   <- t_new
        self$tau <- vectorize_2state(
          t_new, v$t_, v$alpha, v$beta, v$gamma
        )$tau
        self$o   <- as.integer(t_new < v$t_)

        new_pars <- c(v$alpha, v$beta, v$gamma, v$t_, v$scaling)
        self$pars <- cbind(self$pars, new_pars)
        self$loss <- c(self$loss, loss_new)

        self$u0_ <- u_solution(v$t_, 0, v$alpha, v$beta)
        self$s0_ <- s_solution(v$t_, 0, 0, v$alpha, v$beta, v$gamma)
      }

      perform_update
    }
  ),

  private = list(

    # ----------------pretraining phase----------------
    # Provides sufficiently good initial parameters before the main EM loop
    # and reduces the risk of poor local optima.
    .pretrain = function() {
      # Step 1: jointly optimise the shared transcription switch and rate
      self$fit_t_and_alpha()

      # Step 2: optimise scaling if requested
      if (isTRUE(self$fit_scaling)) self$fit_scaling_()

      # Step 3: jointly optimise rates
      self$fit_rates()

      # Step 4: optimise the single transcription switch
      self$fit_t_()

      invisible(self)
    },

    # ----------------main EM loop----------------
    # Two-state: 4-parameter joint optimisation
    # Two-state parameters, plus scaling when enabled
    .run_em = function(assignment_mode) {
      self$assignment_mode <- assignment_mode %||% "projection"
      self$fit_t_and_rates()
      if (isTRUE(self$refit_time)) {
        forced_time <- private$.assign_time(self$get_vars(), force = TRUE)$t
        self$update(t = forced_time, adjust_t_ = FALSE)
      }
      self$fit_t_and_rates()
      invisible(self)
    },

    # ----------------finalisation----------------
    # Store likelihood/variance and perform a simple convergence check.
    .finalize = function() {
      v <- self$get_vars()
      self$tau <- vectorize_2state(
        self$t, v$t_, v$alpha, v$beta, v$gamma
      )$tau
      self$o <- as.integer(self$t < v$t_)
      self$likelihood <- self$get_likelihood()
      self$varx       <- self$get_variance()

      # If relative improvement is below 1%, EM may not have converged
      if (length(self$loss) >= 2) {
        loss_init  <- self$loss[1]
        loss_final <- tail(self$loss, 1)
        rel_impr   <- (loss_init - loss_final) / (abs(loss_init) + 1e-10)
        if (rel_impr < 0.01) {
          warning(sprintf(
            paste0(
              "Gene '%s': EM may not have converged ",
              "(relative loss improvement: %.2f%%). ",
              "Consider increasing max_iter or checking data quality."
            ),
            self$gene, rel_impr * 100
          ))
        }
      }
      invisible(self)
    },

    # ----------------single-parameter optimisation----------------
    .fit_single_param = function(param_name, init_val, sight, fn, ...) {
      # Coarse search
      vals <- init_val + seq(-1, 1, length.out = 4) * init_val * sight
      for (v in vals) {
        args <- list();  args[[param_name]] <- v
        do.call(self$update, args)
      }
      # Nelder-Mead refinement
      res <- private$.run_optim(fn, init_val)
      args <- list();  args[[param_name]] <- res[1]
      do.call(self$update, args)
      invisible(self)
    },

    # ----------------run optimiser----------------
    # One-dimensional optimisation uses optimize();
    # multi-dimensional optimisation uses Nelder-Mead.
    .run_optim = function(fn, x0, tol = NULL) {
      # One-dimensional optimisation: Nelder-Mead is unstable and may warn
      # in the scalar case, so use stats::optimize() instead.
      if (length(x0) == 1L) {
        lower <- max(x0 * 0.01, 1e-8)
        upper <- x0 * 10 + 0.1
        if (lower >= upper) { lower <- 0; upper <- 1 }
        res <- tryCatch(
          stats::optimize(fn, interval = c(lower, upper))$minimum,
          error = function(e) x0
        )
        return(res)
      }

      # Multi-dimensional optimisation: Nelder-Mead
      ctrl <- self$optim_control %||% list(maxit = max(self$max_iter %/% 5, 2L))
      if (!is.null(tol)) ctrl$reltol <- tol

      res <- tryCatch(
        stats::optim(
          par     = x0,
          fn      = fn,
          method  = "Nelder-Mead",
          control = ctrl
        )$par,
        error = function(e) x0
      )
      res
    },

    # ----------------trajectory-projection time assignment----------------
    .assign_time = function(v, force = FALSE) {
      n <- length(self$u)
      if (n == 0L) return(list(t = numeric()))

      mode <- self$assignment_mode %||% "projection"
      if (!identical(mode, "projection")) {
        stop(sprintf(
          "Unsupported assignment_mode '%s'; only 'projection' is available.",
          mode
        ))
      }

      valid_current <- !is.null(self$t) && length(self$t) == n &&
        all(is.finite(self$t)) && all(self$t >= 0)
      if (!isTRUE(self$refit_time)) {
        if (force) {
          stop(
            "Cannot force trajectory projection when refit_time = FALSE."
          )
        }
        if (!valid_current) {
          stop(paste0(
            "refit_time = FALSE requires an existing finite, non-negative ",
            "time vector with one value per cell."
          ))
        }
        return(list(t = self$t))
      }

      current_finite <- if (is.null(self$t)) numeric() else {
        self$t[is.finite(self$t) & self$t >= 0]
      }
      grid_upper <- max(c(2 * v$t_, current_finite), na.rm = TRUE)
      grid_upper <- max(grid_upper, v$t_ + 0.01)
      grid_size <- min(
        800L,
        max(200L, as.integer(ceiling(sqrt(n) * 20)))
      )
      time_grid <- sort(unique(c(
        seq(0, grid_upper, length.out = grid_size),
        v$t_
      )))

      reads <- self$get_reads(v$scaling, weighted = NULL)
      vz <- vectorize_2state(
        time_grid, v$t_, v$alpha, v$beta, v$gamma
      )
      grid_u <- u_solution(vz$tau, vz$u0, vz$alpha, v$beta)
      grid_s <- s_solution(
        vz$tau, vz$s0, vz$u0, vz$alpha, v$beta, v$gamma
      )

      current_u <- current_s <- NULL
      if (valid_current) {
        current_vz <- vectorize_2state(
          self$t, v$t_, v$alpha, v$beta, v$gamma
        )
        current_u <- u_solution(
          current_vz$tau, current_vz$u0, current_vz$alpha, v$beta
        )
        current_s <- s_solution(
          current_vz$tau, current_vz$s0, current_vz$u0,
          current_vz$alpha, v$beta, v$gamma
        )
      }

      safe_scale <- function(value) {
        if (length(value) == 1L && is.finite(value) && value > 0) {
          max(value, 1e-8)
        } else {
          1
        }
      }
      scale_u <- safe_scale(self$std_u / v$scaling)
      scale_s <- safe_scale(self$std_s)

      t_new <- numeric(n)
      for (cell_idx in seq_len(n)) {
        distance_sq <- numeric(length(time_grid))
        current_distance_sq <- 0
        observed_dims <- 0L

        if (is.finite(reads$u[cell_idx])) {
          distance_sq <- distance_sq +
            ((grid_u - reads$u[cell_idx]) / scale_u)^2
          if (valid_current) {
            current_distance_sq <- current_distance_sq +
              ((current_u[cell_idx] - reads$u[cell_idx]) / scale_u)^2
          }
          observed_dims <- observed_dims + 1L
        }
        if (is.finite(reads$s[cell_idx])) {
          distance_sq <- distance_sq +
            ((grid_s - reads$s[cell_idx]) / scale_s)^2
          if (valid_current) {
            current_distance_sq <- current_distance_sq +
              ((current_s[cell_idx] - reads$s[cell_idx]) / scale_s)^2
          }
          observed_dims <- observed_dims + 1L
        }

        if (observed_dims == 0L) {
          t_new[cell_idx] <- if (valid_current) self$t[cell_idx] else 0
        } else {
          best_grid_idx <- which.min(distance_sq)
          if (valid_current &&
              current_distance_sq <= distance_sq[best_grid_idx]) {
            t_new[cell_idx] <- self$t[cell_idx]
          } else {
            t_new[cell_idx] <- time_grid[best_grid_idx]
          }
        }
      }

      list(t = t_new)
    },

    # ----------------compute loss with a given parameter set----------------
    .compute_loss_with = function(v, t) {
      old_t <- self$t
      self$t <- t
      on.exit({ self$t <- old_t }, add = TRUE)
      tryCatch(
        self$get_loss(
          alpha = v$alpha, beta = v$beta, gamma = v$gamma,
          scaling = v$scaling, t_ = v$t_
        ),
        error = function(e) Inf
      )
    },

    # ----------------bimodality test (simplified KDE-based version)----------------
    .test_bimodality = function(x, n_grid = 200) {
      x <- x[is.finite(x) & x > 0]
      if (length(x) < 10) return(list(pval = 1, high_mean = max(x, na.rm = TRUE)))

      dens <- tryCatch(stats::density(x, n = n_grid), error = function(e) NULL)
      if (is.null(dens)) return(list(pval = 1, high_mean = max(x, na.rm = TRUE)))

      # Find local minima of the density (valleys between modes)
      y <- dens$y
      local_mins <- which(diff(sign(diff(y))) > 0) + 1

      if (length(local_mins) == 0) {
        return(list(pval = 1, high_mean = quantile(x, 0.9, na.rm = TRUE)))
      }

      # Use the deepest valley as the bimodal split point
      min_idx   <- local_mins[which.min(y[local_mins])]
      split_val <- dens$x[min_idx]

      low_group  <- x[x <= split_val]
      high_group <- x[x >  split_val]

      if (length(low_group) < 3 || length(high_group) < 3) {
        return(list(pval = 1, high_mean = mean(high_group, na.rm = TRUE)))
      }

      # Approximate Hartigan dip-style significance with a separation proxy
      separation <- (mean(high_group) - mean(low_group)) / (sd(x) + 1e-10)
      pval <- if (separation > 3) 1e-4 else if (separation > 2) 1e-2 else 0.5

      list(pval = pval, high_mean = mean(high_group, na.rm = TRUE))
    }
  )
)


# =============================================================================
# Part 3: recover_dynamics() S4 public API
# =============================================================================

# Default parameter names written to @kinetics$params
# (rows = genes, columns = parameters)
.DYNAMICS_MODEL_VERSION <- "ir_excluded_2state_v1"

.DEFAULT_PARS_NAMES <- c(
  "fit_alpha", "fit_beta", "fit_gamma", "fit_t_", "fit_scaling",
  "fit_std_u", "fit_std_s", "fit_likelihood", "fit_u0", "fit_s0",
  "fit_variance"
)

.invalidate_dynamics_outputs <- function(pv) {
  pv@kinetics[["params"]] <- NULL
  pv@layers[c("fit_t", "fit_tau", "velocity", "velocity_u", "velocity_r")] <- NULL
  pv@graphs[c("velocity_graph", "velocity_graph_neg")] <- NULL
  pv@velocity <- list()
  pv@misc[c(
    "recover_dynamics", "velocity_params", "velocity_gene_stats",
    "ir_gene_stats"
  )] <- NULL
  pv
}

#' @title Recover RNA Velocity Dynamics
#' @description
#' Fit the intron-excluded two-state kinetic model independently for selected
#' genes:
#' \deqn{dU/dt = \alpha o - \beta U,}
#' \deqn{dS/dt = \beta U - \gamma S.}
#'
#' Dynamics use only \code{unspliced}/\code{Mu} and
#' \code{spliced}/\code{Ms}. Optional \code{ir}/\code{Mr} data are preserved
#' but never enter fitting, time assignment, likelihood, or velocity.
#'
#' @section Output:
#' The parameter table contains \code{fit_alpha}, \code{fit_beta},
#' \code{fit_gamma}, \code{fit_t_}, \code{fit_scaling}, \code{fit_std_u},
#' \code{fit_std_s}, \code{fit_likelihood}, \code{fit_u0}, \code{fit_s0},
#' \code{fit_variance}, and \code{fit_model}. The status column
#' \code{fit_model} contains only \code{"fitted"} or \code{"unfitted"}.
#' Run metadata records
#' \code{model_version = "ir_excluded_2state_v1"}.
#'
#' @param pv A \code{plantvelo} object.
#' @param genes Character vector or \code{NULL}. Genes to fit.
#' @param n_top_genes Non-negative integer or \code{NULL}. Number of variable
#'   genes selected when \code{genes = NULL}.
#' @param max_iter Positive integer. Maximum EM iterations.
#' @param fit_scaling Logical or positive numeric scalar controlling U/S
#'   scaling.
#' @param fit_connected_states Logical or \code{NULL}. Use neighbour-smoothed
#'   state assignment.
#' @param fit_basal_transcription Logical. Subtract a basal transcription
#'   baseline.
#' @param use_raw Logical. Use raw spliced and unspliced counts instead of
#'   smoothed moments.
#' @param add_key Reserved compatibility argument; must be \code{"fit"}.
#' @param verbose Logical. Print progress messages.
#'
#' @return A \code{plantvelo} object with updated kinetics, fitted times, and
#'   V0.4.0 schema metadata.
#'
#' @examples
#' \dontrun{
#' pv <- recover_dynamics(pv, n_top_genes = 2000)
#' pv <- recover_dynamics(pv, genes = rownames(pv@layers[["spliced"]]))
#' }
#' @export
recover_dynamics <- function(pv,
                             genes                   = NULL,
                             n_top_genes             = 2000L,
                             max_iter                = 10L,
                             fit_scaling             = TRUE,
                             fit_connected_states    = NULL,
                             fit_basal_transcription = FALSE,
                             use_raw                 = FALSE,
                             add_key                 = "fit",
                             verbose                 = TRUE) {

  # ----------------input validation----------------

  if (!inherits(pv, "plantvelo"))
    stop("`pv` must be a plantvelo object.")

  if (!is.null(genes) &&
      (!is.character(genes) || anyNA(genes))) {
    stop("`genes` must be NULL or a character vector without missing values.")
  }

  if (!is.null(n_top_genes)) {
    valid_n_top_genes <- is.numeric(n_top_genes) &&
      length(n_top_genes) == 1L && is.finite(n_top_genes) &&
      n_top_genes >= 0 && n_top_genes == floor(n_top_genes) &&
      n_top_genes <= .Machine$integer.max
    if (!valid_n_top_genes) {
      stop("`n_top_genes` must be NULL or a finite non-negative integer scalar.")
    }
    n_top_genes <- as.integer(n_top_genes)
  }

  valid_max_iter <- is.numeric(max_iter) && length(max_iter) == 1L &&
    is.finite(max_iter) && max_iter >= 1 && max_iter == floor(max_iter) &&
    max_iter <= .Machine$integer.max
  if (!valid_max_iter) {
    stop("`max_iter` must be a finite positive integer scalar.")
  }
  max_iter <- as.integer(max_iter)

  valid_fit_scaling <- isTRUE(fit_scaling) || identical(fit_scaling, FALSE) ||
    (is.numeric(fit_scaling) && length(fit_scaling) == 1L &&
       is.finite(fit_scaling) && fit_scaling > 0)
  if (!valid_fit_scaling) {
    stop("`fit_scaling` must be TRUE, FALSE, or a positive numeric scalar.")
  }

  .validate_logical_scalar <- function(value, name, allow_null = FALSE) {
    if (allow_null && is.null(value)) return(invisible(NULL))
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be a non-missing logical scalar%s.",
                   name, if (allow_null) " or NULL" else ""))
    }
    invisible(NULL)
  }
  .validate_logical_scalar(fit_connected_states, "fit_connected_states",
                           allow_null = TRUE)
  .validate_logical_scalar(fit_basal_transcription, "fit_basal_transcription")
  .validate_logical_scalar(use_raw, "use_raw")
  .validate_logical_scalar(verbose, "verbose")

  if (!identical(add_key, "fit")) {
    stop(paste(
      "`add_key` is fixed to \"fit\" by the V0.4.0 dynamics schema;",
      "custom parameter prefixes are not supported."
    ))
  }

  if (!use_raw) {
    if (is.null(pv@moments[["Ms"]]) || is.null(pv@moments[["Mu"]]))
      stop("@moments$Ms/Mu not found. Please run compute_moments() first.")
  }

  if (is.null(pv@graphs[["neighbors"]]))
    stop(paste(
      "No neighbor graph found in @graphs$neighbors.",
      "Please run build_neighbor_graph() first."
    ))

  if (is.null(fit_connected_states))
    fit_connected_states <- TRUE

  # ----------------build internal list view for the R6 EM classes----------------

  pv <- .invalidate_dynamics_outputs(pv)
  view <- .make_view(pv, use_raw = use_raw)
  all_genes <- view$var_names

  # ----------------resolve target gene list----------------

  gene_names <- .resolve_genes(view, genes, n_top_genes, use_raw)
  if (length(gene_names) == 0)
    stop("No genes selected for dynamics recovery.")

  if (verbose) {
    message(sprintf("Recovering dynamics for %d genes...", length(gene_names)))
  }

  # ----------------initialise output containers----------------

  n_cells   <- ncol(pv@layers[["spliced"]])
  n_genes   <- nrow(pv@layers[["spliced"]])

  T_mat <- matrix(NA_real_, nrow = n_cells, ncol = n_genes)
  Tau   <- matrix(NA_real_, nrow = n_cells, ncol = n_genes)

  pars_df <- as.data.frame(
    matrix(
      NA_real_,
      nrow     = n_genes,
      ncol     = length(.DEFAULT_PARS_NAMES),
      dimnames = list(all_genes, .DEFAULT_PARS_NAMES)
    )
  )
  pars_df[["fit_model"]] <- rep("unfitted", n_genes)

  fitted_genes <- character(0)
  # ----------------single-gene EM fitting loop----------------

  fit_args <- list(
    adata                   = view,
    use_raw                 = use_raw,
    max_iter                = max_iter,
    fit_scaling             = fit_scaling,
    fit_connected_states    = fit_connected_states,
    fit_basal_transcription = fit_basal_transcription
  )

  results <- lapply(gene_names, .fit_single_gene,
                    fit_args = fit_args, verbose = verbose)

  # ----------------collect results----------------

  pcol <- function(name) {
    paste0("fit_", name)
  }

  for (res in results) {
    if (is.null(res)) next
    dm   <- res$dm
    gene <- res$gene
    ix   <- which(all_genes == gene)
    if (length(ix) == 0) next

    fitted_genes <- c(fitted_genes, gene)

    if (!is.null(dm$t))   T_mat[, ix] <- dm$t
    if (!is.null(dm$tau)) Tau[, ix]   <- dm$tau

    pars_df[ix, pcol("alpha")]      <- dm$alpha      %||% NA
    pars_df[ix, pcol("beta")]       <- dm$beta       %||% NA
    pars_df[ix, pcol("gamma")]      <- dm$gamma      %||% NA
    pars_df[ix, pcol("t_")]         <- dm$t_         %||% NA
    pars_df[ix, pcol("scaling")]    <- dm$scaling    %||% NA
    pars_df[ix, pcol("std_u")]      <- dm$std_u      %||% NA
    pars_df[ix, pcol("std_s")]      <- dm$std_s      %||% NA
    pars_df[ix, pcol("likelihood")] <- dm$likelihood %||% NA
    pars_df[ix, pcol("u0")]         <- dm$u0_        %||% NA
    pars_df[ix, pcol("s0")]         <- dm$s0_        %||% NA
    pars_df[ix, pcol("variance")]   <- dm$varx       %||% NA

    pars_df[ix, pcol("model")]      <- res$fit_model
  }

  # ----------------write results back to the plantvelo object----------------

  pv@kinetics[["params"]] <- pars_df

  # Store layers as genes × cells, consistent with other layers
  fit_t_mat   <- t(T_mat)
  fit_tau_mat <- t(Tau)
  rownames(fit_t_mat) <- rownames(fit_tau_mat) <- all_genes
  colnames(fit_t_mat) <- colnames(fit_tau_mat) <- colnames(pv@layers[["spliced"]])
  pv@layers[["fit_t"]]   <- fit_t_mat
  pv@layers[["fit_tau"]] <- fit_tau_mat

  pv@misc[["recover_dynamics"]] <- list(
    fit_connected_states    = fit_connected_states,
    fit_basal_transcription = fit_basal_transcription,
    use_raw                 = use_raw,
    model_version           = .DYNAMICS_MODEL_VERSION,
    max_iter                = max_iter,
    add_key                 = add_key,
    n_fitted                = length(fitted_genes),
    n_unfitted              = sum(pars_df$fit_model == "unfitted")
  )

  n_fitted <- length(fitted_genes)
  n_total  <- length(gene_names)
  fit_rate <- if (n_total > 0) n_fitted / n_total else 0

  if (verbose) {
    message(sprintf(
      "Finished. Successfully fitted %d / %d genes (%.0f%%).",
      n_fitted, n_total, fit_rate * 100
    ))
    message("Results stored in @kinetics$params, @layers$fit_t, @layers$fit_tau.")
  }

  # Warn if global fit rate is too low
  if (fit_rate < 0.5 && n_total >= 10) {
    warning(sprintf(
      paste0(
        "recover_dynamics: only %.0f%% of genes were successfully fitted (%d / %d).\n",
        "Possible causes:\n",
        "  1. Insufficient cells per gene; try reducing n_top_genes.\n",
        "  2. Moments not computed; ensure compute_moments() was run with include_ir = TRUE.\n",
        "  3. Data quality issues; check for excessive zeros or outlier cells."
      ),
      fit_rate * 100, n_fitted, n_total
    ))
  }

  return(pv)
}


# =============================================================================
# Part 4: Internal helper functions
# =============================================================================

#' Align an expression layer to the internal cells-by-genes layout
#'
#' @param layer Matrix-like expression layer.
#' @param layer_name Character scalar used in diagnostics.
#' @param cell_ids Canonical cell names.
#' @param gene_names Canonical gene names.
#' @param source_orientation Either \code{"cells_by_genes"} or
#'   \code{"genes_by_cells"}.
#' @param required Logical. Whether an absent layer is an error.
#' @return The layer reordered to canonical cells by genes, preserving sparse
#'   matrix classes when possible.
#' @keywords internal
.align_internal_layer <- function(layer, layer_name, cell_ids, gene_names,
                                  source_orientation = c(
                                    "cells_by_genes", "genes_by_cells"
                                  ),
                                  required = FALSE) {
  source_orientation <- match.arg(source_orientation)
  if (is.null(layer)) {
    if (required) stop(sprintf("Required layer '%s' is missing.", layer_name))
    return(NULL)
  }
  if (length(dim(layer)) != 2L) {
    stop(sprintf("Layer '%s' must be a two-dimensional matrix.", layer_name))
  }

  source_row_names <- if (identical(source_orientation, "cells_by_genes")) {
    cell_ids
  } else {
    gene_names
  }
  source_col_names <- if (identical(source_orientation, "cells_by_genes")) {
    gene_names
  } else {
    cell_ids
  }
  expected_dim <- c(length(source_row_names), length(source_col_names))
  if (!identical(as.integer(dim(layer)), as.integer(expected_dim))) {
    stop(sprintf(
      "Layer '%s' has dimensions %d x %d; expected %d x %d (%s).",
      layer_name, nrow(layer), ncol(layer), expected_dim[1], expected_dim[2],
      source_orientation
    ))
  }

  align_axis <- function(observed_names, expected_names, axis_name) {
    if (is.null(observed_names)) return(seq_along(expected_names))
    if (anyNA(observed_names) || any(!nzchar(observed_names)) ||
        anyDuplicated(observed_names)) {
      stop(sprintf(
        "Layer '%s' has missing, empty, or duplicated %s names.",
        layer_name, axis_name
      ))
    }
    axis_index <- match(expected_names, observed_names)
    if (anyNA(axis_index)) {
      missing_names <- expected_names[is.na(axis_index)]
      stop(sprintf(
        "Layer '%s' %s names do not match the canonical object; missing: %s",
        layer_name, axis_name, paste(missing_names, collapse = ", ")
      ))
    }
    axis_index
  }

  row_index <- align_axis(rownames(layer), source_row_names, "row")
  col_index <- align_axis(colnames(layer), source_col_names, "column")
  layer <- layer[row_index, col_index, drop = FALSE]
  if (identical(source_orientation, "genes_by_cells")) layer <- t(layer)
  dimnames(layer) <- list(cell_ids, gene_names)
  layer
}

#' Build an internal list view for the R6 EM classes
#'
#' Converts a \code{plantvelo} S4 object into the list structure expected by
#' the R6 EM classes.
#'
#' All layers are stored as \code{cells × genes} (columns = genes, rows =
#' cells), matching the orientation of \code{@moments} and opposite to the
#' orientation of \code{@layers} (\code{genes × cells}).
#'
#' @param pv A \code{plantvelo} object.
#' @param use_raw Logical. Whether raw expression layers will be used
#'   downstream.
#'
#' @return A named list containing \code{var_names}, \code{obs_names},
#'   \code{layers}, and neighbour information.
#' @keywords internal
.make_view <- function(pv, use_raw = FALSE) {
  spliced <- pv@layers[["spliced"]]
  if (is.null(spliced) || length(dim(spliced)) != 2L) {
    stop("@layers$spliced must be a two-dimensional matrix.")
  }
  all_genes <- rownames(spliced)
  cell_ids  <- colnames(spliced)
  valid_names <- function(values) {
    !is.null(values) && length(values) > 0L && !anyNA(values) &&
      all(nzchar(values)) &&
      !anyDuplicated(values)
  }
  if (!valid_names(all_genes)) {
    stop("@layers$spliced must have unique, non-missing gene row names.")
  }
  if (!valid_names(cell_ids)) {
    stop("@layers$spliced must have unique, non-missing cell column names.")
  }

  layers <- list()

  layers[["Ms"]] <- .align_internal_layer(
    pv@moments[["Ms"]], "Ms", cell_ids, all_genes,
    source_orientation = "cells_by_genes", required = !use_raw
  )
  layers[["Mu"]] <- .align_internal_layer(
    pv@moments[["Mu"]], "Mu", cell_ids, all_genes,
    source_orientation = "cells_by_genes", required = !use_raw
  )
  layers[["spliced"]] <- .align_internal_layer(
    spliced, "spliced", cell_ids, all_genes,
    source_orientation = "genes_by_cells", required = TRUE
  )
  layers[["unspliced"]] <- .align_internal_layer(
    pv@layers[["unspliced"]], "unspliced", cell_ids, all_genes,
    source_orientation = "genes_by_cells", required = TRUE
  )
  list(
    var_names = all_genes,
    obs_names = cell_ids,
    n_obs     = length(cell_ids),
    n_vars    = length(all_genes),
    layers    = layers,
    model_version = .DYNAMICS_MODEL_VERSION,
    uns       = list(
      neighbors = pv@graphs[["neighbors"]],
      recover_dynamics = list(model_version = .DYNAMICS_MODEL_VERSION)
    )
  )
}

#' Extract one gene from an internal expression layer
#'
#' @param adata Internal list view.
#' @param layer_name Layer name.
#' @param gene Gene name.
#' @return Numeric vector, or \code{NULL} when the layer is absent.
#' @keywords internal
.get_view_gene_vector <- function(adata, layer_name, gene) {
  layer <- adata$layers[[layer_name]]
  if (is.null(layer)) return(NULL)

  layer_genes <- colnames(layer)
  gene_idx <- if (!is.null(layer_genes)) {
    match(gene, layer_genes)
  } else {
    match(gene, adata$var_names)
  }
  if (is.na(gene_idx) || gene_idx > ncol(layer)) {
    stop(sprintf("Gene '%s' is not aligned in layer '%s'.", gene, layer_name))
  }

  as.vector(make_dense_matrix(layer[, gene_idx, drop = FALSE]))
}

#' Resolve the target gene list
#'
#' Filters target genes against an aligned internal view. Automatic selection
#' optionally keeps only the top \eqn{N} genes ranked by mean expression.
#'
#' @param adata Internal list view with aligned expression layers.
#' @param genes Character vector of requested genes, or \code{NULL}.
#' @param n_top_genes Integer or \code{NULL}. Number of top genes retained by
#'   mean expression only when \code{genes = NULL}.
#' @param use_raw Logical. Whether to rank genes using raw instead of smoothed
#'   expression.
#'
#' @return Character vector of resolved gene names.
#' @keywords internal
.resolve_genes <- function(adata, genes, n_top_genes, use_raw) {
  all_genes <- adata$var_names

  if (!is.null(genes)) {
    unknown_genes <- unique(genes[!genes %in% all_genes])
    if (length(unknown_genes) > 0L) {
      warning(sprintf(
        "Ignoring unknown `genes`: %s",
        paste(unknown_genes, collapse = ", ")
      ), call. = FALSE)
    }
    return(unique(genes[genes %in% all_genes]))
  }

  genes <- all_genes
  if (!is.null(n_top_genes) && length(genes) > n_top_genes) {
    ranking_layer <- if (!use_raw && !is.null(adata$layers[["Ms"]])) {
      adata$layers[["Ms"]]
    } else {
      adata$layers[["spliced"]]
    }
    X_sub <- ranking_layer[, match(genes, colnames(ranking_layer)), drop = FALSE]
    mean_expr <- colMeans(X_sub, na.rm = TRUE)
    top_idx <- order(mean_expr, decreasing = TRUE)[seq_len(n_top_genes)]
    genes   <- genes[top_idx]
  }

  genes
}

#' Fit a single gene
#'
#' Creates a \code{DynamicsRecovery} R6 object for one gene and runs the full
#' EM fitting workflow.
#'
#' @param gene Character scalar. Gene name.
#' @param fit_args Named list of arguments passed to
#'   \code{DynamicsRecovery$new()}.
#' @param verbose Logical. Print progress messages.
#'
#' @return A list with components \code{dm}, \code{gene}, and
#'   \code{fit_model}, or \code{NULL} if fitting fails or the gene is not
#'   recoverable.
#' @keywords internal
.fit_single_gene <- function(gene, fit_args, verbose = FALSE) {
  adata                   <- fit_args$adata
  use_raw                 <- fit_args$use_raw
  max_iter                <- fit_args$max_iter
  fit_scaling             <- fit_args$fit_scaling
  fit_connected_states    <- fit_args$fit_connected_states
  fit_basal_transcription <- fit_args$fit_basal_transcription
  dm <- tryCatch(
    {
      DynamicsRecovery$new(
        adata                   = adata,
        gene                    = gene,
        use_raw                 = use_raw,
        max_iter                = max_iter,
        fit_scaling             = fit_scaling,
        fit_connected_states    = fit_connected_states,
        fit_basal_transcription = fit_basal_transcription
      )
    },
    error = function(e) {
      if (verbose)
        message(sprintf("  Gene '%s': initialization failed: %s", gene, e$message))
      NULL
    }
  )

  if (is.null(dm) || !dm$recoverable) {
    if (verbose && !is.null(dm))
      message(sprintf("  Gene '%s': not recoverable (insufficient data).", gene))
    return(NULL)
  }

  fit_ok <- tryCatch({
    dm$fit()
    TRUE
  }, error = function(e) {
    if (verbose)
      message(sprintf("  Gene '%s': fit failed: %s", gene, e$message))
    FALSE
  })

  if (!fit_ok) return(NULL)

  if (verbose) {
    loss_val <- tail(dm$loss, 1)
    message(sprintf(
      "  Gene '%s': fitted (loss = %.4f)",
      gene, if (length(loss_val) > 0) loss_val else NA
    ))
  }

  list(
    dm = dm,
    gene = gene,
    fit_model = "fitted"
  )
}
