# =============================================================================
# PlantVelocity: utils.R
# General utility functions + dynamics utility functions
# (the former dynamics_utils.R has been integrated here)
#
# File structure:
#   Part 1 ── General utility functions
#   Part 2 ── Dynamics utility functions
#              (used internally for EM parameter estimation)
# =============================================================================


# =============================================================================
# Part 1: General utility functions
# =============================================================================

#' NULL-coalescing operator
#'
#' Returns \code{a} if \code{a} is not \code{NULL}; otherwise returns
#' \code{b}.
#'
#' @param a Primary value.
#' @param b Fallback value.
#'
#' @return Either \code{a} or \code{b}.
#' @name null_coalesce
#' @aliases %||%
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Safely compute reciprocals
#'
#' Returns 0 for values close to zero, avoiding division-by-zero
#' \code{Inf}/\code{NaN}.
#'
#' Corresponds to scVelo: \code{core/_arithmetic.py::invert()}.
#'
#' @param x Numeric vector or matrix.
#' @param eps Numeric. Threshold below which values are treated as zero.
#'   Default \code{1e-10}.
#'
#' @return Numeric object of the same shape as \code{x}, containing safe
#'   reciprocals.
#' @export
invert_safe <- function(x, eps = 1e-10) {
  result <- ifelse(abs(x) < eps, 0, 1 / x)
  return(result)
}

#' Compute a clipped logarithm
#'
#' Computes \code{log(x)} after clipping \code{x} to the interval
#' \code{[lb + eps, ub]}, preventing \code{log(<= 0)}.
#'
#' Corresponds to scVelo: \code{core/_arithmetic.py::clipped_log()}.
#'
#' @param x Numeric input vector.
#' @param lb Numeric. Lower bound. Default \code{0}.
#' @param ub Numeric. Upper bound. Default \code{1}.
#' @param eps Numeric. Offset added to the lower bound. Default \code{1e-6}.
#'
#' @return Numeric vector of clipped log-values.
#' @export
clipped_log <- function(x, lb = 0, ub = 1, eps = 1e-6) {
  x_clipped <- pmax(pmin(x, ub), lb + eps)
  log(x_clipped)
}

#' Normalise probabilities along a matrix margin
#'
#' Performs L1 normalisation along the specified matrix margin so that each row
#' or column sums to 1.
#'
#' Corresponds to scVelo: \code{tools/_em_model_utils.py::normalize()}.
#'
#' @param X Numeric matrix.
#' @param margin Integer. Normalisation direction: \code{1} for rows,
#'   \code{2} for columns.
#' @param min_confidence Numeric or \code{NULL}. Optional offset added to each
#'   row or column sum to avoid zero totals.
#'
#' @return Normalised numeric matrix.
#' @export
normalize_probs <- function(X, margin = 2, min_confidence = NULL) {
  X_sum <- apply(X, margin, sum)
  if (!is.null(min_confidence)) {
    X_sum <- X_sum + min_confidence
  }
  X_sum[X_sum == 0] <- 1
  if (margin == 1) {
    sweep(X, 1, X_sum, "/")
  } else {
    sweep(X, 2, X_sum, "/")
  }
}

#' Apply KNN-weighted convolution
#'
#' Applies a sparse or dense neighbour-weight matrix to a vector or matrix,
#' yielding a KNN-weighted average.
#'
#' Corresponds to scVelo: \code{tools/_em_model_utils.py::convolve()}.
#'
#' @param x Numeric vector or matrix.
#' @param weights Weight matrix (sparse or dense). If \code{NULL}, \code{x} is
#'   returned unchanged.
#'
#' @return Weighted result with compatible dimensions.
#' @export
knn_convolve <- function(x, weights = NULL) {
  if (is.null(weights)) return(x)
  if (inherits(weights, "sparseMatrix")) {
    as.matrix(weights %*% x)
  } else {
    weights %*% x
  }
}

#' Compute a no-intercept linear-regression slope
#'
#' Computes the slope of the no-intercept regression of \code{u} on \code{s}.
#' Commonly used to initialise \eqn{\gamma/\beta}.
#'
#' Corresponds to scVelo: \code{tools/_em_model_utils.py::linreg()}.
#'
#' @param u Numeric response vector (typically unspliced).
#' @param s Numeric predictor vector (typically spliced).
#'
#' @return Numeric scalar slope:
#'   \eqn{\sum(su) / \sum(s^2)}.
#' @export
linreg_slope <- function(u, s) {
  ss <- sum(s^2)
  us <- sum(s * u)
  if (ss == 0) return(0)
  us / ss
}

#' Scale values to the interval \code{[0, 1]}
#'
#' Linearly rescales a numeric vector to the interval \code{[0, 1]}.
#'
#' Corresponds to scVelo: \code{tools/utils.py::scale()}.
#'
#' @param x Numeric vector.
#' @param min_val Numeric. Target minimum value. Default \code{0}.
#' @param max_val Numeric. Target maximum value. Default \code{1}.
#'
#' @return Rescaled numeric vector.
#' @export
scale_to_01 <- function(x, min_val = 0, max_val = 1) {
  x_range <- range(x, na.rm = TRUE)
  if (x_range[1] == x_range[2]) return(rep(0, length(x)))
  (x - x_range[1]) / (x_range[2] - x_range[1]) * (max_val - min_val) + min_val
}

#' Convert to a dense numeric matrix
#'
#' Converts a sparse matrix or data frame into a standard dense numeric matrix.
#'
#' Corresponds to scVelo: \code{tools/utils.py::make_dense()}.
#'
#' @param X Sparse matrix, data frame, or ordinary matrix.
#'
#' @return Dense numeric matrix.
#' @export
make_dense_matrix <- function(X) {
  if (inherits(X, "sparseMatrix")) {
    as.matrix(X)
  } else if (is.data.frame(X)) {
    as.matrix(X)
  } else {
    as.matrix(X)
  }
}


# =============================================================================
# Part 2: Dynamics utility functions
# These functions are used internally for EM parameter estimation
# (recover_dynamics / base_dynamics / dynamics_recovery).
# End users typically do not need to call them directly.
# Corresponds to the first half of
# scVelo: tools/_em_model_utils.py
# =============================================================================

# -----------------------------------------------------------------------------
# Part 2.1: ODE analytical solutions
# -----------------------------------------------------------------------------

.validate_nonnegative_time <- function(time, name = "tau") {
  if (!is.numeric(time) || !is.null(dim(time))) {
    stop(sprintf("`%s` must be a numeric vector.", name), call. = FALSE)
  }
  if (any(is.infinite(time))) {
    stop(sprintf("`%s` must not contain infinite values.", name),
         call. = FALSE)
  }
  if (any(time < 0, na.rm = TRUE)) {
    stop(sprintf("`%s` must not contain negative values.", name), call. = FALSE)
  }
  invisible(time)
}

.prepare_convolution_args <- function(tau, rates, rate_names) {
  .validate_nonnegative_time(tau)
  if (any(!vapply(rates, is.numeric, logical(1)))) {
    stop("Convolution rates must be numeric.", call. = FALSE)
  }
  if (length(tau) == 0L) {
    if (any(vapply(rates, length, integer(1)) != 1L)) {
      stop("Convolution rates must be scalars when `tau` is empty.",
           call. = FALSE)
    }
    return(list(tau = numeric(), rates = rates))
  }

  lengths <- c(length(tau), vapply(rates, length, integer(1)))
  common_length <- max(lengths)
  if (any(!lengths %in% c(1L, common_length))) {
    stop(
      paste0(
        "`tau` and convolution rates (`",
        paste(rate_names, collapse = "`, `"),
        "`) must each have length 1 or the common length."
      ),
      call. = FALSE
    )
  }

  list(
    tau = rep_len(as.numeric(tau), common_length),
    rates = lapply(rates, function(rate) {
      rep_len(as.numeric(rate), common_length)
    })
  )
}

.validate_convolution_tol <- function(tol) {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a finite positive numeric scalar.", call. = FALSE)
  }
  invisible(tol)
}

.exp_convolution <- function(tau, source_rate, target_rate,
                             tol = sqrt(.Machine$double.eps)) {
  .validate_convolution_tol(tol)
  args <- .prepare_convolution_args(
    tau,
    list(source_rate, target_rate),
    c("source_rate", "target_rate")
  )
  if (length(args$tau) == 0L) return(numeric())
  tau <- args$tau
  source_rate <- args$rates[[1L]]
  target_rate <- args$rates[[2L]]

  rate_gap <- abs(target_rate - source_rate)
  scaled_gap <- rate_gap * abs(tau)
  ratio <- -expm1(-scaled_gap) / scaled_gap
  close <- is.finite(scaled_gap) & scaled_gap <= tol
  x <- scaled_gap[close]
  ratio[close] <- 1 - x / 2 + x^2 / 6 - x^3 / 24 + x^4 / 120

  tau * exp(-pmin(source_rate, target_rate) * tau) * ratio
}

.exp_convolution2 <- function(tau, source_rate, middle_rate, target_rate,
                              tol = sqrt(.Machine$double.eps)) {
  .validate_convolution_tol(tol)
  args <- .prepare_convolution_args(
    tau,
    list(source_rate, middle_rate, target_rate),
    c("source_rate", "middle_rate", "target_rate")
  )
  if (length(args$tau) == 0L) return(numeric())
  tau <- args$tau
  source_rate <- args$rates[[1L]]
  middle_rate <- args$rates[[2L]]
  target_rate <- args$rates[[3L]]

  lowest_rate <- pmin(source_rate, middle_rate, target_rate)
  highest_rate <- pmax(source_rate, middle_rate, target_rate)
  central_rate <- source_rate + middle_rate + target_rate -
    lowest_rate - highest_rate
  rate_span <- highest_rate - lowest_rate
  scaled_span <- abs(rate_span * tau)
  use_series <- is.finite(scaled_span) & scaled_span <= sqrt(tol)

  out <- (
    .exp_convolution(tau, lowest_rate, central_rate, tol) -
      .exp_convolution(tau, highest_rate, central_rate, tol)
  ) / rate_span

  if (any(use_series)) {
    tau_series <- tau[use_series]
    mean_rate <- (
      source_rate[use_series] + middle_rate[use_series] +
        target_rate[use_series]
    ) / 3
    z1 <- (source_rate[use_series] - mean_rate) * tau_series
    z2 <- (middle_rate[use_series] - mean_rate) * tau_series
    z3 <- (target_rate[use_series] - mean_rate) * tau_series
    elementary1 <- z1 + z2 + z3
    elementary2 <- z1 * z2 + z1 * z3 + z2 * z3
    elementary3 <- z1 * z2 * z3

    series_sum <- rep(1 / 2, length(tau_series))
    h_nm3 <- numeric(length(tau_series))
    h_nm2 <- numeric(length(tau_series))
    h_nm1 <- rep(1, length(tau_series))
    for (order in seq_len(12L)) {
      h_n <- elementary1 * h_nm1 - elementary2 * h_nm2 +
        elementary3 * h_nm3
      series_sum <- series_sum + (-1)^order * h_n / factorial(order + 2L)
      h_nm3 <- h_nm2
      h_nm2 <- h_nm1
      h_nm1 <- h_n
    }

    out[use_series] <- tau_series^2 * exp(-mean_rate * tau_series) *
      series_sum
  }

  out
}

#' Analytical solution for unspliced mRNA U(t)
#'
#' Analytical solution for \eqn{U} in the two-state model, where
#' \eqn{\kappa = \beta}.
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::unspliced()}.
#'
#' @param tau Non-negative numeric relative time within the current phase.
#'   Missing values are preserved; infinite values are rejected.
#' @param u0 Initial amount of \eqn{U}.
#' @param alpha Transcription rate into \eqn{U}.
#' @param kappa_u Two-state departure rate from \eqn{U}, \eqn{\beta}.
#'
#' @return Predicted \eqn{U(\tau)}.
#' @export
u_solution <- function(tau, u0, alpha, kappa_u) {
  .validate_nonnegative_time(tau)
  u0 * exp(-kappa_u * tau) +
    alpha * .exp_convolution(tau, 0, kappa_u)
}

#' Analytical solution for spliced mRNA S(t) in the two-state model
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::spliced()}.
#'
#' @param tau Non-negative numeric relative time. Missing values are preserved;
#'   infinite values are rejected.
#' @param s0 Initial amount of \eqn{S}.
#' @param u0 Initial amount of \eqn{U}.
#' @param alpha Transcription rate.
#' @param beta Splicing rate.
#' @param gamma Degradation rate.
#'
#' @return Predicted \eqn{S(\tau)} in the two-state model.
#' @export
s_solution <- function(tau, s0, u0, alpha, beta, gamma) {
  .validate_nonnegative_time(tau)
  s0 * exp(-gamma * tau) + beta * (
    u0 * .exp_convolution(tau, beta, gamma) +
      alpha * .exp_convolution2(tau, 0, beta, gamma)
  )
}

# -----------------------------------------------------------------------------
# Part 2.2: Inverse time mapping
# -----------------------------------------------------------------------------

#' Infer relative time \eqn{\tau} from U values
#'
#' Given observed \eqn{U} values, inverts the analytical solution for
#' \eqn{U(t)} to recover the corresponding relative time \eqn{\tau}.
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::tau_inv()} (U branch).
#'
#' @param u Observed \eqn{U} values.
#' @param u0 Initial amount of \eqn{U}.
#' @param alpha Transcription rate into \eqn{U}.
#' @param kappa Total consumption rate of \eqn{U}.
#'
#' @return Numeric vector of inferred times \eqn{\tau}.
#' @export
tau_inv_u <- function(u, u0, alpha, kappa) {
  u_inf <- alpha / kappa
  -1 / kappa * clipped_log((u - u_inf) / (u0 - u_inf))
}

#' Infer time in the two-state model from U and S
#'
#' Given observed \eqn{(u, s)} values, selects the more stable inversion path
#' (the U-based or S-based branch) to recover \eqn{\tau}.
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::tau_inv()}.
#'
#' @param u Observed U vector.
#' @param s Observed S vector. If \code{NULL}, only the U-based path is used.
#' @param u0 Initial U value.
#' @param s0 Initial S value.
#' @param alpha Transcription rate.
#' @param beta Splicing rate.
#' @param gamma Degradation rate. If \code{NULL}, only the U-based path is
#'   used.
#'
#' @return Numeric vector of inferred times \eqn{\tau}.
#' @export
tau_inv_2state <- function(u, s = NULL, u0 = 0, s0 = 0,
                           alpha = NULL, beta = NULL, gamma = NULL) {
  inv_u <- if (!is.null(gamma)) (gamma >= beta) else TRUE
  inv_us <- !inv_u
  any_invu  <- any(inv_u)  || is.null(s)
  any_invus <- any(inv_us) && !is.null(s)

  tau <- numeric(length(u))

  if (any_invus) {
    tau_s  <- -1 / gamma * clipped_log(
      (s - beta * invert_safe(gamma - beta) * u -
         (alpha / gamma - beta * invert_safe(gamma - beta) * (alpha / beta))) /
        (s0 - beta * invert_safe(gamma - beta) * u0 -
           (alpha / gamma - beta * invert_safe(gamma - beta) * (alpha / beta)))
    )
    tau[inv_us] <- tau_s[inv_us]
  }

  if (any_invu) {
    tau_u <- tau_inv_u(u, u0, alpha, beta)
    if (any_invus) {
      tau[inv_u] <- tau_u[inv_u]
    } else {
      tau <- tau_u
    }
  }

  tau
}


# -----------------------------------------------------------------------------
# Part 2.3: Time vectorisation
# -----------------------------------------------------------------------------

#' Vectorise global time for the two-state model
#'
#' Splits global time \code{t} at the switch time \code{t_} into induction and
#' repression phases, and assigns the appropriate initial conditions and
#' transcription rate to each cell.
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::vectorize()}.
#'
#' @param t Global time vector (\eqn{n\_cells}).
#' @param t_ Switch time from induction to repression.
#' @param alpha Induction-phase transcription rate.
#' @param beta Splicing rate.
#' @param gamma Degradation rate. If \code{NULL}, defaults to \code{beta / 2}.
#' @param alpha_ Repression-phase transcription rate. Default \code{0}.
#' @param u0 Initial U value at induction start. Default \code{0}.
#' @param s0 Initial S value at induction start. Default \code{0}.
#' @param sorted Logical. If \code{TRUE}, return results sorted by time.
#'
#' @return Named list with components \code{tau}, \code{alpha}, \code{u0},
#'   and \code{s0}.
#' @export
vectorize_2state <- function(t, t_, alpha, beta, gamma = NULL,
                             alpha_ = 0, u0 = 0, s0 = 0,
                             sorted = FALSE) {
  if (is.null(gamma)) gamma <- beta / 2

  o <- as.integer(t < t_)

  u0_ <- u_solution(t_, u0, alpha, beta)
  s0_ <- s_solution(t_, s0, u0, alpha, beta, gamma)

  tau <- t * o + (t - t_) * (1 - o)

  u0_vec <- u0 * o + u0_ * (1 - o)
  s0_vec <- s0 * o + s0_ * (1 - o)
  alpha_vec <- alpha * o + alpha_ * (1 - o)

  if (sorted) {
    idx <- order(t)
    tau       <- tau[idx]
    alpha_vec <- alpha_vec[idx]
    u0_vec    <- u0_vec[idx]
    s0_vec    <- s0_vec[idx]
  }

  list(tau = tau, alpha = alpha_vec, u0 = u0_vec, s0 = s0_vec)
}

# -----------------------------------------------------------------------------
# Part 2.4: Time-increment adjustment
# -----------------------------------------------------------------------------

#' Adjust unusually large jumps in a time axis
#'
#' Removes time increments exceeding
#' \code{3 × quantile(0.995)} to avoid discontinuities caused by sparse time
#' assignments.
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::adjust_increments()}.
#'
#' @param tau Within-phase relative-time vector for one transcription phase.
#' @param tau_ Optional within-phase relative-time vector for a second phase.
#'   If \code{NULL}, only \code{tau} is processed.
#'
#' @return If \code{tau_} is \code{NULL}, a corrected \code{tau} vector.
#'   Otherwise, a named list with components \code{tau} and \code{tau_}.
#' @export
adjust_increments <- function(tau, tau_ = NULL) {
  .adjust_single <- function(tv) {
    tv_new <- tv
    tv_sorted <- sort(tv)
    dtv <- diff(c(0, tv_sorted))

    if (!is.null(tau_)) {
      tau_sorted_ <- sort(tau_)
      dtv_ <- diff(c(0, tau_sorted_))
      ub <- 3 * quantile(c(dtv, dtv_), 0.995, na.rm = TRUE)
    } else {
      ub <- 3 * quantile(dtv, 0.995, na.rm = TRUE)
    }

    idx <- which(dtv > ub)
    for (i in idx) {
      ti  <- tv_sorted[i]
      dti <- dtv[i]
      tv_new[tv >= ti] <- tv_new[tv >= ti] - dti
    }
    tv_new
  }

  if (is.null(tau_)) {
    return(.adjust_single(tau))
  }

  tau_new  <- .adjust_single(tau)
  tau_new_ <- .adjust_single(tau_)
  list(tau = tau_new, tau_ = tau_new_)
}


# -----------------------------------------------------------------------------
# Part 2.5: Timeline helper functions
# -----------------------------------------------------------------------------

#' Compute time increments with optional clipping
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::compute_dt()}.
#'
#' @param t Time vector (\eqn{n\_cells}).
#' @param clipped Logical. If \code{TRUE}, clip unusually large increments
#'   using a Poisson-style upper bound. Default \code{TRUE}.
#'
#' @return Numeric vector of time increments.
#' @export
compute_dt <- function(t, clipped = TRUE) {
  t_sorted <- sort(t)
  dt <- diff(c(min(t_sorted), t_sorted))
  m_dt <- max(mean(dt), max(t) / length(t), 0)
  if (clipped) {
    ub <- m_dt + 3 * sqrt(max(m_dt, 0))
    dt <- pmin(pmax(dt, 0), ub)
  }
  dt
}

#' Re-root a timeline at a specified root cell
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::root_time()}.
#'
#' @param t Time matrix of size \eqn{n\_cells \times n\_genes}, or a vector.
#' @param root Integer index of the root cell. If \code{NULL}, uses
#'   \code{t = 0} as the origin.
#'
#' @return Named list with components \code{t_rooted} and \code{t_switch}.
#' @export
root_time <- function(t, root = NULL) {
  t_mat <- if (is.vector(t)) matrix(t, nrow = 1) else t
  nans <- apply(t_mat, 2, function(col) any(is.nan(col) | is.na(col)))
  if (any(nans)) t_mat <- t_mat[, !nans, drop = FALSE]

  t_root <- if (is.null(root)) 0 else t_mat[root, ]
  o <- sweep(t_mat, 2, t_root, ">=") * 1L

  t_after  <- (t_mat - matrix(t_root, nrow(t_mat), ncol(t_mat), byrow = TRUE)) * o
  t_origin <- apply(t_after, 2, max)
  t_before <- (t_mat + matrix(t_origin, nrow(t_mat), ncol(t_mat), byrow = TRUE)) * (1 - o)

  t_switch <- apply(t_before, 2, min)
  t_rooted <- t_after + t_before

  list(t_rooted = t_rooted, t_switch = t_switch)
}

#' Compute a shared time axis across multiple genes
#'
#' Corresponds to scVelo:
#' \code{tools/_em_model_utils.py::compute_shared_time()}.
#'
#' @param t Time matrix of size \eqn{n\_cells \times n\_genes}, or a vector.
#' @param perc Numeric vector of percentile thresholds. Default
#'   \code{c(15, 20, 25, 30, 35)}.
#' @param norm Logical. If \code{TRUE}, normalise the result to \code{[0, 1]}.
#'   Default \code{TRUE}.
#'
#' @return Numeric vector of shared time values.
#' @export
compute_shared_time <- function(t, perc = NULL, norm = TRUE) {
  t_mat <- if (is.vector(t)) matrix(t, nrow = 1) else t

  nans <- apply(t_mat, 2, function(col) any(is.nan(col) | is.na(col)))
  if (any(nans)) t_mat <- t_mat[, !nans, drop = FALSE]

  t_mat <- sweep(t_mat, 2, apply(t_mat, 2, min), "-")

  if (is.null(perc)) perc <- c(15, 20, 25, 30, 35)

  tx_list <- lapply(perc, function(p) {
    tx <- apply(t_mat, 1, quantile, probs = p / 100, na.rm = TRUE)
    tx_max <- max(tx)
    if (tx_max == 0) tx_max <- 1
    tx / tx_max
  })
  tx_mat <- do.call(rbind, tx_list)

  mse_vals <- apply(tx_mat, 1, function(tx) {
    tx_sorted <- sort(tx)
    linx <- seq(0, 1, length.out = length(tx_sorted))
    sum((tx_sorted - linx)^2)
  })

  idx_best <- order(mse_vals)[1:min(2, length(mse_vals))]
  t_shared <- colSums(tx_mat[idx_best, , drop = FALSE])

  if (norm && max(t_shared) > 0) {
    t_shared <- t_shared / max(t_shared)
  }
  t_shared
}
