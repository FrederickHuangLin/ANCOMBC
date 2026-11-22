# =============================================================================
# ANCOM-BC3 second stage and DA inference
# =============================================================================
#
# Second-stage computations of the DA component, in the order of the
# manuscript:
#   Taxon-wise model and penalized likelihood: the penalized zero-inflated
#     Poisson-log-normal log-likelihood l_j^pen of taxon j (cell distribution
#     in ancombc3_pln.R), its Newton-Raphson maximization in omega_j from
#     multiple starting vectors, and the update of kappa_j = log sigma_j.
#   Identifying constraints on the sample-specific bias: the basis Q of the
#     residualized tested covariates, the basis H0 of
#     S = {v : 1'v = 0, Q'v = 0}, and the projection used for the initial bias
#     vector.
#   Joint estimation (fit_ancombc3_da(); the alternation is in
#     ancombc3_joint.R): initial taxon-wise fits, the first pass, pass 2 with
#     the cross-taxon penalty on kappa_j, and verification with the final
#     estimates.
#   Cross-taxon correction (score_ancombc3_da()): the primary estimator
#     (maximal overlap), the secondary estimator (Gaussian mixture), and the
#     selection of the common-shift estimate.
#   DA inference: the covariance of the uncorrected coefficients (the
#     sensitivities g_j and the working covariance V_d^con of the bias
#     estimator, ancombc3_add_bias_covariance()), the working variance of the
#     corrected coefficient, the t reference with the normal-equivalent
#     transformation, the empirical null of the DA statistics, and the
#     P values (score_ancombc3_da()).

.ac3_da_clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# -----------------------------------------------------------------------------
# Taxon-wise model and penalized likelihood
# -----------------------------------------------------------------------------
# For taxon j with regression parameters omega = (alpha, beta, a, xi) and
# sample i, the abundance and presence predictors (log_mu and logit_pi) are
#   log mu_i = alpha + x_i' beta + O^A_i,
#   logit pi_i = a + r_i' xi + O^S_i,
# both truncated to [-30, 30], where the abundance offset O^A_i of a fit
# includes the sample-specific bias d_i and O^S_i is the structural offset.
# Given presence, y_i is Poisson-log-normal with mean mu_i and log-scale
# standard deviation sigma = exp(kappa) (ancombc3_pln.R). The penalized
# log-likelihood is
#   l^pen(omega) = sum_i l_i(omega) - (ridge / 2) (||beta||^2 + ||xi||^2)
#                  + (1 / 2) log det(R' W R),
# where l_i is the cell log-likelihood, R = [1, r], and W is diagonal with
#   w_i = pi_i (1 - pi_i)^2 (1 - f0_i) / {1 - pi_i (1 - f0_i)},
# the Fisher information about logit pi_i carried by the zero indicator. The
# first-stage sampling-zero probability f0_i is fixed, so the Firth-type
# penalty depends on the presence coefficients (a, xi) only.

# Data of one taxon-wise fit: X = [1, x], R = [1, r], the offsets, the
# positions of the abundance block (alpha, beta) and of the presence block
# (a, xi) in omega, the penalized positions (beta, xi) and the positions of
# their diagonal elements in the npar x npar information (pen_diag), one minus
# the fixed sampling-zero probability, sigma, the engine of the evaluation of
# the quadrature (.ac3_pln_obs()), and the positions of the positive and the
# zero counts (index, .ac3_pln_count_index()).
.ac3_da_taxon_data <- function(y, x, r, abundance_offset, structural_offset,
                               ridge, sampling_zero_prob, sigma, engine = "r") {
  y <- as.numeric(y)
  x <- as.matrix(x)
  r <- as.matrix(r)
  n <- length(y)
  qx <- ncol(x)
  qr <- ncol(r)
  pos <- y > 0
  npar <- 2L + qx + qr
  i_pen <- c(if (qx > 0L) 1L + seq_len(qx), if (qr > 0L) 2L + qx + seq_len(qr))
  list(
    y = y, X = cbind(1, x), R = cbind(1, r), n = n, qx = qx, qr = qr,
    npar = npar, abundance_offset = as.numeric(abundance_offset),
    structural_offset = as.numeric(structural_offset), pos = pos,
    index = .ac3_pln_count_index(pos), ridge = ridge,
    one_minus_f0 = 1 - as.numeric(sampling_zero_prob),
    i_abundance = seq_len(1L + qx), i_presence = 1L + qx + seq_len(1L + qr),
    i_pen = i_pen, pen_diag = (i_pen - 1L) * npar + i_pen,
    sigma = as.numeric(sigma), engine = engine
  )
}

# Firth-type penalty (1 / 2) log det M, M = R' W R, its gradient in the
# presence coefficients phi = (a, xi), (1 / 2) sum_i h_i (w'_i / w_i) R_i, with
# h_i = w_i R_i' M^-1 R_i and w'_i = dw_i / d logit pi_i, and, with
# hess = TRUE, its Hessian. With A_k = R' diag(w'_i R_ik) R,
#   d2 / dphi_k dphi_l = (1 / 2) [ -tr(M^-1 A_l M^-1 A_k)
#                                  + sum_i w''_i R_ik R_il R_i' M^-1 R_i ],
# where w' = w dlog_w and w'' = w (dlog_w^2 + d dlog_w / d logit pi); both are
# zero where the presence logit is truncated. presence holds the presence
# quantities at phi (.ac3_pln_presence()); they are computed here when it is
# NULL.
.ac3_da_firth_penalty <- function(phi, dat, grad = TRUE, hess = FALSE, presence = NULL) {
  if (is.null(presence)) presence <- .ac3_pln_presence(as.numeric(dat$R %*% phi) + dat$structural_offset)
  pi_i <- presence$pi_i
  one_minus_pi <- presence$one_minus_pi
  R <- dat$R
  one_minus_f0 <- dat$one_minus_f0
  denom <- 1 - one_minus_f0 * pi_i
  w <- one_minus_f0 * pi_i * one_minus_pi^2 / denom
  M <- crossprod(R * w, R)
  chol_M <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(chol_M)) return(list(value = -Inf, grad = rep(NA_real_, length(phi))))
  val <- sum(log(diag(chol_M)))
  if (!grad) return(list(value = val))
  in_pi <- presence$in_pi
  M_inv <- chol2inv(chol_M)
  R_Minv_R <- rowSums((R %*% M_inv) * R)
  dlog_w <- one_minus_pi - 2 * pi_i + one_minus_f0 * pi_i * one_minus_pi / denom
  out <- list(value = val, grad = 0.5 * as.numeric(crossprod(R, (w * R_Minv_R) * (dlog_w * in_pi))))
  if (!hess) return(out)
  d2log_w <- -3 * pi_i * one_minus_pi + one_minus_f0 * (pi_i * one_minus_pi * (1 - 2 * pi_i) * denom +
                                                        one_minus_f0 * pi_i^2 * one_minus_pi^2) / denom^2
  w_prime <- w * dlog_w * in_pi
  w_second <- w * (dlog_w^2 + d2log_w) * in_pi
  k <- ncol(R)
  # Columns of R, the matrices M^-1 A_k and their transposes, and the products
  # w''_i R_ik.
  R_col <- lapply(seq_len(k), function(j) R[, j])
  Minv_A <- lapply(seq_len(k), function(j) M_inv %*% crossprod(R * (w_prime * R_col[[j]]), R))
  t_Minv_A <- lapply(Minv_A, t)
  ws_R <- lapply(R_col, function(r_j) w_second * r_j)
  H <- matrix(0, k, k)
  for (i in seq_len(k)) for (j in seq_len(i)) {
    H[i, j] <- H[j, i] <- 0.5 * (-sum(Minv_A[[j]] * t_Minv_A[[i]]) +
                                   sum(ws_R[[i]] * R_col[[j]] * R_Minv_R))
  }
  out$hess <- H
  out
}

# Penalized log-likelihood l^pen(omega), its analytic gradient g, and the
# analytic penalized observed information H = -d2 l^pen / d omega2 (cell
# log-likelihood, ridge, and Firth-type terms) at fixed sigma and offsets.
# Returns value (l^pen), grad, hess, and the cell derivatives obs.
.ac3_da_taxon_objective <- function(omega, dat, grad = TRUE, hess = TRUE) {
  log_mu_raw <- as.numeric(dat$X %*% omega[dat$i_abundance]) + dat$abundance_offset
  phi <- omega[dat$i_presence]
  presence <- .ac3_pln_presence(as.numeric(dat$R %*% phi) + dat$structural_offset)
  need <- if (hess) 2L else if (grad) 1L else 0L
  o <- .ac3_pln_observation_derivatives(log_mu_raw, NULL, dat$y, dat$pos, dat$sigma, need = need,
                                         engine = dat$engine, presence = presence, index = dat$index)
  val <- sum(o$ll) - 0.5 * dat$ridge * sum(omega[dat$i_pen]^2)
  fp <- .ac3_da_firth_penalty(phi, dat, grad = grad, hess = grad && hess, presence = presence)
  val <- val + fp$value
  if (!grad || !is.finite(val)) return(list(value = val, obs = o))
  g <- c(as.numeric(crossprod(dat$X, o$dl_mu)), as.numeric(crossprod(dat$R, o$dl_pi)))
  g[dat$i_pen] <- g[dat$i_pen] - dat$ridge * omega[dat$i_pen]
  g[dat$i_presence] <- g[dat$i_presence] + fp$grad
  if (!hess) return(list(value = val, grad = g, obs = o))
  H_aa <- crossprod(dat$X * o$d2l_mumu, dat$X)
  H_ap <- crossprod(dat$X * o$d2l_mupi, dat$R)
  H_pp <- crossprod(dat$R * o$d2l_pipi, dat$R)
  H <- -rbind(cbind(H_aa, H_ap), cbind(t(H_ap), H_pp))
  H[dat$pen_diag] <- H[dat$pen_diag] + dat$ridge
  H[dat$i_presence, dat$i_presence] <- H[dat$i_presence, dat$i_presence] - fp$hess
  list(value = val, grad = g, hess = 0.5 * (H + t(H)), obs = o)
}

# Newton-Raphson maximization of l^pen in omega at fixed kappa and offsets.
# From the current value with gradient g and information H, the step is
# omega + t H_tau^-1 g with H_tau = H + tau I: tau = 0 when H is positive
# definite, and otherwise tau starts at 1e-6 s, s = max(1, max_k |H_kk|), and
# is multiplied by 10 until H_tau is positive definite: the trial values are
# 1e-6 s, 1e-5 s, ..., 1e8 s (fifteen trials, the last set to 1e8 s exactly),
# and the iteration stops without convergence when none of them gives a
# Cholesky factor of H_tau. The step length t starts at 1 and is halved, at
# most 40 times, until the objective increases, or decreases by at most
# 1e-12 (1 + |l|) while max |g| decreases; the iteration stops without
# convergence when no step length is accepted. Converged when
# max |g| <= tol_score or, with H positive definite, when g' H^-1 g <= tol_dec;
# at most max_steps steps. A converged run is regular (pd = TRUE) when H is
# positive definite at its final value. Each trial value is evaluated with the
# observed information, so that an accepted value is not evaluated again. The
# result holds the Cholesky factor chol of H at the final value (NULL when H is
# not positive definite).
.ac3_da_newton <- function(dat, init, max_steps = 100L, tol_dec = 1e-12, tol_score = 1e-8) {
  omega <- init
  ev <- .ac3_da_taxon_objective(omega, dat)
  if (!is.finite(ev$value) || any(!is.finite(ev$grad))) {
    return(list(converged = FALSE, omega = omega, ev = ev, pd = FALSE, chol = NULL))
  }
  converged <- FALSE
  dec <- NA_real_
  pd <- FALSE
  steps <- 0L
  repeat {
    g <- ev$grad
    H <- ev$hess
    chol_H <- if (all(is.finite(H))) tryCatch(chol(H), error = function(e) NULL) else NULL
    pd <- !is.null(chol_H)
    if (pd) {
      direction <- backsolve(chol_H, forwardsolve(t(chol_H), g))
      dec <- sum(g * direction)
    } else {
      dec <- NA_real_
    }
    if (max(abs(g)) <= tol_score || (pd && is.finite(dec) && dec <= tol_dec)) {
      converged <- TRUE
      break
    }
    if (steps >= max_steps) break
    if (!pd) {
      damping_scale <- max(1, max(abs(diag(H)), na.rm = TRUE))
      tau <- 1e-6 * damping_scale
      for (trial in 1:15) {
        chol_H_tau <- tryCatch(chol(H + diag(tau, nrow(H))), error = function(e) NULL)
        if (!is.null(chol_H_tau)) break
        tau <- if (trial < 14L) tau * 10 else 1e8 * damping_scale
      }
      if (is.null(chol_H_tau)) break
      direction <- backsolve(chol_H_tau, forwardsolve(t(chol_H_tau), g))
    }
    step_length <- 1
    accepted <- FALSE
    for (k in 0:40) {
      omega_new <- omega + step_length * direction
      ev_new <- .ac3_da_taxon_objective(omega_new, dat, grad = TRUE, hess = TRUE)
      if (is.finite(ev_new$value) && all(is.finite(ev_new$grad)) &&
          (ev_new$value > ev$value ||
             (ev_new$value >= ev$value - 1e-12 * (1 + abs(ev$value)) &&
                max(abs(ev_new$grad)) < max(abs(g))))) {
        accepted <- TRUE
        break
      }
      step_length <- step_length / 2
    }
    if (!accepted) break
    omega <- omega_new
    ev <- ev_new
    steps <- steps + 1L
  }
  list(converged = converged, omega = omega, ev = ev, pd = pd, chol = chol_H)
}

# Penalized log-likelihood of a regular converged Newton-Raphson result, and
# -Inf otherwise.
.ac3_da_newton_value <- function(res) {
  if (isTRUE(res$converged) && isTRUE(res$pd) && is.finite(res$ev$value)) res$ev$value else -Inf
}

# Default starting vector: zero slopes, the abundance intercept equal to the
# mean log positive count net of the abundance offset, and the presence
# intercept that treats half of the observed zeros as structural (logit
# argument restricted to [0.05, 1 - 1e-4]).
.ac3_da_default_start <- function(y, abundance_offset, structural_offset, qx, qr) {
  pos <- y > 0
  log_pos <- log(pmax(y[pos], 1))
  alpha0 <- if (any(pos)) mean(log_pos - abundance_offset[pos]) else log(0.1) - mean(abundance_offset)
  zero_rate <- mean(y == 0)
  a0 <- qlogis(pmin(pmax(1 - zero_rate * 0.5, 0.05), 1 - 1e-4)) - mean(structural_offset)
  c(alpha0, rep(0, qx), a0, rep(0, qr))
}

# First-stage starting vector: the presence coefficients (a, xi) are the
# quasi-binomial regression coefficients of the first-stage at-risk
# probabilities (restricted to [1e-6, 1 - 1e-6]) on R = [1, r] with the
# structural offset (at most 25 iterations), the abundance slopes are zero, and
# the abundance intercept solves sum_i pi_i exp(alpha + O^A_i) = sum_i y_i.
.ac3_da_first_stage_start <- function(y, x, r, abundance_offset, structural_offset, at_risk_prob) {
  pf <- pmin(pmax(as.numeric(at_risk_prob), 1e-6), 1 - 1e-6)
  R <- cbind(1, r)
  gf <- tryCatch(suppressWarnings(stats::glm.fit(R, pf, family = stats::quasibinomial(),
                                                 offset = structural_offset,
                                                 control = stats::glm.control(maxit = 25))),
                 error = function(e) NULL)
  coef_presence <- if (!is.null(gf) && all(is.finite(gf$coefficients))) {
    unname(gf$coefficients)
  } else {
    c(stats::qlogis(mean(pf)) - mean(structural_offset), rep(0, ncol(r)))
  }
  pi_i <- stats::plogis(pmin(pmax(as.numeric(R %*% coef_presence) + structural_offset, -30), 30))
  alpha <- log(max(sum(y), 0.5)) - log(sum(pi_i * exp(abundance_offset)))
  c(alpha, rep(0, ncol(x)), coef_presence)
}

# Random components of the starting vectors: independent N(0, 1) vectors
# drawn once from a fixed seed, so that every taxon-wise fit uses the same
# starting vectors. The caller's random-number state is restored.
.ac3_da_start_perturbations <- function(npar, n_starts) {
  k <- max(0L, as.integer(n_starts) - 1L)
  if (k == 0L) return(list())
  withr::local_seed(20260924L, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                     .rng_sample_kind = "Rejection")
  lapply(seq_len(k), function(m) stats::rnorm(npar))
}

# Starting vectors of a multiple-start fit, in order: the current estimate
# (warm, when valid), the default starting vector, the first-stage starting
# vector, and the default starting vector plus the random vectors of
# .ac3_da_start_perturbations(), up to n_starts in total.
.ac3_da_taxon_starts <- function(y, x, r, abundance_offset, structural_offset, at_risk_prob,
                                 warm, n_starts, perturbations) {
  npar <- 2L + ncol(x) + ncol(r)
  init0 <- .ac3_da_default_start(y, abundance_offset, structural_offset, ncol(x), ncol(r))
  starts <- list()
  if (!is.null(warm) && length(warm) == npar && all(is.finite(warm))) starts$warm <- warm
  starts$default <- init0
  starts$first_stage <- .ac3_da_first_stage_start(y, x, r, abundance_offset, structural_offset,
                                                  at_risk_prob)
  for (m in seq_len(max(0L, n_starts - length(starts)))) {
    starts[[paste0("perturbed_", m)]] <- init0 + perturbations[[m]]
  }
  starts[seq_len(min(length(starts), n_starts))]
}

# Taxon-wise fit record: ok (the run converged with finite estimates), omega
# (with the attributes sigma and kappa_info), the penalized observed
# information H (hessian) and its inverse (vcov, NULL when H is not positive
# definite), the adjusted objective l^A = l^pen - (1 / 2) log det H at the
# optimum (.ac3_pln_adjusted_from_ev(); l^pen itself when H is not positive
# definite), the coefficients alpha, beta, a, and xi, the fitted mu, pi, and
# f0, sigma and kappa_info, the number of positive counts n_positive, and pd
# (H is positive definite). With fitted_f0 = FALSE the element f0 is NULL and is computed
# afterwards (.ac3_da_complete_fit_summaries()); engine selects the evaluation
# of f0 (.ac3_pln_zero_prob()).
.ac3_da_fit_record <- function(res, y, x, r, abundance_offset, structural_offset,
                               sigma = NA_real_, kappa_info = NA_real_, fitted_f0 = TRUE,
                               engine = "r") {
  qx <- ncol(x)
  qr <- ncol(r)
  ok <- isTRUE(res$converged) && all(is.finite(res$omega))
  omega <- res$omega
  vcov <- NULL
  if (ok && isTRUE(res$pd)) {
    chol_H <- if (is.null(res$chol)) tryCatch(chol(res$ev$hess), error = function(e) NULL) else res$chol
    vcov <- if (is.null(chol_H)) NULL else chol2inv(chol_H)
    if (!is.null(vcov) && any(!is.finite(vcov))) vcov <- NULL
  }
  alpha <- omega[1L]
  beta <- if (qx > 0L) omega[1L + seq_len(qx)] else numeric(0)
  a <- omega[1L + qx + 1L]
  xi <- if (qr > 0L) omega[1L + qx + 1L + seq_len(qr)] else numeric(0)
  names(beta) <- colnames(x)
  names(xi) <- colnames(r)
  log_mu <- pmin(pmax(as.numeric(alpha + x %*% beta + abundance_offset), -30), 30)
  logit_pi <- pmin(pmax(as.numeric(a + r %*% xi + structural_offset), -30), 30)
  mu <- exp(log_mu)
  list(
    ok = ok,
    omega = structure(omega, sigma = sigma, kappa_info = kappa_info),
    hessian = if (ok) res$ev$hess else NULL,
    vcov = vcov,
    objective = if (ok && isTRUE(res$pd)) .ac3_pln_adjusted_from_ev(res$ev, res$chol) else res$ev$value,
    alpha = alpha,
    beta = beta,
    a = a,
    xi = xi,
    mu = mu,
    pi = stats::plogis(logit_pi),
    f0 = if (fitted_f0) .ac3_pln_zero_prob(log_mu, sigma, engine = engine) else NULL,
    sigma = sigma,
    kappa_info = kappa_info,
    n_positive = sum(y > 0),
    pd = isTRUE(res$pd)
  )
}

# Completion of a fit record computed without its summaries
# (.ac3_da_taxon_fit() with summaries = FALSE): kappa_info, when kappa was
# updated in the fit (summaries_pending is TRUE), and f0. dat holds the data
# of the taxon-wise fit at the offsets of that fit (.ac3_da_taxon_data(); its
# sigma is not used) and x the abundance design. The completed record equals
# the record of the same fit computed with summaries = TRUE.
.ac3_da_complete_fit_summaries <- function(f, x, dat) {
  if (is.null(f$summaries_pending)) return(f)
  if (isTRUE(f$summaries_pending)) {
    info <- .ac3_pln_kappa_information(dat, as.numeric(f$omega), f$sigma, f$hessian)
    f$kappa_info <- info
    attr(f$omega, "kappa_info") <- info
  }
  log_mu <- pmin(pmax(as.numeric(f$alpha + x %*% f$beta + dat$abundance_offset), -30), 30)
  f$f0 <- .ac3_pln_zero_prob(log_mu, f$sigma, engine = dat$engine)
  f$summaries_pending <- NULL
  f
}

# Quasi-Newton fallback: BFGS minimization of -l^pen with the analytic
# gradient from the default starting vector (relative tolerance 1e-12, at most
# maxit iterations), followed by Newton-Raphson from its solution. NULL when
# BFGS returns no finite solution.
.ac3_da_bfgs_fallback <- function(dat, init, maxit) {
  fn <- function(p) {
    v <- .ac3_da_taxon_objective(p, dat, grad = FALSE, hess = FALSE)$value
    if (is.finite(v)) -v else 1e300
  }
  gr <- function(p) -.ac3_da_taxon_objective(p, dat, grad = TRUE, hess = FALSE)$grad
  ob <- tryCatch(
    stats::optim(init, fn, gr, method = "BFGS",
                 control = list(maxit = as.integer(maxit), reltol = 1e-12)),
    error = function(e) NULL
  )
  if (is.null(ob) || !all(is.finite(ob$par))) return(NULL)
  .ac3_da_newton(dat, ob$par)
}

# Taxon-wise fit of omega_j and kappa_j at fixed offsets.
#   mode "multistart": Newton-Raphson from every starting vector of
#     .ac3_da_taxon_starts(). The fit retained is the regular converged run
#     with the largest l^pen; values within 1e-8 (1 + |l|) of the largest are
#     treated as equal, and the earliest starting vector is chosen among them.
#     When no run is regular but at least one converged, the converged run
#     with the largest l^pen is retained; its information H is not positive
#     definite, so the record has no vcov. When no run converges, the fit
#     uses the quasi-Newton fallback, whose result is the Newton-Raphson run
#     from the quasi-Newton solution.
#   mode "warm": Newton-Raphson from the current estimate init; when that run
#     is not regular and converged, the taxon is fitted in mode "multistart".
# sigma starts at the sigma of the current estimate (an attribute of init),
# else at sigma_init, else at one, restricted to .ac3_pln_sigma_bounds. For a
# regular retained run, including a regular result of the quasi-Newton
# fallback, kappa = log sigma and omega are then updated by
# .ac3_pln_kappa_update() with the Gaussian prior on kappa (prior =
# list(m, tau2); NULL for none), which also gives kappa_info, the profile
# information about kappa. A retained run that is not regular keeps the
# starting sigma and has kappa_info = NA, which excludes the fit from the
# cross-taxon penalty (.ac3_pln_cross_taxon_penalty()).
# With summaries = FALSE, the quantities of the record that the alternation
# does not use, f0 and kappa_info, are not computed: the record has the element
# summaries_pending (TRUE when kappa was updated in the fit, so that kappa_info
# is to be computed; FALSE otherwise), and .ac3_da_complete_fit_summaries()
# completes it. kappa_info is computed in the fit when the record does not
# hold the information H that its computation uses. engine selects the
# evaluation of the quadrature (.ac3_pln_obs()) in every evaluation of the fit.
.ac3_da_taxon_fit <- function(y,
                              x,
                              r,
                              abundance_offset,
                              structural_offset,
                              sampling_zero_prob,
                              at_risk_prob,
                              ridge = 1e-6,
                              maxit = 500L,
                              init = NULL,
                              mode = c("multistart", "warm"),
                              n_starts = 10L,
                              perturbations = NULL,
                              sigma_init = NULL,
                              prior = NULL,
                              summaries = TRUE,
                              engine = "r") {
  mode <- match.arg(mode)
  y <- as.numeric(y)
  x <- as.matrix(x)
  r <- as.matrix(r)
  if (nrow(x) != length(y) || nrow(r) != length(y)) {
    stop("x, r, and y must have the same sample dimension.")
  }
  abundance_offset <- as.numeric(abundance_offset)
  structural_offset <- as.numeric(structural_offset)
  if (is.null(perturbations)) {
    perturbations <- .ac3_da_start_perturbations(2L + ncol(x) + ncol(r), n_starts)
  }
  ## kappa_info is set only by an update of kappa in this fit; it is not
  ## carried over from the current estimate.
  kappa_info <- NA_real_
  sigma_warm <- if (!is.null(init)) attr(init, "sigma") else NULL
  sigma <- if (!is.null(sigma_warm) && is.finite(sigma_warm) && sigma_warm > 0) {
    sigma_warm
  } else if (!is.null(sigma_init) && is.finite(sigma_init) && sigma_init > 0) {
    sigma_init
  } else {
    1
  }
  sigma <- min(max(sigma, .ac3_pln_sigma_bounds[1L]), .ac3_pln_sigma_bounds[2L])
  dat <- .ac3_da_taxon_data(y, x, r, abundance_offset, structural_offset, ridge,
                            sampling_zero_prob, sigma = sigma, engine = engine)
  record <- function(res, sig = sigma, sig_kappa_info = kappa_info) {
    .ac3_da_fit_record(res, y, x, r, abundance_offset, structural_offset,
                       sigma = sig, kappa_info = sig_kappa_info, fitted_f0 = summaries,
                       engine = engine)
  }
  pending <- function(out, res, kappa_updated) {
    if (summaries) return(out)
    if (kappa_updated && is.null(out$hessian)) {
      info <- .ac3_pln_kappa_information(dat, res$omega, out$sigma, res$ev$hess, chol_H = res$chol)
      out$kappa_info <- info
      attr(out$omega, "kappa_info") <- info
      kappa_updated <- FALSE
    }
    out$summaries_pending <- kappa_updated
    out
  }
  warm_ok <- !is.null(init) && length(init) == dat$npar && all(is.finite(init))
  if (mode == "warm" && warm_ok) {
    res <- .ac3_da_newton(dat, as.numeric(init))
    if (is.finite(.ac3_da_newton_value(res))) {
      upd <- .ac3_pln_kappa_update(dat, res, max_cycles = 2L, prior = prior, information = summaries)
      res <- upd$res
      sigma <- upd$sigma
      kappa_info <- upd$kappa_info
      out <- record(res, sigma, kappa_info)
      return(pending(out, res, TRUE))
    }
  }
  starts <- .ac3_da_taxon_starts(y, x, r, abundance_offset, structural_offset,
                                 at_risk_prob, if (warm_ok) as.numeric(init) else NULL,
                                 n_starts, perturbations)
  results <- lapply(starts, function(s) .ac3_da_newton(dat, s))
  vals <- vapply(results, .ac3_da_newton_value, numeric(1))
  method <- "newton"
  if (all(!is.finite(vals))) {
    conv <- vapply(results, function(rr) isTRUE(rr$converged) && is.finite(rr$ev$value), logical(1))
    if (any(conv)) {
      vals[conv] <- vapply(results[conv], function(rr) rr$ev$value, numeric(1))
      method <- "newton_no_pd"
    }
  }
  if (all(!is.finite(vals))) {
    init0 <- .ac3_da_default_start(y, abundance_offset, structural_offset, ncol(x), ncol(r))
    res <- .ac3_da_bfgs_fallback(dat, init0, maxit)
    if (is.null(res)) {
      res <- list(converged = FALSE, omega = init0, ev = .ac3_da_taxon_objective(init0, dat, hess = FALSE),
                  pd = FALSE)
    }
    kappa_updated <- is.finite(.ac3_da_newton_value(res))
    if (kappa_updated) {
      upd <- .ac3_pln_kappa_update(dat, res, prior = prior, information = summaries)
      res <- upd$res
      sigma <- upd$sigma
      kappa_info <- upd$kappa_info
    }
    out <- record(res, sigma, kappa_info)
    return(pending(out, res, kappa_updated))
  }
  top <- max(vals)
  best <- which(vals >= top - 1e-8 * (1 + abs(top)))[1L]
  res_best <- results[[best]]
  if (identical(method, "newton")) {
    upd <- .ac3_pln_kappa_update(dat, res_best, prior = prior, information = summaries)
    res_best <- upd$res
    sigma <- upd$sigma
    kappa_info <- upd$kappa_info
  }
  out <- record(res_best, sigma, kappa_info)
  pending(out, res_best, identical(method, "newton"))
}

# Fit of taxon j from the shared inputs s of fit_ancombc3_da_taxa().
.ac3_da_fit_one_taxon <- function(j, s) {
  .ac3_da_taxon_fit(
    y = s$y[, j],
    x = s$x,
    r = s$r,
    abundance_offset = s$abundance_offset[, j],
    structural_offset = s$structural_offset[, j],
    sampling_zero_prob = s$sampling_zero_prob[, j],
    at_risk_prob = s$at_risk_prob[, j],
    ridge = s$ridge,
    maxit = s$maxit,
    init = s$init_list[[j]],
    mode = s$mode,
    n_starts = s$n_starts,
    perturbations = s$perturbations,
    sigma_init = if (!is.null(s$sigma_init)) s$sigma_init[j] else NULL,
    prior = s$prior,
    summaries = s$summaries,
    engine = s$engine
  )
}

# Worker side of the distributed taxon-wise fits: the shared inputs are sent
# once per worker and stored in the environment .ac3_worker_state of the
# worker.
# .ac3_da_update_shared_inputs() replaces the elements of the stored inputs
# that are named in s.
.ac3_da_set_shared_inputs <- function(s) {
  assign("da_shared_inputs", s, envir = .ac3_worker_state)
  NULL
}
.ac3_da_update_shared_inputs <- function(s) {
  cur <- get("da_shared_inputs", envir = .ac3_worker_state)
  cur[names(s)] <- s
  assign("da_shared_inputs", cur, envir = .ac3_worker_state)
  NULL
}

# Transfer of the shared inputs of fit_ancombc3_da_taxa() to the workers of
# cl. A cluster made by .ac3_da_make_cluster() holds the inputs last sent (the
# environment in its attribute ac3_sent); only the elements that are not
# identical to those are sent. Without that attribute every element is sent.
.ac3_da_send_shared_inputs <- function(cl, shared) {
  sent <- attr(cl, "ac3_sent")
  if (is.null(sent) || is.null(sent$inputs)) {
    parallel::clusterCall(cl, .ac3_da_set_shared_inputs, shared)
  } else {
    same <- vapply(names(shared), function(nm) identical(shared[[nm]], sent$inputs[[nm]]), logical(1))
    parallel::clusterCall(cl, .ac3_da_update_shared_inputs, shared[!same])
  }
  if (!is.null(sent)) sent$inputs <- shared
  invisible(NULL)
}
.ac3_da_fit_shared_taxon <- function(j) {
  .ac3_da_fit_one_taxon(j, get("da_shared_inputs", envir = .ac3_worker_state))
}

# Distributed taxon-wise fits of taxa 1, ..., p in consecutive chunks of
# max(number of workers, ceiling(p / 20)) taxa, with progress(done, p) called
# before the first chunk and after every chunk. Each fit is that of
# .ac3_da_fit_shared_taxon(), so the fits equal those of a single call of
# parallel::parLapplyLB() over all taxa.
.ac3_da_chunked_fits <- function(cl, p, progress) {
  size <- max(length(cl), ceiling(p / 20))
  chunks <- split(seq_len(p), ceiling(seq_len(p) / size))
  fits <- vector("list", p)
  progress(0L, p)
  for (idx in chunks) {
    fits[idx] <- parallel::parLapplyLB(cl, idx, .ac3_da_fit_shared_taxon, chunk.size = 1L)
    progress(max(idx), p)
  }
  fits
}

# Inputs held by a worker process between calls.
.ac3_worker_state <- new.env(parent = emptyenv())

# Cluster of n_cores worker processes with the namespace of the package
# loaded. One BLAS/OpenMP thread per worker.
.ac3_da_make_cluster <- function(n_cores) {
  cl <- parallel::makeCluster(n_cores)
  ok <- FALSE
  on.exit(if (!ok) parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
    loadNamespace("ANCOMBC")
    NULL
  })
  attr(cl, "ac3_sent") <- new.env(parent = emptyenv())
  ok <- TRUE
  cl
}

# Taxon-wise fits of every taxon at the given offsets (abundance_offset
# includes the sample-specific bias) and the raw coefficient table of the
# tested terms. init_list holds optional current estimates, one per taxon;
# sigma_init holds the starting values of sigma_j for the taxa without a
# current estimate; prior is the Gaussian prior on kappa_j (NULL for none).
# cl is an optional cluster of worker processes (.ac3_da_make_cluster()) over
# which the fits are distributed; each fit is deterministic and independent of
# the others, so the result does not depend on cl. With summaries = FALSE the
# fit records are computed without f0 and kappa_info (.ac3_da_taxon_fit()),
# and the columns of the matrix f0 are NA. The taxa are those of y
# (its column names, or T1, ..., Tp when it has none); the n x p offset and
# probability matrices are matched to them by name when they have column
# names, and by position otherwise. Returns the fit records
# (fits), the fitted mean, sampling-zero probability, and presence probability
# matrices (mu, f0, pi), sigma_hat, and raw, with one row per taxon and tested
# term: the coefficient beta_raw, its model variance var_raw_model (the
# diagonal element of the inverse information at the optimum, with
# se_raw_model), the number of positive counts n_positive, and the number of
# abundance coefficients n_abundance_coef. progress is an optional function called
# as progress(done, total) with the number of completed fits: once with zero
# before the fits, then after every taxon (without cl) or after every chunk of
# taxa (with cl; .ac3_da_chunked_fits()); the fits do not depend on it.
# engine selects the evaluation of the quadrature (.ac3_pln_obs()); it is part
# of the shared inputs sent to the workers.
fit_ancombc3_da_taxa <- function(y,
                                 x,
                                 r,
                                 abundance_offset = NULL,
                                 structural_offset = NULL,
                                 sampling_zero_prob,
                                 at_risk_prob,
                                 test_terms = colnames(x),
                                 ridge = 1e-6,
                                 maxit = 500L,
                                 init_list = NULL,
                                 mode = c("multistart", "warm"),
                                 n_starts = 10L,
                                 perturbations = NULL,
                                 sigma_init = NULL,
                                 prior = NULL,
                                 cl = NULL,
                                 progress = NULL,
                                 summaries = TRUE,
                                 engine = "r") {
  mode <- match.arg(mode)
  y <- as.matrix(y)
  storage.mode(y) <- "double"
  x <- as.matrix(x)
  r <- as.matrix(r)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  if (is.null(colnames(r))) colnames(r) <- paste0("r", seq_len(ncol(r)))
  colnames(y) <- .ac3_taxon_names(y)
  taxa <- colnames(y)
  if (is.null(abundance_offset)) abundance_offset <- matrix(0, nrow(y), ncol(y))
  if (is.null(structural_offset)) structural_offset <- matrix(0, nrow(y), ncol(y))
  abundance_offset <- .ac3_match_taxa(abundance_offset, taxa, "abundance_offset")
  structural_offset <- .ac3_match_taxa(structural_offset, taxa, "structural_offset")
  sampling_zero_prob <- .ac3_match_taxa(sampling_zero_prob, taxa, "sampling_zero_prob")
  at_risk_prob <- .ac3_match_taxa(at_risk_prob, taxa, "at_risk_prob")
  if (!identical(dim(abundance_offset), dim(y))) stop("abundance_offset must be n x p.")
  if (!identical(dim(structural_offset), dim(y))) stop("structural_offset must be n x p.")
  if (!identical(dim(sampling_zero_prob), dim(y))) stop("sampling_zero_prob must be n x p.")
  if (!identical(dim(at_risk_prob), dim(y))) stop("at_risk_prob must be n x p.")
  test_terms <- intersect(test_terms, colnames(x))
  if (length(test_terms) == 0L) stop("No test_terms are present in x.")
  if (is.null(perturbations)) {
    perturbations <- .ac3_da_start_perturbations(2L + ncol(x) + ncol(r), n_starts)
  }

  shared <- list(
    y = y, x = x, r = r, abundance_offset = abundance_offset,
    structural_offset = structural_offset, sampling_zero_prob = sampling_zero_prob,
    at_risk_prob = at_risk_prob, ridge = ridge, maxit = maxit,
    init_list = init_list, mode = mode, n_starts = n_starts, perturbations = perturbations,
    sigma_init = sigma_init, prior = prior, summaries = summaries, engine = engine
  )
  if (!is.null(cl)) {
    .ac3_da_send_shared_inputs(cl, shared)
    if (is.null(progress)) {
      fits <- parallel::parLapplyLB(cl, seq_len(ncol(y)), .ac3_da_fit_shared_taxon, chunk.size = 1L)
    } else {
      fits <- .ac3_da_chunked_fits(cl, ncol(y), progress)
    }
  } else if (is.null(progress)) {
    fits <- lapply(seq_len(ncol(y)), .ac3_da_fit_one_taxon, s = shared)
  } else {
    progress(0L, ncol(y))
    fits <- vector("list", ncol(y))
    for (j in seq_len(ncol(y))) {
      fits[[j]] <- .ac3_da_fit_one_taxon(j, s = shared)
      progress(j, ncol(y))
    }
  }
  names(fits) <- colnames(y)
  mu_hat <- f0_hat <- pi_hat <- matrix(NA_real_, nrow(y), ncol(y), dimnames = dimnames(y))
  # Columns of the raw coefficient table, one row per taxon and tested term in
  # the order of the taxa.
  n_rows <- ncol(y) * length(test_terms)
  raw_taxon <- raw_term <- character(n_rows)
  raw_beta <- raw_se <- raw_var <- rep(NA_real_, n_rows)
  raw_positive <- integer(n_rows)
  row <- 0L

  for (j in seq_len(ncol(y))) {
    fit <- fits[[j]]
    if (isTRUE(fit$ok)) {
      mu_hat[, j] <- fit$mu
      if (!is.null(fit$f0)) f0_hat[, j] <- fit$f0
      pi_hat[, j] <- fit$pi
    }
    for (term in test_terms) {
      idx_beta <- 1L + match(term, colnames(x))
      beta_raw <- var_raw <- se_raw <- NA_real_
      if (isTRUE(fit$ok)) {
        beta_raw <- unname(fit$beta[term])
        if (!is.null(fit$vcov) && idx_beta <= nrow(fit$vcov) &&
            is.finite(fit$vcov[idx_beta, idx_beta]) &&
            fit$vcov[idx_beta, idx_beta] > 0) {
          var_raw <- fit$vcov[idx_beta, idx_beta]
          se_raw <- sqrt(var_raw)
        }
      }
      row <- row + 1L
      raw_taxon[row] <- colnames(y)[j]
      raw_term[row] <- term
      raw_beta[row] <- beta_raw
      raw_se[row] <- se_raw
      raw_var[row] <- var_raw
      raw_positive[row] <- if (!is.null(fit$n_positive)) fit$n_positive else sum(y[, j] > 0)
    }
  }

  # Data frame with automatic row names; its attributes are in the order names,
  # row.names, class.
  raw <- structure(
    list(
      taxon = raw_taxon,
      term = raw_term,
      beta_raw = raw_beta,
      se_raw_model = raw_se,
      var_raw_model = raw_var,
      n_positive = raw_positive,
      n_abundance_coef = rep(1L + ncol(x), n_rows)
    ),
    row.names = .set_row_names(n_rows),
    class = "data.frame"
  )

  list(
    raw = raw,
    fits = fits,
    mu = mu_hat,
    f0 = f0_hat,
    pi = pi_hat,
    sigma_hat = vapply(fits, function(f) if (!is.null(f$sigma)) as.numeric(f$sigma) else NA_real_,
                       numeric(1))
  )
}

# -----------------------------------------------------------------------------
# Identifying constraints on the sample-specific bias
# -----------------------------------------------------------------------------
# The bias d is restricted to S = {v : 1'v = 0, Q'v = 0}, where the columns of
# Q are an orthonormal basis of the span of the residualized tested covariates
# u_T.

# Orthonormal basis Q of the span of the residualized tested covariates. For
# each tested column, u_T is its residual on the intercept and the other
# columns of x, centered; a residual with sum of squares at most the machine
# precision is dropped. Returns an n x 0 matrix when no tested column remains.
.ac3_da_residualized_tested_basis <- function(x, test_terms) {
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  test_terms <- intersect(test_terms, colnames(x))
  if (length(test_terms) == 0L) return(matrix(0, nrow(x), 0L))
  u <- lapply(test_terms, function(term) {
    nuisance <- setdiff(colnames(x), term)
    z <- cbind("(Intercept)" = 1, x[, nuisance, drop = FALSE])
    u_T <- tryCatch(lm.fit(z, x[, term])$residuals, error = function(e) NULL)
    if (is.null(u_T) || !all(is.finite(u_T))) return(NULL)
    u_T <- as.numeric(u_T)
    u_T <- u_T - mean(u_T)
    if (sum(u_T^2) <= .Machine$double.eps) return(NULL)
    u_T
  })
  u <- u[!vapply(u, is.null, logical(1))]
  if (length(u) == 0L) return(matrix(0, nrow(x), 0L))
  U <- do.call(cbind, u)
  colnames(U) <- test_terms[seq_len(ncol(U))]
  qr_U <- qr(U)
  if (qr_U$rank == 0L) return(matrix(0, nrow(x), 0L))
  qr.Q(qr_U)[, seq_len(qr_U$rank), drop = FALSE]
}

# Orthonormal basis H0 of S, the orthogonal complement of span(1, Q), from a
# complete QR decomposition.
.ac3_da_constraint_basis <- function(x, test_terms) {
  Q <- .ac3_da_residualized_tested_basis(x, test_terms)
  n <- nrow(x)
  z <- cbind(rep(1 / sqrt(n), n), Q)
  qz <- qr(z)
  qr.Q(qz, complete = TRUE)[, -seq_len(qz$rank), drop = FALSE]
}

# Projection P of the initialization: the projection of d onto S followed,
# when d_bound is finite and positive, by three cycles of clipping to
# [-d_bound, d_bound] and projection onto S. The final projection enforces the
# identifying constraints exactly and the bound approximately.
.ac3_da_project_bias <- function(d,
                                 x,
                                 test_terms,
                                 d_bound = Inf) {
  n_cycles <- 3L
  d_out <- as.numeric(d)
  Q <- .ac3_da_residualized_tested_basis(x, test_terms)
  project_once <- function(v) {
    v <- as.numeric(v)
    v <- v - mean(v)
    if (ncol(Q) > 0L) {
      v <- v - as.numeric(Q %*% crossprod(Q, v))
    }
    v - mean(v)
  }
  d_out <- project_once(d_out)
  if (is.finite(d_bound) && d_bound > 0) {
    for (cycle in seq_len(n_cycles)) {
      d_out <- pmin(pmax(d_out, -d_bound), d_bound)
      d_out <- project_once(d_out)
    }
  }
  d_out
}

# -----------------------------------------------------------------------------
# Joint estimation of the taxon-wise parameters and the sample-specific bias
# -----------------------------------------------------------------------------

# Taxon-wise fits of every taxon at the bias d in the given mode
# (.ac3_da_taxon_fit()). Without previous fits, sigma starts at
# inp$sigma_start. With previous fits, a taxon whose previous fit succeeded
# starts from that estimate and its sigma, and a taxon whose previous fit
# failed is fitted as without previous fits. In every fit, kappa is updated
# with the Gaussian prior (NULL for none). summaries = FALSE gives fit records
# without f0 and kappa_info (.ac3_da_taxon_fit()).
.ac3_da_refit_taxa <- function(inp, d, prev_fits = NULL, mode = "multistart", prior = NULL,
                               summaries = TRUE) {
  init <- NULL
  if (!is.null(prev_fits)) {
    init <- lapply(seq_along(prev_fits), function(j) {
      f <- prev_fits[[j]]
      if (isTRUE(f$ok)) f$omega else NULL
    })
  }
  fit_ancombc3_da_taxa(
    y = inp$y,
    x = inp$x,
    r = inp$r,
    abundance_offset = inp$abundance_offset + matrix(d, nrow(inp$y), ncol(inp$y)),
    structural_offset = inp$structural_offset,
    sampling_zero_prob = inp$sampling_zero_prob,
    at_risk_prob = inp$at_risk_prob,
    test_terms = inp$test_terms,
    ridge = inp$ridge,
    maxit = inp$maxit,
    init_list = init,
    mode = mode,
    n_starts = inp$n_starts,
    perturbations = inp$perturbations,
    sigma_init = inp$sigma_start,
    prior = prior,
    cl = inp$cl,
    progress = inp$progress,
    summaries = summaries,
    engine = inp$engine
  )
}

# Replaces the fits of taxa idx, with their raw rows and fitted-matrix columns,
# by those of the fit ver (a verification pass, or the previous fit when a refit
# fails).
.ac3_da_replace_taxon_fits <- function(fit, ver, idx) {
  taxa <- names(fit$fits)[idx]
  for (j in idx) fit$fits[[j]] <- ver$fits[[j]]
  keep <- !(fit$raw$taxon %in% taxa)
  raw <- rbind(fit$raw[keep, , drop = FALSE], ver$raw[ver$raw$taxon %in% taxa, , drop = FALSE])
  ord <- order(match(raw$taxon, names(fit$fits)), match(raw$term, unique(fit$raw$term)))
  fit$raw <- raw[ord, , drop = FALSE]
  rownames(fit$raw) <- NULL
  for (nm in c("mu", "f0", "pi")) fit[[nm]][, idx] <- ver[[nm]][, idx]
  fit
}

# Second-stage fit of the taxon-wise parameters psi_j = (omega_j, kappa_j) and
# the sample-specific bias d.
#   y, x, r: counts, abundance design, and prevalence design. The taxa are
#     the column names of y, or T1, ..., Tp when y has none.
#   sigma_vae: first-stage log-scale standard deviations sigma_j^VAE, the
#     starting values of sigma_j (restricted to [0.05, 3]).
#   The per-taxon inputs (sigma_vae and the n x p matrices below) are matched
#   to the taxa of y by name when they carry taxon names (as the outputs of
#   ancombc3_fitted_vae_offsets() do), and by position otherwise.
#   abundance_offset, structural_offset: n x p first-stage offsets O^A and O^S
#     (structural_offset NULL for a zero offset).
#   sampling_zero_prob: n x p first-stage sampling-zero probabilities f0^full,
#     the fixed weights of the Firth-type penalty.
#   at_risk_prob: n x p first-stage at-risk probabilities pi^full, used for the
#     first-stage starting vector of each taxon-wise fit and, by default, for
#     the DA eligibility.
#   test_terms: tested columns of x.
#   fitting_set: the taxa B^(d) used to estimate d (taxon names, or one logical
#     value per taxon). NULL gives the union of eligible_taxa_by_term.
#   eligible_taxa_by_term: DA-eligible taxa of each tested term
#     (ancombc3_da_eligibility()$by_term), stored in the result. When both
#     fitting_set and eligible_taxa_by_term are NULL, eligible_taxa_by_term is
#     computed by ancombc3_da_eligibility() from at_risk_prob, y, and x at its
#     default thresholds, those of the manuscript, with the type of each tested
#     column inferred from its values.
#   d_tol, d_max_iter: stopping rule of each update of d
#     (.ac3_joint_d_update()): the Newton-Raphson steps stop when
#     max_i |Delta_i| / se_i <= d_tol, or after d_max_iter steps (defaults 0.01
#     and 50).
#   n_starts: number of starting vectors of a multiple-start fit.
#   maxit_taxon: iteration limit of the quasi-Newton fallback of a taxon-wise
#     fit.
#   d_bound: numerical bound on max_i |d_i| (default 10), enforced by the step
#     fraction and the clipping and projection of the update of d.
#   n_cores: number of worker processes over which the taxon-wise fits are
#     distributed (default 1, no workers); the result does not depend on it.
#   cluster: with repeated measures, the subject of each sample (NULL for
#     independent samples). The estimates are unchanged. The covariances of
#     the taxon-wise estimates and of the bias estimator are replaced by their
#     cluster-robust versions, the variances of raw are recomputed from them,
#     and raw receives the column df_count (.ac3_cluster_da_covariance(),
#     ancombc3_cluster.R); the default eligibility is computed with cluster.
#   progress: optional function passed to every pass of taxon-wise fits
#     (fit_ancombc3_da_taxa()); the result does not depend on it.
#   engine: evaluation of the Poisson-log-normal quadrature in every
#     evaluation of the second stage: "r" (default) for the R code, "compiled"
#     for the compiled routine (ancombc3_engine.R), which is then loaded, and
#     built when it is not in the cache, in this process and in every worker
#     process.
# Estimation (.ac3_joint_estimate(), ancombc3_joint.R):
#   1. Initial taxon-wise fits: every taxon is fitted from the full set of
#      starting vectors at the two initial bias vectors, zero and the
#      centered log library sizes after projection and clipping
#      (.ac3_da_project_bias()). The fitting set B consists of the taxa of
#      fitting_set whose initial fits succeed at both initial bias vectors, so
#      that both first-pass runs maximize L over the same taxa.
#   2. First pass: from each initial fit, the taxon-wise parameters and d are
#      updated in turn without the penalty on kappa_j; the run with the larger
#      L is retained (ties go to d = 0).
#   3. Pass 2: the alternation continues from the retained first-pass fit with
#      the cross-taxon penalty held fixed.
#   4. Verification and final estimates: every taxon is refitted at the
#      current d from the full set of starting vectors; a refit replaces the
#      current fit when it increases the taxon's objective, and pass 2 resumes
#      when a replaced taxon belongs to B.
# The bias covariance (ancombc3_add_bias_covariance()) is computed at the final
# fits. Besides raw, fits, mu, f0, pi, sigma_hat (the final sigma_j), and
# da_variance_inputs (the input of score_ancombc3_da()), the result holds
# d_hat; fitting_set, the set B; eligible_taxa_by_term (supplied or computed;
# NULL when only fitting_set is supplied); cross_taxon_penalty, the penalty
# list(m, tau2, n_taxa, mad, sampling_sd) (absent when pass 2 is omitted);
# objective_L, L at the solution; start_selected, the retained initial bias
# vector; pass1_iterations; pass2_iterations, the outer iterations of the
# final run of the alternation (pass 2 with its resumptions after
# verification, or the retained first-pass run with its resumptions when
# pass 2 is omitted); converged and alternation_history, the convergence and
# the per-iteration record of that run; verification_passes; and
# verification_replaced, the number of replaced taxon-wise fits.
fit_ancombc3_da <- function(y,
                            x,
                            r,
                            sigma_vae,
                            abundance_offset,
                            structural_offset = NULL,
                            sampling_zero_prob,
                            at_risk_prob,
                            test_terms = colnames(x),
                            fitting_set = NULL,
                            eligible_taxa_by_term = NULL,
                            ridge = 1e-6,
                            maxit_taxon = 500L,
                            d_tol = 0.01,
                            d_max_iter = 50L,
                            n_starts = 10L,
                            d_bound = 10,
                            n_cores = 1L,
                            cluster = NULL,
                            progress = NULL,
                            engine = "r") {
  engine <- .ac3_engine_check(engine)
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) {
    stop("n_cores must be a positive integer.")
  }
  y <- as.matrix(y)
  x <- as.matrix(x)
  r <- as.matrix(r)
  colnames(y) <- .ac3_taxon_names(y)
  taxa <- colnames(y)
  n <- nrow(y)
  p <- ncol(y)
  sigma_vae <- .ac3_match_taxa_vector(sigma_vae, taxa, "sigma_vae")
  if (any(!is.finite(sigma_vae)) || any(sigma_vae <= 0)) {
    stop("sigma_vae must have one positive finite value per taxon.")
  }
  if (is.null(structural_offset)) structural_offset <- matrix(0, n, p)
  abundance_offset <- .ac3_match_taxa(abundance_offset, taxa, "abundance_offset")
  structural_offset <- .ac3_match_taxa(structural_offset, taxa, "structural_offset")
  sampling_zero_prob <- .ac3_match_taxa(sampling_zero_prob, taxa, "sampling_zero_prob")
  at_risk_prob <- .ac3_match_taxa(at_risk_prob, taxa, "at_risk_prob")
  if (!identical(dim(abundance_offset), dim(y))) {
    stop("abundance_offset must be n x p.")
  }
  if (!identical(dim(structural_offset), dim(y))) {
    stop("structural_offset must be n x p.")
  }
  if (!identical(dim(sampling_zero_prob), dim(y))) stop("sampling_zero_prob must be n x p.")
  if (!identical(dim(at_risk_prob), dim(y))) stop("at_risk_prob must be n x p.")
  if (any(!is.finite(sampling_zero_prob)) || any(sampling_zero_prob < 0 | sampling_zero_prob > 1)) {
    stop("sampling_zero_prob must contain probabilities.")
  }
  if (any(!is.finite(at_risk_prob)) ||
      any(at_risk_prob < 0 | at_risk_prob > 1)) {
    stop("at_risk_prob must contain probabilities.")
  }
  if (length(d_tol) != 1L || !is.finite(d_tol) || d_tol <= 0) stop("d_tol must be positive.")
  d_max_iter <- as.integer(d_max_iter)
  if (length(d_max_iter) != 1L || is.na(d_max_iter) || d_max_iter < 0L) {
    stop("d_max_iter must be a non-negative integer.")
  }
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts < 1L) {
    stop("n_starts must be a positive integer.")
  }
  if (is.null(fitting_set) && is.null(eligible_taxa_by_term)) {
    ## DA eligibility at the thresholds of the manuscript.
    eligible_taxa_by_term <- ancombc3_da_eligibility(at_risk_prob = at_risk_prob, y = y, x = x,
                                                     test_terms = test_terms, cluster = cluster)$by_term
  }
  if (is.null(fitting_set) && !is.null(eligible_taxa_by_term)) {
    fitting_set <- unique(unlist(eligible_taxa_by_term, use.names = FALSE))
  }
  fitting_taxa <- if (is.null(fitting_set)) {
    taxa
  } else if (is.character(fitting_set)) {
    taxa[taxa %in% fitting_set]
  } else {
    in_set <- as.logical(fitting_set)
    if (length(in_set) != p) stop("fitting_set must identify one value per taxon.")
    taxa[in_set %in% TRUE]
  }
  inp <- list(
    y = y, x = x, r = r, abundance_offset = abundance_offset,
    structural_offset = structural_offset, sampling_zero_prob = sampling_zero_prob,
    at_risk_prob = at_risk_prob, ridge = ridge, maxit = maxit_taxon,
    test_terms = test_terms, d_bound = d_bound, d_tol = d_tol, d_max_iter = d_max_iter,
    fitting_set = fitting_taxa, n_starts = n_starts,
    perturbations = .ac3_da_start_perturbations(2L + ncol(x) + ncol(r), n_starts),
    sigma_start = .ac3_pln_sigma_start(sigma_vae),
    cl = NULL,
    progress = progress,
    engine = engine
  )
  if (n_cores > 1L) {
    inp$cl <- .ac3_da_make_cluster(n_cores)
    on.exit(parallel::stopCluster(inp$cl), add = TRUE)
    parallel::clusterCall(inp$cl, .ac3_joint_set_worker_inputs,
                          inp[c("y", "x", "r", "abundance_offset", "structural_offset",
                                "sampling_zero_prob", "ridge", "engine")])
  }
  log_lib <- log(pmax(rowSums(y), 1))
  d_library <- .ac3_da_project_bias(log_lib - mean(log_lib), x = x,
                                    test_terms = test_terms, d_bound = d_bound)
  starts <- list(zero = numeric(n), library_size = d_library)
  # Initial taxon-wise fits at both initial bias vectors. The fitting set B
  # holds the taxa of fitting_set whose initial fits succeed at both.
  fit0 <- lapply(starts, function(d0) .ac3_da_refit_taxa(inp, d0, mode = "multistart"))
  ok_both <- Reduce(`&`, lapply(fit0, function(f) {
    vapply(f$fits, function(ff) isTRUE(ff$ok), logical(1))
  }))
  fitting_idx <- which(taxa %in% fitting_taxa & ok_both)
  if (length(fitting_idx) == 0L) {
    stop("No taxon of the fitting set has a successful initial fit at both initial bias vectors.")
  }
  je <- .ac3_joint_estimate(inp, starts, fit0, fitting_idx, x, test_terms)
  final_fit <- ancombc3_add_bias_covariance(
    fit = je$best$fit, y = y, x = x, r = r, test_terms = test_terms, engine = engine
  )
  if (!is.null(cluster)) {
    final_fit <- .ac3_cluster_da_covariance(final_fit, inp = inp, d = je$best$d, cluster = cluster)
  }
  final_fit$sigma_hat <- vapply(final_fit$fits, function(f) {
    if (!is.null(f$sigma)) as.numeric(f$sigma) else NA_real_
  }, numeric(1))
  final_fit$d_hat <- je$best$d
  final_fit$fitting_set <- taxa[fitting_idx]
  final_fit$eligible_taxa_by_term <- eligible_taxa_by_term
  final_fit$cross_taxon_penalty <- je$prior
  final_fit$objective_L <- je$best$L
  final_fit$start_selected <- je$chosen
  final_fit$pass1_iterations <- je$best$pass1_iterations %||% je$best$iterations
  final_fit$pass2_iterations <- je$best$iterations
  final_fit$converged <- je$best$converged
  final_fit$alternation_history <- je$best$history
  final_fit$verification_passes <- je$best$verifications
  final_fit$verification_replaced <- je$best$n_replaced
  final_fit
}

# -----------------------------------------------------------------------------
# Cross-taxon correction: primary estimator, maximal overlap
# -----------------------------------------------------------------------------
# For the eligible coefficients beta_j with a finite working standard error
# s_j > 0 (the set E_T), each coefficient defines the interval
# [beta_j - z s_j, beta_j + z s_j], and N_z(c) counts the intervals that contain
# c. N_z is computed exactly by sorting the interval endpoints, with left
# endpoints before right endpoints at equal values, so that the intervals are
# closed. The maximal-overlap estimate is computed in four steps:
#   1. bandwidth: the smallest z of bandwidths with max_c N_z(c) >= min_count,
#      or the largest z of bandwidths when no value meets the criterion
#      (status "bandwidth_criterion_not_met");
#   2. center: the midpoint c_T of the leftmost interval of locations that
#      attain the maximal overlap at the selected bandwidth;
#   3. estimated null set A^MO: the taxa whose intervals contain c_T;
#   4. delta_mo: the inverse-variance-weighted mean of the coefficients of
#      A^MO, with s_j^2 bounded below by the machine precision.
# Degenerate cases, in which every correction weight is zero: fewer than five
# coefficients in E_T, delta_mo = 0 (status "insufficient_taxa"); an empty
# null set in floating-point arithmetic, delta_mo = c_T
# ("empty_null_set_center_used"); a nonfinite sum of weights, delta_mo = 0
# ("nonfinite_weights"). Returns delta_mo, var_delta_mo (the inverse of the
# sum of weights, the variance of delta_mo with the null set held fixed),
# n_null_set, bandwidth, center, overlap_count (the maximal overlap at the
# selected bandwidth), status, null_set_pos (the positions of the null set
# among the input coefficients), and nu (the correction weights of the null
# set, which sum to one).
.ac3_da_maximal_overlap <- function(beta,
                                    se,
                                    min_count = 25L,
                                    bandwidths = seq(0.5, 3.0, by = 0.1)) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  idx_ok <- which(ok)
  beta <- as.numeric(beta[ok])
  se <- as.numeric(se[ok])
  n_eligible <- length(beta)
  empty <- function(status, delta_mo = 0, var_delta_mo = 0, center = NA_real_,
                    bandwidth = NA_real_, overlap_count = 0L, n_null_set = 0L) {
    list(
      delta_mo = delta_mo,
      var_delta_mo = var_delta_mo,
      n_null_set = as.integer(n_null_set),
      bandwidth = bandwidth,
      center = center,
      overlap_count = as.integer(overlap_count),
      status = status,
      null_set_pos = integer(0),
      nu = numeric(0)
    )
  }
  if (n_eligible < 5L) {
    return(empty("insufficient_taxa"))
  }

  overlap_maximum <- function(z) {
    left <- beta - z * se
    right <- beta + z * se
    event_x <- c(left, right)
    event_d <- c(rep(1L, n_eligible), rep(-1L, n_eligible))
    ord <- order(event_x, -event_d)
    run <- cumsum(event_d[ord])
    mx <- max(run)
    first <- which(run == mx)[1L]
    x0 <- event_x[ord][first]
    x1 <- if (first < length(run)) event_x[ord][first + 1L] else x0
    list(center = (x0 + x1) / 2, overlap_count = mx, left = left, right = right)
  }

  bandwidths <- sort(unique(as.numeric(bandwidths[is.finite(bandwidths) & bandwidths > 0])))
  if (length(bandwidths) == 0L) stop("bandwidths must contain a finite, positive value.")
  overlaps <- vapply(bandwidths, function(z) overlap_maximum(z)$overlap_count, numeric(1))
  hit <- which(overlaps >= as.integer(min_count))[1L]
  if (is.na(hit)) {
    bandwidth <- max(bandwidths)
    status <- "bandwidth_criterion_not_met"
  } else {
    bandwidth <- bandwidths[hit]
    status <- "ok"
  }

  mp <- overlap_maximum(bandwidth)
  in_null_set <- mp$left <= mp$center & mp$right >= mp$center
  if (!any(in_null_set)) {
    return(empty(
      "empty_null_set_center_used",
      delta_mo = mp$center,
      var_delta_mo = 0,
      center = mp$center,
      bandwidth = bandwidth,
      overlap_count = mp$overlap_count
    ))
  }
  w <- 1 / pmax(se[in_null_set]^2, .Machine$double.eps)
  wsum <- sum(w)
  if (!is.finite(wsum) || wsum <= 0) {
    return(empty(
      "nonfinite_weights",
      center = mp$center,
      bandwidth = bandwidth,
      overlap_count = mp$overlap_count,
      n_null_set = sum(in_null_set)
    ))
  }
  delta_mo <- sum(beta[in_null_set] * w) / wsum
  list(
    delta_mo = delta_mo,
    var_delta_mo = 1 / wsum,
    n_null_set = as.integer(sum(in_null_set)),
    bandwidth = bandwidth,
    center = mp$center,
    overlap_count = as.integer(mp$overlap_count),
    status = status,
    null_set_pos = idx_ok[in_null_set],
    ## Correction weights nu of the null set: delta_mo = sum(nu * beta).
    nu = w / wsum
  )
}

# -----------------------------------------------------------------------------
# Cross-taxon correction: secondary estimator, Gaussian mixture
# -----------------------------------------------------------------------------
# The coefficients of E_T (v_j = s_j^2) are modeled by the mixture
#   beta_j ~ pi_0 N(mu_0, v_j) + sum_k pi_k N(mu_k, v_j + tau2_k),  k = 1, ..., K,
# with K = 1 or 2. The null component has the working variance v_j only, and
# each further component has an excess variance tau2_k for the variation of
# non-null effects. With p_null_j the fitted posterior probability that taxon
# j belongs to the null component, the mixture estimate of the common shift is
#   delta_gm = sum_j nu_gm_j beta_j,   nu_gm_j = (p_null_j / v_j) / sum_l (p_null_l / v_l).

# EM iterations for the mixture above from mu0, mu, tau2, and
# pi_all = (pi_0, pi_1, ..., pi_K); mu0 is held fixed when fix_mu0 is TRUE.
# Each iteration updates pi_k = max(rbar_k, 1e-8) / sum_k' max(rbar_k', 1e-8),
# mu0 as the mean of b weighted by r_j0 / v_j, mu_k as the mean of b weighted
# by r_jk / (v_j + tau2_k), and tau2_k by three fixed-point steps, bounded
# below by tau2_floor. A mixture density that evaluates to zero is replaced by
# the smallest positive floating-point number. The iteration stops when the
# log-likelihood changes by less than tol (1 + |l|), or after max_iter
# iterations. Returns mu0, var_mu0 (the inverse of sum_j p_null_j / v_j),
# loglik, pi0, p_null, mu, tau2, pi_alt (pi_1, ..., pi_K), and
# bic = -2 loglik + (1{mu0 estimated} + 3K) log(length(b)).
.ac3_da_gaussian_mixture_em <- function(b, v, mu0, mu, tau2, pi_all, fix_mu0,
                                        max_iter = 400L, tol = 1e-8, tau2_floor = 1e-8) {
  K <- length(mu)
  densities <- function(mu0, mu, tau2, pi_all) {
    sd_alt <- sqrt(outer(v, tau2, `+`))
    cbind(
      pi_all[1L] * stats::dnorm(b, mu0, sqrt(v)),
      sweep(stats::dnorm(outer(b, mu, `-`) / sd_alt) / sd_alt, 2L, pi_all[-1L], `*`)
    )
  }
  ll_old <- -Inf
  for (it in seq_len(max_iter)) {
    dens <- densities(mu0, mu, tau2, pi_all)
    tot <- rowSums(dens)
    tot[tot <= 0 | !is.finite(tot)] <- .Machine$double.xmin
    ll <- sum(log(tot))
    resp <- dens / tot
    pi_all <- pmax(colMeans(resp), 1e-8)
    pi_all <- pi_all / sum(pi_all)
    if (!fix_mu0) {
      w0 <- resp[, 1L] / v
      mu0 <- sum(w0 * b) / sum(w0)
    }
    for (k in seq_len(K)) {
      r_k <- resp[, k + 1L]
      w <- 1 / (v + tau2[k])
      mu[k] <- sum(r_k * w * b) / sum(r_k * w)
      e2 <- (b - mu[k])^2
      for (inner in 1:3) {
        w <- 1 / (v + tau2[k])
        tau2[k] <- max(sum(r_k * w^2 * (e2 - v)) / sum(r_k * w^2), tau2_floor)
      }
    }
    if (abs(ll - ll_old) < tol * (1 + abs(ll))) break
    ll_old <- ll
  }
  dens <- densities(mu0, mu, tau2, pi_all)
  tot <- pmax(rowSums(dens), .Machine$double.xmin)
  resp <- dens / tot
  list(
    mu0 = mu0,
    var_mu0 = 1 / sum(resp[, 1L] / v),
    loglik = sum(log(tot)),
    pi0 = pi_all[1L],
    p_null = resp[, 1L],
    mu = mu,
    tau2 = tau2,
    pi_alt = pi_all[-1L],
    bic = -2 * sum(log(tot)) + ((!fix_mu0) + 3 * K) * log(length(b))
  )
}

# Starting values of mu0, in order: the maximal-overlap center c_T; the
# locations of up to six local maxima of the overlap count N_1.5(c) whose
# height is at least 0.3 max_c N_1.5(c), in decreasing order of height; and the
# 0.1, 0.3, 0.5, 0.7, and 0.9 quantiles of the coefficients. A starting value
# within 0.5 median s_j of an earlier one is removed.
.ac3_da_gaussian_mixture_starts <- function(b, se, center) {
  m <- length(b)
  event_x <- c(b - 1.5 * se, b + 1.5 * se)
  event_d <- c(rep(1L, m), rep(-1L, m))
  ord <- order(event_x, -event_d)
  x_sorted <- event_x[ord]
  run <- cumsum(event_d[ord])
  nx <- length(x_sorted)
  seg_a <- x_sorted[-nx]
  seg_b <- x_sorted[-1L]
  keep <- seg_b > seg_a
  height <- run[-nx][keep]
  mid <- ((seg_a + seg_b) / 2)[keep]
  peaks <- numeric(0)
  if (length(height) > 0L) {
    nh <- length(height)
    prev <- c(-Inf, height[-nh])
    nxt <- c(height[-1L], -Inf)
    is_peak <- which(height >= prev & height >= nxt & height >= 0.3 * max(height))
    peaks <- mid[is_peak][order(-height[is_peak])]
  }
  st <- c(center, utils::head(peaks, 6L),
          as.numeric(stats::quantile(b, c(0.1, 0.3, 0.5, 0.7, 0.9), names = FALSE)))
  ms <- stats::median(se)
  keep_st <- logical(length(st))
  for (i in seq_along(st)) keep_st[i] <- !any(abs(st[i] - st[keep_st]) < 0.5 * ms)
  st[keep_st]
}

# Best fit over the starting values for K further components. With
# fix_at = NULL, mu0 is estimated from each starting value; otherwise mu0 is
# held at fix_at, and fix_at followed by the starting values position the
# further components. For a starting value mu0^(0), the further components
# start at the quantiles (k - 1/2) / K of the coefficients farther than
# 2 median s_j from mu0^(0) (of all coefficients when fewer than five are that
# far), with tau2_k = max(var(b) - median(v), 1e-4) and
# (pi_0, pi_k) = (0.3, 0.7 / K). A fit is retained when its log-likelihood is
# finite and pi0 >= 0.05, and the retained fit with the largest log-likelihood
# is returned (NULL when none is retained).
.ac3_da_gaussian_mixture_fit <- function(b, v, starts, K, fix_at = NULL) {
  ms <- stats::median(sqrt(v))
  vb <- max(stats::var(b) - stats::median(v), 1e-4)
  s_list <- if (is.null(fix_at)) starts else unique(c(fix_at, starts))
  best <- NULL
  for (s0 in s_list) {
    ctr <- if (is.null(fix_at)) s0 else fix_at
    rest <- b[abs(b - s0) > 2 * ms]
    if (length(rest) < 5L) rest <- b
    mu_init <- as.numeric(stats::quantile(rest, probs = (seq_len(K) - 0.5) / K, names = FALSE))
    f <- tryCatch(
      .ac3_da_gaussian_mixture_em(b, v, mu0 = ctr, mu = mu_init, tau2 = rep(vb, K),
                                  pi_all = c(0.3, rep(0.7 / K, K)), fix_mu0 = !is.null(fix_at)),
      error = function(e) NULL
    )
    if (is.null(f) || !is.finite(f$loglik) || f$pi0 < 0.05) next
    if (is.null(best) || f$loglik > best$loglik) best <- f
  }
  best
}

# -----------------------------------------------------------------------------
# Selection of the common-shift estimate
# -----------------------------------------------------------------------------
# The Gaussian mixture is fitted with mu0 estimated (K = 1 and 2, K selected by
# BIC) and, at the selected K, with mu0 fixed at the maximal-overlap estimate
# delta_mo; the remaining parameters are estimated in each fit. The mixture
# estimate delta_gm is selected when the log-likelihood gain
# l_free - l_fixed, the difference of the maximized log-likelihoods of the two
# fits, exceeds min_gain (default 5) and |delta_gm - delta_mo| exceeds
# min_separation times the median working standard error; otherwise the
# maximal-overlap estimate is retained. Both are prespecified selection
# criteria, not a likelihood-ratio test. The comparison requires at least
# min_taxa coefficients and a finite delta_mo. Statuses:
# "gaussian_mixture_selected", "maximal_overlap_retained",
# "insufficient_taxa", "free_fit_failed" (no acceptable fit with mu0
# estimated), and "fixed_fit_failed" (no acceptable fit with mu0 fixed).
# Returns selected, loglik_gain, delta_gm, var_delta_gm, n_components (the
# selected K), p_null and nu_gm (over all input coefficients; NA where a
# coefficient or its variance is not finite), and status.
.ac3_da_shift_selection <- function(beta,
                                    var,
                                    delta_mo,
                                    center,
                                    min_gain = 5,
                                    min_separation = 1,
                                    min_taxa = 10L) {
  ok <- is.finite(beta) & is.finite(var) & var > 0
  b <- as.numeric(beta[ok])
  v <- as.numeric(var[ok])
  out <- list(
    selected = FALSE,
    loglik_gain = NA_real_,
    delta_gm = NA_real_,
    var_delta_gm = NA_real_,
    n_components = NA_integer_,
    p_null = rep(NA_real_, length(beta)),
    nu_gm = rep(NA_real_, length(beta)),
    status = "insufficient_taxa"
  )
  if (length(b) < min_taxa || !is.finite(delta_mo)) return(out)
  se <- sqrt(v)
  start_center <- if (is.finite(center)) center else delta_mo
  starts <- .ac3_da_gaussian_mixture_starts(b, se, start_center)
  fit_free <- lapply(1:2, function(K) .ac3_da_gaussian_mixture_fit(b, v, starts, K))
  ok_fit <- !vapply(fit_free, is.null, logical(1))
  if (!any(ok_fit)) {
    out$status <- "free_fit_failed"
    return(out)
  }
  K <- which(ok_fit)[which.min(vapply(fit_free[ok_fit], function(f) f$bic, numeric(1)))]
  fit_fixed <- .ac3_da_gaussian_mixture_fit(b, v, starts, K, fix_at = delta_mo)
  # Mixture estimate: the weighted mean of the coefficients with the weights
  # nu_gm (posterior null probability divided by v_j, normalized), so that the
  # shift applied on selection equals sum(nu_gm * beta) exactly.
  w0 <- fit_free[[K]]$p_null / v
  nu_gm <- w0 / sum(w0)
  out$delta_gm <- sum(nu_gm * b)
  out$var_delta_gm <- fit_free[[K]]$var_mu0
  out$n_components <- as.integer(K)
  out$p_null[ok] <- fit_free[[K]]$p_null
  out$nu_gm[ok] <- nu_gm
  if (is.null(fit_fixed)) {
    out$status <- "fixed_fit_failed"
    return(out)
  }
  out$loglik_gain <- fit_free[[K]]$loglik - fit_fixed$loglik
  separation <- abs(out$delta_gm - delta_mo) / stats::median(se)
  out$selected <- is.finite(out$loglik_gain) && out$loglik_gain > min_gain &&
    is.finite(separation) && separation > min_separation
  out$status <- if (out$selected) "gaussian_mixture_selected" else "maximal_overlap_retained"
  out
}

# Common-shift estimate of each tested covariate of coefs (columns taxon, term,
# beta_raw, var_raw). For each covariate, the maximal-overlap estimate
# (.ac3_da_maximal_overlap()) is computed over its eligible taxa,
# eligible_taxa_by_term[[term]] (every taxon of the covariate when the list or
# the entry is NULL), and the selection rule (.ac3_da_shift_selection())
# chooses between it and the Gaussian-mixture estimate. Returns one row per covariate. The attributes null_set and
# correction_weights hold, per covariate, the reported estimated null set (A^MO
# under the maximal-overlap estimator; the taxa with p_null >= 0.5 under the
# mixture estimator) and the correction weights nu of the selected estimator,
# named by taxon: delta = sum(nu * beta_raw), and nu is empty when the estimate
# is a constant (no finite coefficient, or a degenerate case of the
# maximal-overlap estimator).
.ac3_da_common_shift <- function(coefs,
                                 eligible_taxa_by_term = NULL,
                                 min_count = 25L,
                                 bandwidths = seq(0.5, 3.0, by = 0.1),
                                 min_gain = 5,
                                 min_separation = 1) {
  coefs <- as.data.frame(coefs)
  min_count <- as.integer(min_count)
  if (length(min_count) != 1L || is.na(min_count)) {
    stop("settings$da_maximal_overlap_min_count must be a single integer.")
  }
  ## The estimated null set and the correction weights are collected alongside
  ## the per-term summary for the corrected variance.
  shift_env <- new.env(parent = emptyenv())
  shift_env$null_set <- list()
  shift_env$weights <- list()
  out <- do.call(rbind, lapply(split(coefs, coefs$term), function(g) {
    term <- g$term[1L]
    if (!is.null(eligible_taxa_by_term)) {
      eligible <- eligible_taxa_by_term[[term]]
      if (!is.null(eligible)) g <- g[g$taxon %in% eligible, , drop = FALSE]
    }
    g <- g[is.finite(g$beta_raw) & is.finite(g$var_raw) & g$var_raw > 0, , drop = FALSE]
    if (nrow(g) == 0L) {
      shift_env$null_set[[term]] <- character(0)
      shift_env$weights[[term]] <- numeric(0)
      return(data.frame(
        term = term,
        delta = 0,
        n_eligible_coefficients = 0L,
        maximal_overlap_n_null_set = 0L,
        maximal_overlap_bandwidth = NA_real_,
        maximal_overlap_center = NA_real_,
        maximal_overlap_min_count = min_count,
        maximal_overlap_status = "no_finite_coefficients",
        maximal_overlap_delta = 0,
        maximal_overlap_var_delta = 0,
        shift_selection_loglik_gain = NA_real_,
        gaussian_mixture_delta = NA_real_,
        gaussian_mixture_selected = FALSE,
        shift_selection_status = "insufficient_taxa",
        row.names = NULL
      ))
    }
    mo <- .ac3_da_maximal_overlap(
      beta = g$beta_raw,
      se = sqrt(g$var_raw),
      min_count = min_count,
      bandwidths = bandwidths
    )
    ## Variance of the maximal-overlap estimate with the null set held fixed.
    mo_var <- if (is.finite(mo$var_delta_mo)) mo$var_delta_mo else 0
    sel <- .ac3_da_shift_selection(
      beta = g$beta_raw,
      var = g$var_raw,
      delta_mo = mo$delta_mo,
      center = mo$center,
      min_gain = min_gain,
      min_separation = min_separation
    )
    if (isTRUE(sel$selected)) {
      ## The mixture estimate is selected. It is sum(nu_gm * beta_raw) with the
      ## weights nu_gm (posterior null probability divided by s_j^2,
      ## normalized), which are the correction weights. The reported null set
      ## holds the taxa with posterior null probability at least 0.5.
      ok <- is.finite(sel$nu_gm)
      shift_env$null_set[[term]] <- as.character(g$taxon[ok & sel$p_null >= 0.5])
      shift_env$weights[[term]] <- stats::setNames(sel$nu_gm[ok], as.character(g$taxon[ok]))
    } else {
      shift_env$null_set[[term]] <- as.character(g$taxon[mo$null_set_pos])
      shift_env$weights[[term]] <- stats::setNames(mo$nu, as.character(g$taxon[mo$null_set_pos]))
    }
    data.frame(
      term = term,
      delta = if (isTRUE(sel$selected)) sel$delta_gm else mo$delta_mo,
      n_eligible_coefficients = nrow(g),
      maximal_overlap_n_null_set = mo$n_null_set,
      maximal_overlap_bandwidth = mo$bandwidth,
      maximal_overlap_center = mo$center,
      maximal_overlap_min_count = min_count,
      maximal_overlap_status = mo$status,
      maximal_overlap_delta = mo$delta_mo,
      maximal_overlap_var_delta = mo_var,
      shift_selection_loglik_gain = sel$loglik_gain,
      gaussian_mixture_delta = sel$delta_gm,
      gaussian_mixture_selected = isTRUE(sel$selected),
      shift_selection_status = sel$status,
      row.names = NULL
    )
  }))
  attr(out, "null_set") <- shift_env$null_set
  attr(out, "correction_weights") <- shift_env$weights
  out
}

# -----------------------------------------------------------------------------
# Covariance of the uncorrected coefficients
# -----------------------------------------------------------------------------
# The errors of the uncorrected coefficients are approximated by
# e_j ~ eps_j + g_j' eps_d, where g_j = d beta_hat_j / d d at d_hat is the
# sensitivity of the tested coefficient to the sample-specific bias and eps_d
# is the estimation error of d_hat. The working covariance is
#   Cov(e_j, e_l) = 1{j = l} V_j,model + g_j' V_d^con g_l,
# where V_j,model is the model variance of the tested coefficient with d held
# fixed at d_hat and V_d^con is the working covariance of the bias estimator
# under both identifying constraints. Its diagonal is the marginal working
# variance s_j^2 = V_j,model + g_j' V_d^con g_j.

# Information blocks of taxon j at the final fit: I_psi_d = -d2 l_j / d psi_j d d'
# (one column per sample); I_dd_diag, the diagonal of I_dd = -d2 l_j / d d d',
# which is diagonal because the contribution of sample i depends only on d_i;
# and C, the fitted conditional covariance of psi_j with d held fixed. The
# second derivatives in the predictors are the analytic cell derivatives
# (.ac3_pln_observation_derivatives()) at the fitted values and sigma_j. When
# the fit holds the covariance vcov_psi of psi_j = (omega_j, kappa_j)
# (.ac3_joint_augment()), I_psi_d has a further row for kappa_j, a central
# difference in kappa_j (step 0.05) of the analytic derivative dl_ij / dnu_ij,
# and C is vcov_psi; otherwise psi_j is omega_j and C is vcov. engine selects
# the evaluation of the quadrature (.ac3_pln_obs()).
.ac3_da_taxon_bias_information <- function(fit_j,
                                           y_j,
                                           x,
                                           r,
                                           engine = "r") {
  if (!isTRUE(fit_j$ok) || is.null(fit_j$vcov) || is.null(fit_j$hessian)) {
    return(NULL)
  }
  if (any(!is.finite(fit_j$vcov)) || any(!is.finite(fit_j$hessian))) {
    return(NULL)
  }
  X <- cbind("(Intercept_mu)" = 1, x)
  R <- cbind("(Intercept_prev)" = 1, r)
  n_par <- ncol(X) + ncol(R)
  if (nrow(fit_j$vcov) != n_par || ncol(fit_j$vcov) != n_par) {
    return(NULL)
  }
  log_mu <- log(pmax(fit_j$mu, .Machine$double.xmin))
  logit_pi <- qlogis(.ac3_da_clip(fit_j$pi, 1e-12, 1 - 1e-12))
  pos <- y_j > 0
  sigma_j <- fit_j$sigma
  if (is.null(sigma_j) || !is.finite(sigma_j) || sigma_j <= 0) return(NULL)
  o <- .ac3_pln_observation_derivatives(log_mu, logit_pi, y_j, pos, sigma_j, need = 2L, engine = engine)
  ok_obs <- is.finite(y_j) & is.finite(log_mu) & is.finite(logit_pi)
  I_dd_diag <- ifelse(ok_obs & is.finite(o$d2l_mumu), -o$d2l_mumu, 0)
  I_logit_pi_d <- ifelse(ok_obs & is.finite(o$d2l_mupi), -o$d2l_mupi, 0)
  I_psi_d <- rbind(t(X * I_dd_diag), t(R * I_logit_pi_d))
  if (!is.null(fit_j$vcov_psi)) {
    # kappa_j is a taxon-wise parameter. Its cross-information with d,
    # -d2 l / d kappa_j d d_i, is a central difference in kappa_j (step 0.05)
    # of the analytic derivative in the log mean.
    hk <- 0.05
    op <- .ac3_pln_observation_derivatives(log_mu, logit_pi, y_j, pos, sigma_j * exp(hk), need = 1L,
                                           engine = engine)
    om <- .ac3_pln_observation_derivatives(log_mu, logit_pi, y_j, pos, sigma_j * exp(-hk), need = 1L,
                                           engine = engine)
    I_kappa_d <- -(op$dl_mu - om$dl_mu) / (2 * hk)
    I_kappa_d[!ok_obs | !is.finite(I_kappa_d)] <- 0
    return(list(I_psi_d = rbind(I_psi_d, I_kappa_d), I_dd_diag = I_dd_diag, C = fit_j$vcov_psi))
  }
  list(I_psi_d = I_psi_d, I_dd_diag = I_dd_diag, C = fit_j$vcov)
}

# Regularized inverse of a symmetric matrix: each positive, finite eigenvalue
# lambda_k is replaced by max(lambda_k, lambda_min) and inverted, with
# lambda_min = max(1e-8, 1e-8 lambda_med) and lambda_med the median of the
# positive, finite eigenvalues (1 when there is none); nonpositive or
# nonfinite eigenvalues receive weight zero.
.ac3_da_symmetric_inverse_with_floor <- function(mat) {
  eigen_floor_rel <- 1e-8
  eigen_floor_abs <- 1e-8
  mat <- 0.5 * (mat + t(mat))
  eig <- eigen(mat, symmetric = TRUE)
  positive <- eig$values[eig$values > 0 & is.finite(eig$values)]
  lambda_med <- if (length(positive) > 0L) stats::median(positive) else 1
  lambda_min <- max(eigen_floor_abs, eigen_floor_rel * lambda_med)
  positive_idx <- is.finite(eig$values) & eig$values > 0
  values_used <- ifelse(positive_idx, pmax(eig$values, lambda_min), Inf)
  inverse_values <- ifelse(positive_idx, 1 / values_used, 0)
  inv <- eig$vectors %*% diag(inverse_values, nrow = length(inverse_values)) %*%
    t(eig$vectors)
  list(
    inverse = 0.5 * (inv + t(inv)),
    eigenvalues = eig$values,
    lambda_min = lambda_min,
    n_floored = sum(eig$values < lambda_min | !is.finite(eig$values)),
    n_nonpos = sum(eig$values <= 0 | !is.finite(eig$values))
  )
}

# Working covariance V_d^con of the bias estimator under both identifying
# constraints. With H0 an orthonormal basis of S (.ac3_da_constraint_basis()),
# the restricted information K_theta = H0' K_d H0 is inverted with the
# eigenvalue floor (.ac3_da_symmetric_inverse_with_floor()) and returned to
# sample coordinates, V_d^con = H0 K_theta^-1 H0', so that V_d^con 1 = 0 and
# V_d^con Q = 0. A zero-dimensional S gives the zero matrix.
.ac3_da_bias_working_covariance <- function(K_d, x, test_terms) {
  n <- nrow(K_d)
  H0 <- .ac3_da_constraint_basis(x, test_terms)
  m <- ncol(H0)
  if (m == 0L) {
    return(list(V_d = matrix(0, n, n), H0 = H0, m = 0L, eigenvalues = numeric(0),
                lambda_min = NA_real_, n_nonpos = 0L, n_floored = 0L))
  }
  K_theta <- crossprod(H0, K_d %*% H0)
  inv <- .ac3_da_symmetric_inverse_with_floor(K_theta)
  V_d <- H0 %*% inv$inverse %*% t(H0)
  list(
    V_d = 0.5 * (V_d + t(V_d)),
    H0 = H0,
    m = m,
    eigenvalues = inv$eigenvalues,
    lambda_min = inv$lambda_min,
    n_nonpos = inv$n_nonpos,
    n_floored = inv$n_floored
  )
}

# Propagation of the uncertainty of the estimated sample-specific bias d to the
# uncorrected DA coefficients. The information about d with the taxon-wise
# parameters maximized out (the negative Hessian of the sum of the taxon-wise
# objectives evaluated at psi_j(d)),
#   K_d = sum_j (I_dd,j - I_psi_d,j' C_j I_psi_d,j),
# is accumulated over the taxa with a valid fit
# (.ac3_da_taxon_bias_information()) and symmetrized, and V_d^con imposes both
# identifying constraints (.ac3_da_bias_working_covariance()). For each tested
# coefficient, g_j = -e_T' C_j I_psi_d,j, the bias component of the variance is
# g_j' V_d^con g_j (var_raw_bias), and the marginal working variance is
# s_j^2 = V_j,model + g_j' V_d^con g_j (var_raw_marginal, with
# se_raw_marginal). The sensitivities G (one matrix per tested term, whose row
# j is g_j') and V_d^con are returned in fit$da_variance_inputs, the input of
# score_ancombc3_da(). engine selects the evaluation of the quadrature
# (.ac3_pln_obs()).
ancombc3_add_bias_covariance <- function(fit,
                                         y,
                                         x,
                                         r,
                                         test_terms,
                                         engine = "r") {
  if (is.null(fit$raw) || is.null(fit$fits)) return(fit)
  y <- as.matrix(y)
  x <- as.matrix(x)
  r <- as.matrix(r)
  n <- nrow(y)
  p <- ncol(y)
  if (p != length(fit$fits)) stop("fit and y have incompatible taxon dimensions.")
  if (is.null(colnames(y))) colnames(y) <- names(fit$fits)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  test_terms <- intersect(test_terms, colnames(x))

  blocks <- vector("list", p)
  names(blocks) <- names(fit$fits)
  K_d <- matrix(0, nrow = n, ncol = n)
  i_diag <- seq_len(n) * (n + 1L) - n
  taxa_used <- 0L
  for (j in seq_len(p)) {
    block <- .ac3_da_taxon_bias_information(
      fit_j = fit$fits[[j]],
      y_j = y[, j],
      x = x,
      r = r,
      engine = engine
    )
    if (!is.null(block)) {
      K_d[i_diag] <- K_d[i_diag] + block$I_dd_diag
      K_d <- K_d - t(block$I_psi_d) %*% block$C %*% block$I_psi_d
      blocks[[j]] <- block
      taxa_used <- taxa_used + 1L
    }
  }
  K_d <- 0.5 * (K_d + t(K_d))
  raw <- fit$raw
  raw$var_raw_bias <- NA_real_
  raw$var_raw_marginal <- NA_real_
  raw$se_raw_marginal <- NA_real_
  if (taxa_used == 0L || length(test_terms) == 0L) {
    fit$raw <- raw
    fit$da_variance_inputs <- NULL
    return(fit)
  }

  con <- .ac3_da_bias_working_covariance(K_d, x = x, test_terms = test_terms)
  V_d <- con$V_d
  ## Row j of G[[term]] is g_j'.
  G <- lapply(stats::setNames(test_terms, test_terms), function(term) {
    matrix(NA_real_, nrow = p, ncol = n, dimnames = list(names(fit$fits), NULL))
  })
  for (rr in seq_len(nrow(raw))) {
    taxon <- raw$taxon[rr]
    term <- raw$term[rr]
    j <- match(taxon, names(fit$fits))
    if (is.na(j) || !(term %in% test_terms) || is.null(blocks[[j]])) next
    idx_beta <- 1L + match(term, colnames(x))
    if (is.na(idx_beta)) next
    e_T <- rep(0, nrow(blocks[[j]]$C))
    e_T[idx_beta] <- 1
    g_j <- -drop(e_T %*% blocks[[j]]$C %*% blocks[[j]]$I_psi_d)
    if (!all(is.finite(g_j))) next
    G[[term]][j, ] <- g_j
    bias_var <- drop(g_j %*% V_d %*% g_j)
    if (!is.finite(bias_var) || !is.finite(raw$var_raw_model[rr]) || raw$var_raw_model[rr] <= 0) next
    raw$var_raw_bias[rr] <- max(bias_var, 0)
    raw$var_raw_marginal[rr] <- raw$var_raw_model[rr] + raw$var_raw_bias[rr]
    raw$se_raw_marginal[rr] <- sqrt(raw$var_raw_marginal[rr])
  }
  fit$raw <- raw
  fit$da_variance_inputs <- list(
    taxa = names(fit$fits),
    G = G,
    V = V_d
  )
  fit
}

# Variance contributed by the latent term that the abundance offset carries
# (settings$da_offset_variance of ancombc3()). The first-stage
# abundance offset O^A estimates the latent term of each taxon after its
# projection onto the orthogonal complement of the residualized tested
# covariates. The component of the latent term along a tested covariate is
# therefore not removed: it remains in the uncorrected coefficient, as
# lambda_j' (zbar_1 - zbar_0) for a binary covariate, while the offset removes
# the latent term from the residuals from which V_j,model is computed. The
# latent scores are treated as random. To first order the coefficient error
# carries -g_j' L_j, where L_j is the latent term of taxon j and
# g_j = d beta_hat_j / d d is the sensitivity of the tested coefficient
# (ancombc3_add_bias_covariance()), which equals the sensitivity to the
# offset of taxon j because d_i and O^A_ij enter the log mean in the same way.
# Its variance is estimated from the offset, the estimate of the projected
# latent term, by the robust form over independent units (the subjects with
# cluster, the samples otherwise): with a_ij = g_ji O^A_ij centered at its
# mean over the samples and the unit totals A_uj = sum_{i in u} a_ij,
#   column j of A = {m / (m - 1)}^(1/2) (A_1j, ..., A_mj)',
# with m the number of units. The covariance of the offset terms of two
# coefficients is a_j' a_l, so the term has rank at most m across taxa, and
# the variance of the corrected coefficient includes it through
# .ac3_da_corrected_variance(). The result holds fit with A (one m x p matrix
# per tested term, columns in the order of da_variance_inputs$taxa) added to
# da_variance_inputs, and with raw receiving var_raw_offset = ||a_j||^2,
# added to var_raw_marginal and se_raw_marginal. A taxon without a finite
# sensitivity has a zero column.
ancombc3_add_offset_covariance <- function(fit, abundance_offset, cluster = NULL) {
  vi <- fit$da_variance_inputs
  if (is.null(vi) || is.null(vi$G)) return(fit)
  O <- .ac3_match_taxa(as.matrix(abundance_offset), vi$taxa, "abundance_offset")
  n <- nrow(O)
  ids <- if (is.null(cluster)) seq_len(n) else .ac3_cluster_index(cluster)
  if (length(ids) != n) stop("cluster must have one value per sample.")
  m <- max(ids)
  if (m < 2L) return(fit)
  scale <- sqrt(m / (m - 1))
  vi$A <- lapply(vi$G, function(G) {
    a <- t(G) * O
    a[!is.finite(a)] <- 0
    a <- sweep(a, 2L, colMeans(a), "-")
    A <- rowsum(a, ids, reorder = TRUE) * scale
    dimnames(A) <- list(NULL, vi$taxa)
    A
  })
  raw <- fit$raw
  raw$var_raw_offset <- NA_real_
  for (rr in seq_len(nrow(raw))) {
    term <- raw$term[rr]
    j <- match(raw$taxon[rr], vi$taxa)
    if (is.na(j) || is.null(vi$A[[term]]) || !is.finite(raw$var_raw_marginal[rr])) next
    v_off <- sum(vi$A[[term]][, j]^2)
    raw$var_raw_offset[rr] <- v_off
    raw$var_raw_marginal[rr] <- raw$var_raw_marginal[rr] + v_off
    raw$se_raw_marginal[rr] <- sqrt(raw$var_raw_marginal[rr])
  }
  fit$raw <- raw
  fit$da_variance_inputs <- vi
  fit
}

# -----------------------------------------------------------------------------
# Working variance of the corrected DA coefficient
# -----------------------------------------------------------------------------
# Working variance V_j,primary of beta_hat_j - delta_hat, where
# delta_hat = sum_l nu_l beta_hat_l with fixed correction weights nu (zero
# outside the taxa that contribute to the estimate; all zero when the
# common-shift estimate is a constant). With
# Sigma = diag(V_model) + G V_d^con G', whose row j of G is g_j',
#   V_j,primary = (1 - nu_j)^2 V_j,model + sum_{l != j} nu_l^2 V_l,model
#                 + (g_j - gbar)' V_d^con (g_j - gbar),   gbar = sum_l nu_l g_l,
# with the last term bounded below by zero. Also returned: the variance of the
# common-shift estimate, V_delta = nu' Sigma nu (var_delta), and its covariance
# with each coefficient, C_j,delta = (Sigma nu)_j = nu_j V_j,model + B_j
# (cov_delta), with B_j = g_j' V_d^con gbar. When the offset term A of
# ancombc3_add_offset_covariance() is supplied, Sigma includes A' A.
.ac3_da_corrected_variance <- function(var_model, G, V, nu, A = NULL) {
  var_model <- as.numeric(var_model)
  nu <- as.numeric(nu)
  G <- as.matrix(G)
  gv <- G %*% V
  gbar <- drop(crossprod(G, nu))
  v_gbar <- drop(V %*% gbar)
  own <- rowSums(gv * G)
  B_j <- drop(G %*% v_gbar)
  shared_delta <- sum(gbar * v_gbar)
  model_delta <- sum(nu^2 * var_model)
  out <- list(
    var_primary = (1 - nu)^2 * var_model + (model_delta - nu^2 * var_model) +
      pmax(own - 2 * B_j + shared_delta, 0),
    var_delta = model_delta + shared_delta,
    cov_delta = nu * var_model + B_j
  )
  if (!is.null(A)) {
    # Offset term (ancombc3_add_offset_covariance()): with the columns a_j of
    # A and abar = sum_l nu_l a_l, the term ||a_j - abar||^2 is added to
    # V_j,primary, ||abar||^2 to V_delta, and a_j' abar to C_j,delta.
    A <- as.matrix(A)
    A[!is.finite(A)] <- 0
    abar <- drop(A %*% nu)
    B_A <- drop(crossprod(A, abar))
    shared_A <- sum(abar^2)
    out$var_primary <- out$var_primary + colSums(A^2) - 2 * B_A + shared_A
    out$var_delta <- out$var_delta + shared_A
    out$cov_delta <- out$cov_delta + B_A
  }
  out
}

# -----------------------------------------------------------------------------
# Empirical null of the DA statistics
# -----------------------------------------------------------------------------
# The normal-equivalent statistics z_tilde_j of one covariate are centered at
# the common-shift estimate (step 1 of score_ancombc3_da()), so the null
# location is zero. They are fitted by maximum likelihood to the two-group
# model
#   z_tilde_j ~ pi_0 N(0, omega_0^2) + sum_k pi_k N(mu_k, omega_0^2 + tau2_k),  k = 1, ..., K,
# with K = 0, 1, or 2. The null component has location zero and standard
# deviation omega_0; each further component has the variance
# omega_0^2 + tau2_k, with tau2_k >= 1e-8. The model with K = 0, in which
# every statistic follows the null component, has the maximum-likelihood
# estimate omega_0^2 = mean(z_tilde_j^2). For K = 1 and 2, fits whose null
# weight pi_0 is below 0.05 are discarded. The BIC of a fit has 1 + 3K
# free parameters. The empirical-null scale is omega_T = max(1, omega_hat),
# where omega_hat is selected from the fits by the rule of
# .ac3_da_empirical_null(); omega_T = 1 when fewer than
# min_n = 10 statistics are finite.

# EM iterations for the model above with K >= 1 from omega2 = omega_0^2, mu,
# tau2, and pi_all = (pi_0, pi_1, ..., pi_K). Each iteration updates the
# weights as in .ac3_da_gaussian_mixture_em(), omega2 as the
# null-responsibility-weighted mean of z_tilde_j^2 (bounded below by
# omega2_floor = 1), mu_k as the responsibility-weighted mean of z_tilde_j, and
# tau2_k as the corresponding weighted variance minus omega2 (bounded below by
# tau2_floor). With omega2_fixed = TRUE, omega2 is held at its starting value,
# and the null density is evaluated once. The iteration stops when the
# log-likelihood changes by less than tol (1 + |l|), or after max_iter
# iterations.
.ac3_da_empirical_null_em <- function(z_tilde, omega2, mu, tau2, pi_all,
                                      max_iter = 500L, tol = 1e-8,
                                      omega2_floor = 1, tau2_floor = 1e-8, omega2_fixed = FALSE) {
  K <- length(mu)
  dnorm <- stats::dnorm
  z2 <- z_tilde^2
  null_fixed <- if (isTRUE(omega2_fixed)) dnorm(z_tilde, 0, sqrt(omega2)) else NULL
  densities <- function(omega2, mu, tau2, pi_all) {
    dens <- matrix(0, length(z_tilde), K + 1L)
    null_dens <- if (is.null(null_fixed)) dnorm(z_tilde, 0, sqrt(omega2)) else null_fixed
    dens[, 1L] <- pi_all[1L] * null_dens
    for (k in seq_len(K)) {
      sd_k <- sqrt(omega2 + tau2[k])
      dens[, k + 1L] <- dnorm((z_tilde - mu[k]) / sd_k) / sd_k * pi_all[k + 1L]
    }
    dens
  }
  ll_old <- -Inf
  for (it in seq_len(max_iter)) {
    dens <- densities(omega2, mu, tau2, pi_all)
    tot <- .ac3_pln_clip(rowSums(dens), .Machine$double.xmin, Inf)
    ll <- sum(log(tot))
    resp <- dens / tot
    pi_all <- .ac3_pln_clip(colMeans(resp), 1e-8, Inf)
    pi_all <- pi_all / sum(pi_all)
    if (!isTRUE(omega2_fixed)) omega2 <- max(sum(resp[, 1L] * z2) / sum(resp[, 1L]), omega2_floor)
    for (k in seq_len(K)) {
      r_k <- resp[, k + 1L]
      sum_r_k <- sum(r_k)
      mu[k] <- sum(r_k * z_tilde) / sum_r_k
      e2 <- (z_tilde - mu[k])^2
      tau2[k] <- max(sum(r_k * e2) / sum_r_k - omega2, tau2_floor)
    }
    if (abs(ll - ll_old) < tol * (1 + abs(ll))) break
    ll_old <- ll
  }
  dens <- densities(omega2, mu, tau2, pi_all)
  tot <- pmax(rowSums(dens), .Machine$double.xmin)
  list(
    omega2 = omega2,
    mu = mu,
    tau2 = tau2,
    pi0 = pi_all[1L],
    pi_all = pi_all,
    loglik = sum(log(tot)),
    bic = -2 * sum(log(tot)) + (1 + 3 * K) * log(length(z_tilde))
  )
}

# Empirical-null scale of one covariate. K = 0 is fitted in closed form. For
# K = 1 and 2, six starting values: the further components start at the
# quantiles (k - 1/2) / K of the statistics with |z_tilde_j| > 2 (of all
# statistics when fewer than two exceed 2), at the same quantiles of all
# statistics, and at the median (K = 1) or the minimum and maximum (K = 2) of
# the former set, each with omega_0 = 1 and omega_0 = 1.5; tau2_k starts at
# max(var(z_tilde) - 1, 1), and the weights at pi_0 = 0.5 and pi_k = 0.5 / K.
# A fit with K >= 1 is retained when its null weight pi_0 is at least 0.05,
# and omega_0 >= 1 in every fit (for K = 0, omega_0^2 = max{mean(z_tilde_j^2),
# 1}). Let BIC_min be the smallest BIC over the fit with K = 0 and the
# retained fits of all starting values, and let BIC_p(omega) be the smallest
# BIC over K = 0, 1, 2 of the fits with the null standard deviation held at
# omega (same penalty 1 + 3K). omega_hat is the largest omega >= 1 on the grid
# with BIC_p(omega) <= BIC_min + 3.84 (the 0.95 quantile of the chi-square
# distribution with one degree of freedom): the upper limit of the profile
# region of omega_0. The grid has step 0.05 from 1 to 1.25 times max{1, root
# mean square of z_tilde}, refined with step 0.005 above the last grid value
# inside the region. The fits at a grid value start from the three starting
# vectors above and from the fit at the preceding grid value.
# Returns factor (omega_T = max(1, omega_hat)), omega_hat, and n_components
# (the K of the fit that supplies omega_hat; NA when fewer than min_n
# statistics are finite).
.ac3_da_empirical_null <- function(z_tilde, min_n = 10L) {
  window <- 3.84
  z_tilde <- z_tilde[is.finite(z_tilde)]
  out <- list(factor = 1, omega_hat = NA_real_, n_components = NA_integer_)
  n <- length(z_tilde)
  if (n < min_n) return(out)
  rest <- z_tilde[abs(z_tilde) > 2]
  if (length(rest) < 2L) rest <- z_tilde
  tau2_init <- max(stats::var(z_tilde) - 1, 1)
  mu_starts <- lapply(1:2, function(K) {
    probs <- (seq_len(K) - 0.5) / K
    list(
      stats::quantile(rest, probs = probs, names = FALSE),
      stats::quantile(z_tilde, probs = probs, names = FALSE),
      if (K == 2L) c(min(rest), max(rest)) else stats::median(rest)
    )
  })
  null_only <- function(omega2) {
    ll <- sum(stats::dnorm(z_tilde, 0, sqrt(omega2), log = TRUE))
    list(omega2 = omega2, loglik = ll, bic = -2 * ll + log(n), K = 0L)
  }
  all_fits <- list(null_only(max(mean(z_tilde^2), 1)))
  for (K in 1:2) {
    for (mu_init in mu_starts[[K]]) for (omega_init in c(1, 1.5)) {
      f <- tryCatch(
        .ac3_da_empirical_null_em(z_tilde, omega2 = omega_init^2, mu = as.numeric(mu_init),
                                  tau2 = rep(tau2_init, K), pi_all = c(0.5, rep(0.5 / K, K))),
        error = function(e) NULL
      )
      if (is.null(f) || !is.finite(f$loglik) || f$pi0 < 0.05) next
      f$K <- K
      all_fits[[length(all_fits) + 1L]] <- f
    }
  }
  bic_of <- function(fits) vapply(fits, function(f) f$bic, numeric(1))
  finish <- function(omega_hat, K) {
    list(factor = max(1, omega_hat), omega_hat = omega_hat, n_components = as.integer(K))
  }
  bic_min <- min(bic_of(all_fits))
  # Profile of the BIC in omega_0.
  warm <- list(NULL, NULL)
  profile_at <- function(omega) {
    best <- null_only(omega^2)
    for (K in 1:2) {
      starts <- lapply(mu_starts[[K]], function(m) {
        list(mu = as.numeric(m), tau2 = rep(tau2_init, K), pi_all = c(0.5, rep(0.5 / K, K)))
      })
      if (!is.null(warm[[K]])) starts <- c(list(warm[[K]]), starts)
      best_K <- NULL
      for (st in starts) {
        f <- tryCatch(
          .ac3_da_empirical_null_em(z_tilde, omega2 = omega^2, mu = st$mu, tau2 = st$tau2,
                                    pi_all = st$pi_all, omega2_fixed = TRUE),
          error = function(e) NULL
        )
        if (is.null(f) || !is.finite(f$loglik) || f$pi0 < 0.05) next
        if (is.null(best_K) || f$loglik > best_K$loglik) best_K <- f
      }
      if (is.null(best_K)) next
      warm[[K]] <<- list(mu = best_K$mu, tau2 = best_K$tau2, pi_all = best_K$pi_all)
      best_K$K <- K
      if (best_K$bic < best$bic) best <- best_K
    }
    best
  }
  upper <- 1.25 * max(1, sqrt(mean(z_tilde^2)))
  grid <- seq(1, max(upper, 1.05), by = 0.05)
  prof <- lapply(grid, profile_at)
  inside <- which(bic_of(prof) <= bic_min + window)
  if (length(inside) == 0L) {
    f <- all_fits[[which.min(bic_of(all_fits))]]
    return(finish(sqrt(f$omega2), f$K))
  }
  last <- max(inside)
  omega_hat <- grid[last]
  K_hat <- prof[[last]]$K
  if (last < length(grid)) {
    warm <- list(NULL, NULL)
    for (omega in seq(grid[last] + 0.005, grid[last] + 0.045, by = 0.005)) {
      f <- profile_at(omega)
      if (f$bic > bic_min + window) break
      omega_hat <- omega
      K_hat <- f$K
    }
  }
  finish(omega_hat, K_hat)
}

# -----------------------------------------------------------------------------
# DA inference: score_ancombc3_da()
# -----------------------------------------------------------------------------
# DA statistics and P values of the coefficients in raw (fit$raw of
# fit_ancombc3_da()); variance_inputs is fit$da_variance_inputs.
# eligible_taxa_by_term holds the DA-eligible taxa of each tested covariate
# (ancombc3_da_eligibility()$by_term, or fit$eligible_taxa_by_term); a NULL
# list or a NULL entry selects every scored taxon of the covariate. This
# default, for callers that do not apply the eligibility criteria, departs
# from the manuscript, which restricts the common-shift estimate and the scale
# adjustment to the DA-eligible taxa. var_raw and se_raw are set to the
# marginal working variance s_j^2 and its square root (var_raw_marginal and
# se_raw_marginal). For each tested covariate:
#   1. Common-shift estimate: the maximal-overlap estimate over the eligible
#      taxa (.ac3_da_maximal_overlap(), with the bandwidth grid
#      maximal_overlap_bandwidths and the overlap criterion
#      maximal_overlap_min_count), or the Gaussian-mixture estimate when the
#      selection rule (.ac3_da_shift_selection(), with the thresholds
#      shift_selection_min_gain of the log-likelihood gain l_free - l_fixed
#      and shift_selection_min_separation) selects it. The correction weights nu of the selected estimator
#      are zero outside the taxa that contribute to it,
#      delta = sum(nu * beta_raw), and nu is all zero when the estimate is a
#      constant.
#   2. DA statistic: Z_j = (beta_raw_j - delta) / V_j,primary^(1/2), with the
#      working variance V_j,primary of the corrected coefficient
#      (.ac3_da_corrected_variance(), with the model variances var_raw_model,
#      the sensitivities G and the bias covariance V of variance_inputs, and
#      nu). Taxa whose V_j,primary is not finite and positive have no
#      statistic or P value.
#   3. Degrees of freedom: df_j = max(n_j - 1 - q_x, 1), where n_j is the
#      number of positive counts of taxon j (column n_positive) and 1 + q_x is
#      the number of abundance coefficients, the intercept and the q_x
#      abundance covariate columns (column n_abundance_coef). When raw has the
#      column df_count (repeated measures: the number of subjects with a
#      positive count of the taxon), n_j is df_count.
#   4. Normal-equivalent statistic: z_tilde_j = sign(Z_j) Phi^-1{F_df_j(|Z_j|)},
#      where F_df is the t distribution function with df degrees of freedom;
#      it is computed from the log tail probability. It is the input of step 5.
#   5. Empirical null: the scale omega_T = max(1, omega_hat) of the null
#      component, located at zero, of the two-group model
#      (.ac3_da_empirical_null()) fitted to the normal-equivalent statistics
#      of the eligible taxa with a finite statistic (the set F_T). omega_T^2
#      multiplies the working variance: se_scaled = omega_T se_corrected.
#   6. Statistic and P value: T_j = beta_corrected / se_scaled = Z_j / omega_T
#      and p_j = 2 {1 - F_df_j(|T_j|)}, adjusted by p_adj_method
#      (BH by default) over the scored taxa of the covariate, eligible or not
#      (q_value). ancombc3() recomputes the q-values over the DA-eligible
#      taxa.
# The P value is computed on the log scale. The result has one row
# per taxon and covariate with a finite coefficient and a positive variance,
# eligible or not. It holds delta, var_delta (V_delta), cov_delta
# (C_j,delta), correction_weight (nu_j), beta_corrected, var_primary
# (V_j,primary), se_corrected (its square root), se_scaled (omega_T
# se_corrected), statistic (T_j), statistic_normal_equivalent (z_tilde_j),
# df_t (df_j), adjustment_factor (omega_T),
# adjustment_factor_hat (omega_hat) and
# adjustment_factor_components (K; NA when fewer than ten statistics are
# finite), in_null_set (membership of the reported estimated null set),
# method = "ancombc3_da", and the per-covariate summaries of the
# maximal-overlap estimator (n_eligible_coefficients and the columns
# maximal_overlap_*), of the Gaussian mixture (gaussian_mixture_*), and of the
# selection rule (shift_selection_*). raw must hold the columns n_positive
# and n_abundance_coef, which fit_ancombc3_da() returns.
score_ancombc3_da <- function(raw,
                              variance_inputs,
                              eligible_taxa_by_term = NULL,
                              p_adj_method = "BH",
                              maximal_overlap_min_count = 25L,
                              maximal_overlap_bandwidths = seq(0.5, 3.0, by = 0.1),
                              shift_selection_min_gain = 5,
                              shift_selection_min_separation = 1) {
  raw <- as.data.frame(raw)
  if (!all(c("n_positive", "n_abundance_coef") %in% names(raw))) {
    stop("score_ancombc3_da requires the columns n_positive and n_abundance_coef of raw.")
  }
  raw$var_raw <- NA_real_
  raw$se_raw <- NA_real_
  if ("var_raw_marginal" %in% names(raw)) {
    use <- is.finite(raw$var_raw_marginal) & raw$var_raw_marginal > 0
    raw$var_raw[use] <- raw$var_raw_marginal[use]
    raw$se_raw[use] <- raw$se_raw_marginal[use]
  }
  raw_fit <- raw[is.finite(raw$beta_raw) & is.finite(raw$var_raw) & raw$var_raw > 0, , drop = FALSE]
  if (nrow(raw_fit) > 0L && (is.null(variance_inputs) || is.null(variance_inputs$V) ||
                             is.null(variance_inputs$G))) {
    stop("score_ancombc3_da requires variance_inputs (fit$da_variance_inputs).")
  }
  if (nrow(raw_fit) == 0L) {
    raw$p_value <- NA_real_
    raw$q_value <- 1
    raw$beta_corrected <- NA_real_
    raw$se_corrected <- NA_real_
    raw$se_scaled <- NA_real_
    raw$df_t <- NA_real_
    raw$delta <- NA_real_
    raw$var_delta <- NA_real_
    raw$in_null_set <- FALSE
    raw$var_primary <- NA_real_
    raw$statistic <- NA_real_
    raw$statistic_normal_equivalent <- NA_real_
    raw$adjustment_factor <- NA_real_
    raw$adjustment_factor_hat <- NA_real_
    raw$adjustment_factor_components <- NA_integer_
    return(raw)
  }
  shift <- .ac3_da_common_shift(
    raw_fit,
    eligible_taxa_by_term = eligible_taxa_by_term,
    min_count = maximal_overlap_min_count,
    bandwidths = maximal_overlap_bandwidths,
    min_gain = shift_selection_min_gain,
    min_separation = shift_selection_min_separation
  )
  null_set <- attr(shift, "null_set")
  correction_weights <- attr(shift, "correction_weights")
  scored <- do.call(rbind, lapply(split(raw_fit, raw_fit$term), function(g) {
    term <- g$term[1L]
    shift_term <- shift[shift$term == term, , drop = FALSE]
    in_null_set <- g$taxon %in% (null_set[[term]] %||% character(0))
    nu_term <- correction_weights[[term]] %||% numeric(0)
    nu <- unname(nu_term[match(g$taxon, names(nu_term))])
    nu[is.na(nu)] <- 0
    cv <- .ac3_da_corrected_variance(
      var_model = g$var_raw_model,
      G = variance_inputs$G[[term]][match(g$taxon, variance_inputs$taxa), , drop = FALSE],
      V = variance_inputs$V,
      nu = nu,
      A = if (is.null(variance_inputs$A[[term]])) NULL else
        variance_inputs$A[[term]][, match(g$taxon, variance_inputs$taxa), drop = FALSE]
    )
    # Steps 1 and 2: corrected coefficient, its working variance, and Z_j.
    z <- g
    z$delta <- as.numeric(shift_term$delta)
    z$var_delta <- as.numeric(cv$var_delta)
    z$beta_corrected <- z$beta_raw - z$delta
    z$var_primary <- as.numeric(cv$var_primary)
    z$se_corrected <- sqrt(pmax(z$var_primary, 0))
    z_da <- z$beta_corrected / pmax(z$se_corrected, .Machine$double.eps)
    z_da[!(is.finite(z$var_primary) & z$var_primary > 0)] <- NA_real_
    z$statistic <- z_da
    z$p_value <- NA_real_
    z$q_value <- NA_real_
    z$method <- "ancombc3_da"
    z$cov_delta <- cv$cov_delta
    z$correction_weight <- nu
    # Steps 3 and 4: degrees of freedom max(n_positive - n_abundance_coef, 1)
    # and the normal-equivalent statistic, from the log tail probability so
    # that small tail probabilities are kept.
    n_df <- if ("df_count" %in% names(g)) g$df_count else g$n_positive
    df_t <- pmax(as.numeric(n_df) - as.numeric(g$n_abundance_coef), 1)
    log_tail <- stats::pt(-abs(z_da), df = df_t, log.p = TRUE)
    z_tilde <- sign(z_da) * stats::qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)
    z$df_t <- df_t
    # Step 5: empirical null N(0, omega_T^2) of the two-group model fitted to
    # the finite normal-equivalent statistics of the eligible taxa of the
    # covariate, eligible_taxa_by_term[[term]] (every taxon of g when no set
    # is given); omega_T = 1 leaves the working variance unchanged.
    z$statistic_normal_equivalent <- z_tilde
    F_taxa <- if (is.null(eligible_taxa_by_term)) NULL else eligible_taxa_by_term[[term]]
    in_F_T <- if (is.null(F_taxa)) rep(TRUE, nrow(z)) else z$taxon %in% F_taxa
    en <- .ac3_da_empirical_null(z_tilde[in_F_T & is.finite(z_tilde)])
    z$adjustment_factor <- en$factor
    z$adjustment_factor_hat <- en$omega_hat
    z$adjustment_factor_components <- en$n_components
    # Step 6: scaled statistic and P values, adjusted within the covariate.
    z$se_scaled <- z$adjustment_factor * z$se_corrected
    z$statistic <- z_da / z$adjustment_factor
    z$p_value <- pmin(2 * exp(stats::pt(-abs(z$statistic), df = df_t, log.p = TRUE)), 1)
    z$q_value <- stats::p.adjust(z$p_value, method = p_adj_method)
    z$in_null_set <- in_null_set
    z$maximal_overlap_center <- shift_term$maximal_overlap_center
    z$n_eligible_coefficients <- shift_term$n_eligible_coefficients
    z$maximal_overlap_n_null_set <- shift_term$maximal_overlap_n_null_set
    z$maximal_overlap_bandwidth <- shift_term$maximal_overlap_bandwidth
    z$maximal_overlap_min_count <- shift_term$maximal_overlap_min_count
    z$maximal_overlap_status <- shift_term$maximal_overlap_status
    z$maximal_overlap_delta <- shift_term$maximal_overlap_delta
    z$maximal_overlap_var_delta <- shift_term$maximal_overlap_var_delta
    z$shift_selection_loglik_gain <- shift_term$shift_selection_loglik_gain
    z$gaussian_mixture_delta <- shift_term$gaussian_mixture_delta
    z$gaussian_mixture_selected <- shift_term$gaussian_mixture_selected
    z$shift_selection_status <- shift_term$shift_selection_status
    z
  }))
  rownames(scored) <- NULL
  scored
}
