# =============================================================================
# ANCOM-BC3 Poisson-log-normal count distribution of the second stage
# =============================================================================
#
# Code names of the symbols of the manuscript: nu_ij = log mu_ij is `log_mu`;
# zeta_ij = logit pi_ij is `logit_pi`; derivative fields dl_mu, dl_pi,
# d2l_mumu, d2l_mupi, d2l_pipi.
# kappa_j = log sigma_j is `kappa`; the regression block
# omega_j = (alpha_j, beta_j, a_j, xi_j) is `omega`; pi_i, p_i, and f_0,i are
# `pi_i`, `p0`, and `f0`; the normal mean c = log mu - sigma^2 / 2 of the log
# intensity is `center`.
#
# At-risk count distribution of the taxon-wise fits of fit_ancombc3_da(),
# parameterized by its mean mu and the taxon log-scale standard deviation sigma:
#   Y | Lambda ~ Poisson(Lambda),  log Lambda ~ N(c, sigma^2),  c = log mu - sigma^2 / 2,
# so that E(Y) = mu, with mass function
#   f(y; mu, sigma) = int exp(y u - e^u - log y!) phi((u - c) / sigma) / sigma du.
# For fixed sigma, dc = d nu with nu = log mu, and the derivatives in nu are
# moments under the conditional distribution of u given y. For y > 0 they are
# moments of the standardized log intensity z = (u - c) / sigma,
#   d log f / d nu    = E(z | y) / sigma,
#   d2 log f / d nu2  = {Var(z | y) - 1} / sigma^2,
# which are of the size of the derivatives for all counts; they equal
# y - E(Lambda | y) and -E(Lambda | y) + Var(Lambda | y), Lambda = e^u, whose
# terms are of size y and y^2. For y = 0 they are moments of Lambda,
#   d f0 / d nu       = -f0 E(Lambda | 0),
#   d2 f0 / d nu2     =  f0 {E(Lambda^2 | 0) - E(Lambda | 0)},
# where f0 = f(0; mu, sigma) is the sampling-zero probability. The integration
# routines below take the normal mean c (center).
#
# Numerical integration. Let A = y for y > 0 and A = e^c for y = 0, and let
# v = u - log A, so that Lambda = e^u = A e^v, and d = c - log A (d = 0 for
# y = 0). With z = (v - d) / sigma, the log integrand is, up to a constant,
#   L(v) = y {v - expm1(v)} - B expm1(v) - (v - d)^2 / (2 sigma^2),
# where B = e^c for y = 0 and B = 0 for y > 0, and the constant is
#   y log y - y - log y! - B,
# evaluated as dpois(y, y, log = TRUE) - B, so that every term of L is of the
# size of L itself for all counts. L is strictly concave with a unique mode v_hat. The interval
# [v_lo, v_hi] = {v : L(v) >= L(v_hat) - 30} is found by Newton-Raphson on each
# side of the mode, and the integral over it is evaluated by a composite
# Gauss-Legendre rule (12 panels of 12 nodes) in v; the integral outside the
# interval is at most e^{-30} times the total.

.ac3_pln_level_nats <- 30

# Composite Gauss-Legendre rule on [0, 1]: P panels of m nodes each, as positions
# in [0, 1] and log weights (summing to one).
.ac3_pln_gl_nodes <- local({
  cache <- new.env(parent = emptyenv())
  function(m = 12L, P = 12L) {
    key <- paste(m, P)
    hit <- cache[[key]]
    if (!is.null(hit)) return(hit)
    i <- seq_len(m - 1L)
    jm <- matrix(0, m, m)
    b <- i / sqrt(4 * i^2 - 1)
    jm[cbind(i, i + 1L)] <- b
    jm[cbind(i + 1L, i)] <- b
    e <- eigen(jm, symmetric = TRUE)
    nodes <- e$values
    w <- 2 * e$vectors[1L, ]^2
    pos <- as.numeric(outer((nodes + 1) / 2, seq_len(P) - 1L, `+`)) / P
    logw <- rep(log(w / (2 * P)), times = P)
    out <- list(pos = pos, logw = logw)
    cache[[key]] <- out
    out
  }
})

# Reference intensity of the integration variable v = u - log A: B = e^c for
# y = 0 and B = 0 for y > 0, A = y + B, d = c - log A, and the constant of the
# log integrand.
.ac3_pln_reference <- function(center, y) {
  zero <- !(y > 0)
  B <- ifelse(zero, exp(center), 0)
  list(B = B, A = y + B, d = ifelse(zero, 0, center - log(pmax(y, 1e-300))),
       const = stats::dpois(y, y, log = TRUE) - B)
}

# Derivative of L in v, L'(v) = -y expm1(v) - B e^v - (v - d) / sigma^2, from
# em1 = expm1(v).
.ac3_pln_score_v <- function(v, em1, y, ref, sigma) {
  -y * em1 - ref$B * (1 + em1) - (v - ref$d) / sigma^2
}

# Mode of L and the curvature in z, -d2 L / dz2 = sigma^2 Lambda_hat + 1 with
# Lambda_hat = A exp(v_hat). L' is decreasing and concave, so every
# Newton-Raphson iterate is an upper bound of the root of L', and the iterates
# decrease monotonically to it. For y > 0 the start is also an upper bound:
# z = min{sigma y, (log y - c) / sigma} when log y > c, and z = 0 otherwise.
# For y = 0 the start z = -log(1 + sigma^2 exp(c)) / sigma is a lower bound,
# and the first iterate moves past the root to an upper bound. The iteration
# stops when the largest |dL / dz| = sigma |L'(v)| is below 1e-11. Returns the
# mode as v and as z = (v - d) / sigma, and the curvature.
.ac3_pln_mode <- function(center, y, sigma, max_iter = 60L) {
  ref <- .ac3_pln_reference(center, y)
  v <- ifelse(y > 0, ifelse(ref$d < 0, pmin(0, ref$d + sigma^2 * y), ref$d),
              -log1p(sigma^2 * ref$B))
  for (it in seq_len(max_iter)) {
    em1 <- expm1(v)
    g <- .ac3_pln_score_v(v, em1, y, ref, sigma)
    if (max(abs(sigma * g)) < 1e-11) break
    v <- v + g / (ref$A * (1 + em1) + 1 / sigma^2)
  }
  list(z = (v - ref$d) / sigma, curv = sigma^2 * ref$A * (1 + expm1(v)) + 1, v = v)
}

# L(v) without its constant.
.ac3_pln_L <- function(v, em1, y, ref, sigma) {
  y * (v - em1) - ref$B * em1 - 0.5 * ((v - ref$d) / sigma)^2
}

# Endpoint of the integration interval on one side of the mode: the root of
# L(v) = L(v_hat) - depth, by Newton-Raphson from the normal approximation of
# the endpoint. side = -1 gives v_lo and side = +1 gives v_hi.
.ac3_pln_level_boundary <- function(y, ref, sigma, vhat, Lhat, curv, depth, side, max_iter = 60L) {
  target <- Lhat - depth
  v <- vhat + side * sigma * sqrt(2 * depth / curv)
  for (it in seq_len(max_iter)) {
    em1 <- expm1(v)
    h <- .ac3_pln_L(v, em1, y, ref, sigma) - target
    if (max(abs(h)) < 1e-9) break
    v_new <- v - h / .ac3_pln_score_v(v, em1, y, ref, sigma)
    # The iterate stays strictly on the chosen side of the mode.
    v <- if (side < 0) pmin(v_new, vhat - sigma * 1e-6 / sqrt(curv)) else pmax(v_new, vhat + sigma * 1e-6 / sqrt(curv))
  }
  v
}

# log f(y; mu, sigma) and the conditional moments given y, for vectors center
# and y of a common length and a scalar sigma: E1 = E(Lambda | y),
# E2 = E(Lambda^2 | y), Ez = E(z | y), and Vz = Var(z | y), where
# z = (log Lambda - c) / sigma is the standardized log intensity. With
# moments = FALSE only log f is computed (the moments are NULL). With
# by_count = TRUE, the R code computes E1 and E2 for the zero counts and Ez
# and Vz for the positive counts only (NA elsewhere), the moments used by the
# cell derivatives; the compiled routine returns every moment.
# engine = "compiled" evaluates the same quantities with the compiled routine
# (.ac3_engine_pln_obs(), ancombc3_engine.R); any other value, including NULL,
# uses the R code below. The log standard normal density of the integrand is
# written out as -(log sqrt(2 pi) + z^2 / 2), and the row maxima are taken with
# max.col(); both give the values of dnorm(log = TRUE) and pmax() exactly.
.ac3_pln_log_sqrt_2pi <- 0.918938533204672741780329736406

.ac3_pln_obs <- function(center, y, sigma, moments = TRUE, engine = "r", by_count = FALSE) {
  if (identical(engine, "compiled")) return(.ac3_engine_pln_obs(center, y, sigma, moments))
  gl <- .ac3_pln_gl_nodes()
  n <- length(center)
  ref <- .ac3_pln_reference(center, y)
  md <- .ac3_pln_mode(center, y, sigma)
  Lhat <- .ac3_pln_L(md$v, expm1(md$v), y, ref, sigma)
  depth <- .ac3_pln_level_nats
  v_lo <- .ac3_pln_level_boundary(y, ref, sigma, md$v, Lhat, md$curv, depth, -1)
  v_hi <- .ac3_pln_level_boundary(y, ref, sigma, md$v, Lhat, md$curv, depth, +1)
  width <- v_hi - v_lo
  v <- v_lo + outer(width, gl$pos)
  em1 <- expm1(v)
  zu <- (v - ref$d) / sigma
  lg <- y * (v - em1) - ref$B * em1 + ref$const + (-(.ac3_pln_log_sqrt_2pi + 0.5 * zu * zu)) - log(sigma) +
    log(width) + rep(gl$logw, each = n)
  m <- lg[cbind(seq_len(n), max.col(lg, ties.method = "first"))]
  e <- exp(lg - m)
  s0 <- rowSums(e)
  logf <- m + log(s0)
  if (!moments) return(list(logf = logf, E1 = NULL, E2 = NULL, Ez = NULL, Vz = NULL))
  if (!by_count) {
    lam <- ref$A * (1 + em1)
    E1 <- rowSums(e * lam) / s0
    E2 <- rowSums(e * lam^2) / s0
    Ez <- rowSums(e * zu) / s0
    Vz <- rowSums(e * (zu - Ez)^2) / s0
    return(list(logf = logf, E1 = E1, E2 = E2, Ez = Ez, Vz = Vz))
  }
  # The cell derivatives use E1 and E2 on the zero counts and Ez and Vz on the
  # positive counts.
  E1 <- E2 <- Ez <- Vz <- rep(NA_real_, n)
  zero <- which(!(y > 0))
  if (length(zero) > 0L) {
    e_zero <- e[zero, , drop = FALSE]
    lam <- ref$A[zero] * (1 + em1[zero, , drop = FALSE])
    E1[zero] <- rowSums(e_zero * lam) / s0[zero]
    E2[zero] <- rowSums(e_zero * lam^2) / s0[zero]
  }
  positive <- which(y > 0)
  if (length(positive) > 0L) {
    e_pos <- e[positive, , drop = FALSE]
    zu_pos <- zu[positive, , drop = FALSE]
    Ez[positive] <- rowSums(e_pos * zu_pos) / s0[positive]
    Vz[positive] <- rowSums(e_pos * (zu_pos - Ez[positive])^2) / s0[positive]
  }
  list(logf = logf, E1 = E1, E2 = E2, Ez = Ez, Vz = Vz)
}

# Sampling-zero probability f0 = f(0; mu, sigma) from log_mu, truncated to
# [-30, 30]; engine as in .ac3_pln_obs().
.ac3_pln_zero_prob <- function(log_mu, sigma, engine = "r") {
  log_mu <- pmin(pmax(log_mu, -30), 30)
  exp(.ac3_pln_obs(log_mu - sigma^2 / 2, rep(0, length(log_mu)), sigma, moments = FALSE,
                   engine = engine)$logf)
}

# Elementwise restriction of v to [lo, hi] for scalar bounds; equal to
# pmin(pmax(v, lo), hi), with missing values kept.
.ac3_pln_clip <- function(v, lo, hi) {
  v[v < lo] <- lo
  v[v > hi] <- hi
  v
}

# Presence quantities of a vector of presence predictors before truncation
# (logit_pi_raw): the truncated predictor logit_pi (to [-30, 30]), pi_i and
# 1 - pi_i (one_minus_pi), and the indicator in_pi that the predictor is not
# truncated.
.ac3_pln_presence <- function(logit_pi_raw) {
  logit_pi <- .ac3_pln_clip(logit_pi_raw, -30, 30)
  list(logit_pi = logit_pi, pi_i = plogis(logit_pi), one_minus_pi = plogis(-logit_pi),
       in_pi = (logit_pi_raw > -30 & logit_pi_raw < 30) * 1)
}

# Positions of the positive counts (pos) and of the zero counts (zero) for the
# indicator pos = (y > 0).
.ac3_pln_count_index <- function(pos) {
  list(pos = which(pos), zero = which(!pos))
}

# Cell log-likelihood l_i of the zero-inflated Poisson-log-normal model and its
# first and second derivatives in the predictors log mu and logit pi, for vectors
# of the predictors before truncation (log_mu_raw, logit_pi_raw), the counts y,
# the indicator pos = (y > 0), and the taxon's sigma. The predictors are
# truncated to [-30, 30], and a derivative is set to zero where a predictor is
# truncated or where the probability of the cell is restricted to
# [1e-12, 1 - 1e-12]. need (0, 1, or 2) is the highest derivative order
# computed, and engine selects the evaluation of the quadrature
# (.ac3_pln_obs()). presence (.ac3_pln_presence() of logit_pi_raw) and index
# (.ac3_pln_count_index() of pos) are computed here when they are NULL.
# Returns ll, dl_mu, dl_pi, d2l_mumu, d2l_pipi, d2l_mupi, mu, pi_i, and f0 (the
# sampling-zero probability on the zero cells, 0 elsewhere).
.ac3_pln_observation_derivatives <- function(log_mu_raw, logit_pi_raw, y, pos, sigma, need = 2L,
                                             engine = "r", presence = NULL, index = NULL) {
  if (is.null(presence)) presence <- .ac3_pln_presence(logit_pi_raw)
  if (is.null(index)) index <- .ac3_pln_count_index(pos)
  ip <- index$pos
  iz <- index$zero
  log_mu <- .ac3_pln_clip(log_mu_raw, -30, 30)
  mu <- exp(log_mu)
  pi_i <- presence$pi_i
  pl <- .ac3_pln_obs(log_mu - sigma^2 / 2, y, sigma, moments = need >= 1L, engine = engine,
                     by_count = TRUE)
  zeros <- y * 0
  f0 <- zeros
  fz <- exp(pl$logf[iz])
  f0[iz] <- fz
  ll <- zeros
  pi_pos <- pi_i[ip]
  ll[ip] <- log(.ac3_pln_clip(pi_pos, 1e-12, 1 - 1e-12)) + pl$logf[ip]
  # p_i = Pr(Y = 0) = (1 - pi_i) + pi_i f0: a structural zero or a sampling
  # zero of an at-risk taxon.
  pz <- pi_i[iz]
  p0_raw <- (1 - pz) + pz * fz
  p0 <- .ac3_pln_clip(p0_raw, 1e-12, 1 - 1e-12)
  ll[iz] <- log(p0)
  if (need < 1L) return(list(ll = ll, mu = mu, pi_i = pi_i, f0 = f0))
  in_mu <- (log_mu_raw > -30 & log_mu_raw < 30) * 1
  in_pi <- presence$in_pi
  one_minus_pi <- presence$one_minus_pi
  dl_mu <- dl_pi <- d2l_mumu <- d2l_pipi <- d2l_mupi <- zeros
  # Positive counts.
  inside_pos <- (pi_pos > 1e-12 & pi_pos < 1 - 1e-12) * 1
  om_pos <- one_minus_pi[ip]
  dl_mu[ip] <- pl$Ez[ip] / sigma
  dl_pi[ip] <- om_pos * inside_pos
  # Zeros.
  oz <- one_minus_pi[iz]
  inside <- (p0_raw > 1e-12 & p0_raw < 1 - 1e-12) * 1
  E1z <- pl$E1[iz]
  df0_mu <- -fz * E1z
  dp0_mu <- pz * df0_mu
  dp0_pi <- -pz * oz * (1 - fz)
  dl_mu_zero <- dp0_mu / p0
  dl_pi_zero <- dp0_pi / p0
  dl_mu[iz] <- dl_mu_zero * inside
  dl_pi[iz] <- dl_pi_zero * inside
  if (need >= 2L) {
    d2l_mumu[ip] <- (pl$Vz[ip] - 1) / sigma^2
    d2l_pipi[ip] <- -pi_pos * om_pos * inside_pos
    d2f0_mumu <- fz * (pl$E2[iz] - E1z)
    d2p0_mumu <- pz * d2f0_mumu
    d2p0_pipi <- -(1 - fz) * pz * oz * (1 - 2 * pz)
    d2p0_mupi <- pz * oz * df0_mu
    d2l_mumu[iz] <- (d2p0_mumu / p0 - dl_mu_zero^2) * inside
    d2l_pipi[iz] <- (d2p0_pipi / p0 - dl_pi_zero^2) * inside
    d2l_mupi[iz] <- (d2p0_mupi / p0 - dl_mu_zero * dl_pi_zero) * inside
  }
  list(ll = ll, dl_mu = dl_mu * in_mu, dl_pi = dl_pi * in_pi, d2l_mumu = d2l_mumu * in_mu,
       d2l_pipi = d2l_pipi * in_pi, d2l_mupi = d2l_mupi * in_mu * in_pi, mu = mu, pi_i = pi_i,
       f0 = f0)
}

# Bounds of sigma_j in the taxon-wise fits; kappa_j = log sigma_j lies in
# [log 0.02, log 4].
.ac3_pln_sigma_bounds <- c(0.02, 4)

# Starting value of sigma_j: the first-stage estimate sigma_j^VAE, restricted
# to [0.05, 3].
.ac3_pln_sigma_start <- function(sigma_vae) {
  pmin(pmax(as.numeric(sigma_vae), 0.05), 3)
}

# Adjusted taxon objective for the estimation of kappa_j = log sigma_j: the
# Cox-Reid adjusted profile log-likelihood
#   l^A(omega, kappa) = l^pen(omega, kappa) - (1 / 2) log det H(omega, kappa),
# where H = -d2 l^pen / d omega2 is the penalized observed information of the
# regression block omega (.ac3_da_taxon_objective()). The adjustment removes
# the downward bias of the maximum-likelihood estimate of sigma_j that comes
# from estimating omega_j from the same counts; in the normal linear model it
# reduces to the REML estimate of the variance. .ac3_pln_kappa_update()
# maximizes it over kappa with omega re-estimated at every kappa (the adjusted
# profile likelihood), and the same quantity at the fitted omega_j is the
# objective of a taxon-wise fit (.ac3_da_fit_record()) and of L
# (.ac3_joint_L()). -Inf when H is not positive definite. chol_H is the
# Cholesky factor of H when it is available (the element chol of a
# Newton-Raphson result, .ac3_da_newton()), and is computed here otherwise.
.ac3_pln_adjusted_from_ev <- function(ev, chol_H = NULL) {
  if (is.null(ev$hess) || !is.finite(ev$value) || any(!is.finite(ev$hess))) return(-Inf)
  if (is.null(chol_H)) chol_H <- tryCatch(chol(ev$hess), error = function(e) NULL)
  if (is.null(chol_H)) return(-Inf)
  ev$value - sum(log(diag(chol_H)))
}
.ac3_pln_adjusted_value <- function(omega, dat) {
  .ac3_pln_adjusted_from_ev(.ac3_da_taxon_objective(omega, dat, grad = TRUE, hess = TRUE))
}

# Update of kappa_j = log sigma_j in a taxon-wise fit: safeguarded
# Newton-Raphson steps in kappa on the adjusted profile objective
#   l^A(kappa) = l^pen(omega_hat(kappa), kappa) - (1 / 2) log det H(omega_hat(kappa), kappa),
# where omega_hat(kappa) is the Newton-Raphson maximizer of l^pen in omega_j at
# that kappa, started from the current omega_j (.ac3_pln_adjusted_from_ev());
# a kappa at which the Newton-Raphson run is not regular and converged is
# excluded (value -Inf). Each cycle evaluates l^A at the current kappa and at
# kappa -/+ 0.1, forms the Newton-Raphson step from these central differences
# (a step of 0.5 in the ascent direction when the curvature is not negative),
# bounds the step by 1 and the result by [log 0.02, log 4], and accepts the
# first of the full step, the half step, and the better of kappa -/+ 0.1 that
# increases l^A over its value at the current kappa (ascent rule); when none
# does, the current kappa and omega are kept. An accepted kappa has a regular
# converged Newton-Raphson update of omega_j, which becomes the current omega_j.
# Cycles are repeated until the accepted change in kappa is at most tol_kappa
# or after max_cycles cycles: 2 in the warm refits of the outer alternation,
# which supplies the remaining iterations, and 8 in the multistart fits.
# prior = list(m, tau2) adds the Gaussian log prior -(kappa - m)^2 / (2 tau2)
# of pass 2 (NULL in the first pass). res is a regular converged
# Newton-Raphson result at dat$sigma. Returns the final result (res),
# sigma = exp(kappa), the number of cycles, and kappa_info, the observed
# information about kappa with omega maximized out, I_sigma of the manuscript,
# at the optimum (NA when not finite).
# Evaluations at the current kappa. res holds the evaluation ev of l^pen, its
# gradient, and H at (omega, sigma_ev), where sigma_ev is the sigma of the
# Newton-Raphson run that produced res. A Newton-Raphson run at the current
# kappa would start at its own converged value, evaluate ev again, and stop
# without a step; its result and the value of l^pen at the current kappa are
# therefore taken from res without a new evaluation.
# With information = FALSE, kappa_info is not computed (NA); the caller
# computes it afterwards for the fits that are retained
# (.ac3_pln_kappa_information()).
.ac3_pln_kappa_update <- function(dat, res, max_cycles = 8L, tol_kappa = 1e-3, prior = NULL,
                                  information = TRUE) {
  omega <- res$omega
  sigma <- dat$sigma
  sigma_ev <- dat$sigma
  bounds <- log(.ac3_pln_sigma_bounds)
  cycles <- 0L
  penalty <- function(kappa) if (is.null(prior)) 0 else -(kappa - prior$m)^2 / (2 * prior$tau2)
  h <- 0.1
  clip <- function(k) min(max(k, bounds[1L]), bounds[2L])
  for (cyc in seq_len(max_cycles)) {
    cache <- new.env(parent = emptyenv())
    f <- function(kappa) {
      key <- format(kappa, digits = 15L)
      hit <- cache[[key]]
      if (!is.null(hit)) return(hit$value)
      dat$sigma <- exp(kappa)
      r <- .ac3_da_newton(dat, omega)
      v <- if (isTRUE(r$converged) && isTRUE(r$pd)) .ac3_pln_adjusted_from_ev(r$ev, r$chol) + penalty(kappa) else -Inf
      if (!is.finite(v)) v <- -1e300
      cache[[key]] <- list(value = v, res = r, sigma = dat$sigma)
      v
    }
    kappa_cur <- log(sigma)
    if (identical(sigma, sigma_ev) && isTRUE(res$converged) && isTRUE(res$pd)) {
      v_cur <- .ac3_pln_adjusted_from_ev(res$ev, res$chol) + penalty(kappa_cur)
      if (!is.finite(v_cur)) v_cur <- -1e300
      cache[[format(kappa_cur, digits = 15L)]] <- list(value = v_cur, res = res, sigma = sigma_ev)
    }
    f_cur <- f(kappa_cur)
    if (f_cur <= -1e300) break
    k_up <- clip(kappa_cur + h)
    k_dn <- clip(kappa_cur - h)
    f_up <- f(k_up)
    f_dn <- f(k_dn)
    if (f_up <= -1e300 && f_dn <= -1e300) break
    gradient <- (f_up - f_dn) / (k_up - k_dn)
    curvature <- (f_up - 2 * f_cur + f_dn) / h^2
    step <- if (is.finite(curvature) && curvature < 0) -gradient / curvature else sign(gradient) * 0.5
    if (!is.finite(step) || step == 0) break
    step <- sign(step) * min(abs(step), 1)
    candidates <- unique(c(clip(kappa_cur + step), clip(kappa_cur + step / 2),
                           if (f_up > f_cur) k_up, if (f_dn > f_cur) k_dn))
    accepted <- NA_real_
    for (kappa_try in candidates) {
      if (f(kappa_try) > f_cur + 1e-10 * (1 + abs(f_cur))) {
        accepted <- kappa_try
        break
      }
    }
    if (is.na(accepted)) {
      cycles <- cyc
      break
    }
    res_new <- cache[[format(accepted, digits = 15L)]]$res
    cycles <- cyc
    converged <- abs(accepted - kappa_cur) <= tol_kappa
    res <- res_new
    omega <- res$omega
    sigma <- exp(accepted)
    sigma_ev <- cache[[format(accepted, digits = 15L)]]$sigma
    if (converged) break
  }
  mid <- if (identical(sigma, sigma_ev)) res$ev else NULL
  info <- if (information) .ac3_pln_kappa_information(dat, omega, sigma, res$ev$hess, mid, res$chol) else NA_real_
  list(res = res, sigma = sigma, cycles = cycles, kappa_info = info)
}

# Observed information about kappa at (omega, sigma), with omega maximized out
# (a Schur complement): I_sigma = h_kk - c' H^-1 c, where h_kk = -d2 l / d kappa^2
# at fixed omega and c = d2 l / d omega d kappa are central differences
# (step 0.05 in kappa) of the penalized log-likelihood l^pen without the
# prior (the curvature of the adjustment term is of smaller order), and H is
# the analytic penalized observed information in omega. The sampling
# variance of kappa_hat is its reciprocal. dat holds the data of the taxon-wise
# fit (its sigma is not used), and mid is an evaluation of l^pen at
# (omega, sigma) when one is available (NULL otherwise). chol_H is the
# Cholesky factor of H when it is available, and is computed here otherwise.
# NA when the information is not finite.
.ac3_pln_kappa_information <- function(dat, omega, sigma, H, mid = NULL, chol_H = NULL) {
  at <- function(kappa, grad) {
    dat$sigma <- exp(kappa)
    .ac3_da_taxon_objective(omega, dat, grad = grad, hess = FALSE)
  }
  kappa <- log(sigma)
  h <- 0.05
  up <- at(kappa + h, TRUE)
  if (is.null(mid)) mid <- at(kappa, FALSE)
  dn <- at(kappa - h, TRUE)
  info <- NA_real_
  if (all(is.finite(c(up$value, mid$value, dn$value, up$grad, dn$grad)))) {
    h_kk <- -(up$value - 2 * mid$value + dn$value) / h^2
    c_j <- (up$grad - dn$grad) / (2 * h)
    if (is.null(chol_H)) {
      chol_H <- if (!is.null(H) && all(is.finite(H))) tryCatch(chol(H), error = function(e) NULL) else NULL
    }
    if (!is.null(chol_H)) info <- h_kk - sum(backsolve(chol_H, c_j, transpose = TRUE)^2)
  }
  if (is.finite(info)) info else NA_real_
}

# Cross-taxon penalty on kappa_j, computed from the taxon-wise fits retained by
# the first pass. Eligible taxa: a successful fit,
# 1.001 x 0.02 < sigma_hat_j < 0.999 x 4, at least min_positive positive
# counts, and a positive, finite information kappa_info. The center m
# is the median of log sigma_hat_j over the eligible taxa, and the variance
# tau2 is the squared median absolute deviation (constant 1.4826) of
# log sigma_hat_j minus the mean sampling variance 1 / kappa_info, bounded
# below by tau_min^2. Returns list(m, tau2, n_taxa, mad, sampling_sd), or NULL
# when fewer than min_taxa taxa are eligible.
.ac3_pln_cross_taxon_penalty <- function(fits, min_positive = 20L, min_taxa = 10L, tau_min = 0.10) {
  ok <- vapply(fits, function(f) isTRUE(f$ok) && !is.null(f$sigma) && is.finite(f$sigma), logical(1))
  if (!any(ok)) return(NULL)
  s <- vapply(fits[ok], function(f) as.numeric(f$sigma), numeric(1))
  npos <- vapply(fits[ok], function(f) as.numeric(f$n_positive), numeric(1))
  info <- vapply(fits[ok], function(f) if (is.null(f$kappa_info)) NA_real_ else as.numeric(f$kappa_info), numeric(1))
  b <- .ac3_pln_sigma_bounds
  use <- s > b[1L] * 1.001 & s < b[2L] * 0.999 & npos >= min_positive & is.finite(info) & info > 0
  if (sum(use) < min_taxa) return(NULL)
  kappa <- log(s[use])
  m <- stats::median(kappa)
  v_total <- stats::mad(kappa, constant = 1.4826)^2
  v_sampling <- mean(1 / info[use])
  tau2 <- max(v_total - v_sampling, tau_min^2)
  list(m = m, tau2 = tau2, n_taxa = sum(use), mad = sqrt(v_total), sampling_sd = sqrt(v_sampling))
}
