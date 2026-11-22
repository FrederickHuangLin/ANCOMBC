# =============================================================================
# ANCOM-BC3 second stage: alternating estimation of the taxon-wise parameters
# and the sample-specific bias (fit_ancombc3_da())
# =============================================================================
#
# The taxon-wise parameter vector is psi_j = (omega_j, kappa_j), with
# kappa_j = log sigma_j. For the fitting set B, the objective is
#   L({psi_j}, d) = sum_{j in B} l_j(psi_j, d),
# where l_j is the adjusted taxon objective l_j^A = l_j^pen - (1 / 2) log det H_j
# (.ac3_pln_adjusted_value(); l_j^pen is the penalized log-likelihood of
# ancombc3_da.R and H_j its penalized observed information in omega_j) in the
# first pass and l_j^(2) = l_j^A - (kappa_j - m)^2 / (2 tau2) in pass 2, with
# the cross-taxon penalty (m, tau2) held fixed. The updates of d and omega_j
# maximize the l_j^pen terms, whose maximizers are those of l_j^A to first
# order, and the update of kappa_j maximizes l_j^A with omega_j re-estimated
# at every kappa_j (the adjusted profile likelihood). Each outer iteration of a
# run of the alternation
#   1. updates d with the taxon-wise parameters fixed: Newton-Raphson on the
#      d-dependent terms of L within S = {v : 1'v = 0, Q'v = 0}, with the
#      numerical bound max_i |d_i| <= d_bound (.ac3_joint_d_update());
#   2. updates every psi_j at the new d: Newton-Raphson in omega_j from its
#      previous estimate, alternated with a one-dimensional maximization in
#      kappa_j (.ac3_da_taxon_fit() in mode "warm", .ac3_pln_kappa_update()).
# A run stops when max_i |d_i^(t+1) - d_i^(t)| <= tol_d and either
# |L^(t+1) - L^(t)| <= tol_L (1 + |L^(t+1)|) with a finite difference
# L^(t+1) - L^(t) or max_{j in B} |kappa_j^(t+1) - kappa_j^(t)| <= tol_kappa
# (the taxon-wise fits are then at a fixed point of the updates), or after
# max_outer iterations.
# The updates of d and omega_j are nondecreasing in l_j^pen and the update of
# kappa_j in l_j^A up to their tolerances, so L can change by small amounts in
# either direction; convergence to the global maximum is not implied.
# Estimation (.ac3_joint_estimate()):
#   first pass: one run from each initial fit (d = 0 and the projected log
#     library sizes) without the penalty on kappa_j; the run with the larger L
#     is retained;
#   pass 2: the cross-taxon penalty is computed from the retained first-pass
#     fits (.ac3_pln_cross_taxon_penalty()), and one run continues from the
#     retained fit with the penalty;
#   verification and final estimates: every taxon is refitted at the current d
#     from the full set of starting vectors; a refit replaces the current fit
#     when it increases the taxon's objective, and pass 2 resumes when a
#     replaced taxon belongs to B. The information matrix J_j of psi_j gives
#     the covariance of the final estimates (.ac3_joint_augment()).

# Stopping rule of the outer alternation and maximum number of verification
# passes, as in the Supplementary Table of numerical settings of the
# second-stage estimation: relative change in L 1e-8, change in d 1e-3,
# change in kappa_j 1e-3, at most 200 outer iterations, and at most 10
# verification passes.
.ac3_joint_tol_L <- 1e-8
.ac3_joint_tol_d <- 1e-3
.ac3_joint_tol_kappa <- 1e-3
.ac3_joint_max_outer <- 200L
.ac3_joint_max_verifications <- 10L

# Data of the taxon-wise objective of taxon j at the bias d and sigma, with
# the engine of the inputs.
.ac3_joint_taxon_dat <- function(inp, j, d, sigma) {
  .ac3_da_taxon_data(inp$y[, j], inp$x, inp$r, inp$abundance_offset[, j] + d,
                     inp$structural_offset[, j], inp$ridge, inp$sampling_zero_prob[, j],
                     sigma = sigma, engine = inp$engine)
}

# Gaussian log prior of kappa = log sigma, -(kappa - m)^2 / (2 tau2), up to a
# constant (0 without a prior).
.ac3_joint_log_prior <- function(sigma, prior) {
  if (is.null(prior)) return(0)
  -(log(sigma) - prior$m)^2 / (2 * prior$tau2)
}

# Largest change in kappa_j = log sigma_j over the fitting set between two
# sets of taxon-wise fits (Inf when no taxon has a finite sigma in both).
.ac3_joint_kappa_change <- function(old, new, fitting_idx) {
  ch <- vapply(fitting_idx, function(j) {
    a <- old[[j]]$sigma
    b <- new[[j]]$sigma
    if (isTRUE(old[[j]]$ok) && isTRUE(new[[j]]$ok) && is.finite(a) && is.finite(b)) abs(log(b) - log(a)) else NA_real_
  }, numeric(1))
  if (all(is.na(ch))) Inf else max(ch, na.rm = TRUE)
}

# Distributed evaluations of the alternation. With worker processes
# (inp$cl), the taxon terms of L and of the d-dependent terms are evaluated by
# the workers and accumulated here in the order of the fitting set, so that
# the sums are those of the evaluation without workers. The inputs that are
# constant during the fit (the counts, the designs, the first-stage offsets,
# the sampling-zero probabilities, the ridge, and the engine) are sent to the
# workers once
# (.ac3_joint_set_worker_inputs(), called by fit_ancombc3_da()), and the
# taxon-wise parameters of the fitting set before each evaluation of L and
# each update of d (.ac3_joint_export_fits()); both are stored in the
# environment .ac3_worker_state of the worker.
.ac3_joint_set_worker_inputs <- function(s) {
  assign("joint_worker_inputs", s, envir = .ac3_worker_state)
  NULL
}
.ac3_joint_set_worker_fits <- function(s) {
  assign("joint_worker_fits", s, envir = .ac3_worker_state)
  NULL
}
.ac3_joint_export_fits <- function(cl, fits, fitting_idx) {
  parallel::clusterCall(cl, .ac3_joint_set_worker_fits, list(
    idx = fitting_idx,
    omega = lapply(fits[fitting_idx], function(f) as.numeric(f$omega)),
    sigma = lapply(fits[fitting_idx], function(f) f$sigma)
  ))
  invisible(NULL)
}
.ac3_joint_worker_L <- function(k, d) {
  s <- get("joint_worker_fits", envir = .ac3_worker_state)
  .ac3_joint_L_taxon(get("joint_worker_inputs", envir = .ac3_worker_state), d, s$idx[k],
                     s$omega[[k]], s$sigma[[k]])
}
.ac3_joint_worker_d_terms <- function(k, d, need) {
  s <- get("joint_worker_fits", envir = .ac3_worker_state)
  .ac3_joint_d_terms_taxon(get("joint_worker_inputs", envir = .ac3_worker_state), d, s$idx[k],
                           s$omega[[k]], s$sigma[[k]], need)
}

# Adjusted objective l_j^A of taxon j at the bias d, omega, and sigma
# (.ac3_pln_adjusted_value()).
.ac3_joint_L_taxon <- function(inp, d, j, omega, sigma) {
  .ac3_pln_adjusted_value(omega, .ac3_joint_taxon_dat(inp, j, d, sigma))
}

# Objective L at the bias d for the taxon-wise parameters stored in fits: the
# sum over the fitting set of the adjusted taxon objectives l_j^A
# (.ac3_joint_L_taxon()) and the log prior of kappa_j; -Inf when a taxon of
# the fitting set has no successful fit or a nonfinite objective. current
# holds the taxa whose fit records were computed at d; for those with a
# positive definite information, the objective of the record is l_j^A at d
# (.ac3_da_fit_record()) and is used without a new evaluation.
.ac3_joint_L <- function(inp, d, fits, fitting_idx, prior, current = integer(0)) {
  if (!all(vapply(fits[fitting_idx], function(f) isTRUE(f$ok), logical(1)))) return(-Inf)
  stored <- fitting_idx %in% current & vapply(fits[fitting_idx], function(f) isTRUE(f$pd), logical(1))
  values <- vector("list", length(fitting_idx))
  values[stored] <- lapply(fits[fitting_idx[stored]], function(f) f$objective)
  todo <- which(!stored)
  if (length(todo) > 0L) {
    if (is.null(inp$cl)) {
      values[todo] <- lapply(todo, function(k) {
        f <- fits[[fitting_idx[k]]]
        .ac3_joint_L_taxon(inp, d, fitting_idx[k], as.numeric(f$omega), f$sigma)
      })
    } else {
      .ac3_joint_export_fits(inp$cl, fits, fitting_idx)
      values[todo] <- parallel::parLapply(inp$cl, todo, .ac3_joint_worker_L, d = d)
    }
  }
  tot <- 0
  for (k in seq_along(fitting_idx)) {
    f <- fits[[fitting_idx[k]]]
    v <- values[[k]]
    if (!is.finite(v)) return(-Inf)
    tot <- tot + v + .ac3_joint_log_prior(f$sigma, prior)
  }
  tot
}

# Terms of taxon j that depend on d at omega and sigma: the summed cell
# log-likelihood (value), dl_ij / dnu_ij (grad, for need >= 1), and
# d2l_ij / dnu_ij^2 (D, for need >= 2).
.ac3_joint_d_terms_taxon <- function(inp, d, j, omega, sigma, need) {
  qx <- ncol(inp$x)
  qr <- ncol(inp$r)
  log_mu <- omega[1L] + drop(inp$x %*% omega[1L + seq_len(qx)]) + inp$abundance_offset[, j] + d
  logit_pi <- omega[2L + qx] + drop(inp$r %*% omega[2L + qx + seq_len(qr)]) + inp$structural_offset[, j]
  y <- inp$y[, j]
  o <- .ac3_pln_observation_derivatives(log_mu, logit_pi, y, y > 0, sigma, need = need, engine = inp$engine)
  list(value = sum(o$ll), grad = if (need >= 1L) o$dl_mu, D = if (need >= 2L) o$d2l_mumu)
}

# Terms of L that depend on d at fixed taxon-wise parameters: the summed cell
# log-likelihood of the fitting set, its gradient g_i = sum_j dl_ij / dnu_ij,
# and the diagonal D_ii = sum_j d2l_ij / dnu_ij^2 of its Hessian (the
# penalties and the prior do not depend on d). With distributed = TRUE the
# taxon terms are evaluated by the workers, which hold the taxon-wise
# parameters (.ac3_joint_export_fits()).
.ac3_joint_d_terms <- function(inp, d, fits, fitting_idx, need = 2L, distributed = FALSE) {
  n <- nrow(inp$y)
  val <- 0
  g <- D <- numeric(n)
  terms <- NULL
  if (distributed) {
    terms <- parallel::parLapply(inp$cl, seq_along(fitting_idx), .ac3_joint_worker_d_terms,
                                 d = d, need = need)
  }
  for (k in seq_along(fitting_idx)) {
    j <- fitting_idx[k]
    o <- if (is.null(terms)) {
      .ac3_joint_d_terms_taxon(inp, d, j, as.numeric(fits[[j]]$omega), fits[[j]]$sigma, need)
    } else {
      terms[[k]]
    }
    val <- val + o$value
    if (need >= 1L) g <- g + o$grad
    if (need >= 2L) D <- D + o$D
  }
  list(value = val, grad = g, D = D)
}

# Update of d at fixed taxon-wise parameters: Newton-Raphson on the
# d-dependent terms of L restricted to S (orthonormal basis H0). The step is
# Delta = H0 A^-1 H0' g with A = -H0' diag(D) H0, whose eigenvalues are raised
# to at least 1e-8 times the largest before inversion, and se_i is the square
# root of the i-th diagonal element of H0 A^-1 H0'. d is not updated when the
# largest eigenvalue of A is not positive. When d + Delta violates the bound,
# the largest feasible step fraction gamma_step is computed by bisection, and
# when gamma_step < 1e-6 the full step is clipped and projected
# (.ac3_da_project_bias()). The step is halved, at most max_halvings times,
# until L does not decrease by more than 1e-10 (1 + |L|); after the last
# halving the candidate is used. L is the full objective, the d-dependent terms
# plus the penalties and the prior; the latter do not depend on d and are
# computed once as L_total, the value of L at the entry value of d, minus the
# d-dependent terms there (0 when L_total is not finite). Stops when
# max_i |Delta_i| / se_i <= tol, or after max_iter steps (d_tol and d_max_iter
# of fit_ancombc3_da(), defaults 0.01 and 50). With worker processes, the
# d-dependent terms are evaluated by the workers (.ac3_joint_d_terms()).
.ac3_joint_d_update <- function(inp, d, fits, fitting_idx, H0, L_total, tol = 0.01,
                                max_iter = 50L, max_halvings = 20L) {
  distributed <- !is.null(inp$cl)
  if (distributed) .ac3_joint_export_fits(inp$cl, fits, fitting_idx)
  cur <- .ac3_joint_d_terms(inp, d, fits, fitting_idx, distributed = distributed)
  L_offset <- L_total - cur$value
  if (!is.finite(L_offset)) L_offset <- 0
  it <- 0L
  for (it in seq_len(max_iter)) {
    A <- crossprod(H0 * (-cur$D), H0)
    eig_A <- eigen(0.5 * (A + t(A)), symmetric = TRUE)
    lambda_top <- max(eig_A$values)
    if (!is.finite(lambda_top) || lambda_top <= 0) break
    lambda <- pmax(eig_A$values, 1e-8 * lambda_top)
    H0U <- H0 %*% eig_A$vectors
    Delta <- drop(H0U %*% (crossprod(H0U, cur$grad) / lambda))
    se <- sqrt(pmax(rowSums(sweep(H0U^2, 2L, lambda, "/")), 1e-300))
    if (max(abs(Delta) / se) <= tol) break
    gamma_step <- 1
    clip <- FALSE
    if (max(abs(d + Delta)) > inp$d_bound) {
      lo <- 0
      hi <- 1
      for (k in 1:50) {
        mid <- (lo + hi) / 2
        if (max(abs(d + mid * Delta)) > inp$d_bound) hi <- mid else lo <- mid
      }
      gamma_step <- lo
      if (gamma_step < 1e-6) {
        gamma_step <- 1
        clip <- TRUE
      }
    }
    for (k in 0:max_halvings) {
      d_try <- d + gamma_step * Delta
      if (clip) {
        d_try <- .ac3_da_project_bias(d_try, x = inp$x, test_terms = inp$test_terms,
                                      d_bound = inp$d_bound)
      }
      try_terms <- .ac3_joint_d_terms(inp, d_try, fits, fitting_idx, distributed = distributed)
      if (is.finite(try_terms$value) &&
          try_terms$value >= cur$value - 1e-10 * (1 + abs(cur$value + L_offset))) break
      gamma_step <- gamma_step / 2
    }
    d <- d_try
    cur <- try_terms
  }
  list(d = d, iterations = it)
}

# Completion of the fit records of the taxa idx that were computed without
# their summaries (.ac3_da_complete_fit_summaries()), and of the corresponding
# columns of the matrix f0. d is the bias at which these records were fitted.
# With worker processes the records are completed by the workers.
.ac3_joint_worker_complete <- function(item, d) {
  inp <- get("joint_worker_inputs", envir = .ac3_worker_state)
  .ac3_da_complete_fit_summaries(item$f, inp$x, .ac3_joint_taxon_dat(inp, item$j, d, item$f$sigma))
}
.ac3_joint_complete_summaries <- function(inp, d, fit, idx = seq_along(fit$fits)) {
  idx <- idx[!vapply(fit$fits[idx], function(f) is.null(f$summaries_pending), logical(1))]
  if (length(idx) == 0L) return(fit)
  items <- lapply(idx, function(j) list(j = j, f = fit$fits[[j]]))
  done <- if (is.null(inp$cl)) {
    lapply(items, function(item) {
      .ac3_da_complete_fit_summaries(item$f, inp$x, .ac3_joint_taxon_dat(inp, item$j, d, item$f$sigma))
    })
  } else {
    parallel::parLapply(inp$cl, items, .ac3_joint_worker_complete, d = d)
  }
  for (k in seq_along(idx)) {
    fit$fits[[idx[k]]] <- done[[k]]
    if (isTRUE(done[[k]]$ok)) fit$f0[, idx[k]] <- done[[k]]$f0
  }
  fit
}

# Taxon-wise fits at d, each from its previous estimate with kappa updated
# under the prior; a taxon whose refit fails keeps its previous fit. The fit
# records are computed without their summaries, which the alternation does not
# use; .ac3_joint_path() completes the records of its final fits. d_fit is the
# bias at which the records of fit that lack their summaries were fitted; a
# previous fit that is kept is completed first. Returns the fits (fit) and the
# taxa whose records were computed at d (current).
.ac3_joint_refit <- function(inp, d, fit, prior, d_fit = NULL) {
  new <- .ac3_da_refit_taxa(inp, d, prev_fits = fit$fits, mode = "warm", prior = prior,
                            summaries = FALSE)
  failed <- which(!vapply(new$fits, function(f) isTRUE(f$ok), logical(1)) &
                    vapply(fit$fits, function(f) isTRUE(f$ok), logical(1)))
  if (length(failed)) {
    fit <- .ac3_joint_complete_summaries(inp, d_fit, fit, failed)
    new <- .ac3_da_replace_taxon_fits(new, fit, failed)
  }
  list(fit = new, current = setdiff(seq_along(new$fits), failed))
}

# One run of the alternation from the bias d0, with the taxon-wise fits fit0
# at d0 as starting values. The update of d stops by the rule of inp$d_tol and
# inp$d_max_iter. The fit records of the iterations are computed without f0
# and kappa_info (.ac3_joint_refit()); those of the returned fits are completed
# at the final d.
.ac3_joint_path <- function(inp, d0, fit0, fitting_idx, prior, H0, tol_L = .ac3_joint_tol_L,
                            tol_d = .ac3_joint_tol_d, tol_kappa = .ac3_joint_tol_kappa,
                            max_outer = .ac3_joint_max_outer) {
  d <- d0
  rf <- .ac3_joint_refit(inp, d, fit0, prior, d_fit = d)
  fit <- rf$fit
  L_cur <- .ac3_joint_L(inp, d, fit$fits, fitting_idx, prior, rf$current)
  rows <- list(data.frame(iteration = 0L, L = L_cur, change_L = NA_real_, max_change_d = NA_real_,
                          max_change_kappa = NA_real_, d_newton_steps = NA_integer_))
  converged <- FALSE
  it <- 0L
  for (it in seq_len(max_outer)) {
    du <- .ac3_joint_d_update(inp, d, fit$fits, fitting_idx, H0, L_total = L_cur,
                              tol = inp$d_tol, max_iter = inp$d_max_iter)
    rf <- .ac3_joint_refit(inp, du$d, fit, prior, d_fit = d)
    fit_new <- rf$fit
    L_new <- .ac3_joint_L(inp, du$d, fit_new$fits, fitting_idx, prior, rf$current)
    change_d <- max(abs(du$d - d))
    change_L <- L_new - L_cur
    change_kappa <- .ac3_joint_kappa_change(fit$fits, fit_new$fits, fitting_idx)
    rows[[length(rows) + 1L]] <- data.frame(iteration = it, L = L_new, change_L = change_L,
                                            max_change_d = change_d, max_change_kappa = change_kappa,
                                            d_newton_steps = du$iterations)
    d <- du$d
    fit <- fit_new
    L_cur <- L_new
    # The criterion on L requires a finite change of L: L is -Inf when a taxon
    # of the fitting set has no successful fit (.ac3_joint_L()), and a change
    # of L from or to -Inf, or between two values -Inf, does not satisfy it.
    L_settled <- is.finite(change_L) && abs(change_L) <= tol_L * (1 + abs(L_cur))
    if (change_d <= tol_d && (L_settled || change_kappa <= tol_kappa)) {
      converged <- TRUE
      break
    }
  }
  fit <- .ac3_joint_complete_summaries(inp, d, fit)
  list(d = d, fit = fit, L = L_cur, history = do.call(rbind, rows), converged = converged,
       iterations = it)
}

# Information matrix J_j of psi_j = (omega_j, kappa_j) at the fitted values,
# with the prior curvature 1 / tau2 in pass 2:
#   J_j = [ H_j      -c_j
#           -c_j'    h_kk + 1 / tau2 ],
# where H_j is the analytic penalized observed information in omega_j, and
# c_j = d2 l_j / d omega_j d kappa_j and h_kk = -d2 l_j / d kappa_j^2 are
# central differences with step h = 0.05 in kappa. When J_j is positive
# definite, its inverse is the covariance of the final estimates (vcov_psi),
# and the omega_j block of the inverse replaces vcov.
.ac3_joint_augment <- function(f, inp, j, d, prior, h = 0.05) {
  if (!isTRUE(f$ok) || is.null(f$hessian) || !is.finite(f$sigma)) return(f)
  dat <- .ac3_joint_taxon_dat(inp, j, d, f$sigma)
  omega <- as.numeric(f$omega)
  at <- function(kappa, grad) {
    dat$sigma <- exp(kappa)
    .ac3_da_taxon_objective(omega, dat, grad = grad, hess = FALSE)
  }
  kappa <- log(f$sigma)
  up <- at(kappa + h, TRUE)
  mid <- at(kappa, FALSE)
  dn <- at(kappa - h, TRUE)
  j_kk <- -(up$value - 2 * mid$value + dn$value) / h^2 + (if (is.null(prior)) 0 else 1 / prior$tau2)
  c_j <- (up$grad - dn$grad) / (2 * h)
  J <- rbind(cbind(f$hessian, -c_j), c(-c_j, j_kk))
  chol_J <- tryCatch(chol(0.5 * (J + t(J))), error = function(e) NULL)
  if (is.null(chol_J)) return(f)
  vcov_psi <- chol2inv(chol_J)
  k <- length(omega)
  f$vcov <- vcov_psi[seq_len(k), seq_len(k), drop = FALSE]
  f$vcov_psi <- vcov_psi
  f
}

# Estimation for fit_ancombc3_da(): the first pass from both initial fits, the
# cross-taxon penalty from the retained first-pass fits, pass 2 from the
# retained first-pass solution, and verification. Returns the final solution
# (best), the penalty (prior), the first-pass runs, their L, and the retained
# initial bias vector (chosen).
.ac3_joint_estimate <- function(inp, starts, fit0, fitting_idx, x, test_terms,
                                max_verifications = .ac3_joint_max_verifications) {
  H0 <- .ac3_da_constraint_basis(x, test_terms)
  runs <- lapply(names(starts), function(nm) {
    .ac3_joint_path(inp, starts[[nm]], fit0[[nm]], fitting_idx, prior = NULL, H0)
  })
  names(runs) <- names(starts)
  Ls <- vapply(runs, function(r) r$L, numeric(1))
  chosen <- if (!is.finite(Ls[["library_size"]]) ||
                (is.finite(Ls[["zero"]]) && Ls[["zero"]] >= Ls[["library_size"]])) {
    "zero"
  } else {
    "library_size"
  }
  prior <- .ac3_pln_cross_taxon_penalty(runs[[chosen]]$fit$fits)
  best <- if (is.null(prior)) {
    runs[[chosen]]
  } else {
    .ac3_joint_path(inp, runs[[chosen]]$d, runs[[chosen]]$fit, fitting_idx, prior, H0)
  }
  best$pass1_iterations <- runs[[chosen]]$iterations
  ## Verification: every taxon is refitted at d from the full set of starting
  ## vectors (including its current estimate), with kappa updated under the
  ## prior. A refit replaces the current fit when its objective (the adjusted
  ## objective l_j^A of the fit record) is larger by more than
  ## 1e-8 (1 + |current value|). When a
  ## replaced taxon belongs to B, the alternation resumes from the updated
  ## fits. At most max_verifications passes are made.
  n_ver <- 0L
  n_replaced <- 0L
  repeat {
    if (n_ver >= max_verifications) break
    ver <- .ac3_da_refit_taxa(inp, best$d, prev_fits = best$fit$fits, mode = "multistart",
                              prior = prior)
    n_ver <- n_ver + 1L
    obj <- function(f) if (isTRUE(f$ok)) f$objective + .ac3_joint_log_prior(f$sigma, prior) else -Inf
    cur_v <- vapply(best$fit$fits, obj, numeric(1))
    new_v <- vapply(ver$fits, obj, numeric(1))
    better <- which(is.finite(new_v) & (!is.finite(cur_v) | (new_v - cur_v) > 1e-8 * (1 + abs(cur_v))))
    if (length(better) == 0L) break
    n_replaced <- n_replaced + length(better)
    best$fit <- .ac3_da_replace_taxon_fits(best$fit, ver, better)
    if (!any(better %in% fitting_idx)) break
    resumed <- .ac3_joint_path(inp, best$d, best$fit, fitting_idx, prior, H0)
    resumed$pass1_iterations <- best$pass1_iterations
    resumed$iterations <- best$iterations + resumed$iterations
    best <- resumed
  }
  best$verifications <- n_ver
  best$n_replaced <- n_replaced
  fits <- lapply(seq_along(best$fit$fits), function(j) {
    .ac3_joint_augment(best$fit$fits[[j]], inp, j, best$d, prior)
  })
  names(fits) <- names(best$fit$fits)
  best$fit$fits <- fits
  ## Model variances of the raw coefficients from the covariance J_j^-1.
  raw <- best$fit$raw
  for (rr in seq_len(nrow(raw))) {
    f <- fits[[match(raw$taxon[rr], names(fits))]]
    idx <- 1L + match(raw$term[rr], colnames(x))
    v <- if (isTRUE(f$ok) && !is.null(f$vcov)) f$vcov[idx, idx] else NA_real_
    ok <- is.finite(v) && v > 0
    raw$var_raw_model[rr] <- if (ok) v else NA_real_
    raw$se_raw_model[rr] <- if (ok) sqrt(v) else NA_real_
  }
  best$fit$raw <- raw
  list(chosen = chosen, best = best, runs = runs, L = Ls, prior = prior)
}
