# DP test with the second-stage sampling-zero probability and the two-step
# variance, on a small simulated data set with a second-stage fit.

# Simulated data: n samples, p taxa, a continuous and a binary covariate in
# both designs, zero-inflated Poisson-log-normal counts with abundance effects
# on the binary covariate (taxa 1-4) and prevalence effects (taxa 5-8).
set.seed(20261021L)
n <- 80L
p <- 16L
x <- cbind(cont = as.numeric(scale(stats::rnorm(n))), binary = rep(0:1, length.out = n))
r <- x
alpha <- stats::rnorm(p, mean = 1, sd = 0.5)
beta <- matrix(0, p, 2L)
beta[1:4, 2L] <- -1.5
gamma <- matrix(0, p, 2L)
gamma[5:8, 2L] <- -1.5
log_mu <- matrix(alpha, n, p, byrow = TRUE) + x %*% t(beta)
lambda <- exp(log_mu + matrix(stats::rnorm(n * p, sd = 0.6), n, p))
present <- matrix(stats::rbinom(n * p, 1L, stats::plogis(2 + r %*% t(gamma))), n, p)
y <- matrix(stats::rpois(n * p, lambda), n, p) * present
colnames(y) <- paste0("T", seq_len(p))
rownames(y) <- paste0("S", seq_len(n))
mu_first_stage <- matrix(pmax(colMeans(y), 0.1), n, p, byrow = TRUE)
f0_first_stage <- (2 / (2 + mu_first_stage))^2
dimnames(f0_first_stage) <- dimnames(y)
test_terms <- c("cont", "binary")
X <- cbind(1, x)
D <- cbind("(Intercept)" = 1, r)
ridge <- 1e-3

second_stage_fit <- ac3_cached(function() {
  fit_ancombc3_da(
    y = y, x = x, r = r, sigma_vae = rep(0.6, p),
    abundance_offset = matrix(0, n, p), structural_offset = matrix(0, n, p),
    sampling_zero_prob = f0_first_stage, at_risk_prob = matrix(0.9, n, p),
    test_terms = test_terms, fitting_set = colnames(y),
    eligible_taxa_by_term = setNames(rep(list(colnames(y)), 2L), test_terms), n_starts = 2L
  )
})

# Taxa with a usable second-stage fit and at least 10 zero counts.
usable_taxa <- function(fit) {
  vapply(colnames(y), function(tx) {
    !is.null(.dp_second_stage_taxon(fit$fits[[tx]], X)) && sum(y[, tx] == 0) >= 10L
  }, logical(1))
}

# Fit record of a taxon with theta = (alpha, beta, kappa) changed by delta.
shift_record <- function(f, delta) {
  f$mu <- f$mu * exp(as.numeric(X %*% delta[seq_len(ncol(X))]))
  f$sigma <- f$sigma * exp(delta[ncol(X) + 1L])
  f
}

test_that("analytic derivatives of the second-stage f0 agree with central differences", {
  fit <- second_stage_fit()
  usable <- usable_taxa(fit)
  expect_true(sum(usable) >= 4L, label = "at least four taxa with a usable second-stage fit and zeros")
  check_taxa <- head(colnames(y)[usable], 4L)
  expect_true(all(vapply(check_taxa, function(tx) .dp_second_stage_taxon(fit$fits[[tx]], X)$kappa_included,
                         logical(1))),
              label = "theta includes kappa when vcov_psi is available")
  h <- 1e-4
  err_jac <- max(vapply(check_taxa, function(tx) {
    f <- fit$fits[[tx]]
    ss <- .dp_second_stage_taxon(f, X)
    fd <- vapply(seq_len(ncol(ss$jacobian)), function(k) {
      e <- replace(numeric(ncol(ss$jacobian)), k, h)
      up <- shift_record(f, e)
      dn <- shift_record(f, -e)
      (.ac3_pln_zero_prob(log(up$mu), up$sigma) - .ac3_pln_zero_prob(log(dn$mu), dn$sigma)) / (2 * h)
    }, numeric(n))
    max(abs(ss$jacobian - fd)) / max(abs(ss$jacobian))
  }, numeric(1)))
  expect_true(is.finite(err_jac) && err_jac <= 1e-6, label = "analytic d f0 / d theta matches finite differences")
  expect_true(isTRUE(all.equal(.dp_second_stage_taxon(fit$fits[[check_taxa[1L]]], X)$f0,
                               fit$fits[[check_taxa[1L]]]$f0, tolerance = 1e-12)),
              label = "second-stage f0 equals the f0 of the fit record")
})

test_that("score and sensitivity G of the zero-pattern model agree with central differences", {
  fit <- second_stage_fit()
  check_taxa <- head(colnames(y)[usable_taxa(fit)], 4L)
  # Score of the zero-pattern model with the ridge penalty, and the root of the
  # score equations by Newton-Raphson.
  score_dp <- function(Z, f0, coef) {
    eta <- as.numeric(D %*% coef)
    pi_dp <- stats::plogis(eta)
    q <- (1 - pi_dp) + pi_dp * f0
    as.numeric(crossprod(D, (1 - pi_dp) * (q - Z) / q)) - c(0, rep(ridge, ncol(D) - 1L)) * coef
  }
  root_dp <- function(Z, f0, start) {
    coef <- start
    for (it in 1:100) {
      s <- score_dp(Z, f0, coef)
      H <- .dp_zero_pattern_information(Z, f0, D, coef, ridge)$H
      step <- solve(H, s)
      coef <- coef + step
      if (max(abs(step)) < 1e-12) break
    }
    coef
  }
  num_grad <- function(fn, par, hh = 1e-5) {
    vapply(seq_along(par), function(k) {
      e <- replace(numeric(length(par)), k, hh)
      (fn(par + e) - fn(par - e)) / (2 * hh)
    }, numeric(1))
  }
  hg <- 1e-4
  res_G <- do.call(rbind, lapply(check_taxa, function(tx) {
    f <- fit$fits[[tx]]
    Z <- as.numeric(y[, tx] == 0)
    ss <- .dp_second_stage_taxon(f, X)
    nll0 <- .dp_zero_pattern_nll(Z, ss$f0, D, ridge = ridge, firth = FALSE)
    coef_test <- c(0.5, 0.2, -0.3)
    err_score <- max(abs(score_dp(Z, ss$f0, coef_test) + num_grad(nll0, coef_test))) /
      max(abs(score_dp(Z, ss$f0, coef_test)))
    coef_hat <- root_dp(Z, ss$f0, rep(0, ncol(D)))
    tv <- .dp_two_step_variance(Z, D, coef_hat, ridge, ss)
    G_fd <- vapply(seq_len(ncol(ss$jacobian)), function(k) {
      e <- replace(numeric(ncol(ss$jacobian)), k, hg)
      up <- .dp_second_stage_taxon(shift_record(f, e), X)$f0
      dn <- .dp_second_stage_taxon(shift_record(f, -e), X)$f0
      (root_dp(Z, up, coef_hat) - root_dp(Z, dn, coef_hat)) / (2 * hg)
    }, numeric(ncol(D)))
    # Sensitivity of the Firth-penalized estimator, -F_cc^-1 F_ct, from central
    # second differences (step 1e-3) of the penalized objective F(coef, theta)
    # at its minimizer.
    f0_at <- function(delta) {
      pmin(pmax(.dp_second_stage_taxon(shift_record(f, delta), X)$f0, 1e-8), 1 - 1e-8)
    }
    F <- function(cf, f0v) .dp_zero_pattern_nll(Z, f0v, D, ridge = ridge)(cf)
    f0_hat <- f0_at(rep(0, ncol(ss$jacobian)))
    cf <- optim(coef_hat, function(b) F(b, f0_hat), method = "BFGS",
                control = list(maxit = 2000, reltol = 1e-15))$par
    hf <- 1e-3
    kc <- ncol(D)
    kt <- ncol(ss$jacobian)
    F_cc <- matrix(0, kc, kc)
    F_ct <- matrix(0, kc, kt)
    for (a in seq_len(kc)) {
      ea <- replace(numeric(kc), a, hf)
      for (b in seq_len(kc)) {
        eb <- replace(numeric(kc), b, hf)
        F_cc[a, b] <- (F(cf + ea + eb, f0_hat) - F(cf + ea - eb, f0_hat) -
                         F(cf - ea + eb, f0_hat) + F(cf - ea - eb, f0_hat)) / (4 * hf^2)
      }
      for (b in seq_len(kt)) {
        up <- f0_at(replace(numeric(kt), b, hf))
        dn <- f0_at(replace(numeric(kt), b, -hf))
        F_ct[a, b] <- (F(cf + ea, up) - F(cf + ea, dn) - F(cf - ea, up) + F(cf - ea, dn)) / (4 * hf^2)
      }
    }
    G_firth <- -solve(F_cc, F_ct)
    G_at_firth <- .dp_two_step_variance(Z, D, cf, ridge, ss)$G
    data.frame(
      taxon = tx,
      score = err_score,
      G = max(abs(tv$G - G_fd)) / max(abs(G_fd)),
      G_firth = max(abs(G_at_firth - G_firth)) / max(abs(G_firth))
    )
  }))
  expect_true(all(res_G$score <= 1e-6), label = "analytic DP score matches finite differences of the objective")
  expect_true(all(res_G$G <= 1e-4), label = "sensitivity G matches finite differences of the refitted root")
  expect_true(all(res_G$G_firth <= 0.25), label = "sensitivity G agrees with the Firth-penalized estimator")
})

test_that("DP test with the second-stage fit: total variance, statistic, and fallback to the supplied f0", {
  fit <- second_stage_fit()
  usable <- usable_taxa(fit)
  check_taxa <- head(colnames(y)[usable], 4L)
  # The fit of the first usable taxon is marked as not converged and the
  # second usable taxon is removed from the fitting set.
  fit_fb <- fit
  fallback_taxon <- check_taxa[1L]
  fit_fb$fits[[fallback_taxon]]$ok <- FALSE
  outside_taxon <- check_taxa[2L]
  fit_fb$fitting_set <- setdiff(colnames(y), outside_taxon)
  dp2 <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0_first_stage, test_terms = test_terms,
                                      ridge = ridge, da_fit = fit_fb, x = x)
  fit_none <- fit
  fit_none$fitting_set <- character(0)
  dp1 <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0_first_stage, test_terms = test_terms,
                                      ridge = ridge, da_fit = fit_none, x = x)
  ss_rows <- dp2$f0_used == "second_stage"
  expect_true(sum(ss_rows) >= 2L * (sum(usable) - 2L), label = "second-stage f0 used for the taxa with a usable fit")
  expect_true(all(dp2$var_propagated[ss_rows] >= 0 & dp2$var_dp[ss_rows] > 0) &&
                all(dp2$var_dp[ss_rows] + dp2$var_propagated[ss_rows] >= dp2$var_dp[ss_rows]) &&
                all(dp2$variance_ratio > 0 & dp2$variance_ratio <= 1),
              label = "total variance is at least V^DP; variance ratio in (0, 1]")
  expect_true(any(dp2$variance_ratio[ss_rows] < 0.999), label = "propagated variance is positive for some taxa")
  expect_true(isTRUE(all.equal(dp2$statistic, dp2$variance_ratio * dp2$statistic_lr, tolerance = 1e-12)) &&
                isTRUE(all.equal(dp2$p_value, stats::pchisq(dp2$statistic, 1, lower.tail = FALSE),
                                 tolerance = 1e-12)),
              label = "statistic is the variance ratio times the likelihood-ratio statistic")
  fb2 <- dp2[dp2$taxon == fallback_taxon, ]
  fb1 <- dp1[dp1$taxon == fallback_taxon, ]
  expect_true(all(fb2$f0_used == "supplied_fallback") && all(fb2$variance_ratio == 1) &&
                isTRUE(all.equal(fb2$p_value, fb1$p_value, tolerance = 0)),
              label = "taxon without a converged fit uses the supplied f0")
  ou2 <- dp2[dp2$taxon == outside_taxon, ]
  ou1 <- dp1[dp1$taxon == outside_taxon, ]
  expect_true(all(ou2$f0_used == "supplied_outside_fitting_set") && all(ou2$variance_ratio == 1) &&
                isTRUE(all.equal(ou2$p_value, ou1$p_value, tolerance = 0)),
              label = "taxon outside the fitting set uses the supplied f0")
  expect_true(all(dp1$f0_used == "supplied_outside_fitting_set") && all(dp1$variance_ratio == 1) &&
                identical(dp1$statistic, dp1$statistic_lr),
              label = "empty fitting set: the statistic is the likelihood-ratio statistic")
})
