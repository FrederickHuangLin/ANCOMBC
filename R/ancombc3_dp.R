# =============================================================================
# ANCOM-BC3 differential-prevalence inference
# =============================================================================
#
# Zero-pattern model of the DP test. For taxon j, the observed-zero indicator
# Z_ij = 1(Y_ij = 0) has
#   q_ij = Pr(Z_ij = 1) = (1 - pi_ij^DP) + pi_ij^DP f0_ij,
#   pi_ij^DP = expit(eta_ij^DP),  eta_ij^DP = a_j^DP + r_i' gamma_j,
# where f0_ij, the probability that an at-risk taxon has a zero count, is held
# fixed in the DP fit. The DP coefficients coef = (a_j^DP, gamma_j) are
# estimated by maximizing the Bernoulli log-likelihood of the indicators with a
# Firth penalty and a ridge penalty on the non-intercept coefficients. A
# positive gamma_j,T means higher prevalence. The DP eligibility
# (dp_eligibility()) uses the structural-absence weights s_ij^el of the first
# stage.
#
# Sampling-zero probability. The DP test of the method uses the second-stage
# sampling-zero probability
#   f0_ij = f(0; mu_hat_ij, sigma_hat_j),
#   nu_hat_ij = log mu_hat_ij = alpha_hat_j + x_i' beta_hat_j + O^A_ij + d_hat_i,
# the zero probability of the Poisson-log-normal distribution at the at-risk
# mean and the log-scale standard deviation sigma_hat_j = exp(kappa_hat_j) of
# the second-stage fit (fit_ancombc3_da()). This f0_ij changes with the fitted
# abundance effects of the taxon, so that the sampling zeros created by an
# abundance effect are not attributed to structural absence. For a taxon whose
# second-stage fit did not converge or has no covariance matrix, and for a
# taxon outside the fitting set of the second stage (the DA-eligible taxa,
# whose positive counts suffice to estimate the abundance coefficients), f0_ij
# is the first-stage probability f0_ij^full.
#
# Two-step variance. With the second-stage f0, coef_hat_j depends on the
# second-stage estimate theta_hat_j = (alpha_hat_j, beta_hat_j, kappa_hat_j)
# through f0. The two-step (Murphy-Topel) variance of coef_hat_j is
#   V_j^tot = V_j^DP + G_j V_j^theta G_j',
# where V_j^DP = H_j^-1 is the inverse observed information of the DP fit at
# fixed f0, V_j^theta is the second-stage covariance of theta_hat_j (the
# (alpha_j, beta_j, kappa_j) block of vcov_psi; the (alpha_j, beta_j) block of
# vcov, with kappa_j fixed, when vcov_psi is not available), and, by the
# implicit function theorem applied to the DP score equations S_j = 0,
#   G_j = d coef_hat_j / d theta_j' = H_j^-1 sum_i (d S_j / d f0_ij) (d f0_ij / d theta_j'),
#   d f0_ij / d nu_ij = -f0_ij E(Lambda_ij | y_ij = 0),
#   d nu_ij / d (alpha_j, beta_j')' = (1, x_i')',
#   d f0_ij / d kappa_j = sigma_j^2 f0_ij E(Lambda_ij^2 | y_ij = 0),
# with the posterior moments of the Poisson-log-normal intensity Lambda_ij
# (ancombc3_pln.R). The sample-specific bias d_hat is held fixed: it is common
# to all taxa and restricted to be orthogonal to the tested covariates. The
# covariance between the DP and the second-stage estimating functions, which
# share the zero counts, is not included.
#
# Test. The penalized likelihood-ratio statistic T_j^LR of gamma_j,T = 0 at
# fixed f0 is multiplied by the variance ratio
# rho_j,T = V_j,TT^DP / V_j,TT^tot, which is at most one:
#   T_j^DP = rho_j,T T_j^LR,  P value 1 - F_chi2_1(T_j^DP).
# V^DP, G, and the propagated term are evaluated at the restricted fit
# (gamma_j,T = 0, the other DP coefficients re-estimated). Under the null
# hypothesis, T_j^LR is asymptotically equivalent to the score statistic for
# gamma_j,T at the restricted fit, U_T^2 V_j,TT^DP with U_T the efficient
# score, and the two-step variance of U_T is the variance at fixed f0
# multiplied by 1 / rho_j,T; T_j^DP is therefore the score statistic with the
# two-step variance (a first-order Rao-Scott type correction). The rescaled
# likelihood-ratio statistic retains the behavior of the Firth-penalized
# likelihood-ratio test when the zero indicators are nearly separated by the
# covariate, where the Wald statistic loses power, and the evaluation at the
# restricted fit avoids the flat likelihood of the unrestricted fit in that
# case. For a taxon with the supplied f0, rho_j,T = 1 and T_j^DP = T_j^LR.
#
# Dimension convention (n samples, p taxa, q prevalence covariates):
#   y, structural_absence_weight, f0: n x p (samples by taxa).
#   r: n x q prevalence design, without an intercept.
#   x: n x q_x abundance design of the second stage, without an intercept.
#   Returned data frames: one row per taxon and tested term.

.as_numeric_matrix_dp <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (anyNA(x)) {
    stop(name, " contains missing values.")
  }
  x
}

# Taxon names of a sample-by-taxon matrix: its column names, or T1, ..., Tp.
.dp_taxon_names <- function(mat) .ac3_taxon_names(mat)

# DP eligibility and q-values (p_adj_method) within one tested term. The eligibility table
# of dp_eligibility() is matched by taxon and term (every taxon is eligible
# when eligibility is NULL), and the q-values are computed over the eligible
# taxa with a finite P value.
.dp_apply_eligibility <- function(results, eligibility = NULL, p_adj_method = "BH", alpha = 0.1) {
  results$q_value <- 1
  results$discovered <- FALSE
  if (is.null(eligibility)) {
    results$dp_eligible <- TRUE
    results$dp_eligibility_reason <- "eligible_no_criteria"
  } else {
    key <- paste(results$taxon, results$test_term, sep = "\r")
    elig_key <- paste(eligibility$taxon, eligibility$test_term, sep = "\r")
    idx <- match(key, elig_key)
    results$dp_eligible <- ifelse(is.na(idx), FALSE, eligibility$dp_eligible[idx])
    results$dp_eligibility_reason <- ifelse(
      is.na(idx),
      "missing_eligibility",
      eligibility$dp_eligibility_reason[idx]
    )
    if ("absence_mass" %in% names(eligibility)) {
      results$absence_mass <- ifelse(is.na(idx), NA_real_, eligibility$absence_mass[idx])
      results$at_risk_mass <- ifelse(is.na(idx), NA_real_, eligibility$at_risk_mass[idx])
      results$min_side_absence_mass <- ifelse(is.na(idx), NA_real_, eligibility$min_side_absence_mass[idx])
      results$min_side_at_risk_mass <- ifelse(is.na(idx), NA_real_, eligibility$min_side_at_risk_mass[idx])
    }
  }
  ok <- is.finite(results$p_value) & results$dp_eligible
  if (any(ok)) {
    results$q_value[ok] <- p.adjust(results$p_value[ok], method = p_adj_method)
    results$discovered[ok] <- results$q_value[ok] <= alpha
  }
  results
}

# DP eligibility from the structural-absence weights
# s_ij^el = 1{y_ij = 0} w_ij of the first stage (ancombc3_fitted_vae_offsets()),
# restricted to [0, 1], and the at-risk weights 1 - s_ij^el. The expected
# structural-zero count M_j^0 = sum_i s_ij^el and the expected at-risk count
# M_j^+ = sum_i (1 - s_ij^el) must be at least min_total_absence and
# min_total_at_risk. On each side of the tested covariate (the two groups of a
# binary covariate; for a continuous covariate, the samples at or below the
# median and the samples above it), the expected structural-zero and at-risk
# counts must be at least min_side_absence and min_side_at_risk. For a
# continuous covariate, the covariate variances weighted by s^el and by
# 1 - s^el must both be at least min_weighted_variance. The default thresholds
# are those of the manuscript (Methods, table of eligibility criteria). The
# taxa are the column names of structural_absence_weight (those of the counts
# for the output of ancombc3_fitted_vae_offsets()), or T1, ..., Tp when it has
# none; fit_ancombc3_dp_zero_pattern() matches the returned table by taxon
# name. With repeated measures, cluster is the subject of each sample. For a
# tested covariate that is constant within every subject, the expected counts
# of each side are then computed over subjects (.ac3_cluster_side_mass()): the
# sum over the subjects of the side of the largest weight of the subject. The
# total expected counts and the weighted variances are unchanged.
dp_eligibility <- function(structural_absence_weight,
                           r,
                           test_terms = colnames(r),
                           min_total_absence = 1,
                           min_total_at_risk = 5,
                           min_side_absence = 0.5,
                           min_side_at_risk = 2,
                           min_weighted_variance = 0.01,
                           cluster = NULL) {
  structural_absence_weight <- .as_numeric_matrix_dp(structural_absence_weight, "structural_absence_weight")
  r <- .as_numeric_matrix_dp(r, "r")
  if (nrow(r) != nrow(structural_absence_weight)) {
    stop("r and structural_absence_weight must have the same sample dimension.")
  }
  if (is.null(colnames(r))) colnames(r) <- paste0("r", seq_len(ncol(r)))
  colnames(structural_absence_weight) <- .dp_taxon_names(structural_absence_weight)
  test_terms <- intersect(test_terms, colnames(r))
  if (length(test_terms) == 0L) stop("No test_terms are present in r.")

  weighted_var <- function(x, w) {
    sw <- sum(w)
    if (sw <= 0) return(0)
    mu <- sum(w * x) / sw
    sum(w * (x - mu)^2) / sw
  }

  rows <- list()
  for (term in test_terms) {
    x <- as.numeric(r[, term])
    ux <- sort(unique(x[is.finite(x)]))
    is_binary <- length(ux) == 2L && all(ux %in% c(0, 1))
    side <- if (is_binary) x else as.numeric(x > stats::median(x, na.rm = TRUE))
    by_subject <- !is.null(cluster) && .ac3_cluster_constant_within(x, .ac3_cluster_index(cluster))

    for (j in seq_len(ncol(structural_absence_weight))) {
      s <- pmin(pmax(structural_absence_weight[, j], 0), 1)
      a <- 1 - s
      absence_mass <- sum(s)
      at_risk_mass <- sum(a)
      side_absence <- tapply(s, side, sum)
      side_at_risk <- tapply(a, side, sum)
      if (by_subject) {
        side_absence <- .ac3_cluster_side_mass(s, side, .ac3_cluster_index(cluster))
        side_at_risk <- .ac3_cluster_side_mass(a, side, .ac3_cluster_index(cluster))
      }
      min_side_s <- if (length(side_absence) >= 2L) min(side_absence) else 0
      min_side_a <- if (length(side_at_risk) >= 2L) min(side_at_risk) else 0
      s_var <- if (is_binary) min_side_s else weighted_var(x, s)
      a_var <- if (is_binary) min_side_a else weighted_var(x, a)

      reasons <- character(0)
      if (absence_mass < min_total_absence) reasons <- c(reasons, "low_total_absence")
      if (at_risk_mass < min_total_at_risk) reasons <- c(reasons, "low_total_at_risk")
      if (min_side_s < min_side_absence) reasons <- c(reasons, "low_side_absence")
      if (min_side_a < min_side_at_risk) reasons <- c(reasons, "low_side_at_risk")
      if (!is_binary && (s_var < min_weighted_variance || a_var < min_weighted_variance)) {
        reasons <- c(reasons, "low_weighted_variance")
      }
      rows[[length(rows) + 1L]] <- data.frame(
        taxon = colnames(structural_absence_weight)[j],
        test_term = term,
        dp_eligible = length(reasons) == 0L,
        dp_eligibility_reason = if (length(reasons) == 0L) "eligible" else paste(reasons, collapse = ";"),
        absence_mass = absence_mass,
        at_risk_mass = at_risk_mass,
        min_side_absence_mass = min_side_s,
        min_side_at_risk_mass = min_side_a,
        term_type = if (is_binary) "binary" else "continuous",
        row.names = NULL
      )
    }
  }
  do.call(rbind, rows)
}

# Negative penalized log-likelihood of the zero-pattern model of one taxon as a
# function of the free DP coefficients coef_free (coef = (a^DP, gamma), with
# the positions fix_idx held at fix_vals):
#   -sum_i {Z_i log q_i + (1 - Z_i) log pi_i^DP}
#     + (lambda_DP / 2) ||coef_free (non-intercept)||^2 - (1 / 2) log det I,
# with eta_i^DP = D_i' coef truncated to [-35, 35] and
# q_i = (1 - pi_i^DP) + pi_i^DP f0_i restricted to [1e-12, 1 - 1e-12]. For a
# positive count, log(1 - q_i) = log pi_i^DP + log(1 - f0_i), and the term
# log(1 - f0_i), which does not depend on coef, is omitted. The Firth penalty
# uses the expected information I = D' W D, with
#   W_ii = [pi_i (1 - pi_i) (A_i + B_i)]^2 q_i (1 - q_i),
#   A_i = (1 - f0_i) / q_i,  B_i = 1 / pi_i,
# (pi_i restricted to [1e-9, 1 - 1e-9]), which equals
# (1 - pi_i)^2 (1 - q_i) / q_i; the penalty is omitted when det I is not
# positive and finite. firth = FALSE omits the Firth penalty (used by the
# numerical checks of the two-step sensitivity in the tests of the package).
.dp_zero_pattern_nll <- function(Z, f0, D, ridge = 1e-3,
                                 fix_idx = integer(0), fix_vals = numeric(0),
                                 firth = TRUE) {
  pcol <- ncol(D)
  free_idx <- setdiff(seq_len(pcol), fix_idx)
  noninter_idx <- setdiff(seq_len(pcol), 1L)
  pen_free <- intersect(free_idx, noninter_idx)
  pen_loc <- match(pen_free, free_idx)
  function(coef_free) {
    coef <- numeric(pcol)
    coef[free_idx] <- coef_free
    if (length(fix_idx)) coef[fix_idx] <- fix_vals
    eta_dp <- pmax(pmin(as.numeric(D %*% coef), 35), -35)
    pi_dp <- 1 / (1 + exp(-eta_dp))
    q <- pmin(pmax((1 - pi_dp) + pi_dp * f0, 1e-12), 1 - 1e-12)
    penalty <- if (length(pen_loc)) {
      0.5 * ridge * sum(coef_free[pen_loc]^2)
    } else {
      0
    }
    firth_pen <- 0
    if (firth) {
      pi_c <- pmin(pmax(pi_dp, 1e-9), 1 - 1e-9)
      A <- (1 - f0) / q
      B <- 1 / pi_c
      W <- (pi_c * (1 - pi_c) * (A + B))^2 * q * (1 - q)
      info_dp <- crossprod(D, W * D)
      ld <- tryCatch(determinant(info_dp, logarithm = TRUE), error = function(e) NULL)
      if (!is.null(ld) && is.finite(as.numeric(ld$modulus)) && ld$sign > 0) {
        firth_pen <- 0.5 * as.numeric(ld$modulus)
      }
    }
    -sum(Z * log(q) + (1 - Z) * log(pi_dp)) + penalty - firth_pen
  }
}

# Second-stage sampling-zero probability of one taxon and its derivatives in
# the second-stage parameters theta_j, from the taxon-wise fit record fit_j of
# fit_ancombc3_da() and the abundance design X = [1, x]. With
# nu_i = log mu_hat_i (the fitted at-risk mean, which includes the offsets and
# the sample-specific bias) and the posterior moments of Lambda_i given a zero
# count (.ac3_pln_obs()),
#   f0_i = f(0; mu_hat_i, sigma_hat),
#   d f0_i / d nu_i = -f0_i E(Lambda_i | 0)   (zero where nu_i is truncated),
#   d f0_i / d kappa = sigma^2 f0_i E(Lambda_i^2 | 0).
# theta_j is (alpha_j, beta_j, kappa_j) with covariance the corresponding block
# of vcov_psi when vcov_psi is available, and (alpha_j, beta_j) with the
# corresponding block of vcov otherwise. engine selects the evaluation of the
# quadrature. Returns f0, the n x k matrix
# jacobian = d f0 / d theta', the k x k covariance vcov_theta, and
# kappa_included; NULL when the fit did not converge or its covariance is not
# available or not finite.
.dp_second_stage_taxon <- function(fit_j, X, engine = "r") {
  if (is.null(fit_j) || !isTRUE(fit_j$ok) || is.null(fit_j$vcov) || is.null(fit_j$mu)) return(NULL)
  sigma <- as.numeric(fit_j$sigma)
  log_mu <- log(as.numeric(fit_j$mu))
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0 ||
      length(log_mu) != nrow(X) || any(!is.finite(log_mu))) {
    return(NULL)
  }
  i_ab <- seq_len(ncol(X))
  k_omega <- length(fit_j$omega)
  vcov <- fit_j$vcov
  if (nrow(vcov) != k_omega || any(!is.finite(vcov[i_ab, i_ab]))) return(NULL)
  pl <- .ac3_pln_obs(log_mu - sigma^2 / 2, rep(0, length(log_mu)), sigma, moments = TRUE, engine = engine)
  f0 <- exp(pl$logf)
  in_range <- (log_mu > -30 & log_mu < 30) * 1
  df0_dnu <- -f0 * pl$E1 * in_range
  vpsi <- fit_j$vcov_psi
  kappa_included <- !is.null(vpsi) && nrow(vpsi) == k_omega + 1L &&
    all(is.finite(vpsi[c(i_ab, k_omega + 1L), c(i_ab, k_omega + 1L)]))
  if (kappa_included) {
    jacobian <- cbind(df0_dnu * X, sigma^2 * f0 * pl$E2)
    vcov_theta <- vpsi[c(i_ab, k_omega + 1L), c(i_ab, k_omega + 1L), drop = FALSE]
  } else {
    jacobian <- df0_dnu * X
    vcov_theta <- vcov[i_ab, i_ab, drop = FALSE]
  }
  if (any(!is.finite(f0)) || any(!is.finite(jacobian))) return(NULL)
  list(f0 = f0, jacobian = jacobian, vcov_theta = 0.5 * (vcov_theta + t(vcov_theta)),
       kappa_included = kappa_included)
}

# Observed information of the zero-pattern model at coef with the ridge
# penalty, and the derivative of the score in f0. With pi_i = pi_i^DP and q_i
# as in .dp_zero_pattern_nll(), the score of coef is
# S = sum_i (1 - pi_i) (q_i - Z_i) / q_i D_i, the observed information is
#   H = sum_i w_i D_i D_i' + lambda_DP diag(0, 1, ..., 1),
#   w_i = pi_i (1 - pi_i) (1 - Z_i / q_i) + pi_i (1 - pi_i)^2 (1 - f0_i) Z_i / q_i^2,
# and d S / d f0_i = pi_i (1 - pi_i) Z_i / q_i^2 D_i (the rows of dscore_df0).
# The Firth penalty, whose contribution to both is of order one, is not
# included. When H is not positive definite, the expected information
# sum_i pi_i (1 - pi_i)^2 (1 - f0_i) / q_i D_i D_i' with the ridge penalty is
# used. Returns H (NULL when neither is positive definite) and dscore_df0.
.dp_zero_pattern_information <- function(Z, f0, D, coef, ridge) {
  eta <- pmax(pmin(as.numeric(D %*% coef), 35), -35)
  pi_dp <- 1 / (1 + exp(-eta))
  q <- pmin(pmax((1 - pi_dp) + pi_dp * f0, 1e-12), 1 - 1e-12)
  pen <- c(0, rep(ridge, ncol(D) - 1L))
  w_obs <- pi_dp * (1 - pi_dp) * (1 - Z / q) + pi_dp * (1 - pi_dp)^2 * (1 - f0) * Z / q^2
  H <- crossprod(D * w_obs, D) + diag(pen, ncol(D))
  if (is.null(tryCatch(chol(H), error = function(e) NULL))) {
    w_exp <- pi_dp * (1 - pi_dp)^2 * (1 - f0) / q
    H <- crossprod(D * w_exp, D) + diag(pen, ncol(D))
    if (is.null(tryCatch(chol(H), error = function(e) NULL))) H <- NULL
  }
  list(H = H, dscore_df0 = D * (pi_dp * (1 - pi_dp) * Z / q^2))
}

# Two-step variance components of coef_hat at the unrestricted DP fit coef with
# the second-stage f0 of one taxon (ss = .dp_second_stage_taxon()):
# V^DP = H^-1, the sensitivity G = H^-1 (dscore_df0' jacobian), and the
# propagated covariance G V^theta G'. NULL when H is not available.
.dp_two_step_variance <- function(Z, D, coef, ridge, ss) {
  info <- .dp_zero_pattern_information(Z, ss$f0, D, coef, ridge)
  if (is.null(info$H)) return(NULL)
  V_dp <- chol2inv(chol(info$H))
  G <- V_dp %*% crossprod(info$dscore_df0, ss$jacobian)
  V_prop <- G %*% ss$vcov_theta %*% t(G)
  list(V_dp = V_dp, G = G, V_prop = 0.5 * (V_prop + t(V_prop)))
}

# DP test of every taxon and tested term by the zero-pattern model
# (.dp_zero_pattern_nll()), with the prevalence design D = [1, r] and f0
# bounded to [1e-8, 1 - 1e-8]. The taxa are the column names of y, or
# T1, ..., Tp when y has none. f0 is matched to them by name when it has
# column names, and by position otherwise; the second-stage fits
# da_fit$fits, the fitting set, and the eligibility table are matched by
# taxon name.
#   f0: n x p first-stage sampling-zero probabilities f0^full, used for the
#     taxa outside the fitting set of the second stage and for the taxa
#     without a usable second-stage fit.
#   da_fit, x: the second-stage fit (fit_ancombc3_da()) and its abundance
#     design x (without an intercept). f0 of each taxon of the fitting set
#     da_fit$fitting_set (every taxon when it is NULL) is the second-stage
#     sampling-zero probability, and the statistic is rescaled by the
#     two-step variance ratio (see the header of this file).
# The unrestricted fit and the null fit with the tested coefficient fixed at
# zero are computed by BFGS (at most 500 iterations; the null fit starts from
# the unrestricted estimate without the tested coefficient). The penalized
# likelihood-ratio statistic is T^LR = max{2 (l_unrestricted - l_null), 0}; the
# test statistic is T^DP = rho T^LR, with rho = V_TT^DP / V_TT^tot at the null
# fit for a taxon with the second-stage f0 and rho = 1 otherwise, and the P
# value is 1 - F_chi2_1(T^DP). When either fit fails, the statistic is 0 and
# the P value is 1; when the two-step variance is not available for a taxon
# with the second-stage f0, rho = 1. estimate is the tested coefficient gamma_j,T of the
# unrestricted fit (positive = higher prevalence). The eligibility table of
# dp_eligibility() selects the taxa of the q-values (p_adj_method) within each
# tested term (every taxon when eligibility is NULL). The fitted latent-prevalence
# offset of the first stage is not part of the DP predictor; the response is
# the observed zero indicator. Besides estimate, statistic (T^DP), and
# p_value, the result holds f0_used ("second_stage",
# "supplied_outside_fitting_set", or "supplied_fallback" for a taxon of the
# fitting set without a usable second-stage fit), statistic_lr (T^LR),
# var_dp (V_TT^DP) and var_propagated ((G V^theta G')_TT) at the null fit,
# variance_ratio (rho), and kappa_included.
#   joint_terms: optional columns of r (the columns of a group factor). When
#     given, the joint test that their coefficients are all zero and the
#     covariance of their estimates are computed for every taxon from the same
#     unrestricted fit (.ac3_mg_dp_joint_test(), ancombc3_multigroup.R) and
#     returned, named by taxon, in the attribute joint of the result. The rows
#     of the result do not depend on joint_terms. With cluster, the joint
#     quantities use the cluster-robust covariance of the DP coefficients.
# With repeated measures, cluster is the subject of each sample. The variance
# ratio of every taxon is then rho = V_TT^DP / (V_TT^DP,CR + (G V^theta G')_TT)
# at the null fit, with the cluster-robust variance V^DP,CR of the tested
# coefficient at fixed f0 (.ac3_cluster_dp_variance(); the propagated term is
# zero for a taxon with the supplied f0), and the result holds the further
# column var_dp_robust (V_TT^DP,CR). The ratio used is the smaller of this
# value and the ratio without cluster, so that the statistic with cluster does
# not exceed the statistic without cluster. When the variance components are
# not available, the ratio without cluster is used. engine
# selects the evaluation of the quadrature of the second-stage f0
# (.dp_second_stage_taxon()).
fit_ancombc3_dp_zero_pattern <- function(y,
                                         r,
                                         f0,
                                         test_terms = colnames(r),
                                         p_adj_method = "BH",
                                         alpha = 0.1,
                                         ridge = 1e-3,
                                         eligibility = NULL,
                                         da_fit,
                                         x,
                                         joint_terms = NULL,
                                         cluster = NULL,
                                         engine = "r") {
  y <- .as_numeric_matrix_dp(y, "y")
  r <- .as_numeric_matrix_dp(r, "r")
  taxa <- .dp_taxon_names(y)
  f0 <- .ac3_match_taxa(.as_numeric_matrix_dp(f0, "f0"), taxa, "f0")
  if (!identical(dim(y), dim(f0))) stop("y and f0 must have the same n x p dimensions.")
  if (nrow(r) != nrow(y)) stop("r and y must have the same sample dimension.")
  if (is.null(colnames(r))) colnames(r) <- paste0("r", seq_len(ncol(r)))
  test_terms <- intersect(test_terms, colnames(r))
  if (length(test_terms) == 0L) stop("No test_terms are present in r.")
  D <- cbind("(Intercept)" = 1, r)
  pcol <- ncol(D)
  X <- cbind(1, .as_numeric_matrix_dp(x, "x"))
  if (nrow(X) != nrow(y)) stop("x and y must have the same sample dimension.")
  if (is.null(da_fit$fits)) stop("da_fit must hold the taxon-wise fits (fits).")

  rows <- list()
  joint <- list()
  ids <- if (is.null(cluster)) NULL else .ac3_cluster_index(cluster)
  if (!is.null(joint_terms) && !all(joint_terms %in% colnames(r))) {
    stop("joint_terms must be columns of r.")
  }
  for (j in seq_len(ncol(y))) {
    Z <- as.numeric(y[, j] == 0)
    in_set <- is.null(da_fit$fitting_set) || taxa[j] %in% da_fit$fitting_set
    ss <- if (in_set) .dp_second_stage_taxon(da_fit$fits[[taxa[j]]], X, engine = engine) else NULL
    f0_used <- if (!in_set) {
      "supplied_outside_fitting_set"
    } else if (is.null(ss)) {
      "supplied_fallback"
    } else {
      "second_stage"
    }
    f0_j <- pmin(pmax(if (is.null(ss)) f0[, j] else ss$f0, 1e-8), 1 - 1e-8)
    nll <- .dp_zero_pattern_nll(Z, f0_j, D, ridge = ridge)
    fit <- tryCatch(
      optim(rep(0, pcol), nll, method = "BFGS", control = list(maxit = 500)),
      error = function(e) NULL
    )
    ok <- !is.null(fit) && fit$convergence == 0 && all(is.finite(fit$par))
    if (!is.null(ss)) ss$f0 <- f0_j
    if (!is.null(joint_terms)) {
      joint[[taxa[j]]] <- .ac3_mg_dp_joint_test(Z, f0_j, D, ridge, fit, ok, ss,
                                                match(joint_terms, colnames(D)),
                                                ids = ids)
    }
    for (term in test_terms) {
      idx <- match(term, colnames(D))
      est <- stat <- pval <- stat_lr <- NA_real_
      var_dp <- var_prop <- NA_real_
      rho <- 1
      if (ok) {
        est <- fit$par[idx]
      }
      null_nll <- .dp_zero_pattern_nll(
        Z, f0_j, D, ridge = ridge, fix_idx = idx, fix_vals = 0
      )
      init <- if (ok) fit$par[-idx] else rep(0, pcol - 1L)
      null_fit <- tryCatch(
        optim(init, null_nll, method = "BFGS", control = list(maxit = 500)),
        error = function(e) NULL
      )
      null_ok <- !is.null(null_fit) && null_fit$convergence == 0 &&
        all(is.finite(null_fit$par))
      if (null_ok && !is.null(ss)) {
        coef0 <- numeric(pcol)
        coef0[-idx] <- null_fit$par
        tv <- .dp_two_step_variance(Z, D, coef0, ridge, ss)
        if (!is.null(tv)) {
          var_dp <- tv$V_dp[idx, idx]
          var_prop <- max(tv$V_prop[idx, idx], 0)
          if (is.finite(var_dp) && var_dp > 0 && is.finite(var_prop)) rho <- var_dp / (var_dp + var_prop)
        }
      }
      var_dp_robust <- NA_real_
      if (null_ok && !is.null(cluster)) {
        coef0 <- numeric(pcol)
        coef0[-idx] <- null_fit$par
        cr <- .ac3_cluster_dp_variance(Z, f0_j, D, coef0, ridge, idx, ss, ids)
        if (!is.null(cr)) {
          var_dp <- cr$var_dp
          var_prop <- if (is.null(ss)) NA_real_ else cr$var_propagated
          var_dp_robust <- cr$var_dp_robust
          if (is.finite(cr$rho)) rho <- min(rho, cr$rho)
        }
      }
      if (ok && null_ok && is.finite(fit$value) && is.finite(null_fit$value)) {
        stat_lr <- max(2 * ((-fit$value) - (-null_fit$value)), 0)
        stat <- rho * stat_lr
        pval <- pchisq(stat, df = 1, lower.tail = FALSE)
      } else {
        stat_lr <- stat <- 0
        pval <- 1
      }
      rows[[length(rows) + 1L]] <- data.frame(
        taxon = taxa[j], test_term = term,
        estimate = est, statistic = stat, p_value = pval,
        f0_used = f0_used, statistic_lr = stat_lr,
        var_dp = var_dp, var_propagated = var_prop, variance_ratio = rho,
        kappa_included = if (is.null(ss)) NA else ss$kappa_included,
        row.names = NULL
      )
      if (!is.null(cluster)) rows[[length(rows)]]$var_dp_robust <- var_dp_robust
    }
  }
  out <- do.call(rbind, rows)
  out <- do.call(rbind, lapply(
    split(out, out$test_term),
    .dp_apply_eligibility, eligibility = eligibility,
    p_adj_method = p_adj_method, alpha = alpha
  ))
  rownames(out) <- NULL
  if (!is.null(joint_terms)) attr(out, "joint") <- joint
  out
}
