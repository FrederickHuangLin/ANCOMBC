# =============================================================================
# ANCOM-BC3 repeated measures (cluster-robust inference) and bias-corrected
# abundances
# =============================================================================
#
# Repeated measures. The argument cluster of ancombc3() names the subject of
# each sample. The samples of one subject are dependent, and the subjects are
# independent. The estimates of the taxon-wise parameters, of the
# sample-specific bias, and of the DP coefficients are those of the analysis
# without cluster; the covariance estimates below replace the model-based
# ones. The common-shift estimate is computed by the same estimator
# (score_ancombc3_da()) from the cluster-robust standard errors, on which its
# intervals and weights depend. With s = 1, ..., m indexing the subjects:
#
#   1. Taxon-wise parameters (.ac3_cluster_da_covariance()). For taxon j with
#      psi_j = (omega_j, kappa_j), let u_ij = d l_ij / d psi_j be the score of
#      the cell log-likelihood l_ij of sample i at the final estimates: the
#      analytic derivatives in omega_j, and a central difference in kappa_j
#      with step 0.05, the step of the information matrix J_j
#      (.ac3_joint_augment()). The penalties of the taxon objective (ridge,
#      Firth-type penalty, adjustment term, and the cross-taxon penalty on
#      kappa_j) are not sums over samples. Their curvature is part of J_j, and
#      their gradient, which equals minus the sum of the scores at the
#      estimate, is attributed in equal parts to the samples: the scores are
#      centered at their mean over the samples, u_ij - ubar_j. With the subject
#      totals U_sj = sum_{i in s} (u_ij - ubar_j),
#        B_j = sum_s U_sj U_sj',
#        V_j^CR = c J_j^-1 B_j J_j^-1,   c = {m / (m - 1)} {(N - 1) / (N - k)},
#      with the finite-sample factor c of .ac3_cluster_factor(), N the number
#      of samples of the taxon likelihood, and k the length of psi_j.
#      V_j^CR replaces vcov_psi and its omega_j block replaces vcov in the fit
#      record; the model-based matrices are kept as vcov_psi_model and
#      vcov_model. When J_j is not positive definite (no vcov_psi), the same
#      construction is applied to omega_j with the inverse of H_j and k the
#      length of omega_j.
#   2. Sample-specific bias. The sensitivities g_j = d beta_hat_j / d d and the
#      information K_d about d are derivatives of the estimating equations, so
#      they are formed with the inverse information J_j^-1
#      (ancombc3_add_bias_covariance(), before the replacement of item 1). The
#      model-based working covariance V = H0 (H0' K_d H0)^-1 H0' is replaced by
#        V^CR = V M V m / (m - 1),
#      which satisfies both identifying constraints because V does. M is the
#      estimated covariance of the estimating function of d with psi_j
#      maximized out, the quantity whose model-based expectation is K_d. The
#      contribution of sample i and taxon j to that estimating function is the
#      n-vector
#        a_ij = e_i s_ij - I_d,psi,j J_j^-1 (u_ij - ubar_j),
#      where s_ij = d l_ij / d d_i is the score in the log mean, e_i is the
#      i-th unit vector, and I_d,psi,j = -d2 l_j / d d d psi_j'
#      (.ac3_da_taxon_bias_information()). With independent samples and a
#      correct model, sum_i E(a_ij a_ij') = I_dd,j - I_d,psi,j J_j^-1 I_psi,d,j,
#      the term of taxon j in K_d. Taxa are treated as independent, as in K_d,
#      and the samples of one subject as dependent:
#        M = sum_j sum_s A_sj A_sj',   A_sj = sum_{i in s} a_ij,
#      summed over the taxa that contribute to K_d. Without the second term of
#      a_ij, M has the elements sum_j s_ij s_i'j for two samples of one subject
#      and zero otherwise. The marginal variances var_raw_bias and
#      var_raw_marginal are recomputed with V^CR.
#   3. DP test (.ac3_cluster_dp_variance()). The statistic is
#      T^DP = rho T^LR. With the per-sample scores S_i of the zero-pattern
#      model at the restricted fit (.ac3_cluster_dp_scores()) and H its
#      information at fixed f0 (.dp_zero_pattern_information()),
#        V^DP,CR = c H^-1 {sum_s (sum_{i in s} S_i)(sum_{i in s} S_i)'} H^-1,
#        rho = V_TT^DP / (V_TT^DP,CR + (G V^theta G')_TT),
#      where c is the factor of item 1 with N the number of samples and k the
#      number of DP coefficients, V^DP = H^-1 is the model-based variance, and
#      the propagated term
#      uses the second-stage covariance V^theta of item 1. The scores are not
#      centered: at the restricted fit the tested component of their sum is the
#      numerator of the score statistic, and the second moment of the subject
#      totals at the restricted fit is the variance estimate of the generalized
#      score test. T^DP is then the score statistic with the cluster-robust
#      two-step variance. A taxon with the supplied f0 has no propagated term.
#   4. Degrees of freedom. The column df_count of the raw table is the number
#      of distinct subjects with a positive count of the taxon;
#      score_ancombc3_da() uses it in place of n_positive.
#   5. Eligibility (.ac3_cluster_binary_support(), .ac3_cluster_side_mass()).
#      For a tested covariate that is constant within every subject
#      (.ac3_cluster_constant_within()), a subject contributes to a group or
#      side the largest value over its samples: the per-group number of
#      positive counts is the number of subjects with at least one positive
#      count, and a per-group or per-side expected count is the sum over the
#      subjects of the largest probability over the samples of the subject.
#      The totals over all samples and the criteria of a covariate that varies
#      within a subject are those of the analysis without cluster.
#   6. First stage (.ac3_cluster_validation_split()). The withheld samples of
#      the ZIPLN-VAE are those of whole subjects.
#
# Bias-corrected abundances (.ac3_bias_corrected_abundance()): the total
# sample bias b_i = d_hat_i + sum_T delta_hat_T (x_iT - xbar_T) over the
# tested columns T; the bias-corrected log abundances, log y_ij - b_i for a
# positive count and E(log Lambda_ij | Y_ij = 0, S_ij = 1) - b_i for a zero;
# and the posterior probabilities of structural absence (struc_absence_prob
# of ancombc3()), 0 for a positive count and
# (1 - pi_ij) / {(1 - pi_ij) + pi_ij f0_ij} for a zero.

# -----------------------------------------------------------------------------
# Subjects
# -----------------------------------------------------------------------------

# Subject index of each sample: the integers 1, ..., m in the order of first
# appearance of the subjects.
.ac3_cluster_index <- function(cluster) {
  if (anyNA(cluster)) stop("cluster contains missing values.")
  match(cluster, unique(cluster))
}

# TRUE when the covariate v takes one value within every subject (differences
# up to 1e-8 max(1, max |v|) are treated as zero).
.ac3_cluster_constant_within <- function(v, ids) {
  v <- as.numeric(v)
  first <- v[match(ids, ids)]
  all(abs(v - first) <= 1e-8 * max(1, max(abs(v))))
}

# Largest value over the samples of each subject, for every column of an
# n x p matrix: an m x p matrix with the subjects in the order of their index.
.ac3_cluster_subject_max <- function(mat, ids) {
  mat <- as.matrix(mat)
  rows <- split(seq_along(ids), ids)
  out <- matrix(0, length(rows), ncol(mat), dimnames = list(NULL, colnames(mat)))
  for (s in seq_along(rows)) {
    out[s, ] <- do.call(pmax, lapply(rows[[s]], function(i) mat[i, ]))
  }
  out
}

# Subjects of ancombc3(): cluster must name one column of metadata without
# missing values among the analyzed samples. Fewer than 10 subjects give an
# error and fewer than 30 a warning. The number of subjects and the range of
# the number of samples per subject are reported. Returns the subject of each
# analyzed sample as a character vector.
.ac3_cluster_subjects <- function(metadata, cluster, verbose = TRUE) {
  if (!is.character(cluster) || length(cluster) != 1L || is.na(cluster) || !nzchar(cluster)) {
    stop("cluster must be the name of one column of meta_data.")
  }
  if (!cluster %in% names(metadata)) stop("cluster '", cluster, "' is not a column of meta_data.")
  subject <- metadata[[cluster]]
  if (anyNA(subject)) {
    stop("cluster '", cluster, "' has missing values for ", sum(is.na(subject)),
         " analyzed sample(s): ", .ac3_main_list(rownames(metadata)[is.na(subject)]), ".")
  }
  subject <- as.character(subject)
  size <- table(subject)
  m <- length(size)
  if (m < 10L) {
    stop("The analysis has ", m, " subject(s) (cluster '", cluster,
         "'); cluster-robust inference requires at least 10 subjects.")
  }
  if (m < 30L) {
    warning("The analysis has ", m, " subjects; cluster-robust standard errors are unreliable ",
            "with fewer than 30 subjects.", call. = FALSE)
  }
  .ac3_main_say(verbose, sprintf("Repeated measures: %d subjects (cluster '%s'), %d to %d samples per subject.",
                                 m, cluster, min(size), max(size)))
  subject
}

# Group sizes of the binary tested covariates of ancombc3() in subjects: for a
# binary tested covariate that is constant within every subject, each group
# must hold at least min_per_group subjects.
.ac3_cluster_check_groups <- function(metadata, tested, coding, subject, min_per_group) {
  ids <- .ac3_cluster_index(subject)
  first <- match(seq_len(max(ids)), ids)
  for (tn in tested) {
    if (!identical(substr(coding[[tn]], 1L, 6L), "binary")) next
    v <- as.numeric(metadata[[tn]])
    if (!.ac3_cluster_constant_within(v, ids)) next
    sizes <- table(factor(v[first], levels = c(0, 1)))
    if (min(sizes) < min_per_group) {
      stop("Tested covariate '", tn, "' has a group with ", min(sizes), " subject(s); the DA ",
           "eligibility requires at least ", min_per_group,
           " subjects with a positive count per group.")
    }
  }
  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Eligibility
# -----------------------------------------------------------------------------

# DA eligibility quantities of a binary tested covariate that is constant
# within every subject, for every taxon: the smaller of the two group sums of
# the subject-level at-risk probabilities (the largest at-risk probability over
# the samples of a subject), and the smaller of the two group numbers of
# subjects with at least one positive count. group1 is the sample-level
# indicator of the second group.
.ac3_cluster_binary_support <- function(q_at_risk, observed_positive, group1, ids) {
  first <- match(seq_len(max(ids)), ids)
  g1 <- group1[first]
  q_sub <- .ac3_cluster_subject_max(q_at_risk, ids)
  pos_sub <- .ac3_cluster_subject_max(observed_positive * 1, ids)
  list(
    at_risk_support = unname(pmin(colSums(q_sub[!g1, , drop = FALSE]), colSums(q_sub[g1, , drop = FALSE]))),
    min_positive_by_group = unname(pmin(colSums(pos_sub[!g1, , drop = FALSE]),
                                        colSums(pos_sub[g1, , drop = FALSE])))
  )
}

# Per-side expected counts of the DP eligibility for a tested covariate that is
# constant within every subject: the sum, over the subjects of each side, of
# the largest weight over the samples of the subject. side is the sample-level
# side indicator; the result has one element per side, as tapply(w, side, sum).
.ac3_cluster_side_mass <- function(w, side, ids) {
  first <- match(seq_len(max(ids)), ids)
  sub <- vapply(split(w, ids), max, numeric(1))
  tapply(sub, side[first], sum)
}

# -----------------------------------------------------------------------------
# First stage: withheld samples of the ZIPLN-VAE
# -----------------------------------------------------------------------------

# Withheld samples of train_zipln_vae() with repeated measures: whole subjects.
# The subjects are taken in the order of a random permutation, and the number
# k of withheld subjects is the one whose cumulative number of samples is
# closest to n_val (the smaller k when two are equally close; at least one
# subject is kept for training). The permutation is drawn from the current
# random-number state, the state is restored, and the k withheld subjects are
# drawn as sample.int(m, k), the first k subjects of that permutation. With one
# sample per subject, the withheld samples and the random-number consumption
# are those of sample.int(n, n_val), the split without cluster. Returns the
# sorted indices of the withheld samples.
.ac3_cluster_validation_split <- function(cluster, n_val) {
  ids <- .ac3_cluster_index(cluster)
  m <- max(ids)
  size <- tabulate(ids, nbins = m)
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1L)
  perm <- withr::with_preserve_seed(sample.int(m))
  k <- min(which.min(abs(cumsum(size[perm]) - n_val)), m - 1L)
  withheld <- sample.int(m, k)
  sort(which(ids %in% withheld))
}

# -----------------------------------------------------------------------------
# Second stage: per-sample scores and cluster-robust covariances
# -----------------------------------------------------------------------------

# Cell log-likelihoods l_i of one taxon at omega and kappa = log sigma, for the
# data dat of the taxon-wise objective (.ac3_da_taxon_data()): the sample terms
# of the penalized log-likelihood l^pen (.ac3_da_taxon_objective()).
.ac3_cluster_cell_loglik <- function(omega, kappa, dat) {
  log_mu_raw <- as.numeric(dat$X %*% omega[dat$i_abundance]) + dat$abundance_offset
  logit_pi_raw <- as.numeric(dat$R %*% omega[dat$i_presence]) + dat$structural_offset
  .ac3_pln_observation_derivatives(log_mu_raw, logit_pi_raw, dat$y, dat$pos, exp(kappa), need = 0L,
                                   engine = dat$engine)$ll
}

# Per-sample scores of one taxon at omega and sigma = dat$sigma. u is the
# n x k matrix whose row i is d l_i / d psi': the analytic derivatives in
# omega, (dl_i / d log mu_i) X_i and (dl_i / d logit pi_i) R_i, and, when kappa
# is TRUE, a last column for kappa = log sigma, the central difference
# {l_i(kappa + h) - l_i(kappa - h)} / (2 h) with h = 0.05. s is the score in
# the sample-specific bias, s_i = dl_i / d log mu_i.
.ac3_cluster_taxon_scores <- function(omega, dat, kappa = TRUE, h = 0.05) {
  log_mu_raw <- as.numeric(dat$X %*% omega[dat$i_abundance]) + dat$abundance_offset
  logit_pi_raw <- as.numeric(dat$R %*% omega[dat$i_presence]) + dat$structural_offset
  o <- .ac3_pln_observation_derivatives(log_mu_raw, logit_pi_raw, dat$y, dat$pos, dat$sigma, need = 1L,
                                         engine = dat$engine)
  u <- cbind(dat$X * o$dl_mu, dat$R * o$dl_pi)
  if (isTRUE(kappa)) {
    k0 <- log(dat$sigma)
    u <- cbind(u, (.ac3_cluster_cell_loglik(omega, k0 + h, dat) -
                     .ac3_cluster_cell_loglik(omega, k0 - h, dat)) / (2 * h))
  }
  dimnames(u) <- NULL
  list(u = u, s = o$dl_mu)
}

# Finite-sample factor of a cluster-robust covariance with m subjects, N
# samples, and k estimated parameters: {m / (m - 1)} {(N - 1) / (N - k)}. The
# second ratio is omitted when N or k is not given and when N <= k.
.ac3_cluster_factor <- function(m, N = NULL, k = NULL) {
  out <- m / (m - 1)
  if (!is.null(N) && !is.null(k) && N > k) out <- out * ((N - 1) / (N - k))
  out
}

# Cluster-robust covariance c C B C of the estimate with inverse information C
# and per-sample scores u (one row per sample, one column per parameter),
# where B = sum_s U_s U_s', U_s is the total over the samples of subject s of
# the scores centered at their mean over the samples, and c is the factor of
# .ac3_cluster_factor() with N the number of rows and k the number of columns
# of u. Returns the covariance (vcov) and the m x k matrix of the subject
# totals (U), with the subjects in the order of their index.
.ac3_cluster_sandwich <- function(C, u, ids) {
  m <- max(ids)
  U <- rowsum(sweep(u, 2L, colMeans(u), "-"), ids, reorder = TRUE)
  V <- C %*% crossprod(U) %*% C * .ac3_cluster_factor(m, nrow(u), ncol(u))
  list(vcov = 0.5 * (V + t(V)), U = U)
}

# Cluster-robust covariances of the second stage (items 1, 2, and 4 of the
# header). fit is the result of ancombc3_add_bias_covariance() at the final
# estimates, inp the inputs of the taxon-wise fits (fit_ancombc3_da()), d the
# final sample-specific bias, and cluster the subject of each sample. Returns
# fit with
#   fits[[j]]: vcov and vcov_psi replaced by their cluster-robust versions,
#     and the model-based matrices in vcov_model and vcov_psi_model;
#   raw: var_raw_model and se_raw_model from the cluster-robust covariance,
#     var_raw_bias, var_raw_marginal, and se_raw_marginal recomputed with the
#     cluster-robust bias covariance, the model-based values in
#     var_raw_model_based, se_raw_model_based, and var_raw_bias_model_based,
#     and df_count, the number of subjects with a positive count of the taxon;
#   da_variance_inputs: V replaced by V^CR = V M V m / (m - 1), with the
#     model-based matrix in V_model.
# A taxon whose scores or cluster-robust covariance are not finite keeps its
# model-based matrices.
.ac3_cluster_da_covariance <- function(fit, inp, d, cluster) {
  n <- nrow(inp$y)
  p <- ncol(inp$y)
  if (length(cluster) != n) stop("cluster must have one value per sample.")
  ids <- .ac3_cluster_index(cluster)
  m <- max(ids)
  if (m < 2L) stop("cluster-robust inference requires at least two subjects.")
  fits <- fit$fits
  S <- matrix(0, n, p)
  M <- matrix(0, n, n)
  for (j in seq_len(p)) {
    f <- fits[[j]]
    if (!isTRUE(f$ok) || is.null(f$vcov) || is.null(f$sigma) || !is.finite(f$sigma)) next
    has_kappa <- !is.null(f$vcov_psi)
    C <- if (has_kappa) f$vcov_psi else f$vcov
    if (any(!is.finite(C))) next
    omega <- as.numeric(f$omega)
    sc <- .ac3_cluster_taxon_scores(omega, .ac3_joint_taxon_dat(inp, j, d, f$sigma), kappa = has_kappa)
    if (ncol(sc$u) != nrow(C) || any(!is.finite(sc$u))) next
    sw <- .ac3_cluster_sandwich(C, sc$u, ids)
    if (any(!is.finite(sw$vcov))) next
    ## Contribution of taxon j to M, with a_ij = e_i s_ij - P_j (u_ij - ubar_j)
    ## and P_j = I_d,psi,j J_j^-1:
    ##   sum_s A_sj A_sj' = (s_j s_j') o Z Z' - E_j P_j' - P_j E_j' + P_j B_j P_j',
    ## where Z Z' indicates two samples of one subject and row i of E_j is
    ## s_ij U_sj' for the subject s of sample i. The first term is accumulated
    ## over the taxa after the loop.
    block <- .ac3_da_taxon_bias_information(f, inp$y[, j], inp$x, inp$r, engine = inp$engine)
    if (!is.null(block) && all(is.finite(sc$s))) {
      P <- crossprod(block$I_psi_d, block$C)
      E <- sc$s * sw$U[ids, , drop = FALSE]
      M <- M + tcrossprod(P %*% crossprod(sw$U) - E, P) - tcrossprod(P, E)
      S[, j] <- sc$s
    }
    k <- length(omega)
    f$vcov_model <- f$vcov
    f$vcov_psi_model <- f$vcov_psi
    if (has_kappa) f$vcov_psi <- sw$vcov
    f$vcov <- sw$vcov[seq_len(k), seq_len(k), drop = FALSE]
    fits[[j]] <- f
  }
  fit$fits <- fits
  M <- M + tcrossprod(S) * outer(ids, ids, "==")
  M <- 0.5 * (M + t(M))

  vi <- fit$da_variance_inputs
  if (!is.null(vi)) {
    V_cr <- vi$V %*% M %*% vi$V * .ac3_cluster_factor(m)
    vi$V_model <- vi$V
    vi$V <- 0.5 * (V_cr + t(V_cr))
    fit$da_variance_inputs <- vi
  }

  raw <- fit$raw
  taxa <- names(fits)
  subjects_positive <- colSums(rowsum((inp$y > 0) * 1, ids, reorder = TRUE) > 0)
  raw$var_raw_model_based <- raw$var_raw_model
  raw$se_raw_model_based <- raw$se_raw_model
  raw$var_raw_bias_model_based <- raw$var_raw_bias
  raw$df_count <- unname(subjects_positive[match(raw$taxon, taxa)])
  for (rr in seq_len(nrow(raw))) {
    j <- match(raw$taxon[rr], taxa)
    f <- fits[[j]]
    if (is.null(f$vcov_model)) next
    term <- raw$term[rr]
    idx <- 1L + match(term, colnames(inp$x))
    v <- f$vcov[idx, idx]
    ok <- is.finite(v) && v > 0
    raw$var_raw_model[rr] <- if (ok) v else NA_real_
    raw$se_raw_model[rr] <- if (ok) sqrt(v) else NA_real_
  }
  if (!is.null(vi)) {
    raw$var_raw_bias <- NA_real_
    raw$var_raw_marginal <- NA_real_
    raw$se_raw_marginal <- NA_real_
    for (rr in seq_len(nrow(raw))) {
      term <- raw$term[rr]
      j <- match(raw$taxon[rr], vi$taxa)
      if (is.na(j) || is.null(vi$G[[term]])) next
      g_j <- vi$G[[term]][j, ]
      if (!all(is.finite(g_j))) next
      bias_var <- drop(g_j %*% vi$V %*% g_j)
      if (!is.finite(bias_var) || !is.finite(raw$var_raw_model[rr]) || raw$var_raw_model[rr] <= 0) next
      raw$var_raw_bias[rr] <- max(bias_var, 0)
      raw$var_raw_marginal[rr] <- raw$var_raw_model[rr] + raw$var_raw_bias[rr]
      raw$se_raw_marginal[rr] <- sqrt(raw$var_raw_marginal[rr])
    }
  }
  fit$raw <- raw
  fit
}

# -----------------------------------------------------------------------------
# DP test
# -----------------------------------------------------------------------------

# Per-sample scores of the zero-pattern model at coef and fixed f0 (one row per
# sample), S_i = (1 - pi_i) (q_i - Z_i) / q_i D_i, with pi_i and q_i as in
# .dp_zero_pattern_nll(). The ridge and Firth penalties are not sums over
# samples and are not included.
.ac3_cluster_dp_scores <- function(Z, f0, D, coef) {
  eta <- pmax(pmin(as.numeric(D %*% coef), 35), -35)
  pi_dp <- 1 / (1 + exp(-eta))
  q <- pmin(pmax((1 - pi_dp) + pi_dp * f0, 1e-12), 1 - 1e-12)
  D * ((1 - pi_dp) * (q - Z) / q)
}

# Variance components of the tested DP coefficient (position idx) of one taxon
# at the restricted fit coef0 with repeated measures: the model-based variance
# var_dp = (H^-1)_TT, the cluster-robust variance var_dp_robust =
# c {H^-1 B H^-1}_TT with B the sum over the subjects of the outer products of
# the subject totals of the scores and c the factor of .ac3_cluster_factor()
# (N the number of samples, k the number of DP coefficients), the propagated term
# var_propagated = (G V^theta G')_TT (zero without the second-stage
# quantities ss of .dp_second_stage_taxon()), and the variance ratio
# rho = var_dp / (var_dp_robust + var_propagated). NULL when H is not
# available or a component is not finite and positive.
.ac3_cluster_dp_variance <- function(Z, f0, D, coef0, ridge, idx, ss, ids) {
  info <- .dp_zero_pattern_information(Z, f0, D, coef0, ridge)
  if (is.null(info$H)) return(NULL)
  m <- max(ids)
  V_dp <- chol2inv(chol(info$H))
  U <- rowsum(.ac3_cluster_dp_scores(Z, f0, D, coef0), ids, reorder = TRUE)
  V_cr <- V_dp %*% crossprod(U) %*% V_dp * .ac3_cluster_factor(m, nrow(D), ncol(D))
  var_prop <- 0
  if (!is.null(ss)) {
    G <- V_dp %*% crossprod(info$dscore_df0, ss$jacobian)
    var_prop <- max((G %*% ss$vcov_theta %*% t(G))[idx, idx], 0)
  }
  var_dp <- V_dp[idx, idx]
  var_cr <- V_cr[idx, idx]
  if (!is.finite(var_dp) || var_dp <= 0 || !is.finite(var_cr) || var_cr <= 0 || !is.finite(var_prop)) {
    return(NULL)
  }
  list(var_dp = var_dp, var_dp_robust = var_cr, var_propagated = var_prop,
       rho = var_dp / (var_cr + var_prop))
}

# Cluster-robust covariance at fixed f0 of all DP coefficients of one taxon at
# coef, the matrix form of .ac3_cluster_dp_variance() used by the joint DP
# quantities of the multigroup analysis (.ac3_mg_dp_joint_test()):
# V_dp = H^-1 and V_cr = c H^-1 {sum_s (sum_{i in s} S_i)(sum_{i in s} S_i)'} H^-1,
# with the uncentered per-sample scores S_i (.ac3_cluster_dp_scores()) and the
# factor c of .ac3_cluster_factor() (N the number of samples, k the number of DP
# coefficients). Also returns dscore_df0 of .dp_zero_pattern_information() for
# the propagated term. NULL when H is not available or V_cr is not finite.
.ac3_cluster_dp_covariance <- function(Z, f0, D, coef, ridge, ids) {
  info <- .dp_zero_pattern_information(Z, f0, D, coef, ridge)
  if (is.null(info$H)) return(NULL)
  V_dp <- chol2inv(chol(info$H))
  U <- rowsum(.ac3_cluster_dp_scores(Z, f0, D, coef), ids, reorder = TRUE)
  V_cr <- V_dp %*% crossprod(U) %*% V_dp * .ac3_cluster_factor(max(ids), nrow(D), ncol(D))
  if (any(!is.finite(V_cr))) return(NULL)
  list(V_dp = V_dp, V_cr = 0.5 * (V_cr + t(V_cr)), dscore_df0 = info$dscore_df0)
}

# -----------------------------------------------------------------------------
# Bias-corrected abundances
# -----------------------------------------------------------------------------

# Total sample bias, bias-corrected log abundances, and posterior
# probabilities of structural absence of ancombc3(). y is the n x p matrix of
# analyzed counts, d_hat the estimated sample-specific bias, x the abundance
# design, tested_cols its tested columns, and scored the table of
# score_ancombc3_da(). mu, pi, and f0 are the n x p matrices of the fitted
# at-risk means, at-risk probabilities, and sampling-zero probabilities of the
# final second-stage model, and sigma the fitted sigma_j (fit$mu, fit$pi,
# fit$f0, and fit$sigma_hat of fit_ancombc3_da()). The total sample bias is
#   b_i = d_hat_i + sum_T delta_hat_T (x_iT - xbar_T),
# with delta_hat_T the selected common-shift estimate of the tested column T
# (zero when the scored table has no finite estimate for T) and xbar_T the mean
# of the column over the analyzed samples; it is on the natural-log scale. The
# bias-corrected log abundance is log y_ij - b_i for a positive count and
# E(log Lambda_ij | Y_ij = 0, S_ij = 1) - b_i for a zero, the posterior mean
# log intensity of an at-risk taxon under the fitted Poisson-log-normal model:
# with c_ij = log mu_ij - sigma_j^2 / 2, it equals c_ij + sigma_j E(z | y = 0),
# where E(z | y) is the moment Ez of .ac3_pln_obs(). The posterior probability
# of structural absence is 0 for a positive count and
#   (1 - pi_ij) / {(1 - pi_ij) + pi_ij f0_ij}
# for a zero. The zero cells of a taxon without a converged second-stage fit are
# NA in both tables. Both tables have taxa in rows when taxa_in_rows is TRUE
# and samples in rows otherwise.
.ac3_bias_corrected_abundance <- function(y, d_hat, x, tested_cols, scored, taxa_in_rows = TRUE,
                                          mu = NULL, pi = NULL, f0 = NULL, sigma = NULL, engine = "r") {
  b <- as.numeric(d_hat)
  delta <- stats::setNames(numeric(length(tested_cols)), tested_cols)
  for (tn in tested_cols) {
    dl <- scored$delta[scored$term == tn]
    dl <- dl[is.finite(dl)]
    if (length(dl) > 0L) delta[[tn]] <- dl[1L]
    b <- b + delta[[tn]] * (x[, tn] - mean(x[, tn]))
  }
  b <- stats::setNames(as.numeric(b), rownames(y))
  corrected <- matrix(NA_real_, nrow(y), ncol(y), dimnames = dimnames(y))
  pos <- y > 0
  corrected[pos] <- (log(y) - b)[pos]
  absence <- matrix(NA_real_, nrow(y), ncol(y), dimnames = dimnames(y))
  absence[pos] <- 0
  if (!is.null(mu) && !is.null(pi) && !is.null(f0) && !is.null(sigma)) {
    for (j in seq_len(ncol(y))) {
      zero <- which(!pos[, j])
      if (length(zero) == 0L) next
      s <- sigma[j]
      m <- mu[zero, j]
      if (!is.finite(s) || s <= 0 || any(!is.finite(m)) || any(m <= 0)) next
      center <- log(m) - s^2 / 2
      Ez <- .ac3_pln_obs(center, rep(0, length(zero)), s, moments = TRUE, engine = engine)$Ez
      corrected[zero, j] <- center + s * Ez - b[zero]
      absent <- 1 - pi[zero, j]
      absence[zero, j] <- absent / (absent + pi[zero, j] * f0[zero, j])
    }
  }
  list(
    total_sample_bias = b,
    bias_corrected_log_abundance = if (isTRUE(taxa_in_rows)) t(corrected) else corrected,
    structural_absence_probability = if (isTRUE(taxa_in_rows)) t(absence) else absence,
    common_shift = delta
  )
}
