# Repeated-measures analysis: the argument cluster of ancombc3() and the
# cluster-robust covariances of the second stage and of the DP test.

rel_err <- function(a, b) max(abs(a - b) / pmax(1, abs(a)))

# Second-stage fit on simulated independent samples.
set.seed(20261110L)
n <- 60L
p <- 30L
x <- cbind(cont = as.numeric(scale(stats::rnorm(n))), binary = rep(0:1, length.out = n))
r <- x
alpha <- stats::rnorm(p, mean = 2, sd = 0.8)
bias <- stats::rnorm(n, sd = 0.3)
log_mu <- matrix(alpha, n, p, byrow = TRUE) + bias
lambda <- exp(log_mu + matrix(stats::rnorm(n * p, sd = 0.5), n, p))
present <- matrix(stats::rbinom(n * p, 1L, 0.85), n, p)
y <- matrix(stats::rpois(n * p, lambda), n, p) * present
colnames(y) <- paste0("T", seq_len(p))
rownames(y) <- paste0("S", seq_len(n))
mu_first_stage <- matrix(pmax(colMeans(y), 0.1), n, p, byrow = TRUE)
f0_first_stage <- (2 / (2 + mu_first_stage))^2
dimnames(f0_first_stage) <- dimnames(y)
test_terms <- c("cont", "binary")
ridge <- 1e-6
inp <- list(y = y, x = x, r = r, abundance_offset = matrix(0, n, p, dimnames = dimnames(y)),
            structural_offset = matrix(0, n, p, dimnames = dimnames(y)),
            sampling_zero_prob = f0_first_stage, ridge = ridge)
bounds <- .ac3_pln_sigma_bounds

# Second-stage fit, the taxa with a regular fit, and the cluster-robust
# covariances with every sample its own subject.
second_stage <- ac3_cached(function() {
  fit0 <- fit_ancombc3_da(
    y = y, x = x, r = r, sigma_vae = rep(0.5, p),
    abundance_offset = matrix(0, n, p), structural_offset = matrix(0, n, p),
    sampling_zero_prob = f0_first_stage, at_risk_prob = matrix(0.9, n, p),
    test_terms = test_terms, ridge = ridge, n_starts = 2L
  )
  regular <- which(vapply(fit0$fits, function(f) {
    isTRUE(f$ok) && !is.null(f$vcov_psi) && f$sigma > 1.01 * bounds[1L] && f$sigma < 0.99 * bounds[2L]
  }, logical(1)))
  fc1 <- .ac3_cluster_da_covariance(fit0, inp, fit0$d_hat, cluster = rownames(y))
  list(fit0 = fit0, regular = regular, fc1 = fc1)
})
taxon_dat <- function(fit0, j, sigma) .ac3_joint_taxon_dat(inp, j, fit0$d_hat, sigma)

# Repeated measures without covariate effects: m subjects with a fixed number
# of visits, a between-subject binary covariate (group), a within-subject visit
# covariate, Poisson-log-normal counts with a subject random intercept per
# taxon (standard deviation sd_subject), sample-specific depth factors, and
# structural zeros with a subject random intercept in the presence logit
# (standard deviation sd_presence, mean mean_presence over taxa).
simulate_repeated_measures <- function(seed, m = 40L, visits = 3L, p = 80L, sd_subject = 0.7,
                                       sd_presence = 2, mean_presence = 1) {
  set.seed(seed)
  n <- m * visits
  subject <- rep(seq_len(m), each = visits)
  visit <- rep(seq_len(visits) - 1, times = m)
  group <- rep(rep(0:1, length.out = m), each = visits)
  alpha <- stats::rnorm(p, 3, 1)
  a <- stats::rnorm(p, mean_presence, 1)
  b <- matrix(stats::rnorm(m * p, sd = sd_subject), m, p)
  u <- matrix(stats::rnorm(m * p, sd = sd_presence), m, p)
  depth <- stats::rnorm(n, sd = 0.5)
  log_lambda <- matrix(alpha, n, p, byrow = TRUE) + b[subject, ] + depth +
    matrix(stats::rnorm(n * p, sd = 0.5), n, p)
  present <- matrix(stats::rbinom(n * p, 1L, stats::plogis(matrix(a, n, p, byrow = TRUE) + u[subject, ])), n, p)
  y <- matrix(stats::rpois(n * p, exp(log_lambda)), n, p) * present
  dimnames(y) <- list(paste0("S", seq_len(n)), paste0("T", seq_len(p)))
  list(counts = t(y),
       metadata = data.frame(subject = paste0("P", subject), visit = visit, group = group,
                             row.names = rownames(y)))
}

# Steps of ancombc3() at its default settings and its default engine for
# numeric covariates, with the DA table scored with and without the empirical
# null, both without the offset term, and, with cluster, with the empirical
# null and the offset term (ancombc3_add_offset_covariance(), the default of
# ancombc3() with cluster).
run_steps <- function(dat, tested, formula_rhs, cluster = NULL, epochs = 20L, seed = 1L, alpha = 0.10) {
  cfg <- ancombc3_default_settings()
  y <- t(dat$counts)
  storage.mode(y) <- "double"
  subject <- if (is.null(cluster)) NULL else as.character(dat$metadata[[cluster]])
  design <- ancombc3_make_design(dat$metadata, x_formula = formula_rhs, test_terms = tested,
                                 r_formula = formula_rhs)
  offsets <- ancombc3_fitted_vae_offsets(y = y, x = design$x, r = design$r, tested_terms = tested,
                                         epochs = epochs, device = "cpu", seed = seed, cluster = subject)
  da_sets <- ancombc3_da_eligibility(
    at_risk_prob = offsets$at_risk_prob, y = y, x = design$x, test_terms = design$tests[tested],
    min_observed_positive = cfg$min_observed_positive, min_positive_per_group = cfg$min_positive_per_group,
    min_continuous_var_ratio = cfg$min_continuous_var_ratio,
    min_observed_prevalence = cfg$min_observed_prevalence, min_at_risk_eff_n = cfg$min_at_risk_eff_n,
    min_at_risk_design_support = cfg$min_at_risk_design_support, cluster = subject
  )
  fit <- fit_ancombc3_da(
    y = y, x = design$x, r = design$r, sigma_vae = offsets$sigma_vae,
    abundance_offset = offsets$abundance_offset, structural_offset = offsets$structural_offset,
    sampling_zero_prob = offsets$sampling_zero_prob, at_risk_prob = offsets$at_risk_prob,
    test_terms = tested, fitting_set = da_sets$fitting_set, eligible_taxa_by_term = da_sets$by_term,
    ridge = cfg$ridge_direct, maxit_taxon = cfg$maxit_taxon, d_tol = cfg$d_tol, d_max_iter = cfg$d_max_iter,
    n_starts = cfg$n_starts, d_bound = cfg$d_bound,
    cluster = subject, engine = "compiled"
  )
  fit_offset <- if (is.null(subject)) NULL else
    ancombc3_add_offset_covariance(fit, offsets$abundance_offset, cluster = subject)
  score <- function(f = fit) {
    score_ancombc3_da(
      f$raw, variance_inputs = f$da_variance_inputs, p_adj_method = "BH",
      eligible_taxa_by_term = da_sets$by_term, maximal_overlap_min_count = cfg$da_maximal_overlap_min_count,
      maximal_overlap_bandwidths = cfg$da_maximal_overlap_bandwidths
    )
  }
  # P values of the DA statistics before the empirical-null scale:
  # Z_j = T_j omega_T with the t reference of the taxon.
  unscaled <- function(s) {
    s$p_value <- pmin(2 * exp(stats::pt(-abs(s$statistic * s$adjustment_factor), df = s$df_t, log.p = TRUE)), 1)
    s
  }
  dp_elig <- dp_eligibility(
    structural_absence_weight = offsets$structural_absence_weight, r = design$r, test_terms = tested,
    min_total_absence = cfg$dp_min_total_absence, min_total_at_risk = cfg$dp_min_total_at_risk,
    min_side_absence = cfg$dp_min_side_absence, min_side_at_risk = cfg$dp_min_side_at_risk,
    min_weighted_variance = cfg$dp_min_weighted_variance, cluster = subject
  )
  dp <- fit_ancombc3_dp_zero_pattern(
    y = y, r = design$r, f0 = offsets$sampling_zero_prob, test_terms = tested, alpha = alpha,
    p_adj_method = "BH", ridge = cfg$ridge_dp, eligibility = dp_elig,
    da_fit = fit, x = design$x, cluster = subject, engine = "compiled"
  )
  da <- score()
  list(da = da, da_no_empirical_null = unscaled(da),
       da_offset = if (is.null(fit_offset)) NULL else score(fit_offset),
       dp = dp, by_term = da_sets$by_term, d_hat = fit$d_hat)
}

# Independent samples: one visit per subject.
independent <- simulate_repeated_measures(501L, m = 60L, visits = 1L, sd_subject = 0, sd_presence = 0)
independent$metadata$cont <- as.numeric(scale(stats::rnorm(60L)))
n_data <- 4L
repeated <- lapply(seq_len(n_data), function(k) simulate_repeated_measures(600L + k))
jobs <- c(
  list(omitted = list(kind = "api", data = independent, args = list()),
       null = list(kind = "api", data = independent, args = list(cluster = NULL)),
       own = list(kind = "api", data = independent, args = list(cluster = "subject")),
       never = list(kind = "api", data = independent, args = list(settings = list(da_offset_variance = "never"))),
       rm1_api = list(kind = "api", data = repeated[[1L]], args = list(cluster = "subject")),
       rm1_api_never = list(kind = "api", data = repeated[[1L]],
                            args = list(cluster = "subject", settings = list(da_offset_variance = "never")))),
  stats::setNames(lapply(repeated, function(dat) list(kind = "steps", data = dat, cluster = NULL)),
                  paste0("rm", seq_len(n_data), "_without")),
  stats::setNames(lapply(repeated, function(dat) list(kind = "steps", data = dat, cluster = "subject")),
                  paste0("rm", seq_len(n_data), "_with"))
)
run_job <- function(jb) {
  tested <- if (length(unique(jb$data$metadata$visit)) > 1L) c("group", "visit") else c("group", "cont")
  warns <- character(0)
  value <- withCallingHandlers(
    if (identical(jb$kind, "api")) {
      do.call(ancombc3, c(list(data = jb$data$counts, meta_data = jb$data$metadata,
                               fix_formula = stats::reformulate(tested), tested = tested,
                               epochs = 20L, seed = 1L, verbose = FALSE), jb$args))
    } else {
      run_steps(jb$data, tested, paste(tested, collapse = " + "), cluster = jb$cluster)
    },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  list(value = value, warnings = warns)
}
# Result of the job of the given name, computed at its first request.
cluster_job <- local({
  store <- list()
  function(name) {
    if (is.null(store[[name]])) store[[name]] <<- run_job(jobs[[name]])
    store[[name]]
  }
})

test_that("per-sample scores of the taxon-wise objective at a second-stage fit", {
  ss <- second_stage()
  fit0 <- ss$fit0
  regular <- ss$regular
  expect_true(length(regular) >= 20L, label = "at least 20 taxa with a regular fit and the covariance of psi")
  prior <- fit0$cross_taxon_penalty
  k_omega <- 2L + ncol(x) + ncol(r)
  score_checks <- do.call(rbind, lapply(utils::head(regular, 8L), function(j) {
    f <- fit0$fits[[j]]
    omega <- as.numeric(f$omega)
    kappa <- log(f$sigma)
    dat <- taxon_dat(fit0, j, f$sigma)
    sc <- .ac3_cluster_taxon_scores(omega, dat)
    # Central differences of the cell log-likelihoods: step 1e-5 in omega, and
    # the step 0.05 of the kappa score against the step 1e-3.
    fd_omega <- vapply(seq_len(k_omega), function(k) {
      e <- replace(numeric(k_omega), k, 1e-5)
      (.ac3_cluster_cell_loglik(omega + e, kappa, dat) - .ac3_cluster_cell_loglik(omega - e, kappa, dat)) / 2e-5
    }, numeric(n))
    fd_kappa <- (.ac3_cluster_cell_loglik(omega, kappa + 1e-3, dat) -
                   .ac3_cluster_cell_loglik(omega, kappa - 1e-3, dat)) / 2e-3
    # Gradient of l^pen in omega: the sum of the scores plus the gradient of the
    # ridge and Firth-type penalties.
    ev <- .ac3_da_taxon_objective(omega, dat, grad = TRUE, hess = FALSE)
    pen <- numeric(k_omega)
    pen[dat$i_pen] <- -dat$ridge * omega[dat$i_pen]
    pen[dat$i_presence] <- pen[dat$i_presence] + .ac3_da_firth_penalty(omega[dat$i_presence], dat)$grad
    total_omega <- colSums(sc$u[, seq_len(k_omega)]) + pen
    # Derivative in kappa of the adjusted profile objective: the sum of the
    # kappa scores, the derivative of -(1 / 2) log det H along omega_hat(kappa)
    # (central difference with the step 0.1 of the kappa update), and the
    # derivative of the cross-taxon penalty.
    adjustment <- function(kp) {
      dk <- dat
      dk$sigma <- exp(kp)
      nr <- .ac3_da_newton(dk, omega)
      -sum(log(diag(chol(nr$ev$hess))))
    }
    d_adjustment <- (adjustment(kappa + 0.1) - adjustment(kappa - 0.1)) / 0.2
    d_prior <- if (is.null(prior)) 0 else -(kappa - prior$m) / prior$tau2
    data.frame(
      taxon = names(fit0$fits)[j],
      fd_omega = rel_err(sc$u[, seq_len(k_omega)], fd_omega),
      fd_kappa = rel_err(sc$u[, k_omega + 1L], fd_kappa),
      objective_gradient = max(abs(total_omega - ev$grad)),
      omega_equation = max(abs(total_omega)),
      penalty_gradient = max(abs(pen)),
      kappa_equation = abs(sum(sc$u[, k_omega + 1L]) + d_adjustment + d_prior) *
        sqrt(f$vcov_psi[k_omega + 1L, k_omega + 1L])
    )
  }))
  expect_true(all(score_checks$fd_omega <= 1e-6),
              label = "omega scores match central differences of the cell log-likelihoods")
  expect_true(all(score_checks$fd_kappa <= 1e-2), label = "kappa scores match central differences (step 1e-3) within 1e-2")
  expect_true(all(score_checks$objective_gradient <= 1e-10),
              label = "sum of omega scores plus penalty gradient is the objective gradient")
  expect_true(all(score_checks$omega_equation <= 1e-5),
              label = "sum of omega scores plus penalty gradient is zero at the fit")
  expect_true(all(score_checks$kappa_equation <= 0.1), label = "kappa equation within 0.1 standard deviations of kappa_hat")
})

test_that("cluster-robust covariances of the second stage", {
  ss <- second_stage()
  fit0 <- ss$fit0
  fc1 <- ss$fc1
  ids3 <- rep(seq_len(n / 3L), each = 3L)
  m3 <- max(ids3)
  fc3 <- .ac3_cluster_da_covariance(fit0, inp, fit0$d_hat, cluster = paste0("P", ids3))
  # Definitions evaluated subject by subject.
  M_def <- matrix(0, n, n)
  err_vcov <- 0
  for (j in seq_len(p)) {
    f <- fit0$fits[[j]]
    if (is.null(fc3$fits[[j]]$vcov_model)) next
    has_kappa <- !is.null(f$vcov_psi)
    C <- if (has_kappa) f$vcov_psi else f$vcov
    sc <- .ac3_cluster_taxon_scores(as.numeric(f$omega), taxon_dat(fit0, j, f$sigma), kappa = has_kappa)
    uc <- sweep(sc$u, 2L, colMeans(sc$u), "-")
    B <- matrix(0, ncol(uc), ncol(uc))
    for (s in seq_len(m3)) {
      U_s <- colSums(uc[ids3 == s, , drop = FALSE])
      B <- B + tcrossprod(U_s)
    }
    V_def <- C %*% B %*% C * (m3 / (m3 - 1)) * ((n - 1) / (n - ncol(uc)))
    V_new <- if (has_kappa) fc3$fits[[j]]$vcov_psi else fc3$fits[[j]]$vcov
    err_vcov <- max(err_vcov, max(abs(V_new - V_def)) / max(abs(V_def)))
    block <- .ac3_da_taxon_bias_information(f, y[, j], x, r)
    if (is.null(block)) next
    P <- t(block$I_psi_d) %*% block$C
    A <- matrix(0, n, m3)
    for (i in seq_len(n)) {
      a_i <- -drop(P %*% uc[i, ])
      a_i[i] <- a_i[i] + sc$s[i]
      A[, ids3[i]] <- A[, ids3[i]] + a_i
    }
    M_def <- M_def + tcrossprod(A)
  }
  V_model <- fit0$da_variance_inputs$V
  V_def <- V_model %*% M_def %*% V_model * (m3 / (m3 - 1))
  V_cr <- fc3$da_variance_inputs$V
  expect_true(err_vcov <= 1e-10, label = "covariance of psi equals c J^-1 B J^-1 by subject")
  expect_true(max(abs(V_cr - V_def)) / max(abs(V_def)) <= 1e-8,
              label = "bias covariance equals V M V m / (m - 1) with M by subject")
  expect_true(isTRUE(all.equal(.ac3_cluster_factor(20, 60, 8), (20 / 19) * (59 / 52), tolerance = 1e-15)) &&
                identical(.ac3_cluster_factor(20, 8, 8), 20 / 19) && identical(.ac3_cluster_factor(20), 20 / 19),
              label = "finite-sample factor {m / (m - 1)} {(N - 1) / (N - k)}; m / (m - 1) when N <= k")
  Q <- .ac3_da_residualized_tested_basis(x, test_terms)
  expect_true(max(abs(V_cr %*% cbind(1, Q))) <= 1e-8 * max(abs(V_cr)),
              label = "bias covariance satisfies both identifying constraints")
  expect_true(min(eigen(V_cr, symmetric = TRUE, only.values = TRUE)$values) >= -1e-8 * max(abs(V_cr)),
              label = "bias covariance is positive semidefinite")
  expect_true(identical(fc3$raw$df_count,
                        unname(colSums(rowsum((y > 0) * 1, ids3) > 0)[match(fc3$raw$taxon, colnames(y))])),
              label = "df_count is the number of subjects with a positive count")
  idx_T <- 1L + match(fc3$raw$term, colnames(x))
  v_T <- vapply(seq_len(nrow(fc3$raw)), function(rr) {
    fc3$fits[[fc3$raw$taxon[rr]]]$vcov[idx_T[rr], idx_T[rr]]
  }, numeric(1))
  expect_true(isTRUE(all.equal(fc3$raw$var_raw_model, v_T, tolerance = 0)) &&
                isTRUE(all.equal(fc3$raw$var_raw_marginal, fc3$raw$var_raw_model + fc3$raw$var_raw_bias,
                                 tolerance = 1e-14)),
              label = "raw variances are those of the cluster-robust covariances")

  # Every sample its own subject.
  expect_true(identical(fc1$d_hat, fit0$d_hat) && identical(fc1$raw$beta_raw, fit0$raw$beta_raw) &&
                identical(fc1$sigma_hat, fit0$sigma_hat),
              label = "own subjects: estimates unchanged")
  expect_true(identical(fc1$raw$var_raw_model_based, fit0$raw$var_raw_model) &&
                identical(fc1$da_variance_inputs$V_model, fit0$da_variance_inputs$V) &&
                identical(fc1$da_variance_inputs$G, fit0$da_variance_inputs$G) &&
                all(vapply(seq_len(p), function(j) {
                  is.null(fc1$fits[[j]]$vcov_model) ||
                    (identical(fc1$fits[[j]]$vcov_model, fit0$fits[[j]]$vcov) &&
                       identical(fc1$fits[[j]]$vcov_psi_model, fit0$fits[[j]]$vcov_psi))
                }, logical(1))),
              label = "own subjects: model-based matrices kept, sensitivities unchanged")
  expect_true(identical(as.numeric(fc1$raw$df_count), as.numeric(fc1$raw$n_positive)),
              label = "own subjects: df_count equals n_positive")
  ratio_model <- fc1$raw$se_raw_model / fc1$raw$se_raw_model_based
  ratio_bias <- sqrt(fc1$raw$var_raw_bias / fc1$raw$var_raw_bias_model_based)
  sc0 <- score_ancombc3_da(fit0$raw, fit0$da_variance_inputs, eligible_taxa_by_term = fit0$eligible_taxa_by_term)
  sc1 <- score_ancombc3_da(fc1$raw, fc1$da_variance_inputs, eligible_taxa_by_term = fc1$eligible_taxa_by_term)
  ratio_corrected <- sc1$se_corrected / sc0$se_corrected
  expect_true(stats::median(ratio_model) >= 0.85 && stats::median(ratio_model) <= 1.10,
              label = "own subjects: median standard-error ratio at fixed bias in [0.85, 1.10]")
  expect_true(stats::median(ratio_bias) >= 0.80 && stats::median(ratio_bias) <= 1.20,
              label = "own subjects: median ratio of the bias component in [0.80, 1.20]")
  expect_true(stats::median(ratio_corrected) >= 0.85 && stats::median(ratio_corrected) <= 1.10,
              label = "own subjects: median ratio of the corrected standard error in [0.85, 1.10]")
  expect_true(identical(sc1$beta_corrected, sc0$beta_corrected) || max(abs(sc1$delta - sc0$delta)) <= 0.1,
              label = "own subjects: common shift within 0.1 of the model-based analysis")
  expect_true(stats::cor(sc0$statistic_normal_equivalent, sc1$statistic_normal_equivalent) >= 0.95,
              label = "own subjects: DA statistics correlate (>= 0.95) with the model-based ones")
})

test_that("subject-level helpers: withheld samples, eligibility, DP scores, and DP test", {
  ss <- second_stage()
  fit0 <- ss$fit0
  regular <- ss$regular
  fc1 <- ss$fc1
  # Withheld samples of the first stage.
  set.seed(11L)
  v0 <- sort(sample.int(50L, 5L))
  state0 <- .Random.seed
  set.seed(11L)
  v1 <- .ac3_cluster_validation_split(paste0("s", 1:50), 5L)
  expect_true(identical(v0, v1) && identical(state0, .Random.seed),
              label = "own subjects: withheld samples and random-number state of sample.int(n, n_val)")
  size_s <- c(rep(1L, 10L), rep(3L, 20L), rep(6L, 5L))
  subj <- sample(rep(seq_along(size_s), times = size_s))
  set.seed(12L)
  v2 <- .ac3_cluster_validation_split(subj, 10L)
  expect_true(all(subj[-v2] %in% setdiff(subj, subj[v2])) && length(v2) >= 5L && length(v2) <= 15L &&
                !is.unsorted(v2),
              label = "withheld samples are whole subjects, about n_val samples")

  # DA eligibility: subject counts for a subject-constant binary covariate.
  ids_e <- rep(1:20, each = 3L)
  grp_e <- as.numeric(ids_e > 10L)
  vis_e <- rep(0:2, times = 20L)
  x_e <- cbind(grp = grp_e, vis = vis_e, alt = rep(c(0, 1, 0, 1, 1, 0), times = 10L))
  set.seed(13L)
  y_e <- matrix(stats::rpois(60L * 12L, 3) *
                  stats::rbinom(60L * 12L, 1L, rep(seq(0.15, 0.9, length.out = 12L), each = 60L)),
                60L, 12L, dimnames = list(NULL, paste0("T", 1:12)))
  q_e <- matrix(stats::runif(60L * 12L), 60L, 12L)
  term_grp <- list(term = "grp", type = "binary", columns = "grp")
  term_alt <- list(term = "alt", type = "binary", columns = "alt")
  el0 <- ancombc3_da_eligibility_term(q_e, y_e, x_e, term_grp)
  el1 <- ancombc3_da_eligibility_term(q_e, y_e, x_e, term_grp, cluster = ids_e)
  pos_by_subject <- rowsum((y_e > 0) * 1, ids_e) > 0
  q_by_subject <- apply(q_e, 2L, function(v) tapply(v, ids_e, max))
  expect_true(identical(ancombc3_da_eligibility_term(q_e, y_e, x_e, term_grp, cluster = NULL), el0),
              label = "DA eligibility: cluster = NULL gives the result without the argument")
  expect_true(isTRUE(all.equal(el1$min_positive_by_group,
                               unname(pmin(colSums(pos_by_subject[1:10, ]), colSums(pos_by_subject[11:20, ]))))) &&
                isTRUE(all.equal(el1$at_risk_design_support,
                                 unname(pmin(colSums(q_by_subject[1:10, ]), colSums(q_by_subject[11:20, ]))))) &&
                all(el1$min_positive_by_group <= el0$min_positive_by_group) &&
                identical(el1$observed_positive_n, el0$observed_positive_n) &&
                identical(el1$observed_prevalence, el0$observed_prevalence) &&
                identical(el1$at_risk_eff_n, el0$at_risk_eff_n),
              label = "DA eligibility: group quantities by subject, totals by sample")
  expect_true(identical(ancombc3_da_eligibility_term(q_e, y_e, x_e, term_alt, cluster = ids_e),
                        ancombc3_da_eligibility_term(q_e, y_e, x_e, term_alt)),
              label = "DA eligibility: covariate varying within subject is unchanged")
  expect_true(identical(ancombc3_da_eligibility(q_e, y_e, x_e, c("grp", "alt"), cluster = ids_e)$by_term$grp,
                        el1$taxon[el1$da_eligible]),
              label = "DA eligibility: cluster is passed by ancombc3_da_eligibility()")
  # DP eligibility: subject-level side masses.
  s_e <- ifelse(y_e == 0, q_e, 0)
  dp0 <- dp_eligibility(s_e, x_e, test_terms = c("grp", "alt"))
  dp1 <- dp_eligibility(s_e, x_e, test_terms = c("grp", "alt"), cluster = ids_e)
  s_by_subject <- apply(s_e, 2L, function(v) tapply(v, ids_e, max))
  a_by_subject <- apply(1 - s_e, 2L, function(v) tapply(v, ids_e, max))
  g1 <- dp1[dp1$test_term == "grp", ]
  expect_true(identical(dp_eligibility(s_e, x_e, test_terms = c("grp", "alt"), cluster = NULL), dp0),
              label = "DP eligibility: cluster = NULL gives the result without the argument")
  expect_true(isTRUE(all.equal(g1$min_side_absence_mass,
                               unname(pmin(colSums(s_by_subject[1:10, ]), colSums(s_by_subject[11:20, ]))))) &&
                isTRUE(all.equal(g1$min_side_at_risk_mass,
                                 unname(pmin(colSums(a_by_subject[1:10, ]), colSums(a_by_subject[11:20, ]))))) &&
                identical(g1$absence_mass, dp0$absence_mass[dp0$test_term == "grp"]) &&
                identical(dp1[dp1$test_term == "alt", ], dp0[dp0$test_term == "alt", ]),
              label = "DP eligibility: side masses by subject, totals and varying covariate unchanged")

  # DP scores and DP test.
  D <- cbind("(Intercept)" = 1, r)
  j_dp <- regular[which.max(colSums(y[, regular] == 0))]
  Z <- as.numeric(y[, j_dp] == 0)
  f0_j <- pmin(pmax(fit0$fits[[j_dp]]$f0, 1e-8), 1 - 1e-8)
  coef_test <- c(0.5, 0.2, -0.3)
  nll0 <- .dp_zero_pattern_nll(Z, f0_j, D, ridge = 0, firth = FALSE)
  grad_fd <- vapply(seq_along(coef_test), function(k) {
    e <- replace(numeric(3L), k, 1e-5)
    (nll0(coef_test + e) - nll0(coef_test - e)) / 2e-5
  }, numeric(1))
  expect_true(max(abs(colSums(.ac3_cluster_dp_scores(Z, f0_j, D, coef_test)) + grad_fd)) <= 1e-6 * max(abs(grad_fd)),
              label = "DP scores sum to the gradient of the zero-pattern log-likelihood")
  dp_a <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0_first_stage, test_terms = test_terms,
                                       da_fit = fit0, x = x)
  dp_n <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0_first_stage, test_terms = test_terms,
                                       da_fit = fit0, x = x, cluster = NULL)
  dp_c <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0_first_stage, test_terms = test_terms,
                                       da_fit = fc1, x = x, cluster = rownames(y))
  expect_true(identical(dp_a, dp_n), label = "DP test: cluster = NULL gives the result without the argument")
  expect_true(identical(dp_c$estimate, dp_a$estimate) && identical(dp_c$statistic_lr, dp_a$statistic_lr) &&
                identical(dp_c$var_dp[is.finite(dp_a$var_dp)], dp_a$var_dp[is.finite(dp_a$var_dp)]),
              label = "DP test, own subjects: estimates, T^LR, and model-based variance unchanged")
  ok_c <- is.finite(dp_c$var_dp_robust)
  total_c <- dp_c$var_dp_robust + ifelse(is.finite(dp_c$var_propagated), dp_c$var_propagated, 0)
  prop_c <- ifelse(is.finite(dp_c$var_propagated), dp_c$var_propagated, 0)
  ratio_c <- pmin(dp_c$var_dp / (dp_c$var_dp + prop_c), dp_c$var_dp / total_c)
  expect_true(all(ok_c) && isTRUE(all.equal(dp_c$variance_ratio, ratio_c, tolerance = 1e-12)) &&
                isTRUE(all.equal(dp_c$statistic, dp_c$variance_ratio * dp_c$statistic_lr, tolerance = 1e-12)),
              label = "DP test: statistic is T^LR times the smaller of the two variance ratios")
  ratio_dp <- sqrt(dp_c$var_dp_robust / dp_c$var_dp)
  expect_true(stats::median(ratio_dp) >= 0.85 && stats::median(ratio_dp) <= 1.10,
              label = "DP test, own subjects: median standard-error ratio in [0.85, 1.10]")
})

test_that("ancombc3(): cluster = NULL, own subjects, equality with the steps, and input checks", {
  skip_if_no_torch()
  skip_if_not_extended()
  torch::torch_set_num_threads(1L)
  strip <- function(res) {
    res$settings$elapsed_seconds <- NULL
    res$settings$stage_seconds <- NULL
    res$settings$call <- NULL
    res
  }
  r_om <- strip(cluster_job("omitted")$value)
  r_null <- strip(cluster_job("null")$value)
  r_own <- cluster_job("own")$value
  expect_true(identical(r_om, r_null), label = "cluster = NULL gives the output of a call without the argument")
  r_never <- strip(cluster_job("never")$value)
  r_never$settings$method$da_offset_variance <- r_om$settings$method$da_offset_variance
  expect_true(identical(r_never, r_om) && identical(cluster_job("never")$value$settings$method$da_offset_variance, "never"),
              label = "without cluster: the default da_offset_variance gives the output of \"never\"")
  expect_true(identical(r_own$sample_bias, r_om$sample_bias) && identical(r_own$sigma, r_om$sigma) &&
                identical(r_own$res$dp_estimate, r_om$res$dp_estimate) &&
                identical(r_own$res$da_eligible, r_om$res$da_eligible) &&
                identical(r_own$res$dp_eligible, r_om$res$dp_eligible),
              label = "own subjects: bias, sigma, DP estimates, and eligibility unchanged")
  # The DA estimate is the uncorrected coefficient minus the common-shift
  # estimate, a function of the standard errors.
  shift_difference <- tapply(r_own$res$da_estimate - r_om$res$da_estimate, r_om$res$covariate,
                             function(v) range(v, na.rm = TRUE))
  expect_true(all(vapply(shift_difference, function(v) diff(v) <= 1e-12 && max(abs(v)) <= 0.05, logical(1))),
              label = "own subjects: DA estimates differ by one constant per covariate (<= 0.05)")
  expect_true(identical(r_own$settings$cluster, "subject") && identical(r_own$settings$n_subjects, 60L) &&
                is.null(r_om$settings$cluster) && identical(names(r_own), names(r_om)),
              label = "settings record cluster and the number of subjects")
  da_ok <- is.finite(r_om$res$da_p)
  dp_ok <- is.finite(r_om$res$dp_p)
  se_ratio <- stats::median(r_own$res$da_se[da_ok] / r_om$res$da_se[da_ok])
  cor_da <- stats::cor(stats::qnorm(r_own$res$da_p[da_ok] / 2), stats::qnorm(r_om$res$da_p[da_ok] / 2))
  cor_dp <- stats::cor(sqrt(r_own$res$dp_statistic[dp_ok]), sqrt(r_om$res$dp_statistic[dp_ok]))
  expect_true(se_ratio >= 0.85 && se_ratio <= 1.10, label = "own subjects: median ratio of da_se in [0.85, 1.10]")
  expect_true(cor_da >= 0.90 && cor_dp >= 0.90, label = "own subjects: DA and DP statistics correlate (>= 0.90)")

  # ancombc3() and its steps on the first repeated-measures data set, at the
  # default settings (offset term included) and with da_offset_variance = "never".
  st <- cluster_job("rm1_with")$value
  equals_steps <- function(res, da_tab) {
    api <- res$res
    da_st <- da_tab[match(paste(api$taxon, api$covariate), paste(da_tab$taxon, da_tab$term)), ]
    dp_st <- st$dp[match(paste(api$taxon, api$covariate), paste(st$dp$taxon, st$dp$test_term)), ]
    el <- api$da_eligible & is.finite(api$da_p)
    de <- api$dp_eligible
    identical(unname(res$sample_bias), unname(st$d_hat)) &&
      identical(api$da_estimate[el], da_st$beta_corrected[el]) &&
      identical(api$da_se[el], da_st$se_scaled[el]) && identical(api$da_p[el], da_st$p_value[el]) &&
      identical(de, dp_st$dp_eligible %in% TRUE) && identical(api$dp_p[de], dp_st$p_value[de]) &&
      identical(api$da_eligible, unname(mapply(function(tx, cv) tx %in% st$by_term[[cv]], api$taxon,
                                               api$covariate)))
  }
  expect_true(equals_steps(cluster_job("rm1_api")$value, st$da_offset),
              label = "repeated measures: ancombc3() equals its steps with cluster and the offset term")
  expect_true(equals_steps(cluster_job("rm1_api_never")$value, st$da),
              label = "repeated measures: \"never\" equals the steps with cluster without the offset term")
  se_off <- cluster_job("rm1_api")$value$res$da_se
  se_never <- cluster_job("rm1_api_never")$value$res$da_se
  both_se <- is.finite(se_off) & is.finite(se_never)
  expect_true(sum(both_se) > 0L && !isTRUE(all.equal(se_off[both_se], se_never[both_se])),
              label = "repeated measures: the default with cluster includes the offset term")

  rm1 <- repeated[[1L]]
  quick <- function(metadata = rm1$metadata, cluster = "subject") {
    ancombc3(data = rm1$counts, meta_data = metadata, fix_formula = ~ group + visit, tested = c("group", "visit"),
             epochs = 1L, verbose = FALSE, cluster = cluster)
  }
  expect_true(ac3_errs(quick(cluster = "patient"), "not a column of meta_data"),
              label = "cluster is not a column of metadata")
  expect_true(ac3_errs(quick(cluster = c("subject", "visit")), "name of one column"),
              label = "cluster is not a single name")
  meta_na <- rm1$metadata
  meta_na$subject[5L] <- NA
  expect_true(ac3_errs(quick(metadata = meta_na), "missing values for 1 analyzed sample"), label = "missing subject")
  meta_few <- rm1$metadata
  meta_few$subject <- paste0("P", rep(1:8, each = 15L))
  meta_few$group <- as.numeric(rep(1:8, each = 15L) > 4L)
  expect_true(ac3_errs(quick(metadata = meta_few), "at least 10 subjects"), label = "fewer than 10 subjects")
  meta_grp <- rm1$metadata
  meta_grp$subject <- paste0("P", rep(1:12, each = 10L))
  meta_grp$group <- as.numeric(rep(1:12, each = 10L) > 8L)
  expect_true(ac3_errs(suppressWarnings(quick(metadata = meta_grp)), "group with 4 subject"),
              label = "binary group with fewer than five subjects")
  msgs <- character(0)
  warns <- character(0)
  subj_24 <- withCallingHandlers(
    .ac3_cluster_subjects(data.frame(id = rep(1:24, times = c(rep(2L, 12L), rep(4L, 12L)))), "id", verbose = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  expect_true(length(subj_24) == 72L && any(grepl("24 subjects", warns)) &&
                any(grepl("24 subjects [(]cluster 'id'[)], 2 to 4 samples per subject", msgs)),
              label = "fewer than 30 subjects: warning; message with the number of subjects")
  expect_true(!any(grepl("subjects", cluster_job("rm1_api")$warnings)) &&
                !any(grepl("subjects", cluster_job("own")$warnings)),
              label = "no warning about subjects with 40 and with 60 subjects")
})

test_that("repeated measures without covariate effects: rejection rates without and with cluster", {
  skip_if_no_torch()
  skip_if_not_extended()
  torch::torch_set_num_threads(1L)
  pooled <- function(suffix, covariate) {
    rr <- lapply(paste0("rm", seq_len(n_data), suffix), function(nm) cluster_job(nm)$value)
    da_rows <- function(z, tab) tab$term == covariate & tab$taxon %in% z$by_term[[covariate]]
    list(
      da = unlist(lapply(rr, function(z) z$da$p_value[da_rows(z, z$da)])),
      da_no_en = unlist(lapply(rr, function(z) {
        z$da_no_empirical_null$p_value[da_rows(z, z$da_no_empirical_null)]
      })),
      omega = vapply(rr, function(z) z$da$adjustment_factor[da_rows(z, z$da)][1L], numeric(1)),
      dp = unlist(lapply(rr, function(z) z$dp$p_value[z$dp$test_term == covariate & z$dp$dp_eligible %in% TRUE]))
    )
  }
  rate <- function(pv) mean(pv <= 0.05, na.rm = TRUE)
  rates <- do.call(rbind, lapply(c("group", "visit"), function(cv) {
    do.call(rbind, lapply(c("_without", "_with"), function(sf) {
      pl <- pooled(sf, cv)
      data.frame(covariate = cv, cluster = sub("^_", "", sf), da_tests = sum(is.finite(pl$da)),
                 da_rate = rate(pl$da), da_rate_no_empirical_null = rate(pl$da_no_en),
                 mean_empirical_null_scale = mean(pl$omega), dp_tests = sum(is.finite(pl$dp)),
                 dp_rate = rate(pl$dp))
    }))
  }))
  rate_of <- function(covariate, cluster, column) {
    rates[[column]][rates$covariate == covariate & rates$cluster == cluster]
  }
  # Group is the between-subject covariate and visit the within-subject
  # covariate; each rate pools about 320 tests (binomial standard error 0.012
  # at a rate of 0.05).
  expect_true(rate_of("group", "without", "da_rate_no_empirical_null") > 0.07 &&
                rate_of("group", "without", "dp_rate") > 0.07,
              label = "without cluster: between-subject DA (no empirical null) and DP rates > 0.07")
  with_cluster <- rates[rates$cluster == "with", c("da_rate", "da_rate_no_empirical_null", "dp_rate")]
  expect_true(all(as.matrix(with_cluster) <= 0.10), label = "with cluster: every rejection rate is at most 0.10")
  expect_true(rate_of("group", "with", "da_rate_no_empirical_null") <
                rate_of("group", "without", "da_rate_no_empirical_null") &&
                rate_of("group", "with", "dp_rate") < rate_of("group", "without", "dp_rate"),
              label = "between-subject covariate: rates with cluster below the rates without")
})
