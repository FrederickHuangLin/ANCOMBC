# Bias-corrected abundances of ancombc3() (elements samp_frac,
# bias_correct_log_table, and struc_absence_prob), on one simulated data set
# with a sample and a taxon without a positive count.

# Simulated counts: n samples, p taxa, a binary and a continuous covariate, a
# sample-specific bias correlated with the binary covariate, abundance effects
# of the binary covariate on the first eight taxa, and structural zeros.
# Sample S3 and taxon T5 have no positive count.
set.seed(20261111L)
n <- 60L
p <- 60L
group <- rep(c("a", "b"), length.out = n)
age <- as.numeric(scale(stats::rnorm(n)))
alpha <- stats::rnorm(p, mean = 3, sd = 1)
beta <- c(rep(1, 8L), rep(0, p - 8L))
bias <- stats::rnorm(n, sd = 0.4) + 0.5 * (group == "b")
log_mu <- matrix(alpha, n, p, byrow = TRUE) + outer(as.numeric(group == "b"), beta) + bias
y <- matrix(stats::rpois(n * p, exp(log_mu + matrix(stats::rnorm(n * p, sd = 0.5), n, p))), n, p) *
  matrix(stats::rbinom(n * p, 1L, 0.8), n, p)
y[3L, ] <- 0
y[, 5L] <- 0
dimnames(y) <- list(paste0("S", seq_len(n)), paste0("T", seq_len(p)))
metadata <- data.frame(group = group, age = age, row.names = rownames(y), stringsAsFactors = FALSE)
kept_samples <- rownames(y)[-3L]
kept_taxa <- colnames(y)[-5L]
y_kept <- y[kept_samples, kept_taxa]

# ancombc3() with taxa in rows and with samples in rows.
orientation_runs <- ac3_cached(function() {
  run_job <- function(counts, taxa_are_rows) {
    suppressWarnings(ancombc3(data = counts, taxa_are_rows = taxa_are_rows, meta_data = metadata,
                              fix_formula = ~ group + age, tested = c("group", "age"), epochs = 20L, seed = 1L,
                              verbose = FALSE))
  }
  list(taxa_rows = run_job(t(y), TRUE), sample_rows = run_job(y, FALSE))
})

test_that("elements are returned and do not depend on the orientation of the count table", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  runs <- orientation_runs()
  tr <- runs$taxa_rows
  sr <- runs$sample_rows
  expect_true(all(c("res", "sample_bias", "samp_frac", "bias_correct_log_table", "struc_absence_prob", "sigma",
                    "settings") %in% names(tr)),
              label = "samp_frac, bias_correct_log_table, and struc_absence_prob are returned")
  expect_true(identical(tr$res, sr$res) && identical(tr$sample_bias, sr$sample_bias) &&
                identical(tr$sigma, sr$sigma) && identical(tr$samp_frac, sr$samp_frac),
              label = "res, sample_bias, sigma, and samp_frac do not depend on the orientation")
})

test_that("bias-corrected log abundances: orientation, dimnames, and the constant samp_frac per sample", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  runs <- orientation_runs()
  tr <- runs$taxa_rows
  sr <- runs$sample_rows
  M <- tr$bias_correct_log_table
  b <- tr$samp_frac
  expect_true(identical(dimnames(M), list(kept_taxa, kept_samples)),
              label = "taxa in rows: analyzed taxa by analyzed samples")
  expect_true(identical(dimnames(sr$bias_correct_log_table), list(kept_samples, kept_taxa)) &&
                identical(sr$bias_correct_log_table, t(M)),
              label = "samples in rows: the transposed matrix")
  expect_true(all(is.finite(M)), label = "every cell has a value")
  P <- tr$struc_absence_prob
  expect_true(identical(dimnames(P), dimnames(M)) && identical(sr$struc_absence_prob, t(P)) &&
                all(P[t(y_kept > 0)] == 0) && all(P[t(y_kept == 0)] >= 0 & P[t(y_kept == 0)] <= 1) &&
                any(P > 0.5) && any(P[t(y_kept == 0)] < 0.5),
              label = "struc_absence_prob: zero for positive counts, in [0, 1] for zero counts")
  diff <- log(t(y_kept)) - M
  diff[t(y_kept == 0)] <- NA
  spread <- apply(diff, 2L, function(v) diff(range(v[is.finite(v)])))
  level <- apply(diff, 2L, function(v) mean(v[is.finite(v)]))
  expect_true(max(spread) <= 1e-12, label = "log(y) minus the corrected value is constant within a sample")
  expect_true(identical(names(b), kept_samples) && max(abs(level - b)) <= 1e-12,
              label = "that constant is samp_frac")
})

test_that("samp_frac equals sample_bias plus the common shifts times the centered tested columns", {
  skip_if_no_torch()
  skip_if_not_extended()
  torch::torch_set_num_threads(1L)
  runs <- orientation_runs()
  tr <- runs$taxa_rows
  b <- tr$samp_frac
  # Steps of ancombc3() up to the scored DA table, at its settings and its
  # default engine.
  cfg <- ancombc3_default_settings()
  meta_kept <- metadata[kept_samples, , drop = FALSE]
  meta_kept$group <- as.numeric(meta_kept$group == "b")
  design <- ancombc3_make_design(meta_kept, x_formula = "group + age", test_terms = c("group", "age"),
                                 r_formula = "group + age")
  y_steps <- y_kept
  storage.mode(y_steps) <- "double"
  offsets <- ancombc3_fitted_vae_offsets(y = y_steps, x = design$x, r = design$r, tested_terms = c("group", "age"),
                                         epochs = 20L, device = "cpu", seed = 1L)
  da_sets <- ancombc3_da_eligibility(
    at_risk_prob = offsets$at_risk_prob, y = y_steps, x = design$x, test_terms = design$tests,
    min_observed_positive = cfg$min_observed_positive, min_positive_per_group = cfg$min_positive_per_group,
    min_continuous_var_ratio = cfg$min_continuous_var_ratio,
    min_observed_prevalence = cfg$min_observed_prevalence, min_at_risk_eff_n = cfg$min_at_risk_eff_n,
    min_at_risk_design_support = cfg$min_at_risk_design_support
  )
  fit <- fit_ancombc3_da(
    y = y_steps, x = design$x, r = design$r, sigma_vae = offsets$sigma_vae,
    abundance_offset = offsets$abundance_offset, structural_offset = offsets$structural_offset,
    sampling_zero_prob = offsets$sampling_zero_prob, at_risk_prob = offsets$at_risk_prob,
    test_terms = c("group", "age"), fitting_set = da_sets$fitting_set, eligible_taxa_by_term = da_sets$by_term,
    ridge = cfg$ridge_direct, maxit_taxon = cfg$maxit_taxon, d_tol = cfg$d_tol, d_max_iter = cfg$d_max_iter,
    n_starts = cfg$n_starts, d_bound = cfg$d_bound,
    n_cores = 2L, engine = "compiled"
  )
  scored <- score_ancombc3_da(
    fit$raw, variance_inputs = fit$da_variance_inputs, p_adj_method = "BH",
    eligible_taxa_by_term = da_sets$by_term, maximal_overlap_min_count = cfg$da_maximal_overlap_min_count,
    maximal_overlap_bandwidths = cfg$da_maximal_overlap_bandwidths
  )
  delta <- vapply(c("group", "age"), function(tn) unique(scored$delta[scored$term == tn]), numeric(1))
  expect_true(identical(unname(fit$d_hat), unname(tr$sample_bias)), label = "recomputed sample bias equals sample_bias")
  b_expected <- tr$sample_bias + delta[["group"]] * (design$x[, "group"] - mean(design$x[, "group"])) +
    delta[["age"]] * (design$x[, "age"] - mean(design$x[, "age"]))
  expect_true(max(abs(b - b_expected)) <= 1e-12, label = "samp_frac = sample_bias + sum_T delta_T (x_T - mean(x_T))")
  expect_true(abs(delta[["group"]]) > 0.05 && max(abs(b - tr$sample_bias)) > 0.01,
              label = "the common shift of the binary covariate contributes to the total bias")
  da_group <- tr$res[tr$res$covariate == "group" & tr$res$da_eligible, ]
  raw_group <- fit$raw[fit$raw$term == "group", ]
  expect_true(max(abs(da_group$da_estimate -
                        (raw_group$beta_raw[match(da_group$taxon, raw_group$taxon)] - delta[["group"]]))) <= 1e-12,
              label = "da_estimate is the uncorrected coefficient minus the same common shift")
})

test_that(".ac3_bias_corrected_abundance() on a small input with known values", {
  y4 <- matrix(c(1, 0, 7, 20, 3, 0), 3L, 2L, dimnames = list(c("s1", "s2", "s3"), c("t1", "t2")))
  x4 <- cbind(g = c(0, 1, 1), z = c(-1, 0, 4))
  scored4 <- data.frame(term = c("g", "g"), delta = c(0.3, 0.3))
  bc4 <- .ac3_bias_corrected_abundance(y4, d_hat = c(0.1, -0.2, 0.1), x = x4, tested_cols = c("g", "z"),
                                       scored = scored4, taxa_in_rows = TRUE)
  b4 <- c(s1 = 0.1 + 0.3 * (0 - 2 / 3), s2 = -0.2 + 0.3 / 3, s3 = 0.1 + 0.3 / 3)
  expect_true(isTRUE(all.equal(bc4$total_sample_bias, b4, tolerance = 1e-14)),
              label = "total bias; a tested column without a scored estimate has shift zero")
  expected4 <- rbind(t1 = c(s1 = log(1) - b4[["s1"]], s2 = NA, s3 = log(7) - b4[["s3"]]),
                     t2 = c(s1 = log(20) - b4[["s1"]], s2 = log(3) - b4[["s2"]], s3 = NA))
  expect_true(isTRUE(all.equal(bc4$bias_corrected_log_abundance, expected4, tolerance = 1e-14)) &&
                identical(is.na(bc4$structural_absence_probability), is.na(expected4)),
              label = "without the fitted quantities: corrected values and missing values for the zero counts")
  mu4 <- matrix(c(2, 0.7, 5, 12, 4, 0.2), 3L, 2L); pi4 <- matrix(c(0.9, 0.4, 0.8, 0.95, 0.7, 0.3), 3L, 2L)
  sigma4 <- c(0.6, 1.4)
  f04 <- sapply(1:2, function(j) .ac3_pln_zero_prob(log(mu4[, j]), sigma4[j]))
  bc4f <- .ac3_bias_corrected_abundance(y4, d_hat = c(0.1, -0.2, 0.1), x = x4, tested_cols = c("g", "z"),
                                        scored = scored4, taxa_in_rows = TRUE, mu = mu4, pi = pi4, f0 = f04,
                                        sigma = sigma4)
  # Posterior mean of u = log Lambda given a zero count: density proportional
  # to exp(-e^u) N(u; log mu - sigma^2 / 2, sigma^2).
  post_mean <- function(mu, sigma) {
    ctr <- log(mu) - sigma^2 / 2
    w <- function(u) exp(-exp(u)) * stats::dnorm(u, ctr, sigma)
    stats::integrate(function(u) u * w(u), -Inf, Inf, rel.tol = 1e-12)$value /
      stats::integrate(w, -Inf, Inf, rel.tol = 1e-12)$value
  }
  zero4 <- rbind(c(2L, 1L), c(3L, 2L))
  imputed <- apply(zero4, 1L, function(k) post_mean(mu4[k[1L], k[2L]], sigma4[k[2L]]) - b4[[k[1L]]])
  absence <- apply(zero4, 1L, function(k) {
    (1 - pi4[k[1L], k[2L]]) / ((1 - pi4[k[1L], k[2L]]) + pi4[k[1L], k[2L]] * f04[k[1L], k[2L]])
  })
  got <- t(bc4f$bias_corrected_log_abundance)
  got_p <- t(bc4f$structural_absence_probability)
  expect_true(max(abs(got[zero4] - imputed)) <= 1e-7 &&
                isTRUE(all.equal(got[y4 > 0], t(expected4)[y4 > 0], tolerance = 1e-14)),
              label = "zero counts: posterior mean log intensity minus the total bias; positive counts unchanged")
  expect_true(max(abs(got_p[zero4] - absence)) <= 1e-14 && all(got_p[y4 > 0] == 0),
              label = "probability of structural absence: formula for zero counts, zero for positive counts")
})
