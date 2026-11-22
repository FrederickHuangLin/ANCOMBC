# Multigroup analysis: level support and eligibility, the DA and DP channels,
# the global statistic, Dunnett-type p-values, the mixed-directional-FDR
# procedures, and ancombc3() with a three-level group.

rel_diff <- function(a, b) {
  ok <- is.finite(a) | is.finite(b)
  if (!any(ok)) return(0)
  if (any(is.finite(a[ok]) != is.finite(b[ok]))) return(Inf)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-8))
}
cfg <- ancombc3_default_settings()

# Small simulated data set with a three-level group for the DA and DP
# channels: n samples, p taxa, a continuous covariate, abundance effects of
# both group columns on taxa 1-4 and of level C on taxa 5-8.
set.seed(20261110L)
n <- 45L
p <- 30L
grp <- factor(rep(c("A", "B", "C"), each = n / 3L))
cont <- as.numeric(scale(stats::rnorm(n)))
x <- cbind(cont = cont, grpB = as.numeric(grp == "B"), grpC = as.numeric(grp == "C"))
r <- x
alpha_j <- stats::rnorm(p, mean = 2.5, sd = 0.8)
beta <- matrix(0, p, 3L)
beta[1:4, 2L] <- 1
beta[1:4, 3L] <- 2
beta[5:8, 3L] <- -1.5
bias <- stats::rnorm(n, sd = 0.3)
log_mu <- matrix(alpha_j, n, p, byrow = TRUE) + x %*% t(beta) + bias
lambda <- exp(log_mu + matrix(stats::rnorm(n * p, sd = 0.5), n, p))
present <- matrix(stats::rbinom(n * p, 1L, 0.9), n, p)
y <- matrix(stats::rpois(n * p, lambda), n, p) * present
colnames(y) <- paste0("T", seq_len(p))
rownames(y) <- paste0("S", seq_len(n))
taxa <- colnames(y)
sigma_vae <- rep(0.5, p)
mu_first_stage <- matrix(pmax(colMeans(y), 0.1), n, p, byrow = TRUE)
sampling_zero_prob <- (2 / (2 + mu_first_stage))^2
at_risk_prob <- matrix(0.9, n, p)
absence_weight <- ifelse(y == 0, 0.6, 0)
cols <- c("grpB", "grpC")
# Nonzero abundance offset and 15 units of three samples within a level for
# the offset term of ancombc3_add_offset_covariance().
O <- matrix(stats::rnorm(n * p, sd = 0.5), n, p, dimnames = list(NULL, taxa))
units <- rep(seq_len(n / 3L), each = 3L)
# Design with level B as the reference: columns grpA and grpC.
x2 <- cbind(cont = cont, grpA = as.numeric(grp == "A"), grpC = as.numeric(grp == "C"))

# Second-stage fits on the design with reference A (fit1) and with reference
# B (fit2), with the same offsets.
mg_fits <- ac3_cached(function() {
  el <- ancombc3_group_eligibility(grp, y, x, r, cols, at_risk_prob, absence_weight)
  mg <- list(name = "grp", levels = levels(grp), columns = cols, reference_column = "grpA", eligibility = el)
  fitting_set <- unique(unlist(el$da_by_column, use.names = FALSE))
  fit_args <- list(y = y, sigma_vae = sigma_vae, abundance_offset = matrix(0, n, p),
                   structural_offset = matrix(0, n, p), sampling_zero_prob = sampling_zero_prob,
                   at_risk_prob = at_risk_prob, fitting_set = fitting_set, n_starts = 4L, n_cores = 2L)
  fit1 <- do.call(fit_ancombc3_da, c(fit_args, list(x = x, r = r, test_terms = cols,
                                                    eligible_taxa_by_term = el$da_by_column)))
  by2 <- list(grpA = el$da_by_column$grpB, grpC = taxa[el$da_pair[, "C_vs_B"]])
  fit2 <- do.call(fit_ancombc3_da, c(fit_args, list(x = x2, r = x2, test_terms = c("grpA", "grpC"),
                                                    eligible_taxa_by_term = by2)))
  list(el = el, mg = mg, by2 = by2, fit1 = fit1, fit2 = fit2)
})

# Simulated data set with three groups for ancombc3(): 90 samples, 80 taxa,
# abundance effects on taxa 1-20 and prevalence effects on taxa 21-28.
set.seed(20261111L)
n5 <- 90L
p5 <- 80L
grp5 <- factor(rep(c("A", "B", "C"), each = n5 / 3L))
age <- as.numeric(scale(stats::rnorm(n5)))
gB <- as.numeric(grp5 == "B")
gC <- as.numeric(grp5 == "C")
alpha5 <- stats::rnorm(p5, 3, 1)
bB <- bC <- numeric(p5)
bB[1:8] <- 1
bC[1:8] <- 2
bB[9:14] <- -1.5
bC[9:14] <- -1.5
bC[15:20] <- 1.5
bias5 <- stats::rnorm(n5, sd = 0.5)
log_mu5 <- matrix(alpha5, n5, p5, byrow = TRUE) + outer(gB, bB) + outer(gC, bC) + bias5
lambda5 <- exp(log_mu5 + matrix(stats::rnorm(n5 * p5, sd = 0.5), n5, p5))
gamma_B <- gamma_C <- numeric(p5)
gamma_B[21:28] <- -1.5
gamma_C[21:28] <- -3
a_dp <- rep(2, p5)
a_dp[21:40] <- 1
a_dp[60:80] <- 0.3
presence <- stats::plogis(matrix(a_dp, n5, p5, byrow = TRUE) + outer(gB, gamma_B) + outer(gC, gamma_C))
counts <- matrix(stats::rpois(n5 * p5, lambda5), n5, p5) * matrix(stats::rbinom(n5 * p5, 1L, presence), n5, p5)
colnames(counts) <- paste0("taxon", seq_len(p5))
rownames(counts) <- paste0("S", seq_len(n5))
meta <- data.frame(grp = as.character(grp5), age = age, row.names = rownames(counts), stringsAsFactors = FALSE)

test_that("level support and eligibility of comparisons, with and without cluster", {
  # Twelve samples in three groups of four and three taxa. The at-risk
  # probability is 0.5 in every cell and the structural-absence weight of a zero
  # count is 0.8. Thresholds: two positive counts and an at-risk mass of 1 per
  # level (DA); a structural-absence mass of 0.5 and an at-risk mass of 2 per
  # level (DP).
  g <- factor(rep(c("A", "B", "C"), each = 4L))
  y_s <- cbind(
    t1 = rep(5, 12L),
    t2 = c(3, 3, 3, 3, 2, 0, 0, 0, 4, 4, 0, 0),
    t3 = c(0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2)
  )
  rownames(y_s) <- paste0("s", 1:12)
  at_risk_s <- matrix(0.5, 12L, 3L, dimnames = dimnames(y_s))
  absence_s <- ifelse(y_s == 0, 0.8, 0)
  sup <- ancombc3_group_level_support(g, y_s, at_risk_s, absence_s, min_positive_per_group = 2,
                                      min_at_risk_design_support = 1, dp_min_side_absence = 0.5,
                                      dp_min_side_at_risk = 2)
  expect_true(identical(unname(sup$n_positive), rbind(c(4, 4, 4), c(4, 1, 2), c(3, 3, 3))),
              label = "numbers of positive counts by level")
  expect_true(identical(unname(sup$da), rbind(c(TRUE, TRUE, TRUE), c(TRUE, FALSE, TRUE), c(TRUE, TRUE, TRUE))),
              label = "DA level support")
  expect_true(isTRUE(all.equal(unname(sup$absence_mass), rbind(c(0, 0, 0), c(0, 2.4, 1.6), c(0.8, 0.8, 0.8)))) &&
                isTRUE(all.equal(unname(sup$dp_at_risk_mass), 4 - unname(sup$absence_mass))),
              label = "DP structural-absence and at-risk masses by level")
  expect_true(identical(unname(sup$dp), rbind(c(FALSE, FALSE, FALSE), c(FALSE, FALSE, TRUE), c(TRUE, TRUE, TRUE))),
              label = "DP level support")
  # Two samples per subject: the two positive counts of t2 in level C belong to
  # one subject, so that the level has one subject with a positive count. The
  # group is constant within subject, so the masses of a level are sums over its
  # subjects of the largest weight of the subject.
  subject <- rep(1:6, each = 2L)
  sup_c <- ancombc3_group_level_support(g, y_s, at_risk_s, absence_s, min_positive_per_group = 2,
                                        min_at_risk_design_support = 1, dp_min_side_absence = 0.5,
                                        dp_min_side_at_risk = 2, cluster = subject)
  expect_true(identical(unname(sup_c$n_positive), rbind(c(2, 2, 2), c(2, 1, 1), c(2, 2, 2))) &&
                identical(unname(sup_c$da[2L, ]), c(TRUE, FALSE, FALSE)),
              label = "cluster: positive counts are counted by subject")
  expect_true(isTRUE(all.equal(unname(sup_c$absence_mass), rbind(c(0, 0, 0), c(0, 1.6, 0.8), c(0.8, 0.8, 0.8)))) &&
                isTRUE(all.equal(unname(sup_c$dp_at_risk_mass), rbind(c(2, 2, 2), c(2, 1.2, 1.2), c(2, 2, 2)))) &&
                identical(unname(sup_c$dp), rbind(c(FALSE, FALSE, FALSE), c(FALSE, FALSE, FALSE), c(TRUE, TRUE, TRUE))),
              label = "cluster: level masses are sums of subject maxima")
  x_s <- cbind(gB = as.numeric(g == "B"), gC = as.numeric(g == "C"))
  el_s <- ancombc3_group_eligibility(g, y_s, x_s, x_s, c("gB", "gC"), at_risk_s, absence_s,
                                     min_positive_per_group = 2, min_at_risk_eff_n = 5,
                                     min_observed_prevalence = 0.1, min_at_risk_design_support = 1,
                                     dp_min_total_absence = 1, dp_min_total_at_risk = 5,
                                     dp_min_side_absence = 0.5, dp_min_side_at_risk = 2)
  expect_true(identical(el_s$pairs$comparison, c("B_vs_A", "C_vs_A", "C_vs_B")) &&
                identical(el_s$pairs$ref, c(1L, 1L, 2L)) && identical(el_s$pairs$num, c(2L, 3L, 3L)),
              label = "comparisons: reference comparisons first, then C_vs_B")
  expect_true(identical(unname(el_s$da_pair), rbind(c(TRUE, TRUE, TRUE), c(FALSE, TRUE, FALSE), c(TRUE, TRUE, TRUE))) &&
                identical(unname(el_s$da_global), c(TRUE, FALSE, TRUE)),
              label = "DA eligibility of comparisons and of the global test")
  expect_true(identical(unname(el_s$dp_pair), rbind(c(FALSE, FALSE, FALSE), c(FALSE, FALSE, FALSE), c(TRUE, TRUE, TRUE))) &&
                identical(unname(el_s$dp_global), c(FALSE, FALSE, TRUE)) &&
                identical(unname(el_s$dp_common), c(FALSE, TRUE, TRUE)),
              label = "DP eligibility of comparisons and of the global test")
  expect_true(identical(el_s$da_by_column, list(gB = c("t1", "t3"), gC = c("t1", "t2", "t3"))),
              label = "eligible taxa of the group columns (level versus reference)")
  el_strict <- ancombc3_group_eligibility(g, y_s, x_s, x_s, c("gB", "gC"), at_risk_s, absence_s,
                                          min_positive_per_group = 2, min_at_risk_eff_n = 7,
                                          min_observed_prevalence = 0.1, min_at_risk_design_support = 1)
  expect_true(!any(el_strict$da_pair) && !any(el_strict$da_global),
              label = "at-risk effective sample size below its threshold: no DA comparison")
})

test_that("DA channel: contrasts, scoring, joint covariance, global statistic, and offset term", {
  fits <- mg_fits()
  el <- fits$el
  mg <- fits$mg
  by2 <- fits$by2
  fit1 <- fits$fit1
  fit2 <- fits$fit2
  con <- .ac3_mg_da_contrast(fit1, colnames(x), col_num = "grpC", col_ref = "grpB", term = "C_vs_B")
  raw2 <- fit2$raw[fit2$raw$term == "grpC", , drop = FALSE]
  raw2 <- raw2[match(taxa, raw2$taxon), , drop = FALSE]
  d_beta <- max(abs(con$raw$beta_raw - raw2$beta_raw), na.rm = TRUE)
  d_model <- rel_diff(con$raw$var_raw_model, raw2$var_raw_model)
  d_marginal <- rel_diff(con$raw$var_raw_marginal, raw2$var_raw_marginal)
  d_G <- max(abs(con$G - fit2$da_variance_inputs$G$grpC), na.rm = TRUE)
  expect_true(sum(is.finite(con$raw$beta_raw)) >= p - 2L && d_beta <= 1e-4,
              label = "contrast: estimates equal those of the releveled fit")
  expect_true(d_model <= 1e-4, label = "contrast: model variances equal those of the releveled fit")
  expect_true(d_marginal <= 1e-4 && d_G <= 1e-4,
              label = "contrast: bias sensitivities and marginal variances equal those of the releveled fit")
  expect_true(identical(setdiff(names(con$raw), names(fit1$raw)), character(0)) &&
                identical(con$raw$n_positive, fit1$raw$n_positive[fit1$raw$term == "grpC"]) &&
                all(con$raw$term == "C_vs_B"),
              label = "contrast: the other columns are copied from the taxon's raw row")
  raw1_b <- fit1$raw[fit1$raw$term == "grpB", , drop = FALSE]
  raw2_a <- fit2$raw[fit2$raw$term == "grpA", , drop = FALSE]
  expect_true(max(abs(raw1_b$beta_raw + raw2_a$beta_raw), na.rm = TRUE) <= 1e-4 &&
                rel_diff(raw1_b$var_raw_model, raw2_a$var_raw_model) <= 1e-4,
              label = "releveled fit: coefficient of A versus B is minus that of B versus A")

  scored1 <- .ac3_mg_score(fit1$raw, fit1$da_variance_inputs, el$da_by_column, cfg, p_adj_method = "BH")
  scored2 <- .ac3_mg_score(fit2$raw, fit2$da_variance_inputs, by2, cfg, p_adj_method = "BH")
  da <- .ac3_mg_da_channel(fit1, scored1, colnames(x), mg, cfg, p_adj_method = "BH", pairs = TRUE)
  s2 <- scored2[scored2$term == "grpC", , drop = FALSE]
  s2 <- s2[match(taxa, s2$taxon), , drop = FALSE]
  d_stat <- max(abs(da$stat[, "C_vs_B"] - s2$statistic), na.rm = TRUE)
  d_est <- max(abs(da$est[, "C_vs_B"] - s2$beta_corrected), na.rm = TRUE)
  expect_true(d_est <= 1e-3 && d_stat <= 1e-3 && rel_diff(da$se[, "C_vs_B"], s2$se_scaled) <= 1e-3,
              label = "scored comparison C_vs_B equals the scoring of the releveled fit")
  s1 <- lapply(cols, function(tn) {
    s <- scored1[scored1$term == tn, , drop = FALSE]
    s[match(taxa, s$taxon), , drop = FALSE]
  })
  expect_true(identical(unname(da$est[, 1L]), s1[[1L]]$beta_corrected) &&
                identical(unname(da$p[, 2L]), s1[[2L]]$p_value),
              label = "comparisons with the reference are the scored group columns")

  jc <- .ac3_mg_da_joint_covariance(fit1, scored1, cols, colnames(x))
  diag_err <- max(vapply(1:2, function(k) {
    rel_diff(jc$Sigma[, k, k], s1[[k]]$var_primary)
  }, numeric(1)))
  expect_true(diag_err <= 1e-10, label = "joint covariance: diagonal equals var_primary")
  # Direct computation from the covariance of all uncorrected group
  # coefficients: Cov(beta_jk, beta_ml) = 1{j = m} V_j[k, l] + g_jk' V_d g_ml,
  # and the corrected coefficient of taxon j and column k is
  # beta_jk - sum_m nu_m^k beta_mk.
  idx <- 1L + match(cols, colnames(x))
  V_d <- fit1$da_variance_inputs$V
  G_all <- rbind(fit1$da_variance_inputs$G[[1L]], fit1$da_variance_inputs$G[[2L]])
  ok_taxa <- which(vapply(fit1$fits, function(f) isTRUE(f$ok) && !is.null(f$vcov), logical(1)) &
                     is.finite(rowSums(G_all[seq_len(p), , drop = FALSE])) &
                     is.finite(rowSums(G_all[p + seq_len(p), , drop = FALSE])))
  G_all[!is.finite(G_all)] <- 0
  full_cov <- G_all %*% V_d %*% t(G_all)
  for (j in ok_taxa) {
    pos <- c(j, p + j)
    full_cov[pos, pos] <- full_cov[pos, pos] + fit1$fits[[j]]$vcov[idx, idx]
  }
  nu <- cbind(s1[[1L]]$correction_weight, s1[[2L]]$correction_weight)
  nu[!is.finite(nu)] <- 0
  direct_err <- max(vapply(ok_taxa, function(j) {
    a <- matrix(0, 2L * p, 2L)
    a[seq_len(p), 1L] <- -nu[, 1L]
    a[p + seq_len(p), 2L] <- -nu[, 2L]
    a[j, 1L] <- a[j, 1L] + 1
    a[p + j, 2L] <- a[p + j, 2L] + 1
    direct <- crossprod(a, full_cov %*% a)
    max(abs(direct - jc$Sigma[j, , ]) / sqrt(outer(diag(direct), diag(direct))))
  }, numeric(1)))
  expect_true(length(ok_taxa) >= p - 2L && direct_err <= 1e-8, label = "joint covariance equals the direct computation")
  expect_true(any(nu[, 1L] != 0) && any(nu[, 2L] != 0), label = "correction weights are not all zero")

  expect_true(rel_diff(abs(da$b) / da$se_eff, abs(da$stat[, 1:2])) <= 1e-10,
              label = "standard errors reproduce the reported statistics")
  global_err <- max(vapply(which(da$global$has), function(j) {
    Sigma <- jc$Sigma[j, , ]
    S <- diag(da$se_eff[j, ])
    Sigma_tilde <- S %*% stats::cov2cor(Sigma) %*% S
    direct <- drop(da$b[j, ] %*% solve(Sigma_tilde, da$b[j, ]))
    abs(direct - da$global$stat[j]) / max(direct, 1e-8)
  }, numeric(1)))
  expect_true(sum(da$global$has) >= 20L && global_err <= 1e-8 &&
                isTRUE(all.equal(da$global$p, stats::pf(da$global$stat / 2, df1 = 2, df2 = da$df, lower.tail = FALSE))),
              label = "DA global statistic equals b' (S R S)^-1 b with an F(2, df) reference")
  expect_true(all(da$global$p[1:4] < 1e-3, na.rm = TRUE), label = "DA global test detects the simulated group effects")

  # Offset term (ancombc3_add_offset_covariance()) added to both fits. The term
  # depends on the fit only through G, so the second-stage fits are not
  # repeated.
  fit1o <- ancombc3_add_offset_covariance(fit1, O, cluster = units)
  fit2o <- ancombc3_add_offset_covariance(fit2, O, cluster = units)
  con_o <- .ac3_mg_da_contrast(fit1o, colnames(x), col_num = "grpC", col_ref = "grpB", term = "C_vs_B")
  raw2o <- fit2o$raw[fit2o$raw$term == "grpC", , drop = FALSE]
  raw2o <- raw2o[match(taxa, raw2o$taxon), , drop = FALSE]
  A_direct <- fit1o$da_variance_inputs$A$grpC - fit1o$da_variance_inputs$A$grpB
  off_share <- stats::median(con_o$raw$var_raw_offset / con_o$raw$var_raw_marginal, na.rm = TRUE)
  d_A <- max(abs(con_o$A - fit2o$da_variance_inputs$A$grpC)) / max(abs(fit2o$da_variance_inputs$A$grpC))
  d_offset <- rel_diff(con_o$raw$var_raw_offset, raw2o$var_raw_offset)
  d_marginal_o <- rel_diff(con_o$raw$var_raw_marginal, raw2o$var_raw_marginal)
  expect_true(off_share >= 0.05, label = "offset term: a material share of the marginal variance of the contrast")
  has_off <- is.finite(con_o$raw$var_raw_offset)
  expect_true(sum(has_off) >= p - 2L && identical(con_o$A[, has_off], A_direct[, has_off]) &&
                identical(unname(con_o$raw$var_raw_offset[has_off]), unname(colSums(A_direct^2)[has_off])) &&
                rel_diff(con_o$raw$var_raw_marginal, con$raw$var_raw_marginal + colSums(A_direct^2)) <= 1e-12,
              label = "offset term: contrast uses A_C - A_B and adds ||a_C - a_B||^2 to the marginal variance")
  expect_true(d_A <= 1e-4 && d_offset <= 1e-4 && d_marginal_o <= 1e-4,
              label = "offset term: contrast equals the offset term of the releveled fit")
  scored1o <- .ac3_mg_score(fit1o$raw, fit1o$da_variance_inputs, el$da_by_column, cfg, p_adj_method = "BH")
  scored2o <- .ac3_mg_score(fit2o$raw, fit2o$da_variance_inputs, by2, cfg, p_adj_method = "BH")
  da_o <- .ac3_mg_da_channel(fit1o, scored1o, colnames(x), mg, cfg, p_adj_method = "BH", pairs = TRUE)
  s2o <- scored2o[scored2o$term == "grpC", , drop = FALSE]
  s2o <- s2o[match(taxa, s2o$taxon), , drop = FALSE]
  d_est_o <- max(abs(da_o$est[, "C_vs_B"] - s2o$beta_corrected), na.rm = TRUE)
  d_stat_o <- max(abs(da_o$stat[, "C_vs_B"] - s2o$statistic), na.rm = TRUE)
  d_se_o <- rel_diff(da_o$se[, "C_vs_B"], s2o$se_scaled)
  expect_true(d_est_o <= 1e-3 && d_stat_o <= 1e-3 && d_se_o <= 1e-3 &&
                rel_diff(da_o$se[, "C_vs_B"], da$se[, "C_vs_B"]) > 1e-3,
              label = "offset term: scored comparison C_vs_B equals the scoring of the releveled fit")
  s1o <- lapply(cols, function(tn) {
    s <- scored1o[scored1o$term == tn, , drop = FALSE]
    s[match(taxa, s$taxon), , drop = FALSE]
  })
  jco <- .ac3_mg_da_joint_covariance(fit1o, scored1o, cols, colnames(x))
  diag_err_o <- max(vapply(1:2, function(k) {
    rel_diff(jco$Sigma[, k, k], s1o[[k]]$var_primary)
  }, numeric(1)))
  nu_o <- cbind(s1o[[1L]]$correction_weight, s1o[[2L]]$correction_weight)
  nu_o[!is.finite(nu_o)] <- 0
  A_all <- cbind(fit1o$da_variance_inputs$A[[1L]], fit1o$da_variance_inputs$A[[2L]])
  full_cov_o <- full_cov + crossprod(A_all)
  direct_err_o <- max(vapply(ok_taxa, function(j) {
    a <- matrix(0, 2L * p, 2L)
    a[seq_len(p), 1L] <- -nu_o[, 1L]
    a[p + seq_len(p), 2L] <- -nu_o[, 2L]
    a[j, 1L] <- a[j, 1L] + 1
    a[p + j, 2L] <- a[p + j, 2L] + 1
    direct <- crossprod(a, full_cov_o %*% a)
    max(abs(direct - jco$Sigma[j, , ]) / sqrt(outer(diag(direct), diag(direct))))
  }, numeric(1)))
  expect_true(diag_err_o <= 1e-10 && rel_diff(jco$Sigma[, 1L, 1L], jc$Sigma[, 1L, 1L]) > 1e-3,
              label = "offset term: joint covariance diagonal equals var_primary")
  expect_true(direct_err_o <= 1e-8, label = "offset term: joint covariance equals the direct computation")
})

test_that("DP channel: joint test and comparison of two non-reference levels", {
  fits <- mg_fits()
  mg <- fits$mg
  fit1 <- fits$fit1
  dp_args <- list(y = y, r = r, f0 = sampling_zero_prob, alpha = 0.1, ridge = cfg$ridge_dp,
                  da_fit = fit1, x = x)
  dp0 <- do.call(fit_ancombc3_dp_zero_pattern, c(dp_args, list(test_terms = cols)))
  dp1 <- do.call(fit_ancombc3_dp_zero_pattern, c(dp_args, list(test_terms = cols, joint_terms = cols)))
  joint <- attr(dp1, "joint")
  attr(dp1, "joint") <- NULL
  expect_true(identical(dp0, dp1) && is.null(attr(dp0, "joint")), label = "rows do not depend on joint_terms")
  expect_true(identical(names(joint), taxa), label = "joint results are named by taxon")
  dp_single <- do.call(fit_ancombc3_dp_zero_pattern, c(dp_args, list(test_terms = "grpB", joint_terms = "grpB")))
  single <- attr(dp_single, "joint")
  expect_true(rel_diff(vapply(single, function(v) v$statistic, numeric(1)), dp_single$statistic) <= 1e-10 &&
                rel_diff(vapply(single, function(v) v$p_value, numeric(1)), dp_single$p_value) <= 1e-10,
              label = "joint test of one column equals the one-degree-of-freedom test")
  # Direct computation of the global test and of the covariance.
  D <- cbind("(Intercept)" = 1, r)
  gi <- match(cols, colnames(D))
  direct_dp <- t(vapply(taxa, function(tx) {
    Z <- as.numeric(y[, tx] == 0)
    ss <- if (tx %in% fit1$fitting_set) .dp_second_stage_taxon(fit1$fits[[tx]], cbind(1, x)) else NULL
    f0_j <- pmin(pmax(if (is.null(ss)) sampling_zero_prob[, tx] else ss$f0, 1e-8), 1 - 1e-8)
    full <- stats::optim(rep(0, ncol(D)), .dp_zero_pattern_nll(Z, f0_j, D, ridge = cfg$ridge_dp),
                         method = "BFGS", control = list(maxit = 500))
    null <- stats::optim(full$par[-gi],
                         .dp_zero_pattern_nll(Z, f0_j, D, ridge = cfg$ridge_dp, fix_idx = gi, fix_vals = c(0, 0)),
                         method = "BFGS", control = list(maxit = 500))
    lr <- max(2 * (null$value - full$value), 0)
    rho <- 1
    cov_u <- matrix(NA_real_, 2L, 2L)
    if (!is.null(ss)) {
      ss$f0 <- f0_j
      coef0 <- numeric(ncol(D))
      coef0[-gi] <- null$par
      tv <- .dp_two_step_variance(Z, D, coef0, cfg$ridge_dp, ss)
      V_dp <- tv$V_dp[gi, gi]
      rho <- 2 / sum(diag(solve(V_dp) %*% (V_dp + tv$V_prop[gi, gi])))
      tu <- .dp_two_step_variance(Z, D, full$par, cfg$ridge_dp, ss)
      cov_u <- (tu$V_dp + tu$V_prop)[gi, gi]
    } else {
      cov_u <- solve(.dp_zero_pattern_information(Z, f0_j, D, full$par, cfg$ridge_dp)$H)[gi, gi]
    }
    c(statistic = rho * lr, rho = rho, c11 = cov_u[1L, 1L], c12 = cov_u[1L, 2L], c22 = cov_u[2L, 2L],
      second_stage = !is.null(ss))
  }, numeric(6)))
  got <- t(vapply(joint, function(v) {
    c(v$statistic, v$variance_ratio, v$covariance[1L, 1L], v$covariance[1L, 2L], v$covariance[2L, 2L])
  }, numeric(5)))
  expect_true(sum(direct_dp[, "second_stage"]) >= 20L && rel_diff(got[, 1L], direct_dp[, "statistic"]) <= 1e-8 &&
                rel_diff(got[, 2L], direct_dp[, "rho"]) <= 1e-8,
              label = "DP global statistic and variance ratio equal the direct computation")
  expect_true(all(got[, 2L] > 0 & got[, 2L] <= 1) && any(got[, 2L] < 1), label = "DP global variance ratio is in (0, 1]")
  expect_true(rel_diff(got[, 3:5], direct_dp[, c("c11", "c12", "c22")]) <= 1e-8 &&
                all(vapply(joint, function(v) identical(v$covariance_source, "unrestricted"), logical(1))),
              label = "covariance of the group coefficients at the unrestricted fit")
  expect_true(isTRUE(all.equal(vapply(joint, function(v) v$p_value, numeric(1)),
                               stats::pchisq(vapply(joint, function(v) v$statistic, numeric(1)), df = 2,
                                             lower.tail = FALSE))),
              label = "DP global p-value from the chi-square(2) reference")

  dp1 <- do.call(fit_ancombc3_dp_zero_pattern, c(dp_args, list(test_terms = cols, joint_terms = cols)))
  dp <- .ac3_mg_dp_channel(y, r, x, sampling_zero_prob, fit1, dp1, mg, cfg, alpha = 0.1, pairs = TRUE)
  # Prevalence design with level B as the reference.
  dp_b <- fit_ancombc3_dp_zero_pattern(y = y, r = x2, f0 = sampling_zero_prob, test_terms = c("grpA", "grpC"),
                                       alpha = 0.1, ridge = cfg$ridge_dp, da_fit = fit1, x = x)
  dp_bc <- dp_b[dp_b$test_term == "grpC", , drop = FALSE]
  dp_bc <- dp_bc[match(taxa, dp_bc$taxon), , drop = FALSE]
  expect_true(identical(unname(dp$est[, "C_vs_B"]), dp_bc$estimate) &&
                identical(unname(dp$stat[, "C_vs_B"]), dp_bc$statistic) &&
                identical(unname(dp$p[, "C_vs_B"]), dp_bc$p_value),
              label = "DP comparison C_vs_B is the test on the design with reference B")
  dp_ba <- dp_b[dp_b$test_term == "grpA", , drop = FALSE]
  dp_ba <- dp_ba[match(taxa, dp_ba$taxon), , drop = FALSE]
  moderate <- abs(dp$est[, "B_vs_A"]) < 3 & is.finite(dp_ba$estimate)
  d_rev <- max(abs(dp$est[moderate, "B_vs_A"] + dp_ba$estimate[moderate]))
  d_rev_stat <- max(abs(dp$stat[moderate, "B_vs_A"] - dp_ba$statistic[moderate]))
  d_tri <- max(abs(dp$est[moderate, "C_vs_B"] - (dp$est[moderate, "C_vs_A"] - dp$est[moderate, "B_vs_A"])))
  expect_true(sum(moderate) >= 15L && d_rev <= 0.05 && d_rev_stat <= 0.05,
              label = "A versus B on the releveled design mirrors B versus A")
  expect_true(d_tri <= 0.05, label = "DP estimate of C_vs_B equals the difference of the reference comparisons")
  z_dp <- sign(dp$est[, 1:2]) * sqrt(dp$stat[, 1:2])
  expect_true(rel_diff(dp$z, z_dp) <= 1e-12 &&
                rel_diff((abs(dp$b) / dp$se_eff)[abs(z_dp) > 0], abs(z_dp)[abs(z_dp) > 0]) <= 1e-10,
              label = "DP signed roots and effective standard errors")
  expect_true(all(vapply(seq_len(p), function(j) {
    isTRUE(all.equal(dp$corr[[j]], stats::cov2cor(joint[[j]]$covariance), check.attributes = FALSE))
  }, logical(1))), label = "DP correlation matrix is that of the returned covariance")
})

test_that("global statistic, Dunnett-type p-values, and mixed directional FDR", {
  R3 <- matrix(c(1, 0.5, 0.3, 0.5, 1, -0.2, 0.3, -0.2, 1), 3L)
  b3 <- c(0.8, -0.3, 1.1)
  se3 <- c(0.4, 0.2, 0.5)
  expect_true(isTRUE(all.equal(.ac3_mg_global_statistic(b3 / se3, R3),
                               drop(b3 %*% solve(diag(se3) %*% R3 %*% diag(se3), b3)))),
              label = "global statistic z' R^-1 z equals b' (S R S)^-1 b")
  near <- matrix(c(4, 5.9999999, 5.9999999, 9), 2L)
  R_near <- .ac3_mg_correlation(near)
  expect_true(isTRUE(all.equal(min(eigen(R_near, symmetric = TRUE)$values), 1e-6, tolerance = 1e-6)) &&
                all(diag(R_near) == 1) && is.null(.ac3_mg_correlation(matrix(c(1, 0, 0, 0), 2L))),
              label = "correlation matrix: smallest eigenvalue raised to 1e-6; NULL without positive variances")
  expect_true(isTRUE(all.equal(.ac3_mg_correlation(diag(se3) %*% R3 %*% diag(se3)), R3)),
              label = "correlation matrix of a well-conditioned covariance is unchanged")

  # Dunnett-type p-values. Two comparisons: one-dimensional numerical
  # integration of P(|Z_1| <= c, |Z_2| <= c) over Z_1.
  rho <- 0.6
  R2 <- matrix(c(1, rho, rho, 1), 2L)
  z2 <- c(2.3, -1.1)
  p2_direct <- vapply(abs(z2), function(bound) {
    inner <- function(u) {
      stats::dnorm(u) * (stats::pnorm((bound - rho * u) / sqrt(1 - rho^2)) -
                           stats::pnorm((-bound - rho * u) / sqrt(1 - rho^2)))
    }
    1 - stats::integrate(inner, -bound, bound, rel.tol = 1e-12, abs.tol = 1e-14)$value
  }, numeric(1))
  expect_true(max(abs(.ac3_mg_dunnett_p(z2, R2) - p2_direct)) <= 1e-9,
              label = "Dunnett p-values, two comparisons: numerical integration")
  expect_true(identical(.ac3_mg_dunnett_p(1.7), 2 * stats::pnorm(-1.7)), label = "Dunnett p-value, one comparison")
  # Multivariate t reference: one comparison; the limit of many degrees of
  # freedom; and, for two comparisons, the one-dimensional integral over the
  # chi-square mixing variable of the bivariate normal rectangle probability.
  expect_true(identical(.ac3_mg_dunnett_p(1.7, df = 12), 2 * stats::pt(-1.7, df = 12)),
              label = "Dunnett p-value with a t reference, one comparison")
  expect_true(max(abs(.ac3_mg_dunnett_p(z2, R2, df = 1e6) - .ac3_mg_dunnett_p(z2, R2))) <= 1e-5,
              label = "Dunnett p-value with a t reference: limit of many degrees of freedom")
  p2_t <- vapply(abs(z2), function(bound) {
    1 - stats::integrate(function(q) vapply(q, function(qq) {
      .ac3_mg_rectangle_exact(bound * sqrt(qq / 9), R2) * stats::dchisq(qq, df = 9)
    }, numeric(1)), 0, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
  expect_true(max(abs(.ac3_mg_dunnett_p(z2, R2, df = 9) - p2_t)) <= 1e-7,
              label = "Dunnett p-value with a t reference, two comparisons, equals the mixture integral")
  # Independence: p_k = 1 - (1 - 2 Phi(-|z_k|))^m.
  z4 <- c(0.4, -2.2, 2.9, 3.6)
  for (m in 2:5) {
    zz <- c(z4, 1.5)[seq_len(m)]
    exact <- 1 - (1 - 2 * stats::pnorm(-abs(zz)))^m
    expect_true(max(abs(.ac3_mg_dunnett_p(zz, diag(m)) - exact) / exact) <= (if (m <= 3L) 1e-9 else 2e-3),
                label = sprintf("Dunnett p-values, %d independent comparisons: closed form", m))
  }
  # Simulation with a fixed seed, three and four comparisons.
  R4 <- 0.5 + 0.5 * diag(4L)
  R4[1L, 4L] <- R4[4L, 1L] <- 0.1
  sim_check <- function(z, R, n_sim = 400000L) {
    set.seed(20261110L)
    draws <- matrix(stats::rnorm(n_sim * length(z)), n_sim) %*% chol(R)
    max_abs <- do.call(pmax, lapply(seq_along(z), function(k) abs(draws[, k])))
    p_sim <- vapply(abs(z), function(bound) mean(max_abs > bound), numeric(1))
    p_fun <- .ac3_mg_dunnett_p(z, R)
    list(p = p_fun, sim = p_sim, z_score = (p_fun - p_sim) / sqrt(pmax(p_sim * (1 - p_sim), 1e-12) / n_sim))
  }
  sim3 <- sim_check(c(0.9, -2.1, 2.8), R3)
  sim4 <- sim_check(z4, R4)
  expect_true(all(abs(sim3$z_score) <= 4), label = "Dunnett p-values, three comparisons: simulation with a fixed seed")
  expect_true(all(abs(sim4$z_score) <= 4), label = "Dunnett p-values, four comparisons: simulation with a fixed seed")
  set.seed(99L)
  state <- .Random.seed
  p_first <- .ac3_mg_dunnett_p(z4, R4)
  expect_true(identical(state, .Random.seed) && identical(p_first, .ac3_mg_dunnett_p(z4, R4)),
              label = "Dunnett p-values are deterministic and leave the random-number stream unchanged")
  expect_true(all(.ac3_mg_dunnett_p(z4, R4) >= 2 * stats::pnorm(-abs(z4))) &&
                all(.ac3_mg_dunnett_p(z4, R4) <= pmin(1, 4 * 2 * stats::pnorm(-abs(z4)))) &&
                identical(.ac3_mg_dunnett_p(z4, NULL), pmin(1, 4 * 2 * stats::pnorm(-abs(z4)))),
              label = "Dunnett p-values lie between the marginal and the Bonferroni p-values")

  # Pairwise mixed directional FDR: six taxa and three comparisons at
  # alpha = 0.1. Screening p-values: taxa 1, 3, and 4 use the global p-value;
  # taxon 2 has one eligible comparison (0.04), taxon 6 two (2 x 0.01 = 0.02);
  # taxon 5 has none. BH over the five screened taxa: sorted p 0.0001, 0.0008,
  # 0.02, 0.04, 0.25 give q 0.0005, 0.002, 0.0333, 0.05, 0.25, so R = 4 taxa are
  # selected and n_tax / R = 1.25. Holm within the selected taxa, times 1.25.
  p_pair <- rbind(c(0.001, 0.02, 0.5), c(0.04, NA, NA), c(0.2, 0.3, 0.9), c(0.0005, 0.004, 0.03),
                  c(NA, NA, NA), c(0.01, 0.6, NA))
  has_pair <- !is.na(p_pair)
  p_glob <- c(0.0008, NA, 0.25, 0.0001, NA, NA)
  has_glob <- c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE)
  md <- .ac3_mg_mdfdr_pairwise(p_pair, has_pair, p_glob, has_glob, alpha = 0.1, fwer_ctrl_method = "holm")
  q_expected <- rbind(c(0.003, 0.04, 0.5) * 1.25, c(0.04 * 1.25, NA, NA), c(1, 1, 1),
                      c(0.0015, 0.008, 0.03) * 1.25, c(NA, NA, NA), c(0.02 * 1.25, 0.6 * 1.25, NA))
  expect_true(isTRUE(all.equal(md$p_screen, c(0.0008, 0.04, 0.25, 0.0001, NA, 0.02))) &&
                isTRUE(all.equal(md$q_screen, c(0.002, 0.05, 0.25, 0.0005, NA, 0.1 / 3))),
              label = "pairwise mdFDR: screening p-values and their BH adjustment")
  expect_true(md$n_tax == 5L && md$n_selected == 4L && identical(md$selected, c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)),
              label = "pairwise mdFDR: screened and selected taxa")
  expect_true(isTRUE(all.equal(md$q, q_expected)), label = "pairwise mdFDR: Holm within taxon times n_tax / R")
  md_bonf <- .ac3_mg_mdfdr_pairwise(p_pair, has_pair, p_glob, has_glob, alpha = 0.1,
                                    fwer_ctrl_method = "bonferroni")
  expect_true(isTRUE(all.equal(md_bonf$q[4L, ], pmin(1, c(0.0015, 0.012, 0.09) * 1.25))),
              label = "pairwise mdFDR: fwer_ctrl_method is applied within taxon")

  # Dunnett-type mixed directional FDR with independent statistics, so that the
  # adjusted p-values have the closed form 1 - (1 - p)^m: five taxa and two
  # comparisons; taxon 3 has one supported comparison, taxon 5 none.
  z_d <- rbind(c(3.5, 0.5), c(1.2, -0.8), c(NA, 2.6), c(-2.9, 2.2), c(NA, NA))
  has_d <- !is.na(z_d)
  corr_d <- rep(list(diag(2L)), 5L)
  dn <- .ac3_mg_mdfdr_dunnett(z_d, has_d, corr_d, alpha = 0.1)
  marg <- 2 * stats::pnorm(-abs(z_d))
  p_adj_expected <- 1 - (1 - marg)^2
  p_adj_expected[3L, 2L] <- marg[3L, 2L]
  screen_expected <- apply(p_adj_expected, 1L, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
  q_screen_expected <- rep(NA_real_, 5L)
  q_screen_expected[1:4] <- stats::p.adjust(screen_expected[1:4], method = "BH")
  selected_expected <- !is.na(q_screen_expected) & q_screen_expected <= 0.1
  q_d_expected <- ifelse(has_d, 1, NA_real_)
  q_d_expected[selected_expected, ] <- pmin(1, p_adj_expected[selected_expected, ] * 4 / sum(selected_expected))
  expect_true(isTRUE(all.equal(dn$p_adjusted, p_adj_expected)) && isTRUE(all.equal(dn$p_screen, screen_expected)),
              label = "Dunnett mdFDR: adjusted and screening p-values")
  expect_true(dn$n_tax == 4L && identical(dn$selected, selected_expected) && sum(selected_expected) == 3L &&
                isTRUE(all.equal(dn$q, q_d_expected)),
              label = "Dunnett mdFDR: q = min(1, adjusted p times n_tax / R) for selected taxa, 1 otherwise")
  dn_null <- .ac3_mg_mdfdr_dunnett(z_d, has_d, vector("list", 5L), alpha = 0.1)
  expect_true(isTRUE(all.equal(dn_null$p_adjusted[1L, ], pmin(1, 2 * marg[1L, ]))),
              label = "Dunnett mdFDR: Bonferroni bound when the correlation is not available")
})

test_that("ancombc3() with a three-level group: tables, consistency, and calls without group or tested", {
  skip_if_no_torch()
  skip_if_not_extended()
  torch::torch_set_num_threads(1L)
  res <- suppressWarnings(ancombc3(
    data = t(counts), meta_data = meta, fix_formula = ~ age + grp, tested = "age", group = "grp",
    global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
    epochs = 20L, seed = 1L, n_cl = 2L, verbose = FALSE
  ))
  kept <- unique(res$res$taxon)
  pair_cols <- c("taxon", "group", "comparison", "level_num", "level_ref", "da_eligible", "da_estimate",
                 "da_se", "da_statistic", "da_p", "da_q", "da_discovery", "da_direction", "dp_eligible",
                 "dp_estimate", "dp_statistic", "dp_p", "dp_q", "dp_discovery", "dp_direction",
                 "omnibus_p", "omnibus_q", "omnibus_discovery")
  global_cols <- c("taxon", "group", "da_eligible", "da_statistic", "da_df", "da_df2", "da_p", "da_q", "da_discovery",
                   "dp_eligible", "dp_statistic", "dp_df", "dp_p", "dp_q", "dp_discovery", "omnibus_p",
                   "omnibus_q", "omnibus_discovery")
  trend_cols <- c("taxon", "group", "da_eligible", "da_statistic", "da_p", "da_q", "da_discovery", "da_pattern",
                  "da_fit_B", "da_fit_C", "dp_eligible", "dp_statistic", "dp_p", "dp_q", "dp_discovery",
                  "dp_pattern", "dp_fit_B", "dp_fit_C", "omnibus_p", "omnibus_q", "omnibus_discovery")
  primary_cols <- c("taxon", "covariate", "covariate_type", "da_eligible", "da_estimate", "da_se",
                    "da_statistic", "da_df", "da_p", "da_q", "da_discovery", "dp_eligible", "dp_estimate",
                    "dp_statistic", "dp_p", "dp_q", "dp_discovery", "omnibus_p", "omnibus_q",
                    "omnibus_discovery", "component_label")
  expect_true(all(c("res", "sample_bias", "samp_frac", "bias_correct_log_table", "sigma", "settings", "res_global",
                    "res_pair", "res_dunn", "res_trend") %in% names(res)),
              label = "returned elements with the four requests")
  expect_true(identical(names(res$res_pair), pair_cols) && identical(names(res$res_dunn), pair_cols),
              label = "res_pair and res_dunn have the specified columns")
  expect_true(identical(names(res$res_global), global_cols), label = "res_global has the specified columns")
  expect_true(identical(names(res$res_trend), trend_cols), label = "res_trend has the specified columns")
  blocks <- unique(res$res[c("covariate", "level")])
  expect_true(identical(names(res$res), append(primary_cols, "level", after = 2L)) &&
                nrow(res$res) == 3L * length(kept) &&
                identical(blocks$covariate, c("age", "grp", "grp")) && identical(blocks$level, c(NA, "B", "C")),
              label = "results: tested covariate and the two group columns, with the level column")
  expect_true(nrow(res$res_pair) == 3L * length(kept) && nrow(res$res_dunn) == 2L * length(kept) &&
                nrow(res$res_global) == length(kept) && nrow(res$res_trend) == length(kept) &&
                identical(unique(res$res_pair$comparison), c("B_vs_A", "C_vs_A", "C_vs_B")) &&
                identical(unique(res$res_dunn$comparison), c("B_vs_A", "C_vs_A")) &&
                identical(res$res_global$taxon, kept) && all(res$res_pair$group == "grp"),
              label = "table dimensions and comparisons")
  expect_true(identical(res$settings$group$levels, c("A", "B", "C")) &&
                identical(res$settings$group$columns, c("grpB", "grpC")) &&
                grepl("reference = A", res$settings$covariate_coding$grp),
              label = "settings record the group, its levels, and its columns")
  # Consistency with the primary table.
  primary_cols_shared <- c("da_eligible", "da_estimate", "da_se", "da_statistic", "dp_eligible", "dp_estimate",
                           "dp_statistic", "omnibus_p", "omnibus_q", "omnibus_discovery")
  same_as_primary <- function(tab, p_cols) {
    all(vapply(c("B", "C"), function(lv) {
      a <- res$res[res$res$level %in% lv, , drop = FALSE]
      b <- tab[tab$level_num == lv & tab$level_ref == "A", , drop = FALSE]
      identical(a$taxon, b$taxon) &&
        all(vapply(c(primary_cols_shared, p_cols), function(nm) identical(a[[nm]], b[[nm]]), logical(1)))
    }, logical(1)))
  }
  expect_true(same_as_primary(res$res_pair, c("da_p", "dp_p")),
              label = "res_pair: comparisons with the reference repeat the primary rows")
  expect_true(same_as_primary(res$res_dunn, character(0)),
              label = "res_dunn: estimates, statistics, and omnibus columns repeat the primary rows")
  pr <- res$res_pair
  dn <- res$res_dunn
  ref_rows <- pr[pr$level_ref == "A", , drop = FALSE]
  expect_true(all(dn$da_p >= ref_rows$da_p - 1e-12, na.rm = TRUE) && all(dn$dp_p >= ref_rows$dp_p - 1e-12, na.rm = TRUE) &&
                all(dn$da_p <= pmin(1, 2 * ref_rows$da_p) + 1e-12, na.rm = TRUE) &&
                identical(is.na(dn$da_p), is.na(ref_rows$da_p)),
              label = "res_dunn: adjusted p-values between the unadjusted p-values and twice them")
  q_rules <- function(tab) {
    all(vapply(c("da", "dp"), function(ch) {
      q <- tab[[paste0(ch, "_q")]]
      pv <- tab[[paste0(ch, "_p")]]
      disc <- tab[[paste0(ch, "_discovery")]]
      dir <- tab[[paste0(ch, "_direction")]]
      est <- tab[[paste0(ch, "_estimate")]]
      identical(is.na(q), is.na(pv)) && all(q >= pv - 1e-12, na.rm = TRUE) && all(q <= 1, na.rm = TRUE) &&
        identical(disc, !is.na(q) & q <= 0.1) && identical(!is.na(dir), disc) &&
        all(dir[disc] == ifelse(est[disc] > 0, "positive", "negative"))
    }, logical(1)))
  }
  expect_true(q_rules(pr), label = "res_pair: q, discoveries, and directions are consistent")
  expect_true(q_rules(dn), label = "res_dunn: q, discoveries, and directions are consistent")
  gl <- res$res_global
  tr <- res$res_trend
  expect_true(all(gl$da_df[gl$da_eligible & !is.na(gl$da_p)] == 2) && all(gl$dp_df[!is.na(gl$dp_p)] == 2) &&
                isTRUE(all.equal(gl$da_p, stats::pf(gl$da_statistic / 2, df1 = 2, df2 = gl$da_df2, lower.tail = FALSE))) &&
                isTRUE(all.equal(gl$da_q[!is.na(gl$da_p)], stats::p.adjust(gl$da_p[!is.na(gl$da_p)], "BH"))),
              label = "res_global: degrees of freedom, p-values, and BH q-values")
  all_levels <- vapply(kept, function(tx) all(pr$da_eligible[pr$taxon == tx]), logical(1))
  expect_true(identical(unname(all_levels), gl$da_eligible) && identical(gl$da_eligible, tr$da_eligible) &&
                identical(gl$dp_eligible, tr$dp_eligible),
              label = "global and trend tests are eligible when all comparisons are eligible")
  expect_true(all(tr$da_p > 0 & tr$da_p <= 1, na.rm = TRUE) && all(tr$da_pattern %in% c("increasing", "decreasing", NA)) &&
                all(abs(tr$da_p[!is.na(tr$da_p)] * 1001 - round(tr$da_p[!is.na(tr$da_p)] * 1001)) < 1e-8) &&
                identical(is.na(tr$da_fit_B), is.na(tr$da_p)),
              label = "res_trend: p-values on the grid of B + 1 = 1001 values, patterns, fits")
  inc <- which(tr$da_pattern == "increasing")
  dec <- which(tr$da_pattern == "decreasing")
  expect_true(all(tr$da_fit_B[inc] >= -1e-8 & tr$da_fit_C[inc] >= tr$da_fit_B[inc] - 1e-8) &&
                all(tr$da_fit_B[dec] <= 1e-8 & tr$da_fit_C[dec] <= tr$da_fit_B[dec] + 1e-8) &&
                isTRUE(all.equal(tr$da_statistic[inc], abs(tr$da_fit_C[inc]))),
              label = "res_trend: constrained estimates satisfy the reported pattern")
  true_da <- paste0("taxon", 1:20)
  true_inc <- paste0("taxon", 1:8)
  power_da <- mean(gl$da_discovery[gl$taxon %in% true_da])
  false_da <- sum(gl$da_discovery[!gl$taxon %in% true_da])
  expect_true(power_da >= 0.8 && false_da <= 5L, label = "global DA test detects the simulated abundance effects")
  expect_true(all(tr$da_discovery[tr$taxon %in% true_inc]) &&
                all(tr$da_pattern[tr$taxon %in% true_inc] == "increasing"),
              label = "trend test detects the simulated increasing effects")
  expect_true(all(pr$da_direction[pr$taxon %in% true_inc & pr$comparison == "C_vs_B"] %in% "positive"),
              label = "pairwise test: C above B for the simulated increasing effects")

  # Call without group: the primary analysis only.
  res0 <- suppressWarnings(ancombc3(
    data = t(counts), meta_data = meta, fix_formula = ~ age + grp, tested = "age",
    epochs = 20L, seed = 1L, n_cl = 2L, verbose = FALSE
  ))
  expect_true(!any(c("res_global", "res_pair", "res_dunn", "res_trend") %in% names(res0)) &&
                identical(names(res0$res), primary_cols) && is.null(res0$settings$group) &&
                identical(unique(res0$res$covariate), "age") &&
                identical(names(res0$settings$covariate_coding), "age"),
              label = "without group: elements and columns of the primary analysis only")
  expect_true(!identical(res0$res$da_p, res$res$da_p[res$res$covariate == "age"]),
              label = "with group, the group columns are tested columns of the fit (age results differ)")

  # Group without tested and without requests, on a subset of the data.
  sub_taxa <- 1:40
  res_g <- suppressWarnings(ancombc3(
    data = t(counts[, sub_taxa]), meta_data = meta, fix_formula = ~ grp, group = "grp",
    epochs = 5L, seed = 1L, n_cl = 2L, verbose = FALSE
  ))
  expect_true(!any(c("res_global", "res_pair", "res_dunn", "res_trend") %in% names(res_g)) &&
                identical(unique(res_g$res$covariate), "grp") &&
                identical(unique(res_g$res$level), c("B", "C")) &&
                identical(res_g$settings$tested, character(0)),
              label = "group without tested and without requests: group columns in the primary table")
})

test_that("argument validation of the multigroup analysis", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  run <- function(metadata = meta, formula = ~ age + grp, tested = "age", group = "grp", ...) {
    ancombc3(data = t(counts), meta_data = metadata, fix_formula = formula, tested = tested, group = group,
             epochs = 1L, verbose = FALSE, ...)
  }
  expect_true(ac3_errs(run(group = NULL, global = TRUE), "require the argument group"), label = "requests without group")
  expect_true(ac3_errs(run(group = "arm"), "not a column of meta_data"), label = "group outside metadata")
  expect_true(ac3_errs(run(tested = c("age", "grp")), "also in tested"), label = "group also in tested")
  expect_true(ac3_errs(run(formula = ~ age), "not a term of fix_formula"), label = "group outside the formula")
  expect_true(ac3_errs(run(formula = ~ age * grp, tested = NULL), "group 'grp' appears in an interaction"),
              label = "group in an interaction")
  m2 <- meta
  m2$grp[m2$grp == "C"] <- "B"
  expect_true(ac3_errs(run(metadata = m2), "requires at least three"), label = "group with two levels")
  m3 <- meta
  m3$grp[5:30] <- "B"
  expect_true(ac3_errs(run(metadata = m3), "each level requires at least"), label = "group level with too few samples")
  m4 <- meta
  m4$grp <- as.numeric(factor(m4$grp))
  expect_true(ac3_errs(run(metadata = m4), "factor or a character"), label = "numeric group")
  expect_true(ac3_errs(run(trend = TRUE, trend_control = list(contrast = list(up = diag(3)), node = list(2))),
                       "2 x 2 matrix"), label = "trend contrast of the wrong dimension")
  expect_true(ac3_errs(run(trend = TRUE, trend_control = list(contrast = list(up = diag(2)), node = list(3))),
                       "single position"), label = "trend node outside 1, ..., K - 1")
  expect_true(ac3_errs(run(trend = TRUE, trend_control = list(B = 0)), "positive integer"), label = "trend B not positive")
  expect_true(ac3_errs(run(trend = TRUE, trend_control = list(draws = 10)), "Unknown element"),
              label = "unknown element of trend_control")
  expect_true(ac3_errs(run(pairwise = TRUE, mdfdr_control = list(fwer_ctrl_method = "tukey")), "must be one of"),
              label = "unknown fwer_ctrl_method")
  expect_true(ac3_errs(run(global = NA), "TRUE or FALSE"), label = "request that is not TRUE or FALSE")
  m5 <- meta
  m5$grp3 <- factor(m5$grp)
  expect_true(ac3_errs(ancombc3(data = t(counts), meta_data = m5, fix_formula = ~ grp3, tested = "grp3", epochs = 1L,
                                verbose = FALSE),
                       "levels"), label = "three-level factor in tested is still an error")
})
