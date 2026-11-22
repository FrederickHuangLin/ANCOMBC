# Multigroup analysis with repeated measures (arguments group and cluster of
# ancombc3() used together).

# Repeated measures with m subjects and a fixed number of visits; the group
# factor (three levels) is constant within subject. Counts are
# Poisson-log-normal with a subject random intercept per taxon, sample-specific
# depth factors, and structural zeros with a subject random intercept in the
# presence logit. Taxa 1 to 8 have an abundance shift of 1.2 in level C, and
# taxa 9 to 14 a presence shift of -2 in level B.
simulate_group_repeated <- function(seed, m = 42L, visits = 3L, p = 60L) {
  set.seed(seed)
  n <- m * visits
  subject <- rep(seq_len(m), each = visits)
  grp <- factor(rep(rep(c("A", "B", "C"), length.out = m), each = visits), levels = c("A", "B", "C"))
  alpha <- stats::rnorm(p, 3, 1)
  a <- stats::rnorm(p, 1.5, 1)
  b <- matrix(stats::rnorm(m * p, sd = 0.7), m, p)
  u <- matrix(stats::rnorm(m * p, sd = 1.5), m, p)
  shift_c <- ifelse(seq_len(p) <= 8L, 1.2, 0)
  shift_b <- ifelse(seq_len(p) >= 9L & seq_len(p) <= 14L, -2, 0)
  depth <- stats::rnorm(n, sd = 0.5)
  log_lambda <- matrix(alpha, n, p, byrow = TRUE) + b[subject, ] + depth +
    outer(as.numeric(grp == "C"), shift_c) + matrix(stats::rnorm(n * p, sd = 0.5), n, p)
  logit_present <- matrix(a, n, p, byrow = TRUE) + u[subject, ] + outer(as.numeric(grp == "B"), shift_b)
  present <- matrix(stats::rbinom(n * p, 1L, stats::plogis(logit_present)), n, p)
  y <- matrix(stats::rpois(n * p, exp(log_lambda)), n, p) * present
  dimnames(y) <- list(paste0("S", seq_len(n)), paste0("T", seq_len(p)))
  list(counts = t(y),
       metadata = data.frame(subject = paste0("P", subject), grp = grp, sample = rownames(y),
                             row.names = rownames(y), stringsAsFactors = FALSE))
}
dat <- simulate_group_repeated(20261112L)
run <- function(cluster, settings = list()) {
  suppressWarnings(ancombc3(
    data = dat$counts, meta_data = dat$metadata, fix_formula = ~ grp, group = "grp", global = TRUE, pairwise = TRUE,
    dunnet = TRUE, trend = TRUE, trend_control = list(B = 200L), cluster = cluster, epochs = 20L,
    seed = 3L, verbose = FALSE, settings = settings
  ))
}
# ancombc3() with cluster = "subject" (c), without cluster (n), with cluster
# and da_offset_variance = "never" (cn), and with one sample per subject (s),
# each computed at its first request.
group_cluster_run <- local({
  store <- list()
  calls <- list(c = function() run("subject"), n = function() run(NULL),
                cn = function() run("subject", list(da_offset_variance = "never")),
                s = function() run("sample"))
  function(name) {
    if (is.null(store[[name]])) store[[name]] <<- calls[[name]]()
    store[[name]]
  }
})

test_that("joint DP quantities with one sample per subject are close to the model-based ones", {
  set.seed(20261111L)
  n <- 300L
  p <- 40L
  lev <- c("A", "B", "C")
  g <- factor(rep(lev, length.out = n), levels = lev)
  r <- cbind(grpB = as.numeric(g == "B"), grpC = as.numeric(g == "C"))
  gamma <- cbind(stats::rnorm(p, 1.5, 0.5), ifelse(seq_len(p) <= 8L, 1, 0), ifelse(seq_len(p) > 32L, -1, 0))
  f0 <- matrix(stats::runif(p, 0.05, 0.4), n, p, byrow = TRUE)
  present <- matrix(stats::rbinom(n * p, 1L, stats::plogis(cbind(1, r) %*% t(gamma))), n, p)
  sampling_zero <- matrix(stats::rbinom(n * p, 1L, f0), n, p)
  y <- present * (1 - sampling_zero) * matrix(1L + stats::rpois(n * p, 5), n, p)
  dimnames(y) <- list(paste0("S", seq_len(n)), paste0("T", seq_len(p)))
  dimnames(f0) <- dimnames(y)
  # A second-stage fit with an empty fitting set: every taxon uses the supplied
  # f0.
  no_fit <- list(fits = list(), fitting_set = character(0))
  dp_model <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0, test_terms = colnames(r),
                                           joint_terms = colnames(r), da_fit = no_fit, x = r)
  dp_robust <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0, test_terms = colnames(r),
                                            joint_terms = colnames(r), cluster = rownames(y),
                                            da_fit = no_fit, x = r)
  jm <- attr(dp_model, "joint")
  jr <- attr(dp_robust, "joint")
  expect_true(identical(dp_model$estimate, dp_robust$estimate) &&
                identical(dp_model$statistic_lr, dp_robust$statistic_lr),
              label = "estimates and likelihood-ratio statistics do not depend on cluster")
  expect_true(all(vapply(jr, function(v) isTRUE(v$robust), logical(1))), label = "joint covariances are cluster-robust")
  stat_ratio <- vapply(names(jm), function(tx) jr[[tx]]$statistic / jm[[tx]]$statistic, numeric(1))
  var_ratio <- unlist(lapply(names(jm), function(tx) diag(jr[[tx]]$covariance) / diag(jm[[tx]]$covariance)))
  cov_model <- vapply(names(jm), function(tx) jm[[tx]]$covariance[1L, 2L], numeric(1))
  cov_robust <- vapply(names(jm), function(tx) jr[[tx]]$covariance[1L, 2L], numeric(1))
  expect_true(stats::median(stat_ratio, na.rm = TRUE) >= 0.85 && stats::median(stat_ratio, na.rm = TRUE) <= 1.15,
              label = "median ratio of the joint DP statistics in [0.85, 1.15]")
  expect_true(stats::median(var_ratio) >= 0.85 && stats::median(var_ratio) <= 1.15,
              label = "median ratio of the joint DP variances in [0.85, 1.15]")
  expect_true(stats::cor(cov_model, cov_robust) > 0.8, label = "covariances of the two group coefficients agree (cor > 0.8)")
  one_col <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = f0, test_terms = "grpB", joint_terms = "grpB",
                                          cluster = rownames(y), da_fit = no_fit, x = r)
  j1 <- attr(one_col, "joint")
  expect_true(isTRUE(all.equal(vapply(j1, function(v) v$statistic, numeric(1)),
                               stats::setNames(one_col$statistic, one_col$taxon), tolerance = 1e-6)),
              label = "joint test of one column equals the 1-df test with cluster")
})

test_that("ancombc3() with group and cluster: tables, settings, diagnostics, and the offset term", {
  skip_if_no_torch()
  skip_if_not_extended()
  torch::torch_set_num_threads(1L)
  res_c <- group_cluster_run("c")
  res_n <- group_cluster_run("n")
  res_cn <- group_cluster_run("cn")
  tables <- c("res_global", "res_pair", "res_dunn", "res_trend")
  expect_true(all(tables %in% names(res_c)) && all(vapply(res_c[tables], nrow, integer(1)) > 0L),
              label = "the four multigroup tables are produced with cluster")
  finite_rows <- function(tab, channel, cols) {
    el <- tab[[paste0(channel, "_eligible")]]
    vals <- as.matrix(tab[el, paste0(channel, "_", cols), drop = FALSE])
    c(n = sum(el), finite = sum(stats::complete.cases(vals) & apply(is.finite(vals), 1L, all)))
  }
  for (tb in tables) {
    for (ch in c("da", "dp")) {
      cols <- if (tb == "res_global") c("statistic", "p", "q") else if (tb == "res_trend") c("statistic", "p", "q")
              else c("estimate", "statistic", "p", "q")
      fr <- finite_rows(res_c[[tb]], ch, cols)
      expect_true(fr[["n"]] > 0L && fr[["finite"]] == fr[["n"]],
                  label = sprintf("%s, %s: %d eligible rows, all with finite values", tb, ch, fr[["n"]]))
    }
  }
  prim <- res_c$res
  expect_true(all(c("grp") %in% prim$covariate) && all(prim$level[prim$covariate == "grp"] %in% c("B", "C")),
              label = "primary table lists the group columns")
  expect_true(identical(res_c$settings$n_subjects, 42L) && identical(res_c$settings$cluster, "subject"),
              label = "settings record the cluster and 42 subjects")
  expect_true(identical(res_c$diagnostics$cluster$n_subjects, 42L) &&
                identical(as.integer(res_c$diagnostics$cluster$samples_per_subject[["3"]]), 42L) &&
                identical(res_c$diagnostics$by_covariate$grpB$da_eligibility$criterion[1L], "multigroup"),
              label = "diagnostics record the subjects and the multigroup eligibility")
  dp_g_c <- res_c$res_global$dp_statistic
  dp_g_n <- res_n$res_global$dp_statistic
  da_g_c <- res_c$res_global$da_statistic
  da_g_n <- res_n$res_global$da_statistic
  expect_true(!isTRUE(all.equal(dp_g_c, dp_g_n)) && !isTRUE(all.equal(da_g_c, da_g_n)),
              label = "DA and DP global statistics differ from the call without cluster")
  expect_true(!isTRUE(all.equal(res_c$res_pair$dp_p, res_n$res_pair$dp_p)) &&
                !isTRUE(all.equal(res_c$res_dunn$da_p, res_n$res_dunn$da_p)),
              label = "pairwise DP and Dunnett-type DA p-values differ from the call without cluster")
  pair_cb <- res_c$res_pair$comparison == "C_vs_B"
  se_ratio <- res_c$res_pair$da_se[pair_cb] / res_cn$res_pair$da_se[pair_cb]
  se_ratio <- se_ratio[is.finite(se_ratio)]
  expect_true(identical(res_c$res_pair$dp_p, res_cn$res_pair$dp_p) &&
                identical(res_c$res_global$dp_statistic, res_cn$res_global$dp_statistic),
              label = "da_offset_variance = \"never\": DP results unchanged")
  expect_true(length(se_ratio) > 0L && any(abs(se_ratio - 1) > 1e-6) &&
                !isTRUE(all.equal(res_c$res_global$da_statistic, res_cn$res_global$da_statistic)),
              label = "default with cluster: offset term in C_vs_B and in the DA global statistics")
})

test_that("ancombc3() with one sample per subject: DP estimates and bias unchanged, DP statistics close", {
  skip_if_no_torch()
  skip_if_not_extended()
  torch::torch_set_num_threads(1L)
  res_n <- group_cluster_run("n")
  res_s <- group_cluster_run("s")
  shift_gap <- tapply(res_s$res$da_estimate - res_n$res$da_estimate, res_n$res$level,
                      function(v) diff(range(v, na.rm = TRUE)))
  expect_true(identical(res_s$res$dp_estimate, res_n$res$dp_estimate) &&
                identical(res_s$sample_bias, res_n$sample_bias) && all(shift_gap <= 1e-12),
              label = "DP estimates and sample bias equal; DA estimates differ by one common shift per level")
  dp_g_s <- res_s$res_global$dp_statistic
  dp_g_n <- res_n$res_global$dp_statistic
  both <- is.finite(dp_g_s) & is.finite(dp_g_n) & dp_g_n > 0
  ratio <- stats::median(dp_g_s[both] / dp_g_n[both])
  expect_true(ratio >= 0.8 && ratio <= 1.25, label = "DP global statistics close to the model-based ones")
})

test_that("a level with fewer subjects than min_positive_per_group is an error", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  few <- dat$metadata
  few$grp[few$grp == "C" & few$subject %in% unique(few$subject[few$grp == "C"])[-(1:4)]] <- "A"
  err <- tryCatch(ancombc3(data = dat$counts, meta_data = few, fix_formula = ~ grp, group = "grp", cluster = "subject",
                           epochs = 1L, verbose = FALSE), error = function(e) conditionMessage(e))
  expect_true(is.character(err) && grepl("4 subject", err, fixed = TRUE), label = "a level with 4 subjects is an error")
})
