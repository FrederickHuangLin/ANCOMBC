# Trend test of the multigroup analysis.

pat3 <- ancombc3_trend_patterns(c("A", "B", "C"))
pat4 <- ancombc3_trend_patterns(c("A", "B", "C", "D"))
umbrella <- rbind(c(1, 0, 0), c(-1, 1, 0), c(0, 1, -1))
pat_u <- ancombc3_trend_patterns(c("A", "B", "C", "D"), contrast = list(peak_C = umbrella), node = list(2))

# For every subset of the constraints taken as equalities, the minimizer of
# (b - b_hat)' Sigma^-1 (b - b_hat) is
#   b_hat - Sigma A_S' (A_S Sigma A_S')^-1 A_S b_hat;
# the constrained estimate is the candidate that satisfies A b >= 0 with the
# smallest value of the objective.
enumerate_active_sets <- function(b_hat, Sigma, A) {
  Sigma_inv <- solve(Sigma)
  best <- NULL
  best_value <- Inf
  for (pattern in 0:(2^nrow(A) - 1L)) {
    active <- which(bitwAnd(pattern, 2^(0:(nrow(A) - 1L))) > 0)
    candidate <- b_hat
    if (length(active) > 0L) {
      A_S <- A[active, , drop = FALSE]
      M <- A_S %*% Sigma %*% t(A_S)
      if (rcond(M) < 1e-12) next
      candidate <- drop(b_hat - Sigma %*% t(A_S) %*% solve(M, A_S %*% b_hat))
    }
    if (all(A %*% candidate >= -1e-9)) {
      value <- drop((candidate - b_hat) %*% Sigma_inv %*% (candidate - b_hat))
      if (value < best_value) {
        best_value <- value
        best <- candidate
      }
    }
  }
  best
}

# Estimates of 2000 taxa drawn from the null distribution with covariance
# Sigma0.
m <- 3L
Sigma0 <- diag(c(0.5, 1, 2)) %*% (0.5 + 0.5 * diag(m)) %*% diag(c(0.5, 1, 2))
n_taxa <- 2000L
b_null <- withr::with_seed(20261112L, matrix(stats::rnorm(n_taxa * m), n_taxa, m) %*% chol(Sigma0))
channel <- list(b = b_null, se = matrix(sqrt(diag(Sigma0)), n_taxa, m, byrow = TRUE),
                corr = rep(list(stats::cov2cor(Sigma0)), n_taxa), eligible = rep(TRUE, n_taxa))

test_that("patterns: default monotone patterns and validation of user-supplied patterns", {
  expect_true(identical(names(pat3$contrast), c("increasing", "decreasing")) &&
                identical(pat3$contrast$increasing, rbind(c(1, 0), c(-1, 1))) &&
                identical(pat3$contrast$decreasing, -rbind(c(1, 0), c(-1, 1))) &&
                identical(pat3$node, list(increasing = 2L, decreasing = 2L)),
              label = "default patterns for three levels")
  expect_true(identical(pat4$contrast$increasing, rbind(c(1, 0, 0), c(-1, 1, 0), c(0, -1, 1))) &&
                identical(pat4$node, list(increasing = 3L, decreasing = 3L)),
              label = "default patterns for four levels")
  expect_true(identical(pat_u$contrast$peak_C, umbrella) && identical(pat_u$node, list(peak_C = 2L)),
              label = "user-supplied umbrella pattern")
  expect_true(ac3_errs(ancombc3_trend_patterns(c("A", "B", "C"), contrast = list(u = umbrella), node = list(2)),
                       "2 x 2 matrix"),
              label = "contrast of the wrong dimension")
  expect_true(ac3_errs(ancombc3_trend_patterns(c("A", "B", "C", "D"), contrast = list(u = umbrella)), "one position"),
              label = "contrast without node")
  expect_true(ac3_errs(ancombc3_trend_patterns(c("A", "B", "C", "D"), contrast = list(u = umbrella), node = list(4)),
                       "single position"),
              label = "node outside 1, ..., K - 1")
  expect_true(ac3_errs(ancombc3_trend_patterns(c("A", "B", "C"), node = list(2)), "only together"),
              label = "node without contrast")
  expect_true(ac3_errs(ancombc3_trend_patterns(c("A", "B")), "at least three levels"), label = "two levels")
})

test_that("constrained estimate agrees with the enumeration of active sets", {
  set.seed(20261110L)
  qp_error <- function(m, constraint_sets, enumerate, n_rep = 150L) {
    worst <- 0
    for (i in seq_len(n_rep)) {
      L <- matrix(stats::rnorm(m * m), m)
      Sigma <- crossprod(L) + diag(0.05, m)
      s <- sqrt(diag(Sigma))
      R <- stats::cov2cor(Sigma)
      b_hat <- drop(stats::rnorm(m) %*% chol(Sigma)) * stats::runif(1, 0.5, 3)
      A <- constraint_sets[[1L + (i - 1L) %% length(constraint_sets)]]
      if (is.null(A)) A <- matrix(stats::rnorm(m * m), m)
      solved <- drop(.ac3_trend_constrained(matrix(b_hat, 1L), s, R, solve(R), A, enumerate = enumerate))
      expected <- enumerate_active_sets(b_hat, Sigma, A)
      worst <- max(worst, max(abs(solved - expected)) / max(s))
    }
    worst
  }
  for (enumerate in c(TRUE, FALSE)) {
    err3 <- qp_error(2L, c(pat3$contrast, list(NULL)), enumerate)
    err4 <- qp_error(3L, c(pat4$contrast, list(umbrella), list(NULL)), enumerate)
    method <- if (enumerate) "Karush-Kuhn-Tucker conditions" else "solve.QP"
    expect_true(err3 <= 1e-8, label = paste0(method, ", K = 3: monotone and arbitrary constraints"))
    expect_true(err4 <= 1e-8, label = paste0(method, ", K = 4: monotone, umbrella, and arbitrary constraints"))
  }
  # Many rows at once, and a constraint matrix with a repeated row.
  R_many <- 0.4 + 0.6 * diag(3L)
  s_many <- c(0.5, 1, 2)
  many <- sweep(matrix(stats::rnorm(3000L), 1000L, 3L) %*% chol(R_many), 2L, s_many, "*")
  repeated <- rbind(c(1, 0, 0), c(1, 0, 0), c(0, -1, 1))
  for (A in list(pat4$contrast$increasing, umbrella, repeated)) {
    by_kkt <- .ac3_trend_constrained(many, s_many, R_many, solve(R_many), A)
    by_qp <- .ac3_trend_constrained(many, s_many, R_many, solve(R_many), A, enumerate = FALSE)
    expect_true(all(is.finite(by_kkt)) && max(abs(by_kkt - by_qp)) <= 1e-8 && all(by_kkt %*% t(A) >= -1e-8),
                label = sprintf("1000 rows: both methods agree (%d distinct constraints)", nrow(unique(A))))
  }
  feasible <- rbind(c(0.5, 1.2, 2), c(-1, 0.3, 0.1))
  solved <- .ac3_trend_constrained(feasible, c(1, 2, 0.5), diag(3L), diag(3L), pat4$contrast$increasing)
  expect_true(identical(solved[1L, ], feasible[1L, ]) && all(pat4$contrast$increasing %*% solved[2L, ] >= -1e-10) &&
                isTRUE(all.equal(solved[2L, ], enumerate_active_sets(feasible[2L, ], diag(c(1, 4, 0.25)),
                                                                     pat4$contrast$increasing))),
              label = "an estimate that satisfies the constraints is unchanged")
})

test_that("statistic, reported pattern, and constrained estimates on constructed examples", {
  expect_true(identical(.ac3_trend_statistic(rbind(c(1, 3, 1), c(-2, -1, 4)), 2L), c(3, 5)) &&
                identical(.ac3_trend_statistic(rbind(c(1, 3, 1), c(-2, -1, 4)), 3L), c(1, 4)),
              label = "statistic max(|b[node]|, |b[node] - b[K - 1]|)")
  E <- matrix(stats::rnorm(999L * 2L), 999L, 2L)
  up <- .ac3_trend_taxon(c(1, 2), c(0.2, 0.2), diag(2L), pat3, E)
  expect_true(identical(up$pattern, "increasing") && identical(up$fit, c(1, 2)) && up$statistic == 2 &&
                up$p_value == 1 / 1000 && up$n_null == 999L,
              label = "increasing estimate: pattern, fit, statistic, smallest p-value 1 / (B + 1)")
  # b_hat = (1, -1) with the identity covariance: the increasing cone gives
  # (0, 0) and the decreasing cone gives (0, -1).
  mixed <- .ac3_trend_taxon(c(1, -1), c(1, 1), diag(2L), pat3, E)
  expect_true(identical(mixed$pattern, "decreasing") && isTRUE(all.equal(mixed$fit, c(0, -1))) &&
                isTRUE(all.equal(mixed$statistic, 1)),
              label = "non-monotone estimate: projection onto the decreasing cone")
  null_stat <- vapply(seq_len(nrow(E)), function(i) {
    max(vapply(pat3$contrast, function(A) {
      bt <- enumerate_active_sets(E[i, ], diag(2L), A)
      abs(bt[2L])
    }, numeric(1)))
  }, numeric(1))
  expect_true(isTRUE(all.equal(mixed$p_value, (1 + sum(null_stat >= mixed$statistic - 1e-12)) / 1000)),
              label = "p-value (1 + #{null >= observed}) / (B + 1) against a direct computation")
  zero <- .ac3_trend_taxon(c(0, 0), c(1, 1), diag(2L), pat3, E)
  expect_true(zero$statistic == 0 && zero$p_value == 1, label = "zero estimate: statistic 0 and p-value 1")
  peak <- .ac3_trend_taxon(c(1, 3, 1), rep(0.1, 3L), diag(3L), pat_u, matrix(stats::rnorm(600L), 200L, 3L))
  expect_true(identical(peak$pattern, "peak_C") && identical(peak$fit, c(1, 3, 1)) && peak$statistic == 3,
              label = "umbrella pattern: statistic at the peak")
})

test_that("p-values are uniform under the null hypothesis within Monte Carlo error", {
  null_run <- ancombc3_trend_test(list(da = channel), pat4, B = 199L, seed = 7L)
  p_null <- null_run$da$p_value
  levels_checked <- c(0.05, 0.10, 0.25, 0.50)
  rates <- vapply(levels_checked, function(a) mean(p_null <= a), numeric(1))
  z_rates <- (rates - levels_checked) / sqrt(levels_checked * (1 - levels_checked) / n_taxa)
  expect_true(all(is.finite(p_null)) && all(abs(z_rates) <= 4),
              label = "null p-values are uniform within Monte Carlo error")
  expect_true(abs(mean(p_null) - 0.5) <= 4 * sqrt(1 / 12 / n_taxa) + 1 / 400, label = "mean null p-value is one half")
  expect_true(abs(mean(null_run$da$pattern == "increasing") - 0.5) <= 4 * sqrt(0.25 / n_taxa),
              label = "the two monotone patterns are reported equally often under the null")
})

test_that("random numbers: reproducible from the seed, independent of other taxa, caller state unchanged", {
  withr::local_preserve_seed()
  small <- list(b = b_null[1:30, ], se = channel$se[1:30, ], corr = channel$corr[1:30], eligible = rep(TRUE, 30L))
  two <- list(da = small, dp = small)
  set.seed(123L)
  state <- .Random.seed
  run1 <- ancombc3_trend_test(two, pat4, B = 99L, seed = 3L)
  expect_true(identical(state, .Random.seed), label = "state of the random-number stream is unchanged")
  run2 <- ancombc3_trend_test(two, pat4, B = 99L, seed = 3L)
  run3 <- ancombc3_trend_test(two, pat4, B = 99L, seed = 4L)
  expect_true(identical(run1, run2), label = "same seed: identical results")
  expect_true(!identical(run1$da$p_value, run3$da$p_value), label = "different seed: different p-values")
  expect_true(!identical(run1$da$p_value, run1$dp$p_value) && identical(run1$da$statistic, run1$dp$statistic),
              label = "channels use separate draws")
  some <- two
  some$da$eligible[c(2, 5, 11)] <- FALSE
  some$dp$eligible[1:10] <- FALSE
  run4 <- ancombc3_trend_test(some, pat4, B = 99L, seed = 3L)
  keep_da <- some$da$eligible
  expect_true(identical(run4$da$p_value[keep_da], run1$da$p_value[keep_da]) && all(is.na(run4$da$p_value[!keep_da])) &&
                identical(run4$dp$p_value[11:30], run1$dp$p_value[11:30]) && all(is.na(run4$dp$statistic[1:10])),
              label = "the result of a taxon does not depend on which other taxa are tested")
  suppressWarnings(rm(".Random.seed", envir = globalenv()))
  # B = 1 is accepted, also with the multivariate t reference.
  with_df <- list(da = c(small, list(df = rep(12, 30L))))
  run_b1 <- tryCatch(ancombc3_trend_test(with_df, pat4, B = 1L, seed = 3L), error = function(e) NULL)
  expect_true(!is.null(run_b1) && all(run_b1$da$p_value %in% c(0.5, 1)), label = "B = 1 is accepted")
  suppressWarnings(rm(".Random.seed", envir = globalenv()))
  invisible(ancombc3_trend_test(two, pat4, B = 9L, seed = 3L))
  expect_true(!exists(".Random.seed", envir = globalenv(), inherits = FALSE),
              label = "no random-number state is created when the caller has none")
})

test_that("near-singular and singular covariance matrices", {
  E <- matrix(stats::rnorm(999L * 2L), 999L, 2L)
  near <- diag(c(1, 2, 3)) %*% matrix(c(1, 1 - 1e-12, 0.5, 1 - 1e-12, 1, 0.5, 0.5, 0.5, 1), 3L) %*% diag(c(1, 2, 3))
  R_near <- .ac3_mg_correlation(near)
  ns <- .ac3_trend_taxon(c(0.5, 1.5, 4), sqrt(diag(near)), R_near, pat4, matrix(stats::rnorm(600L), 200L, 3L))
  expect_true(is.finite(ns$statistic) && is.finite(ns$p_value) && ns$n_null == 200L &&
                all(pat4$contrast[[ns$pattern]] %*% ns$fit >= -1e-8),
              label = "correlation raised to a positive definite matrix: finite statistic and p-value")
  singular <- .ac3_trend_taxon(c(1, 2), c(1, 1), matrix(1, 2L, 2L), pat3, E)
  expect_true(is.na(singular$statistic) && is.na(singular$p_value), label = "singular correlation matrix: no test")
  missing_corr <- ancombc3_trend_test(list(da = list(b = b_null[1:3, ], se = channel$se[1:3, ],
                                                     corr = list(NULL, channel$corr[[1L]], channel$corr[[1L]]),
                                                     eligible = c(TRUE, TRUE, FALSE))), pat4, B = 19L, seed = 1L)
  expect_true(identical(is.na(missing_corr$da$p_value), c(TRUE, FALSE, TRUE)),
              label = "taxa without a correlation matrix or not eligible are not tested")
})
