# =============================================================================
# ANCOM-BC3 trend test over the levels of a group factor
# =============================================================================
#
# For a group factor with K levels, b = (b_2, ..., b_K) denotes the K - 1
# coefficients of the non-reference levels relative to the reference level, in
# factor order (the reference level has coefficient zero). A pattern is a pair
# (A, node): A is a matrix whose rows are linear inequality constraints, and
# the alternative hypothesis of the pattern is A b >= 0 with at least one
# strict inequality; node is the position among the non-reference levels used
# by the test statistic (the last position for a monotone pattern, the
# position of the peak for an umbrella pattern).
#
# Constrained estimate. Given the estimate b_hat and its working covariance
# Sigma, the constrained estimate of a pattern is
#   b_tilde = argmin_b (b - b_hat)' Sigma^-1 (b - b_hat)  subject to  A b >= 0,
# a strictly convex quadratic program with a unique solution
# (.ac3_trend_constrained()).
#
# Statistic. For one pattern, T = max(|b_tilde[node]|, |b_tilde[node] - b_tilde[K - 1]|).
# The test statistic is the maximum of T over the patterns; the pattern that
# attains the maximum is reported with its constrained estimate.
#
# Null distribution. B vectors are drawn from N(0, Sigma), the distribution of
# b_hat under b = 0, and each is passed through the same constrained estimates
# and the same maximum over patterns. With T_obs the observed statistic,
#   p = (1 + #{T_null >= T_obs}) / (B + 1).
#
# Random numbers. The draws are generated from the seed of the analysis with
# the Mersenne-Twister generator, and the state of the random-number stream of
# the caller is restored afterward. Every taxon of the analysis receives its
# own block of standard normal variates, in taxon order and, within a taxon,
# in the order of the channels, whether or not the taxon is tested; the normal
# variates of a taxon therefore depend only on the seed and on its position.
# The chi-square variates of the multivariate t reference are drawn in taxon
# order for the taxa with finite degrees of freedom; those of a taxon depend
# on the seed and on the degrees of freedom of the preceding taxa.

# Patterns of the trend test for a group factor with the given levels.
# contrast is a named list of (K - 1) x (K - 1) matrices A and node a list (or
# vector) of positions in 1, ..., K - 1, one per matrix. With contrast = NULL,
# the patterns are "increasing" (0 <= b_2 <= ... <= b_K) and "decreasing"
# (0 >= b_2 >= ... >= b_K), both with node K - 1. Returns list(contrast, node),
# both named by pattern.
ancombc3_trend_patterns <- function(levels, contrast = NULL, node = NULL) {
  m <- length(levels) - 1L
  if (m < 2L) stop("The trend test requires a group factor with at least three levels.")
  if (is.null(contrast)) {
    if (!is.null(node)) stop("trend_control$node is used only together with trend_control$contrast.")
    increasing <- diag(m)
    increasing[cbind(2:m, 1:(m - 1L))] <- -1
    return(list(contrast = list(increasing = increasing, decreasing = -increasing),
                node = list(increasing = m, decreasing = m)))
  }
  if (is.matrix(contrast)) contrast <- list(contrast)
  if (!is.list(contrast) || length(contrast) == 0L) {
    stop("trend_control$contrast must be a named list of matrices.")
  }
  if (is.null(names(contrast))) names(contrast) <- paste0("pattern_", seq_along(contrast))
  if (any(!nzchar(names(contrast))) || anyDuplicated(names(contrast))) {
    stop("trend_control$contrast must have distinct, non-empty names.")
  }
  if (is.null(node) || length(node) != length(contrast)) {
    stop("trend_control$node must give one position for each matrix of trend_control$contrast.")
  }
  node <- as.list(node)
  for (k in seq_along(contrast)) {
    A <- contrast[[k]]
    if (!is.matrix(A) || !is.numeric(A) || nrow(A) != m || ncol(A) != m || any(!is.finite(A))) {
      stop("trend_control$contrast[[\"", names(contrast)[k], "\"]] must be a finite numeric ", m, " x ", m,
           " matrix (rows: constraints; columns: the non-reference levels ",
           paste(levels[-1L], collapse = ", "), ").")
    }
    nd <- node[[k]]
    if (!is.numeric(nd) || length(nd) != 1L || !is.finite(nd) || nd != round(nd) || nd < 1 || nd > m) {
      stop("trend_control$node[[", k, "]] must be a single position in 1, ..., ", m, ".")
    }
    node[[k]] <- as.integer(nd)
    storage.mode(contrast[[k]]) <- "double"
    dimnames(contrast[[k]]) <- NULL
  }
  names(node) <- names(contrast)
  list(contrast = contrast, node = node)
}

# Constrained estimates of the rows of b (N x m) under A b >= 0 for the working
# covariance Sigma = S R S, with S = diag(s) and R a positive definite
# correlation matrix whose inverse is R_inv. The program is solved in the
# standardized coordinates u = b / s, in which the objective is
# (u - u_hat)' R^-1 (u - u_hat) and the constraints are C u >= 0 with
# C = A S, each row scaled to unit length. A row of b that satisfies the
# constraints is its own constrained estimate. For the other rows, the
# minimizer is unique and is characterized by the Karush-Kuhn-Tucker
# conditions: for the set T of constraints that hold with equality,
#   u = u_hat + R C_T' lambda,   lambda = -(C_T R C_T')^-1 C_T u_hat,
# with lambda >= 0 and C u >= 0. With enumerate = TRUE and at most eight
# constraints, the sets T are examined in turn for all rows at once, and a row
# is assigned the candidate of the first set that satisfies the conditions
# (tolerance 1e-9 in standardized units); a set whose matrix C_T R C_T' is
# numerically singular is skipped. The rows that remain, and every row that
# violates the constraints when enumerate = FALSE or there are more than eight
# constraints, are solved by quadprog::solve.QP(). A row for which the solver
# fails is returned as NA.
.ac3_trend_constrained <- function(b, s, R, R_inv, A, enumerate = TRUE) {
  m <- length(s)
  tol <- 1e-9
  u <- sweep(b, 2L, s, "/")
  C <- sweep(A, 2L, s, "*")
  row_norm <- sqrt(rowSums(C^2))
  C <- C[row_norm > 0, , drop = FALSE] / row_norm[row_norm > 0]
  n_con <- nrow(C)
  if (n_con == 0L) return(b)
  open <- which(rowSums(tcrossprod(u, C) < 0) > 0)
  if (isTRUE(enumerate) && n_con <= 8L) {
    CR <- C %*% R
    for (pattern in seq_len(2^n_con - 1L)) {
      if (length(open) == 0L) break
      active <- which(bitwAnd(pattern, 2^(0:(n_con - 1L))) > 0)
      C_T <- C[active, , drop = FALSE]
      M <- CR[active, , drop = FALSE] %*% t(C_T)
      if (rcond(M) < 1e-10) next
      u_open <- u[open, , drop = FALSE]
      lambda <- -tcrossprod(u_open, C_T) %*% solve(M)
      candidate <- u_open + lambda %*% CR[active, , drop = FALSE]
      solved <- rowSums(lambda < -tol) == 0 & rowSums(tcrossprod(candidate, C) < -tol) == 0
      solved[is.na(solved)] <- FALSE
      u[open[solved], ] <- candidate[solved, , drop = FALSE]
      open <- open[!solved]
    }
  }
  if (length(open) > 0L) {
    d <- u[open, , drop = FALSE] %*% R_inv
    Amat <- t(C)
    bvec <- numeric(n_con)
    for (i in seq_along(open)) {
      u[open[i], ] <- tryCatch(quadprog::solve.QP(R_inv, d[i, ], Amat, bvec)$solution,
                               error = function(e) rep(NA_real_, m))
    }
  }
  sweep(u, 2L, s, "*")
}

# Statistic of one pattern for the rows of the constrained estimates bt.
.ac3_trend_statistic <- function(bt, node) {
  pmax(abs(bt[, node]), abs(bt[, node] - bt[, ncol(bt)]))
}

# Trend test of one taxon: estimate b_hat, standard errors s, positive
# definite correlation matrix R (working covariance S R S), patterns
# (ancombc3_trend_patterns()), the B x m matrix E of standard normal
# variates, and, for a multivariate t reference with df degrees of freedom,
# the B multipliers scale = sqrt(df / Q) with Q chi-square with df degrees of
# freedom (NULL for the multivariate normal reference). The null draws are
# the rows of E U S with R = U'U, each multiplied by its element of scale. Returns statistic, p_value,
# pattern (the name of the maximizing pattern; the first one in case of ties),
# fit (its constrained estimate), and n_null, the number of null draws with a
# finite statistic. When the observed statistic is not finite, statistic and
# p_value are NA.
.ac3_trend_taxon <- function(b_hat, s, R, patterns, E, scale = NULL) {
  m <- length(b_hat)
  out <- list(statistic = NA_real_, p_value = NA_real_, pattern = NA_character_,
              fit = rep(NA_real_, m), n_null = 0L)
  R_inv <- tryCatch(solve(R), error = function(e) NULL)
  U <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(R_inv) || is.null(U)) return(out)
  R_inv <- 0.5 * (R_inv + t(R_inv))
  null_b <- sweep(E %*% U, 2L, s, "*")
  if (!is.null(scale)) null_b <- null_b * scale
  all_b <- rbind(b_hat, null_b)
  n_pattern <- length(patterns$contrast)
  stat <- matrix(NA_real_, nrow(all_b), n_pattern)
  fits <- vector("list", n_pattern)
  for (k in seq_len(n_pattern)) {
    bt <- .ac3_trend_constrained(all_b, s, R, R_inv, patterns$contrast[[k]])
    stat[, k] <- .ac3_trend_statistic(bt, patterns$node[[k]])
    fits[[k]] <- bt[1L, ]
  }
  observed <- stat[1L, ]
  if (any(!is.finite(observed))) return(out)
  best <- which.max(observed)
  null <- do.call(pmax, lapply(seq_len(n_pattern), function(k) stat[-1L, k]))
  ok <- is.finite(null)
  out$statistic <- observed[best]
  out$p_value <- (1 + sum(null[ok] >= observed[best])) / (sum(ok) + 1)
  out$pattern <- names(patterns$contrast)[best]
  out$fit <- fits[[best]]
  out$n_null <- sum(ok)
  out
}

# Trend tests of every taxon in one or several channels. channels is a named
# list; each channel holds b (p x m estimates), se (p x m standard errors),
# corr (a list of p correlation matrices, NULL where not available), eligible
# (p logical values), and optionally df (p degrees of freedom of a
# multivariate t reference; without df, or where df is not finite, the
# reference is multivariate normal). The chi-square variates of a channel
# with df are generated before the normal variates, from seed + 1, so that
# the normal variates are those of a call without df. A taxon is tested in a channel when it is
# eligible and its estimates, standard errors, and correlation matrix are
# available. Returns, per channel, statistic, p_value, pattern (length p), and
# fit (p x m), NA for the taxa that are not tested. The random numbers are
# generated as described in the header of this file.
ancombc3_trend_test <- function(channels, patterns, B = 1000L, seed = 1L) {
  B <- as.integer(B)
  if (length(B) != 1L || is.na(B) || B < 1L) stop("trend_control$B must be a positive integer.")
  p <- nrow(channels[[1L]]$b)
  m <- ncol(channels[[1L]]$b)
  out <- lapply(channels, function(ch) {
    list(statistic = rep(NA_real_, p), p_value = rep(NA_real_, p), pattern = rep(NA_character_, p),
         fit = matrix(NA_real_, p, m))
  })
  withr::local_seed(as.integer(seed) + 1L, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                     .rng_sample_kind = "Rejection")
  scales <- lapply(channels, function(ch) {
    if (is.null(ch$df)) return(NULL)
    matrix(vapply(seq_len(p), function(j) {
      df_j <- ch$df[j]
      if (is.na(df_j) || !is.finite(df_j) || df_j <= 0) return(rep(1, B))
      sqrt(df_j / stats::rchisq(B, df = df_j))
    }, numeric(B)), nrow = B)
  })
  withr::local_seed(as.integer(seed), .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                     .rng_sample_kind = "Rejection")
  for (j in seq_len(p)) {
    for (nm in names(channels)) {
      E <- matrix(stats::rnorm(B * m), B, m)
      ch <- channels[[nm]]
      R <- ch$corr[[j]]
      if (!isTRUE(ch$eligible[j]) || is.null(R) || any(!is.finite(ch$b[j, ])) ||
          any(!is.finite(ch$se[j, ])) || any(ch$se[j, ] <= 0)) {
        next
      }
      res <- .ac3_trend_taxon(ch$b[j, ], ch$se[j, ], R, patterns, E,
                              scale = if (is.null(scales[[nm]])) NULL else scales[[nm]][, j])
      out[[nm]]$statistic[j] <- res$statistic
      out[[nm]]$p_value[j] <- res$p_value
      out[[nm]]$pattern[j] <- res$pattern
      out[[nm]]$fit[j, ] <- res$fit
    }
  }
  out
}
