# Poisson-log-normal quadrature as a function of the normal mean
# c = log mu - sigma^2 / 2 of the log intensity (center), on the grid
# c in {-8, ..., 7}, sigma in {0.03, ..., 3.5}, and y in {0, ..., 2000}.

# Composite Simpson rule on 2 * 70001 + 1 points over the mode +/- max(14, 40 s),
# where s = (-L''(z_hat))^(-1/2) is the scale of the integrand in z (step at
# most 2e-4 times max(1, 40 s / 14)).
ref_pln <- function(center, y, sigma) {
  md <- .ac3_pln_mode(center, y, sigma)
  s <- 1 / sqrt(md$curv)
  hw <- max(14, 40 * s)
  z <- seq(md$z - hw, md$z + hw, length.out = 2L * 70001L + 1L)
  lam <- exp(center + sigma * z)
  lg <- y * (center + sigma * z) - lam - lgamma(y + 1) - 0.5 * z^2 - 0.5 * log(2 * pi)
  m <- max(lg)
  w <- rep(c(2, 4), length.out = length(z))
  w[1L] <- w[length(z)] <- 1
  e <- w * exp(lg - m)
  dz <- z[2L] - z[1L]
  list(logf = m + log(sum(e) * dz / 3), E1 = sum(e * lam) / sum(e), E2 = sum(e * lam^2) / sum(e))
}

# Quadrature, reference, and central finite differences in c (step 1e-4) on
# the grid.
quadrature_grid <- ac3_cached(function() {
  grid <- expand.grid(center = c(-8, -3, -1, 0, 2, 4, 7), sigma = c(0.03, 0.2, 0.7, 1.5, 2.5, 3.5),
                      y = c(0, 1, 2, 5, 20, 100, 2000))
  h <- 1e-4
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    q <- .ac3_pln_obs(g$center, g$y, g$sigma)
    r <- ref_pln(g$center, g$y, g$sigma)
    lp <- .ac3_pln_obs(g$center + h, g$y, g$sigma, moments = FALSE)$logf
    lm <- .ac3_pln_obs(g$center - h, g$y, g$sigma, moments = FALSE)$logf
    d1 <- g$y - q$E1
    d2 <- -q$E1 + (q$E2 - q$E1^2)
    data.frame(g, logf = q$logf,
               logf_err = q$logf - r$logf,
               E1_relerr = (q$E1 - r$E1) / max(abs(r$E1), 1e-300),
               E2_relerr = (q$E2 - r$E2) / max(abs(r$E2), 1e-300),
               d1_err = (d1 - (lp - lm) / (2 * h)) / max(1, abs(d1)),
               d2_err = (d2 - (lp - 2 * q$logf + lm) / h^2) / max(1, abs(d2)))
  }))
})

# The log-density error is measured relative to max(1, |log f|): log f is a
# sum of terms of the size of |log f|, so its rounding error grows with
# |log f|. The finite-difference errors are relative to max(1, |derivative|);
# the second difference loses about 1e-16 |log f| / h^2 to rounding, which
# reaches 1e-4 at y = 2000, so its tolerance for y > 0 is wider than for y = 0.
within_tol <- function(value, tol) is.finite(value) && value <= tol

test_that("log f agrees with the composite Simpson reference", {
  res <- quadrature_grid()
  expect_true(within_tol(max(abs(res$logf_err) / pmax(1, abs(res$logf))), 1e-11),
              label = "log f, error / max(1, |log f|)")
})

test_that("posterior moments E(Lambda | y) and E(Lambda^2 | y) agree with the reference", {
  res <- quadrature_grid()
  expect_true(within_tol(max(abs(res$E1_relerr)), 1e-9), label = "E(Lambda | y), relative error")
  expect_true(within_tol(max(abs(res$E2_relerr)), 1e-7), label = "E(Lambda^2 | y), relative error")
})

test_that("derivative identities in c agree with central finite differences", {
  res <- quadrature_grid()
  expect_true(within_tol(max(abs(res$d1_err)), 1e-6), label = "d log f / dc vs finite difference")
  expect_true(within_tol(max(abs(res$d2_err[res$y > 0])), 1e-3),
              label = "d2 log f / dc2 vs finite difference, y > 0")
  expect_true(within_tol(max(abs(res$d2_err[res$y == 0])), 1e-5),
              label = "d2 log f / dc2 vs finite difference, y = 0")
})

test_that("cell derivatives at counts from 1e6 to 1e10 agree with central finite differences", {
  # For y in {1e6, ..., 1e10}, sigma in {0.05, 1, 3}, and log mu within
  # sigma / 2 of log y + sigma^2 / 2: the cell derivatives of
  # .ac3_pln_observation_derivatives() against central finite differences of
  # log f in c (step 1e-2), relative to max(1, |derivative|).
  large <- expand.grid(y = 10^(6:10), sigma = c(0.05, 1, 3), shift = seq(-0.5, 0.5, length.out = 5))
  h_large <- 1e-2
  large_err <- vapply(seq_len(nrow(large)), function(i) {
    g <- large[i, ]
    nu <- log(g$y) + g$sigma^2 / 2 + g$shift * g$sigma
    o <- .ac3_pln_observation_derivatives(nu, 5, g$y, TRUE, g$sigma)
    lf <- function(d) .ac3_pln_obs(nu - g$sigma^2 / 2 + d, g$y, g$sigma, moments = FALSE)$logf
    c(abs(o$dl_mu - (lf(h_large) - lf(-h_large)) / (2 * h_large)) / max(1, abs(o$dl_mu)),
      abs(o$d2l_mumu - (lf(h_large) - 2 * lf(0) + lf(-h_large)) / h_large^2) / max(1, abs(o$d2l_mumu)))
  }, numeric(2))
  expect_true(within_tol(max(large_err[1L, ]), 1e-7), label = "large counts: d log f / dc vs finite difference")
  expect_true(within_tol(max(large_err[2L, ]), 1e-6), label = "large counts: d2 log f / dc2 vs finite difference")
})
