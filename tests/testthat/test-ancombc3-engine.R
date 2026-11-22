# Compiled engine of ancombc3(): accuracy of the compiled kernel and agreement
# of ancombc3() with engine = "compiled" (the default) and engine = "r".

test_that("compiled kernel agrees with the R code and with the composite Simpson reference", {
  ys <- c(0, 1, 2, 3, 5, 10, 20, 50, 100, 300, 1000, 3000, 1e4, 3e4, 1e5)
  grid <- expand.grid(nu = seq(-30, 30, by = 1), y = ys)
  sigmas <- exp(seq(log(0.02), log(4), length.out = 9))
  scaled <- function(a, b) max(abs(a - b) / pmax(1, abs(a)))
  relative <- function(a, b) max(abs(a - b) / abs(a))
  dev <- matrix(0, 0, 8, dimnames = list(NULL, c("logf", "E1", "E2", "libm", "dl_mu", "d2l_mumu",
                                                  "f0", "other")))
  set.seed(1)
  for (s in sigmas) {
    center <- grid$nu - s^2 / 2
    a <- .ac3_pln_obs(center, grid$y, s)
    b <- .ac3_pln_obs(center, grid$y, s, engine = "compiled")
    l <- .ac3_engine_pln_obs(center, grid$y, s, libm_exp = TRUE)
    logit_pi <- stats::qlogis(runif(nrow(grid), 0.05, 0.95))
    oa <- .ac3_pln_observation_derivatives(grid$nu, logit_pi, grid$y, grid$y > 0, s)
    ob <- .ac3_pln_observation_derivatives(grid$nu, logit_pi, grid$y, grid$y > 0, s, engine = "compiled")
    dev <- rbind(dev, c(
      logf = scaled(a$logf, b$logf), E1 = relative(a$E1, b$E1), E2 = relative(a$E2, b$E2),
      libm = max(relative(a$E1, l$E1), relative(a$E2, l$E2), scaled(a$logf, l$logf)),
      dl_mu = scaled(oa$dl_mu, ob$dl_mu), d2l_mumu = scaled(oa$d2l_mumu, ob$d2l_mumu),
      f0 = max(abs(oa$f0 - ob$f0)),
      other = max(abs(c(oa$ll - ob$ll, oa$dl_pi - ob$dl_pi, oa$d2l_pipi - ob$d2l_pipi,
                        oa$d2l_mupi - ob$d2l_mupi)) / pmax(1, abs(c(oa$ll, oa$dl_pi, oa$d2l_pipi, oa$d2l_mupi))))
    ))
  }
  worst <- apply(dev, 2L, max)
  # The tolerances are about 100 times the largest differences on this grid.
  # The differences of log f, E1, and E2 are of the order of the quadrature
  # error (the integration interval is determined to a tolerance and changes
  # with the last bits of the exponential); the first and second derivatives in
  # the log mean are y - E1 and -E1 + E2 - E1^2, whose cancellation for
  # y = 1e5 amplifies these differences to the size of the error of the R code
  # itself. With the exponential of the C library, the compiled code and the R
  # code differ only by floating-point rounding, which depends on the compiler
  # and the mathematical library of the platform.
  expect_true(worst[["logf"]] <= 1e-9, label = sprintf("log f: max |difference| / max(1, |log f|) %.1e", worst[["logf"]]))
  expect_true(worst[["E1"]] <= 1e-10, label = sprintf("E(Lambda | y): max relative difference %.1e", worst[["E1"]]))
  expect_true(worst[["E2"]] <= 1e-10, label = sprintf("E(Lambda^2 | y): max relative difference %.1e", worst[["E2"]]))
  expect_true(worst[["dl_mu"]] <= 1e-6,
              label = sprintf("dl / dlog mu: max |difference| / max(1, |value|) %.1e", worst[["dl_mu"]]))
  expect_true(worst[["d2l_mumu"]] <= 1e-3,
              label = sprintf("d2l / dlog mu2: max |difference| / max(1, |value|) %.1e", worst[["d2l_mumu"]]))
  expect_true(worst[["f0"]] <= 1e-13, label = sprintf("f0: max |difference| %.1e", worst[["f0"]]))
  expect_true(worst[["other"]] <= 1e-9,
              label = sprintf("ll and the presence derivatives: max scaled difference %.1e", worst[["other"]]))
  expect_true(worst[["libm"]] <= 1e-13,
              label = sprintf("C-library exponential vs R code: max relative difference %.1e", worst[["libm"]]))

  # Composite Simpson reference on 2 * 70001 + 1 points over the mode
  # +/- max(14, 40 s), with s the scale of the integrand in z.
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
  rgrid <- expand.grid(center = c(-8, -3, -1, 0, 2, 4, 7), sigma = c(0.03, 0.2, 0.7, 1.5, 2.5, 3.5),
                       y = c(0, 1, 2, 5, 20, 100, 2000))
  err <- t(vapply(seq_len(nrow(rgrid)), function(i) {
    g <- rgrid[i, ]
    q <- .ac3_pln_obs(g$center, g$y, g$sigma, engine = "compiled")
    r <- ref_pln(g$center, g$y, g$sigma)
    c(abs(q$logf - r$logf) / max(1, abs(r$logf)), abs(q$E1 - r$E1) / abs(r$E1), abs(q$E2 - r$E2) / abs(r$E2))
  }, numeric(3)))
  expect_true(max(err[, 1L]) <= 1e-11, label = sprintf("reference: log f error / max(1, |log f|) %.1e", max(err[, 1L])))
  expect_true(max(err[, 2L]) <= 1e-9, label = sprintf("reference: E(Lambda | y) relative error %.1e", max(err[, 2L])))
  expect_true(max(err[, 3L]) <= 1e-7, label = sprintf("reference: E(Lambda^2 | y) relative error %.1e", max(err[, 3L])))
})

# ancombc3() on the first 30 taxa of the toy data with 5 epochs: a call
# without engine, calls with engine = "compiled" with one and two worker
# processes, a call with engine = "r", and a second call without engine.
engine_runs <- ac3_cached(function() {
  toy <- ac3_toy_data()
  counts <- toy$counts[1:30, ]
  run <- function(...) {
    suppressWarnings(ancombc3(data = counts, meta_data = toy$meta, fix_formula = ~ group + bmi_z,
                              tested = c("group", "bmi_z"), epochs = 5L, seed = 1L, verbose = FALSE, ...))
  }
  list(default_before = run(), compiled_1 = run(engine = "compiled"),
       compiled_2 = run(engine = "compiled", n_cl = 2L), r_after = run(engine = "r"), default_after = run())
})
# Elements that record the call, the computation time, and the engine.
strip <- function(x) {
  x$settings[c("call", "elapsed_seconds", "stage_seconds", "n_cl", "engine")] <- NULL
  x
}

test_that("ancombc3() with both engines: default, worker processes, and agreement within tolerances", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  runs <- engine_runs()
  default_before <- runs$default_before
  compiled_1 <- runs$compiled_1
  compiled_2 <- runs$compiled_2
  r_after <- runs$r_after
  expect_true(identical(strip(default_before), strip(compiled_1)),
              label = "engine = \"compiled\" equals the call without engine")
  expect_true(identical(default_before$settings$engine, "compiled") && identical(r_after$settings$engine, "r"),
              label = "settings$engine records the engine")
  expect_true(identical(strip(compiled_1), strip(compiled_2)),
              label = "compiled engine: identical results with n_cl = 1 and 2")
  num <- c("da_estimate", "da_se", "da_statistic", "da_p", "da_q", "dp_estimate", "dp_statistic", "dp_p",
           "dp_q", "omnibus_p", "omnibus_q")
  diffs <- vapply(num, function(k) {
    a <- r_after$res[[k]]
    b <- compiled_1$res[[k]]
    if (!identical(is.na(a), is.na(b))) return(Inf)
    max(c(0, abs(a - b)), na.rm = TRUE)
  }, numeric(1))
  expect_true(all(diffs <= 1e-7),
              label = sprintf("res columns within 1e-7 (largest %.1e, %s)", max(diffs), names(which.max(diffs))))
  expect_true(max(abs(r_after$sample_bias - compiled_1$sample_bias)) <= 1e-8 &&
                max(abs(r_after$sigma - compiled_1$sigma), na.rm = TRUE) <= 1e-8,
              label = "sample_bias and sigma within 1e-8")
  same_cols <- c("taxon", "covariate", "da_eligible", "dp_eligible", "da_discovery", "dp_discovery",
                 "omnibus_discovery", "component_label")
  expect_true(identical(r_after$res[same_cols], compiled_1$res[same_cols]),
              label = "identical eligibility, discoveries, and component labels")
  estimator <- function(x) vapply(x$diagnostics$by_covariate, function(b) b$common_shift$estimator, character(1))
  expect_true(identical(estimator(r_after), estimator(compiled_1)), label = "identical common-shift estimators")
})

test_that("later calls without engine give the result of the first call", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  runs <- engine_runs()
  expect_true(identical(strip(runs$default_before), strip(runs$default_after)) &&
                identical(runs$default_after$settings$engine, "compiled"),
              label = "a call without engine after calls with both engines equals the call before them")
  center <- seq(-5, 5, length.out = 40)
  y <- rep(c(0, 3, 50, 1000), 10)
  expect_true(identical(.ac3_pln_obs(center, y, 0.8), local({
    r_only <- .ac3_pln_obs
    body(r_only) <- body(r_only)[-2L]
    r_only(center, y, 0.8)
  })), label = ".ac3_pln_obs() without engine evaluates the R code")
})
