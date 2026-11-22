# Common-shift correction and DA eligibility on synthetic coefficient vectors.

test_that("maximal overlap: bandwidth, center, null set, and estimate by direct enumeration", {
  # 30 coefficients near zero with unequal standard errors and 10 coefficients
  # at +/- 2.
  beta <- c(seq(-0.3, 0.3, length.out = 30L), rep(c(-2, 2), each = 5L))
  se <- 0.08 + 0.04 * (seq_along(beta) %% 3L)
  grid_z <- sort(unique(as.numeric(seq(0.5, 3.0, by = 0.1))))
  # Largest number of closed intervals with a common point, attained at a left
  # endpoint.
  coverage <- function(z, at) {
    vapply(at, function(c0) sum(beta - z * se <= c0 & beta + z * se >= c0), numeric(1))
  }
  max_cover <- vapply(grid_z, function(z) max(coverage(z, beta - z * se)), numeric(1))
  z_expected <- grid_z[which(max_cover >= 25)[1L]]
  left <- beta - z_expected * se
  right <- beta + z_expected * se
  n_left <- coverage(z_expected, left)
  a1 <- min(left[n_left == max(n_left)])
  b1 <- min(right[left <= a1 & right >= a1])
  center_expected <- (a1 + b1) / 2
  null_set_expected <- which(left <= center_expected & right >= center_expected)
  w_expected <- 1 / pmax(se[null_set_expected]^2, .Machine$double.eps)
  delta_expected <- sum(beta[null_set_expected] * w_expected) / sum(w_expected)

  mo <- .ac3_da_maximal_overlap(beta, se, min_count = 25L)
  expect_true(isTRUE(z_expected > 0.5) && identical(mo$bandwidth, z_expected) && identical(mo$status, "ok"),
              label = "bandwidth: smallest grid value with 25 overlapping intervals")
  expect_true(identical(mo$overlap_count, as.integer(max_cover[grid_z == z_expected])),
              label = "overlap count at the selected bandwidth")
  expect_true(isTRUE(all.equal(mo$center, center_expected, tolerance = 1e-14)),
              label = "center: midpoint of the leftmost interval of maximal overlap")
  expect_true(identical(mo$null_set_pos, null_set_expected) && mo$n_null_set == length(null_set_expected),
              label = "null set: taxa whose intervals contain the center")
  expect_true(isTRUE(all.equal(mo$delta_mo, delta_expected, tolerance = 1e-14)) &&
                isTRUE(all.equal(mo$var_delta_mo, 1 / sum(w_expected), tolerance = 1e-14)),
              label = "estimate: inverse-variance-weighted mean over the null set")
  expect_true(isTRUE(all.equal(sum(mo$nu), 1, tolerance = 1e-14)) &&
                isTRUE(all.equal(sum(mo$nu * beta[mo$null_set_pos]), mo$delta_mo, tolerance = 1e-14)),
              label = "correction weights of the null set sum to one and reproduce the estimate")
  # Positions refer to the input coefficients when some are not finite.
  beta_na <- c(NA, beta)
  se_na <- c(0.1, se)
  mo_na <- .ac3_da_maximal_overlap(beta_na, se_na, min_count = 25L)
  expect_true(identical(mo_na$null_set_pos, null_set_expected + 1L) && identical(mo_na$delta_mo, mo$delta_mo),
              label = "null-set positions refer to the input coefficients")
  mo_few <- .ac3_da_maximal_overlap(beta[1:10], se[1:10], min_count = 25L)
  expect_true(identical(mo_few$status, "bandwidth_criterion_not_met") && identical(mo_few$bandwidth, 3),
              label = "widest bandwidth when no bandwidth reaches 25 overlapping intervals")
  mo_min <- .ac3_da_maximal_overlap(beta[1:4], se[1:4])
  expect_true(identical(mo_min$status, "insufficient_taxa") && identical(mo_min$delta_mo, 0) &&
                length(mo_min$nu) == 0L,
              label = "fewer than five coefficients: zero estimate and no correction weights")
})

test_that("selection of the common-shift estimate: log-likelihood gain and separation criteria", {
  set.seed(20261015L)
  se_sel <- rep(0.1, 60L)
  beta_sel <- c(stats::rnorm(40L, 0, 0.1), stats::rnorm(20L, 1.5, sqrt(0.3^2 + 0.1^2)))
  var_sel <- se_sel^2
  rule <- function(sel, delta_mo, min_gain, min_separation) {
    separation <- abs(sel$delta_gm - delta_mo) / stats::median(se_sel)
    is.finite(sel$loglik_gain) && sel$loglik_gain > min_gain &&
      is.finite(separation) && separation > min_separation
  }
  delta_right <- .ac3_da_maximal_overlap(beta_sel, se_sel)$delta_mo
  sel_right <- .ac3_da_shift_selection(beta_sel, var_sel, delta_mo = delta_right, center = delta_right)
  sel_wrong <- .ac3_da_shift_selection(beta_sel, var_sel, delta_mo = 1.5, center = 1.5)
  expect_true(!sel_right$selected && identical(sel_right$status, "maximal_overlap_retained") &&
                identical(sel_right$selected, rule(sel_right, delta_right, 5, 1)),
              label = "maximal-overlap estimate retained when the mixture adds little")
  expect_true(sel_wrong$selected && identical(sel_wrong$status, "gaussian_mixture_selected") &&
                identical(sel_wrong$selected, rule(sel_wrong, 1.5, 5, 1)),
              label = "mixture estimate selected when gain and separation exceed the margins")
  # The gain recomputed from the free and fixed fits at the selected K.
  starts_wrong <- .ac3_da_gaussian_mixture_starts(beta_sel, sqrt(var_sel), 1.5)
  K_wrong <- sel_wrong$n_components
  free_wrong <- .ac3_da_gaussian_mixture_fit(beta_sel, var_sel, starts_wrong, K_wrong)
  fixed_wrong <- .ac3_da_gaussian_mixture_fit(beta_sel, var_sel, starts_wrong, K_wrong, fix_at = 1.5)
  expect_true(identical(sel_wrong$loglik_gain, free_wrong$loglik - fixed_wrong$loglik),
              label = "gain is the difference l_free - l_fixed of the maximized log-likelihoods")
  sel_gain <- .ac3_da_shift_selection(beta_sel, var_sel, delta_mo = 1.5, center = 1.5, min_gain = Inf)
  sel_sep <- .ac3_da_shift_selection(beta_sel, var_sel, delta_mo = 1.5, center = 1.5, min_separation = Inf)
  expect_true(!sel_gain$selected && !sel_sep$selected &&
                identical(sel_gain$status, "maximal_overlap_retained") &&
                identical(sel_sep$status, "maximal_overlap_retained"),
              label = "no selection when either the gain or the separation criterion fails")
  sel_open <- .ac3_da_shift_selection(beta_sel, var_sel, delta_mo = delta_right, center = delta_right,
                                      min_gain = -Inf, min_separation = -Inf)
  expect_true(sel_open$selected && identical(sel_open$selected, rule(sel_open, delta_right, -Inf, -Inf)),
              label = "selection when both margins are removed")
  expect_true(isTRUE(all.equal(sum(sel_wrong$nu_gm), 1, tolerance = 1e-14)) &&
                isTRUE(all.equal(sel_wrong$delta_gm, sum(sel_wrong$nu_gm * beta_sel), tolerance = 1e-14)),
              label = "mixture weights sum to one and reproduce the mixture estimate")
  sel_few <- .ac3_da_shift_selection(beta_sel[1:9], var_sel[1:9], delta_mo = 0, center = 0)
  expect_true(!sel_few$selected && identical(sel_few$status, "insufficient_taxa"),
              label = "fewer than ten coefficients: maximal-overlap estimate retained")
})

test_that("empirical null: the scale is the upper limit of the profile region", {
  set.seed(7L)
  z_wide <- stats::rnorm(300L, sd = 1.3)
  en_wide <- .ac3_da_empirical_null(z_wide)
  z_narrow <- stats::rnorm(300L, sd = 0.7)
  en_narrow <- .ac3_da_empirical_null(z_narrow)
  z_mix <- c(stats::rnorm(270L, sd = 1.2), stats::rnorm(30L, mean = 6))
  en_mix <- .ac3_da_empirical_null(z_mix)
  en_few <- .ac3_da_empirical_null(stats::rnorm(9L))
  # Upper limit of the profile region. For K = 0 the profile is
  # 2 n {log(omega / omega_ml) - (m2 / 2) (1 / omega_ml^2 - 1 / omega^2)} with
  # m2 the mean square of the statistics and omega_ml^2 = max(m2, 1).
  profile_k0 <- function(z, omega) {
    m2 <- mean(z^2); o2 <- max(m2, 1)
    2 * length(z) * (log(omega) - 0.5 * log(o2) - 0.5 * m2 * (1 / o2 - 1 / omega^2))
  }
  expect_true(identical(en_wide$factor, en_wide$omega_hat) && en_wide$omega_hat > sqrt(mean(z_wide^2)) &&
                profile_k0(z_wide, en_wide$omega_hat) <= 3.84 &&
                profile_k0(z_wide, en_wide$omega_hat + 0.005) > 3.84,
              label = "normal statistics: the scale is the upper limit of the profile region")
  expect_true(en_narrow$omega_hat >= 1 && identical(en_narrow$factor, en_narrow$omega_hat) &&
                profile_k0(z_narrow, en_narrow$omega_hat) <= 3.84 &&
                profile_k0(z_narrow, en_narrow$omega_hat + 0.005) > 3.84,
              label = "narrow null: the scale is at least one and is the upper limit of the profile region")
  expect_true(en_mix$n_components >= 1L && en_mix$omega_hat > sqrt(mean(z_mix[1:270]^2)) &&
                en_mix$omega_hat < 1.6 && identical(en_mix$factor, en_mix$omega_hat),
              label = "separated further component: upper limit above the root mean square of the central component")
  expect_true(identical(en_few$factor, 1) && is.na(en_few$omega_hat),
              label = "fewer than ten statistics give the scale one")
})

test_that("correction weights of score_ancombc3_da() sum to one and reproduce the common shift", {
  set.seed(3L)
  p <- 60L
  n <- 20L
  taxa <- paste0("T", seq_len(p))
  se_a <- 0.08 + 0.04 * (seq_len(p) %% 3L)
  beta_a <- c(stats::rnorm(45L, 0.3, se_a[1:45]), stats::rnorm(15L, 2, 0.5))
  beta_b <- c(stats::rnorm(3L, 0, 0.1), rep(NA_real_, p - 3L))
  raw <- rbind(
    data.frame(taxon = taxa, term = "a", beta_raw = beta_a, var_raw_model = 0.8 * se_a^2,
               var_raw_marginal = se_a^2, se_raw_marginal = se_a, n_positive = 50L,
               n_abundance_coef = 3L, stringsAsFactors = FALSE),
    data.frame(taxon = taxa, term = "b", beta_raw = beta_b, var_raw_model = 0.8 * se_a^2,
               var_raw_marginal = se_a^2, se_raw_marginal = se_a, n_positive = 50L,
               n_abundance_coef = 3L, stringsAsFactors = FALSE)
  )
  variance_inputs <- list(
    taxa = taxa,
    G = list(a = matrix(stats::rnorm(p * n, sd = 0.02), p, n),
             b = matrix(stats::rnorm(p * n, sd = 0.02), p, n)),
    V = diag(0.05, n)
  )
  weight_sums <- function(scored) {
    vapply(split(scored, scored$term), function(g) sum(g$correction_weight), numeric(1))
  }
  shift_reproduced <- function(scored) {
    all(vapply(split(scored, scored$term), function(g) {
      isTRUE(all.equal(g$delta[1L], sum(g$correction_weight * g$beta_raw), tolerance = 1e-12))
    }, logical(1)))
  }
  scored_mo <- score_ancombc3_da(raw, variance_inputs = variance_inputs)
  scored_gm <- score_ancombc3_da(raw, variance_inputs = variance_inputs,
                                 shift_selection_min_gain = -Inf, shift_selection_min_separation = -Inf)
  term_a_mo <- scored_mo[scored_mo$term == "a", ]
  term_a_gm <- scored_gm[scored_gm$term == "a", ]
  term_b <- scored_mo[scored_mo$term == "b", ]
  expect_true(identical(unique(term_a_mo$shift_selection_status), "maximal_overlap_retained") &&
                isTRUE(all.equal(weight_sums(scored_mo)[["a"]], 1, tolerance = 1e-12)) &&
                sum(term_a_mo$correction_weight > 0) == term_a_mo$maximal_overlap_n_null_set[1L],
              label = "maximal-overlap correction weights sum to one over the null set")
  expect_true(identical(unique(term_a_gm$shift_selection_status), "gaussian_mixture_selected") &&
                isTRUE(unique(term_a_gm$gaussian_mixture_selected)) &&
                isTRUE(all.equal(weight_sums(scored_gm)[["a"]], 1, tolerance = 1e-12)),
              label = "mixture correction weights sum to one")
  expect_true(nrow(term_b) == 3L && all(term_b$correction_weight == 0) && all(term_b$delta == 0) &&
                identical(unique(term_b$maximal_overlap_status), "insufficient_taxa"),
              label = "degenerate case: all correction weights zero")
  expect_true(shift_reproduced(scored_mo) && shift_reproduced(scored_gm),
              label = "common-shift estimate equals sum(correction_weight * beta_raw)")
  expect_true(all(scored_mo$in_null_set == (scored_mo$correction_weight > 0)),
              label = "in_null_set flags the taxa with a positive weight (maximal overlap)")
})

test_that("DA eligibility: the fitting set is the union of the eligible taxa of the tested covariates", {
  set.seed(11L)
  n_e <- 60L
  p_e <- 40L
  x_e <- cbind(cont = stats::rnorm(n_e), binary = rep(0:1, length.out = n_e))
  prevalence <- stats::runif(p_e, 0.05, 0.95)
  y_e <- matrix(stats::rpois(n_e * p_e, 20) * stats::rbinom(n_e * p_e, 1L, rep(prevalence, each = n_e)),
                n_e, p_e, dimnames = list(NULL, paste0("T", seq_len(p_e))))
  at_risk_e <- matrix(rep(pmin(prevalence + 0.05, 1), each = n_e), n_e, p_e)
  elig <- ancombc3_da_eligibility(at_risk_prob = at_risk_e, y = y_e, x = x_e,
                                  test_terms = c("cont", "binary"), min_observed_prevalence = 0.25)
  term_specs <- list(
    cont = list(term = "cont", type = "continuous", coef_name = "cont", columns = "cont"),
    binary = list(term = "binary", type = "binary", coef_name = "binary", columns = "binary")
  )
  elig_list <- ancombc3_da_eligibility(at_risk_prob = at_risk_e, y = y_e, x = x_e,
                                       test_terms = term_specs, min_observed_prevalence = 0.25)
  per_term <- lapply(term_specs, function(term) {
    e <- ancombc3_da_eligibility_term(at_risk_e, y_e, x_e, term, min_observed_prevalence = 0.25)
    e$taxon[e$da_eligible]
  })
  expect_true(identical(elig$fitting_set, unique(unlist(elig$by_term, use.names = FALSE))),
              label = "fitting_set is the union of by_term")
  expect_true(length(elig$fitting_set) > 0L && length(elig$fitting_set) < p_e,
              label = "fitting set is a proper, nonempty subset of the taxa")
  expect_true(identical(elig$by_term, per_term), label = "by_term agrees with ancombc3_da_eligibility_term()")
  expect_true(identical(elig_list, elig), label = "term specifications and column names give the same result")
  expect_true(nrow(elig$diagnostics) == 2L * p_e && identical(unique(elig$diagnostics$term), c("cont", "binary")),
              label = "diagnostics hold one row per taxon and term")
})
