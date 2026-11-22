# Small deterministic second-stage fit on simulated inputs.

# Simulated inputs: n samples, p taxa, a continuous and a binary covariate in
# both designs, zero-inflated Poisson-log-normal counts with a sample-specific
# bias, and simple first-stage quantities.
set.seed(20261005L)
n <- 40L
p <- 30L
x <- cbind(cont = as.numeric(scale(stats::rnorm(n))), binary = rep(0:1, length.out = n))
r <- x
alpha <- stats::rnorm(p, mean = 2, sd = 0.8)
beta <- matrix(0, p, 2L)
beta[1:4, 1L] <- 0.8
beta[5:8, 2L] <- -1
bias <- stats::rnorm(n, sd = 0.3)
log_mu <- matrix(alpha, n, p, byrow = TRUE) + x %*% t(beta) + bias
lambda <- exp(log_mu + matrix(stats::rnorm(n * p, sd = 0.5), n, p))
present <- matrix(stats::rbinom(n * p, 1L, 0.9), n, p)
y <- matrix(stats::rpois(n * p, lambda), n, p) * present
colnames(y) <- paste0("T", seq_len(p))
rownames(y) <- paste0("S", seq_len(n))
sigma_vae <- rep(0.5, p)
mu_first_stage <- matrix(pmax(colMeans(y), 0.1), n, p, byrow = TRUE)
sampling_zero_prob <- (2 / (2 + mu_first_stage))^2
at_risk_prob <- matrix(0.9, n, p)
test_terms <- c("cont", "binary")

# DA eligibility of both tested covariates.
elig <- ancombc3_da_eligibility(at_risk_prob = at_risk_prob, y = y, x = x, test_terms = test_terms)

fit_args <- list(
  y = y, x = x, r = r, sigma_vae = sigma_vae,
  abundance_offset = matrix(0, n, p), structural_offset = matrix(0, n, p),
  sampling_zero_prob = sampling_zero_prob, at_risk_prob = at_risk_prob,
  test_terms = test_terms, fitting_set = elig$fitting_set,
  eligible_taxa_by_term = elig$by_term, n_starts = 4L
)
fit_one_core <- ac3_cached(function() do.call(fit_ancombc3_da, c(fit_args, list(n_cores = 1L))))

test_that("DA eligibility: the fitting set is the union of the eligible taxa", {
  expect_true(identical(names(elig$by_term), test_terms) &&
                identical(elig$fitting_set, unique(unlist(elig$by_term, use.names = FALSE))) &&
                length(elig$fitting_set) >= 10L,
              label = "eligibility: fitting set is the union of the eligible taxa")
  expect_true(identical(unique(elig$diagnostics$term_type), c("continuous", "binary")),
              label = "eligibility: term types inferred from the columns of x")
})

test_that("second-stage fit: identical with one and two workers, converged, default eligibility", {
  fit1 <- fit_one_core()
  fit2 <- do.call(fit_ancombc3_da, c(fit_args, list(n_cores = 2L)))
  same <- vapply(names(fit1), function(nm) identical(fit1[[nm]], fit2[[nm]]), logical(1))
  expect_true(identical(names(fit1), names(fit2)) && all(same), label = "fit identical for n_cores = 1 and n_cores = 2")
  expect_true(sum(vapply(fit1$fits, function(f) isTRUE(f$ok), logical(1))) >= p - 2L,
              label = "taxon-wise fits succeed (at most two failures)")
  expect_true(isTRUE(fit1$converged), label = "alternation converged")
  expect_true(all(fit1$fitting_set %in% elig$fitting_set) && identical(fit1$eligible_taxa_by_term, elig$by_term),
              label = "fit records the fitting set and the eligible taxa")
  expect_true(isTRUE(all.equal(fit1$objective_L,
                               fit1$alternation_history$L[nrow(fit1$alternation_history)], tolerance = 0)),
              label = "objective_L is the last L of the alternation history")
  # Default eligibility: fit without fitting_set and eligible_taxa_by_term.
  default_args <- fit_args[setdiff(names(fit_args), c("fitting_set", "eligible_taxa_by_term"))]
  fit3 <- do.call(fit_ancombc3_da, c(default_args, list(n_cores = 1L)))
  union3 <- unique(unlist(fit3$eligible_taxa_by_term, use.names = FALSE))
  expect_true(identical(fit3$eligible_taxa_by_term, elig$by_term) &&
                identical(fit3$fitting_set, intersect(colnames(y), union3)),
              label = "default eligibility: sets stored, fitting set is their union")
  expect_true(identical(fit3, fit1), label = "default eligibility: fit identical to the fit with the sets")
})

test_that("score_ancombc3_da() on the fit: p-values, empirical-null scale, and working variance", {
  fit1 <- fit_one_core()
  scored <- score_ancombc3_da(fit1$raw, variance_inputs = fit1$da_variance_inputs,
                              eligible_taxa_by_term = elig$by_term)
  n_coef <- sum(is.finite(fit1$raw$beta_raw) & is.finite(fit1$raw$var_raw_marginal))
  expect_true(nrow(scored) == n_coef, label = "scored table has one row per finite coefficient")
  expect_true(all(c("beta_corrected", "se_corrected", "statistic", "statistic_normal_equivalent",
                    "p_value", "q_value", "delta", "correction_weight", "adjustment_factor",
                    "var_primary", "var_raw", "se_raw", "in_null_set", "maximal_overlap_status",
                    "shift_selection_status") %in% names(scored)),
              label = "scored table has the scoring columns")
  expect_true(all(is.finite(scored$p_value) & scored$p_value >= 0 & scored$p_value <= 1),
              label = "P values are finite and in [0, 1]")
  expect_true(all(is.finite(scored$adjustment_factor) & scored$adjustment_factor >= 1),
              label = "empirical-null scale is finite and at least one")
  expect_true(all(is.finite(scored$var_primary) & scored$var_primary > 0) &&
                isTRUE(all.equal(scored$se_corrected, sqrt(scored$var_primary), tolerance = 0)),
              label = "var_primary is positive and se_corrected is its square root")
  expect_true(all(scored$method == "ancombc3_da"), label = "method label of the DA test")
  expect_true(all(vapply(split(scored, scored$term), function(g) {
    isTRUE(all.equal(g$delta[1L], sum(g$correction_weight * g$beta_raw), tolerance = 1e-12))
  }, logical(1))), label = "shift equals sum(correction_weight * beta_raw) per covariate")
})

test_that("zero-inflated Poisson-log-normal cell derivatives agree with central differences", {
  # Analytic first derivatives in the log mean and the presence logit against
  # central differences of the log-likelihood, and analytic second derivatives
  # against central differences of the analytic first derivatives (the cross
  # derivative in both orders), with step 1e-4 in each predictor, at several
  # values of sigma. Errors are relative to max(1, |analytic derivative|). No
  # predictor or probability is truncated.
  set.seed(1L)
  m <- 400L
  log_mu_o <- stats::runif(m, -4, 6)
  logit_pi_o <- stats::runif(m, -3, 4)
  y_o <- ifelse(stats::runif(m) < 0.4, 0, stats::rpois(m, exp(log_mu_o)))
  pos_o <- y_o > 0
  hd <- 1e-4
  rel_err <- function(a, b) max(abs(a - b) / pmax(1, abs(a)))
  deriv_err <- do.call(rbind, lapply(c(0.05, 0.5, 1.5, 2.5), function(sg) {
    deriv <- function(lm, lp) .ac3_pln_observation_derivatives(lm, lp, y_o, pos_o, sg, need = 2L)
    o <- deriv(log_mu_o, logit_pi_o)
    mp <- deriv(log_mu_o + hd, logit_pi_o)
    mm <- deriv(log_mu_o - hd, logit_pi_o)
    pp <- deriv(log_mu_o, logit_pi_o + hd)
    pm <- deriv(log_mu_o, logit_pi_o - hd)
    data.frame(
      sigma = sg,
      dl_mu = rel_err(o$dl_mu, (mp$ll - mm$ll) / (2 * hd)),
      dl_pi = rel_err(o$dl_pi, (pp$ll - pm$ll) / (2 * hd)),
      d2l_mumu = rel_err(o$d2l_mumu, (mp$dl_mu - mm$dl_mu) / (2 * hd)),
      d2l_pipi = rel_err(o$d2l_pipi, (pp$dl_pi - pm$dl_pi) / (2 * hd)),
      d2l_mupi = max(rel_err(o$d2l_mupi, (pp$dl_mu - pm$dl_mu) / (2 * hd)),
                     rel_err(o$d2l_mupi, (mp$dl_pi - mm$dl_pi) / (2 * hd)))
    )
  }))
  err_first <- max(deriv_err$dl_mu, deriv_err$dl_pi)
  err_second <- max(deriv_err$d2l_mumu, deriv_err$d2l_pipi, deriv_err$d2l_mupi)
  expect_true(is.finite(err_first) && err_first <= 1e-6, label = "PLN analytic first derivatives match finite differences")
  expect_true(is.finite(err_second) && err_second <= 1e-6,
              label = "PLN analytic second derivatives match finite differences")
})
