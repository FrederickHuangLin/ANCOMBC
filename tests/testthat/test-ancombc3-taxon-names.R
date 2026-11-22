# Taxon identification by name in the first stage and the downstream steps.

# Simulated counts: n samples, p taxa, a continuous and a binary covariate in
# both designs, Poisson-log-normal counts with abundance effects on taxa 1 to 6
# and structural zeros whose prevalence depends on the binary covariate for
# taxa 7 to 12; taxa 25 to 30 have a low prevalence.
set.seed(20261103L)
n <- 60L
p <- 30L
x <- cbind(cont = as.numeric(scale(stats::rnorm(n))), binary = rep(0:1, length.out = n))
r <- x
beta <- matrix(0, p, 2L)
beta[1:3, 1L] <- 0.8
beta[4:6, 2L] <- -1
gamma_bin <- c(rep(0, 6L), rep(-2.5, 6L), rep(0, p - 12L))
a_prev <- c(rep(2, p - 6L), rep(-2.5, 6L))
log_mu <- matrix(stats::rnorm(p, 3, 0.8), n, p, byrow = TRUE) + x %*% t(beta) +
  stats::rnorm(n, sd = 0.3) + matrix(stats::rnorm(n * p, sd = 0.5), n, p)
presence <- matrix(stats::rbinom(n * p, 1L, stats::plogis(outer(rep(1, n), a_prev) + outer(x[, "binary"], gamma_bin))),
                   n, p)
y <- matrix(stats::rpois(n * p, exp(log_mu)), n, p) * presence
dimnames(y) <- list(paste0("S", seq_len(n)), paste0("T", seq_len(p)))
test_terms <- c("cont", "binary")
new_name <- stats::setNames(sprintf("otu_%s", rev(utils::head(
  as.vector(outer(letters, letters, paste0)), p))), colnames(y))
old_name <- stats::setNames(names(new_name), new_name)
perm <- sample.int(p)
epochs <- 5L
y_renamed <- y
colnames(y_renamed) <- unname(new_name[colnames(y)])

first_stage <- function(y) {
  ancombc3_fitted_vae_offsets(y = y, x = x, r = r, tested_terms = test_terms, epochs = epochs, seed = 1L)
}
first_stage_original <- ac3_cached(function() first_stage(y))
matrix_outputs <- c("abundance_offset", "structural_offset", "completed_zero",
                    "sampling_zero_posterior", "at_risk_prob", "structural_absence_prob",
                    "structural_absence_weight", "sampling_zero_prob")
# First-stage outputs with the taxa renamed by new_name and the columns
# (elements of sigma_vae) in the order cols.
rename_offsets <- function(offsets, cols = seq_len(p)) {
  for (nm in matrix_outputs) {
    colnames(offsets[[nm]]) <- unname(new_name[colnames(offsets[[nm]])])
    offsets[[nm]] <- offsets[[nm]][, cols, drop = FALSE]
  }
  names(offsets$sigma_vae) <- unname(new_name[names(offsets$sigma_vae)])
  offsets$sigma_vae <- offsets$sigma_vae[cols]
  offsets
}

# Downstream steps on counts y and first-stage outputs offsets; one row per
# taxon and tested covariate, with the taxa named by old (the original names).
downstream <- function(y, offsets, old = stats::setNames(colnames(y), colnames(y))) {
  elig <- ancombc3_da_eligibility(at_risk_prob = offsets$at_risk_prob, y = y, x = x,
                                  test_terms = test_terms)
  fit <- fit_ancombc3_da(
    y = y, x = x, r = r, sigma_vae = offsets$sigma_vae,
    abundance_offset = offsets$abundance_offset, structural_offset = offsets$structural_offset,
    sampling_zero_prob = offsets$sampling_zero_prob, at_risk_prob = offsets$at_risk_prob,
    test_terms = test_terms, fitting_set = elig$fitting_set, eligible_taxa_by_term = elig$by_term,
    n_starts = 4L
  )
  scored <- score_ancombc3_da(fit$raw, variance_inputs = fit$da_variance_inputs,
                              eligible_taxa_by_term = elig$by_term)
  dp_elig <- dp_eligibility(structural_absence_weight = offsets$structural_absence_weight,
                            r = r, test_terms = test_terms)
  dp <- fit_ancombc3_dp_zero_pattern(y = y, r = r, f0 = offsets$sampling_zero_prob,
                                     test_terms = test_terms, eligibility = dp_elig,
                                     da_fit = fit, x = x)
  taxa <- colnames(y)
  out <- do.call(rbind, lapply(test_terms, function(tn) {
    blank <- stats::setNames(rep(NA_real_, p), taxa)
    est <- p_da <- blank
    s <- scored[scored$term == tn, , drop = FALSE]
    est[s$taxon] <- s$beta_corrected
    p_da[s$taxon] <- s$p_value
    d <- dp[dp$test_term == tn, , drop = FALSE]
    p_dp <- dp_est <- blank
    p_dp[d$taxon] <- d$p_value
    dp_est[d$taxon] <- d$estimate
    dp_ok <- stats::setNames(taxa %in% d$taxon[d$dp_eligible %in% TRUE], taxa)
    has_da <- stats::setNames(taxa %in% elig$by_term[[tn]], taxa) & is.finite(p_da)
    has_dp <- dp_ok & is.finite(p_dp)
    # p_dp and has_dp are passed in reverse taxon order; they are matched by
    # name to the order of has_da.
    p_omni <- omnibus_pvalue(p_da, rev(p_dp), has_da, rev(has_dp))
    discovered <- (has_da | has_dp) & is.finite(p_omni) & p_omni <= 0.05
    label <- ancombc3_component_labels(discovered, has_da, rev(has_dp), p_da, rev(p_dp))
    data.frame(taxon = unname(old[taxa]), term = tn, da_eligible = unname(has_da),
               da_estimate = unname(est), da_p = unname(p_da), dp_eligible = unname(dp_ok),
               dp_estimate = unname(dp_est), dp_p = unname(p_dp), omnibus_p = p_omni,
               component_label = label, stringsAsFactors = FALSE)
  }))
  out <- out[order(out$term, match(out$taxon, paste0("T", seq_len(p)))), , drop = FALSE]
  rownames(out) <- NULL
  list(table = out, d_hat = fit$d_hat,
       sigma = stats::setNames(unname(fit$sigma_hat), old[names(fit$sigma_hat)])[paste0("T", seq_len(p))])
}
max_diff <- function(a, b) {
  num <- c("da_estimate", "da_p", "dp_estimate", "dp_p", "omnibus_p")
  m <- vapply(num, function(nm) {
    u <- a$table[[nm]]
    v <- b$table[[nm]]
    if (!identical(is.na(u), is.na(v))) return(Inf)
    if (all(is.na(u))) 0 else max(abs(u - v), na.rm = TRUE)
  }, numeric(1))
  c(m, d_hat = max(abs(a$d_hat - b$d_hat)), sigma = max(abs(a$sigma - b$sigma), na.rm = TRUE))
}
same_labels <- function(a, b) {
  identical(a$table[, c("taxon", "term", "da_eligible", "dp_eligible", "component_label")],
            b$table[, c("taxon", "term", "da_eligible", "dp_eligible", "component_label")])
}

test_that("first-stage outputs carry the taxon names; renaming changes the names only", {
  skip_if_no_torch()
  off <- first_stage_original()
  expect_true(all(vapply(matrix_outputs, function(nm) identical(dimnames(off[[nm]]), dimnames(y)), logical(1))),
              label = "every n x p output has the dimnames of the counts")
  expect_true(identical(names(off$sigma_vae), colnames(y)), label = "sigma_vae is named by taxon")
  off_renamed <- first_stage(y_renamed)
  expect_true(identical(rename_offsets(off), off_renamed),
              label = "renamed counts: identical outputs with the new names")
})

test_that("downstream steps match taxa by name", {
  skip_if_no_torch()
  skip_if_not_extended()
  off <- first_stage_original()
  base <- downstream(y, off)
  expect_true(sum(base$table$da_eligible) >= 10L && sum(base$table$dp_eligible) >= 2L,
              label = "the example has DA-eligible and DP-eligible taxa")
  # Counts renamed in their original order; first-stage outputs renamed and
  # permuted.
  res_a <- downstream(y_renamed, rename_offsets(off, perm), old = old_name)
  diff_a <- max_diff(base, res_a)
  expect_true(all(diff_a == 0) && same_labels(base, res_a),
              label = "permuted first-stage outputs: identical results by name")
  # Counts permuted and renamed; first-stage outputs renamed in the original
  # order. The sums over taxa in the estimation of the sample-specific bias are
  # then formed in a different order.
  y_perm <- y_renamed[, perm, drop = FALSE]
  res_b <- downstream(y_perm, rename_offsets(off), old = old_name)
  diff_b <- max_diff(base, res_b)
  expect_true(all(diff_b <= 1e-6) && same_labels(base, res_b),
              label = "permuted counts: same results by name (to 1e-6)")
  # Unnamed first-stage outputs of the original order would be attached to the
  # permuted taxa by position; the named outputs remove this dependence.
  expect_true(ac3_errs(fit_ancombc3_dp_zero_pattern(y = y_perm, r = r, f0 = off$sampling_zero_prob,
                                                    test_terms = test_terms), "no column for"),
              label = "first-stage outputs with other taxon names give an error")
})

test_that("DA eligibility attaches to the right taxa with a subset of the taxa", {
  skip_if_no_torch()
  off <- first_stage_original()
  elig_full <- ancombc3_da_eligibility(off$at_risk_prob, y, x, test_terms)
  sub <- c(12L, 3L, 27L, 8L, 1L, 29L, 19L, 5L, 22L, 14L, 9L, 26L)
  y_sub <- y[, sub, drop = FALSE]
  expected <- elig_full$diagnostics[elig_full$diagnostics$taxon %in% colnames(y_sub), , drop = FALSE]
  expected <- expected[order(expected$term, match(expected$taxon, colnames(y_sub))), , drop = FALSE]
  compare_sub <- function(e) {
    got <- e$diagnostics[order(e$diagnostics$term, match(e$diagnostics$taxon, colnames(y_sub))), ,
                         drop = FALSE]
    rownames(got) <- rownames(expected) <- NULL
    identical(got, expected) &&
      identical(e$by_term, lapply(elig_full$by_term, function(v) colnames(y_sub)[colnames(y_sub) %in% v]))
  }
  elig_named <- ancombc3_da_eligibility(off$at_risk_prob, y_sub, x, test_terms)
  expect_true(compare_sub(elig_named), label = "at-risk matrix of all taxa, matched by name")
  at_risk_sub <- unname(off$at_risk_prob[, sub, drop = FALSE])
  expect_true(compare_sub(ancombc3_da_eligibility(at_risk_sub, y_sub, x, test_terms)),
              label = "unnamed at-risk matrix of the subset, matched by position")
  expect_true(ac3_errs(ancombc3_da_eligibility(off$at_risk_prob[, -12L], y_sub, x, test_terms), "no column for"),
              label = "named at-risk matrix without a taxon of the counts gives an error")
  expect_true(sum(elig_named$diagnostics$da_eligible) > 0L &&
                sum(!elig_named$diagnostics$da_eligible) > 0L,
              label = "the subset has eligible and ineligible taxa")
})

test_that("ancombc3() on the renamed count table gives the results with the new taxon names", {
  skip_if_no_torch()
  meta <- data.frame(cont = x[, "cont"], binary = x[, "binary"], row.names = rownames(y))
  run_wrapper <- function(counts) {
    suppressWarnings(ancombc3(data = counts, meta_data = meta, fix_formula = ~ cont + binary,
                              tested = test_terms, taxa_are_rows = FALSE, epochs = epochs,
                              seed = 1L, verbose = FALSE))
  }
  w0 <- run_wrapper(y)
  w1 <- run_wrapper(y_renamed)
  r1 <- w1$res
  r1$taxon <- unname(old_name[r1$taxon])
  expect_true(identical(r1, w0$res), label = "results identical with the new taxon names")
  expect_true(identical(unname(w1$sigma), unname(w0$sigma)) &&
                identical(names(w1$sigma), unname(new_name[names(w0$sigma)])) &&
                identical(w1$sample_bias, w0$sample_bias),
              label = "sigma and sample_bias identical")
  y_blank <- y
  colnames(y_blank)[3L] <- ""
  expect_true(ac3_errs(run_wrapper(y_blank), "missing or empty taxon names"),
              label = "an empty taxon name gives an error")
})
