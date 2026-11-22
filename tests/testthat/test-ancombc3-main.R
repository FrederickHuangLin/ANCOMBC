# User-level function ancombc3(): input validation, the interface, and the
# multiple-testing adjustment p_adj_method, on the toy data (60 samples, 80
# taxa; covariates group and bmi_z).

toy <- ac3_toy_data()
small <- toy$counts
meta <- toy$meta

# ancombc3() at 20 epochs with the given additional arguments.
run_adj <- function(...) {
  suppressWarnings(ancombc3(
    data = small, meta_data = meta, fix_formula = ~ bmi_z + group,
    tested = c("bmi_z", "group"), alpha = 0.10, epochs = 20L, seed = 1L, verbose = FALSE,
    ...
  ))
}

test_that("input validation: errors for invalid input; all-zero samples and taxa dropped", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  run <- function(counts = small, metadata = meta, tested = c("bmi_z", "group"),
                  formula = ~ bmi_z + group) {
    ancombc3(data = counts, meta_data = metadata, fix_formula = formula, tested = tested, epochs = 1L,
             verbose = FALSE)
  }
  unnamed <- small
  colnames(unnamed) <- NULL
  expect_true(ac3_errs(run(counts = unnamed), "no sample names"), label = "missing sample names")
  expect_true(ac3_errs(run(metadata = meta[-1, ]), "have no row in meta_data"), label = "sample without metadata")
  frac <- small
  frac[1, 1] <- 0.5
  expect_true(ac3_errs(run(counts = frac), "non-negative integers"), label = "non-integer counts")
  neg <- small
  neg[1, 1] <- -1
  expect_true(ac3_errs(run(counts = neg), "negative"), label = "negative counts")
  m2 <- meta
  m2$const <- 1
  expect_true(ac3_errs(run(metadata = m2, tested = "const", formula = ~ const + group), "constant"),
              label = "constant tested covariate")
  m3 <- meta
  m3$grp <- factor(rep(c("a", "b", "c"), length.out = nrow(m3)))
  expect_true(ac3_errs(run(metadata = m3, tested = "grp", formula = ~ grp), "levels"),
              label = "tested factor with three levels")
  expect_true(ac3_errs(run(counts = small[, 1:8], metadata = meta[1:8, ]), "Too few samples"),
              label = "too few samples")
  m1 <- meta
  m1$binary_group <- meta$group
  expect_true(ac3_errs(run(metadata = m1, tested = "binary_group", formula = ~ bmi_z + group), "not a term"),
              label = "tested covariate outside the formula")
  # The first sample and the first taxon have no positive count, and the
  # covariate flag is constant over the other samples, so that the call stops
  # after the samples and taxa are dropped.
  z <- small
  z[, 1] <- 0
  z[1, ] <- 0
  m4 <- meta
  m4$flag <- as.numeric(seq_len(nrow(m4)) == 1L)
  msgs <- character(0)
  err <- withCallingHandlers(
    tryCatch(ancombc3(data = z, meta_data = m4, fix_formula = ~ bmi_z + flag, tested = "flag", verbose = TRUE),
             error = function(e) conditionMessage(e)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  expect_true(any(grepl("1 sample\\(s\\) with no positive count: sample_01", msgs)) &&
                any(grepl("taxon/taxa with no positive count", msgs)) &&
                grepl("constant", err),
              label = "all-zero samples and taxa dropped with a message")
})

test_that("interface: the argument names of v0.3.1 are unused arguments", {
  old_names <- list(
    counts = quote(ancombc3(counts = small, meta_data = meta, fix_formula = ~ bmi_z, tested = "bmi_z")),
    metadata = quote(ancombc3(data = small, metadata = meta, fix_formula = ~ bmi_z, tested = "bmi_z")),
    abundance_formula = quote(ancombc3(data = small, meta_data = meta, abundance_formula = ~ bmi_z,
                                       tested = "bmi_z")),
    prevalence_formula = quote(ancombc3(data = small, meta_data = meta, fix_formula = ~ bmi_z,
                                        prevalence_formula = ~ bmi_z, tested = "bmi_z")),
    n_cores = quote(ancombc3(data = small, meta_data = meta, fix_formula = ~ bmi_z, tested = "bmi_z",
                             n_cores = 2L)),
    dunnett = quote(ancombc3(data = small, meta_data = meta, fix_formula = ~ bmi_z, tested = "bmi_z",
                             dunnett = TRUE)),
    bias_corrected = quote(ancombc3(data = small, meta_data = meta, fix_formula = ~ bmi_z, tested = "bmi_z",
                                    bias_corrected = TRUE))
  )
  for (nm in names(old_names)) {
    msg <- tryCatch({ eval(old_names[[nm]]); "" }, error = function(e) conditionMessage(e))
    expect_true(grepl(paste0("unused argument \\(", nm, " = "), msg), label = paste("old argument", nm, "is unused"))
  }
})

test_that("p_adj_method: BH default, holm and bonferroni q-values and labels, and settings validation", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  res <- suppressWarnings(ancombc3(
    data = small, meta_data = meta, fix_formula = ~ bmi_z + group,
    tested = c("bmi_z", "group"), alpha = 0.10, epochs = 20L, seed = 1L, verbose = FALSE
  ))
  res_bh <- run_adj(p_adj_method = "BH", n_cl = 2L)
  res_holm <- run_adj(p_adj_method = "holm", n_cl = 2L)
  res_bonf <- run_adj(p_adj_method = "bonferroni", n_cl = 2L)
  expect_true(identical(res_bh$res, res$res) && identical(res_bh$sample_bias, res$sample_bias) &&
                identical(res_bh$sigma, res$sigma) && identical(res_bh$samp_frac, res$samp_frac) &&
                identical(res_bh$settings$method, res$settings$method),
              label = "explicit \"BH\" equals the default")
  rh <- res_holm$res
  expect_true(identical(rh$da_p, res$res$da_p) && identical(rh$dp_p, res$res$dp_p) &&
                identical(rh$omnibus_p, res$res$omnibus_p),
              label = "holm: p-values unchanged")
  holm_ok <- vapply(unique(rh$covariate), function(cv) {
    g <- rh[rh$covariate == cv, , drop = FALSE]
    da <- is.finite(g$da_p)
    dp <- is.finite(g$dp_p)
    om <- is.finite(g$omnibus_p)
    isTRUE(all.equal(g$da_q[da], stats::p.adjust(g$da_p[da], "holm"), tolerance = 0)) &&
      isTRUE(all.equal(g$dp_q[dp], stats::p.adjust(g$dp_p[dp], "holm"), tolerance = 0)) &&
      isTRUE(all.equal(g$omnibus_q[om], stats::p.adjust(g$omnibus_p[om], "holm"), tolerance = 0)) &&
      identical(g$da_discovery, (g$da_q <= 0.10) %in% TRUE)
  }, logical(1))
  expect_true(all(holm_ok), label = "holm: DA, DP, and omnibus q-values equal p.adjust(p, \"holm\")")
  expect_true(identical(res_holm$settings$p_adj_method, "holm"), label = "holm: recorded in settings")
  rb <- res_bonf$res
  label_ok <- vapply(unique(rb$covariate), function(cv) {
    g <- rb[rb$covariate == cv, , drop = FALSE]
    has_da <- is.finite(g$da_p)
    has_dp <- is.finite(g$dp_p)
    universe <- has_da | has_dp
    tau <- 0.10 / sum(universe) / 2
    expected <- ifelse(universe, "not_discovered", "not_in_universe")
    disc <- g$omnibus_discovery
    expected[disc & has_da & !has_dp] <- "DA_driven"
    expected[disc & has_dp & !has_da] <- "DP_driven"
    both <- disc & has_da & has_dp
    da_pass <- has_da & g$da_p <= tau
    dp_pass <- has_dp & g$dp_p <= tau
    expected[both] <- ifelse(da_pass[both] & dp_pass[both], "DA_DP",
                             ifelse(da_pass[both], "DA_driven", ifelse(dp_pass[both], "DP_driven", "omnibus_only")))
    identical(g$component_label, expected) &&
      isTRUE(all.equal(g$omnibus_q[universe], stats::p.adjust(g$omnibus_p[universe], "bonferroni"), tolerance = 0))
  }, logical(1))
  expect_true(all(label_ok), label = "bonferroni: component labels agree with a direct computation")
  unit <- ancombc3_component_labels(
    discovered = c(TRUE, TRUE, TRUE, FALSE), has_da = rep(TRUE, 4L), has_dp = rep(TRUE, 4L),
    p_da = c(0.001, 0.03, 0.001, 0.5), p_dp = c(0.03, 0.001, 0.001, 0.5), alpha = 0.1, p_adj_method = "holm"
  )
  expect_true(identical(unit, c("DA_driven", "DP_driven", "DA_DP", "not_discovered")),
              label = "holm: labels at half the cutoff alpha / (|U| - R + 1) = 0.025")
  by_unit <- ancombc3_component_labels(
    discovered = c(TRUE, TRUE), has_da = c(TRUE, TRUE), has_dp = c(TRUE, TRUE),
    p_da = c(0.0329, 0.0335), p_dp = c(0.5, 0.5), alpha = 0.1, p_adj_method = "BY"
  )
  expect_true(identical(by_unit, c("DA_driven", "omnibus_only")), label = "BY: labels at half the cutoff R alpha / (|U| c(|U|))")
  expect_true(ac3_errs(run_adj(p_adj_method = "fdr", n_cl = 2L), "p_adj_method must be one of"),
              label = "unsupported p_adj_method is an error")
  expect_true(ac3_errs(run_adj(settings = list(da_null_scale = "released"), n_cl = 2L), "Unknown element") &&
                ac3_errs(run_adj(settings = list(omnibus_rule = "simes"), n_cl = 2L), "Unknown element"),
              label = "a removed setting is an error")
  expect_true(ac3_errs(run_adj(settings = list(da_p_adjust = "BY"), n_cl = 2L), "Unknown element") &&
                ac3_errs(run_adj(settings = list(omnibus_p_adjust = "BY"), n_cl = 2L), "Unknown element") &&
                ac3_errs(run_adj(settings = list(epochs = 20L), n_cl = 2L), "Unknown element"),
              label = "settings da_p_adjust, omnibus_p_adjust, and epochs are rejected")
  expect_true(ac3_errs(run_adj(settings = list(da_offset_variance = "alw"), n_cl = 2L), "must be one of"),
              label = "an abbreviated value of da_offset_variance is an error")
  expect_true(ac3_errs(run_adj(settings = list(0.02), n_cl = 2L), "must be named") &&
                ac3_errs(run_adj(settings = list(d_tol = 0.02, 50L), n_cl = 2L), "must be named") &&
                ac3_errs(run_adj(settings = list(d_tol = 0.02, d_tol = 0.03), n_cl = 2L), "duplicated names"),
              label = "an unnamed, partially named, or duplicated settings list is an error")
})
