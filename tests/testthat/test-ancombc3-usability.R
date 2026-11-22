# Reporting functions of ancombc3() and the progress arguments of the stage
# functions, on the first 40 taxa of the toy data.

toy <- ac3_toy_data()
counts <- toy$counts[1:40, ]
meta <- toy$meta
y <- t(counts)
storage.mode(y) <- "double"
design <- ancombc3_make_design(meta[rownames(y), ], x_formula = "group + bmi_z", test_terms = c("group", "bmi_z"),
                               r_formula = "group + bmi_z")
tt <- names(design$tests)

# Constructed training history whose validation loss decreases after the
# warm-up.
epochs <- 1:100
kl <- pmin(1, epochs / 40)
decreasing <- data.frame(epoch = epochs, train_loss = 100, val_loss = 100 * 0.997^epochs, kl_weight = kl)

# ancombc3() with verbose = FALSE, and with verbose = TRUE and a progress
# reporter installed, with the reporters installed and the messages issued.
verbose_runs <- ac3_cached(function() {
  run <- function(verbose) {
    suppressWarnings(ancombc3(data = counts, meta_data = meta, fix_formula = ~ group + bmi_z,
                              tested = c("group", "bmi_z"), epochs = 5L, seed = 1L, verbose = verbose))
  }
  quiet <- run(FALSE)
  reporters <- character(0)
  msgs <- character(0)
  loud <- with_mocked_bindings(
    withCallingHandlers(run(TRUE), message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }),
    .ac3_main_progress = function(verbose, kind) {
      reporters <<- c(reporters, kind)
      .ac3_progress_reporter(kind)
    }
  )
  list(quiet = quiet, loud = loud, reporters = reporters, msgs = msgs)
})

test_that("progress callbacks of the stage functions leave the results and the random-number state unchanged", {
  skip_if_no_torch()
  skip_if_not_extended()
  torch::torch_set_num_threads(1L)
  epoch_calls <- list()
  recorder <- function(...) epoch_calls[[length(epoch_calls) + 1L]] <<- c(...)
  set.seed(11L)
  off_plain <- ancombc3_fitted_vae_offsets(y, design$x, design$r, tested_terms = tt, epochs = 6L, seed = 2L)
  seed_plain <- .Random.seed
  set.seed(11L)
  off_cb <- ancombc3_fitted_vae_offsets(y, design$x, design$r, tested_terms = tt, epochs = 6L, seed = 2L,
                                        progress = recorder)
  seed_cb <- .Random.seed
  expect_true(identical(off_plain, off_cb), label = "first stage: identical output with a callback")
  expect_true(identical(seed_plain, seed_cb), label = "first stage: same random-number state after the call")
  expect_true(length(epoch_calls) == 6L && identical(vapply(epoch_calls, `[`, numeric(1), 1L), as.numeric(1:6)) &&
                isTRUE(all.equal(vapply(epoch_calls, `[`, numeric(1), 4L), off_cb$training_history$val_loss)),
              label = "first stage: one callback per epoch with the recorded losses")

  da_args <- list(y = y, x = design$x, r = design$r, sigma_vae = off_plain$sigma_vae,
                  abundance_offset = off_plain$abundance_offset, structural_offset = off_plain$structural_offset,
                  sampling_zero_prob = off_plain$sampling_zero_prob, at_risk_prob = off_plain$at_risk_prob,
                  test_terms = tt, n_starts = 2L)
  fit_calls <- list()
  fit_recorder <- function(done, total) fit_calls[[length(fit_calls) + 1L]] <<- c(done, total)
  fit_plain <- do.call(fit_ancombc3_da, da_args)
  fit_cb <- do.call(fit_ancombc3_da, c(da_args, list(progress = fit_recorder)))
  calls_serial <- do.call(rbind, fit_calls)
  fit_calls <- list()
  fit_par <- do.call(fit_ancombc3_da, c(da_args, list(progress = fit_recorder, n_cores = 2L)))
  calls_par <- do.call(rbind, fit_calls)
  expect_true(identical(fit_plain, fit_cb), label = "second stage: identical fit with a callback")
  expect_true(identical(fit_plain$raw, fit_par$raw) && identical(fit_plain$d_hat, fit_par$d_hat) &&
                identical(fit_plain$da_variance_inputs, fit_par$da_variance_inputs),
              label = "second stage: identical fit with a callback and two workers")
  n_pass <- sum(calls_serial[, 1L] == 0)
  expect_true(n_pass >= 3L && nrow(calls_serial) == n_pass * (ncol(y) + 1L) && all(calls_serial[, 2L] == ncol(y)),
              label = sprintf("second stage: %d passes, one update per taxon", n_pass))
  expect_true(sum(calls_par[, 1L] == 0) == n_pass && nrow(calls_par) == n_pass * 21L &&
                all(calls_par[calls_par[, 1L] > 0, 1L] %in% seq(2, 40, by = 2)),
              label = "second stage with workers: one update per chunk of two taxa")
  set.seed(5L)
  state <- .Random.seed
  pb <- .ac3_progress_reporter("taxon_fits")
  for (k in 0:3) pb$update(k, 3)
  pb$done()
  pe <- .ac3_progress_reporter("epochs")
  for (k in 1:3) pe$update(k, 3, 1, 2, 0.5)
  expect_true(identical(state, .Random.seed), label = "progress reporter leaves the random-number stream unchanged")
})

test_that("training history: one row per epoch, losses passed to the progress function, nothing printed", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  reported <- list()
  printed <- utils::capture.output(
    trained <- train_zipln_vae(y = y, x = design$x, r = design$r, N = rowSums(y), epochs = 12L, seed = 3L,
                               progress = function(...) reported[[length(reported) + 1L]] <<- c(...))
  )
  h <- trained$history
  expect_true(nrow(h) == 12L && identical(h$epoch, 1:12) && all(is.finite(as.matrix(h))),
              label = "history has one row per epoch")
  reported <- do.call(rbind, reported)
  expect_true(length(printed) == 0L && nrow(reported) == 12L &&
                identical(unname(reported[, 3:5]), unname(as.matrix(h[c("train_loss", "val_loss", "kl_weight")]))),
              label = "losses passed to progress are those of the history; nothing is printed")
  expect_true(isTRUE(all.equal(h$kl_weight, pmin(1, (1:12) / 40))), label = "kl_weight follows the warm-up schedule")
  # A last mini-batch of one sample, and the random-number state of the caller.
  n_train <- nrow(y) - floor(0.1 * nrow(y))
  set.seed(99L)
  state <- .Random.seed
  one <- tryCatch(train_zipln_vae(y = y, x = design$x, r = design$r, N = rowSums(y), epochs = 2L, seed = 3L,
                                  batch_size = n_train - 1L),
                  error = function(e) NULL)
  expect_true(!is.null(one) && all(is.finite(as.matrix(one$history))), label = "a mini-batch of one sample is accepted")
  expect_true(identical(state, .Random.seed), label = "train_zipln_vae() leaves the random-number stream unchanged")
})

test_that("convergence check of the first stage on constructed histories", {
  flat <- decreasing
  flat$val_loss <- 80 + 0.001 * sin(epochs)
  short <- decreasing[1:60, ]
  cd <- ancombc3_vae_convergence(decreasing)
  cf <- ancombc3_vae_convergence(flat)
  cs <- ancombc3_vae_convergence(short)
  expect_true(isTRUE(cd$checked) && isTRUE(cd$still_decreasing) && cd$relative_decrease > 0.01 &&
                cd$n_after_warmup == 61L,
              label = sprintf("decreasing history: still decreasing (relative decrease %.4f)", cd$relative_decrease))
  expect_true(isTRUE(cf$checked) && !cf$still_decreasing, label = "flat history: not still decreasing")
  expect_true(!cs$checked && cs$n_after_warmup == 21L, label = "21 epochs after the warm-up: not checked")
  w_dec <- tryCatch(.ac3_main_convergence_warning(cd), warning = function(w) conditionMessage(w))
  w_flat <- tryCatch(.ac3_main_convergence_warning(cf), warning = function(w) conditionMessage(w))
  expect_true(is.character(w_dec) && grepl("still decreasing", w_dec) && isFALSE(w_flat),
              label = "warning for the decreasing history only")
})

test_that("ancombc3() with verbose = FALSE and with verbose = TRUE and a progress reporter", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  runs <- verbose_runs()
  quiet <- runs$quiet
  loud <- runs$loud
  msgs <- runs$msgs
  strip <- function(res) {
    res$settings$call <- NULL
    res$settings$elapsed_seconds <- NULL
    res$settings$stage_seconds <- NULL
    res
  }
  expect_true(identical(runs$reporters, c("epochs", "taxon_fits")),
              label = "progress reporters installed for both long steps")
  expect_true(identical(strip(quiet), strip(loud)), label = "identical results with verbose = TRUE and progress reporters")
  expected_msgs <- c("^Analyzing 60 samples and 40 taxa", "^First stage: training", "^First stage completed",
                     "^DA eligibility:", "^Second stage: taxon-wise fits", "^Second stage completed",
                     "^Common shift of group", "^DP test completed", "^group: ", "^Total time")
  expect_true(all(vapply(expected_msgs, function(pt) any(grepl(pt, msgs)), logical(1))), label = "step messages issued")
  expect_true(identical(quiet$training_history$epoch, 1:5) &&
                identical(names(quiet$training_history), c("epoch", "train_loss", "val_loss", "kl_weight")),
              label = "training_history has one row per epoch")
  expect_true(identical(class(quiet), c("ancombc3", "list")) && is.list(quiet) &&
                identical(quiet$res, quiet[["res"]]) && is.data.frame(quiet$res),
              label = "class c(\"ancombc3\", \"list\") with list access")
  dg <- quiet$diagnostics
  expect_true(identical(names(dg), c("by_covariate", "fit", "first_stage")) &&
                identical(names(dg$by_covariate), c("group", "bmi_z")),
              label = "diagnostics: elements and covariates")
  bc <- dg$by_covariate$group
  expect_true(identical(names(bc), c("n_da_eligible", "da_eligibility", "common_shift", "n_dp_eligible",
                                     "dp_eligibility")) &&
                bc$n_da_eligible == sum(quiet$res$da_eligible[quiet$res$covariate == "group"]) &&
                bc$n_dp_eligible == sum(quiet$res$dp_eligible[quiet$res$covariate == "group"]) &&
                sum(bc$dp_eligibility) == 40L && nrow(bc$da_eligibility) == 40L,
              label = "diagnostics: eligibility counts agree with res")
  shift <- c(group = dg$by_covariate$group$common_shift$delta, bmi_z = dg$by_covariate$bmi_z$common_shift$delta)
  analyzed <- meta[names(quiet$samp_frac), ]
  x_g <- as.numeric(analyzed$group == "treated")
  x_b <- analyzed$bmi_z
  reconstructed <- quiet$sample_bias + shift[["group"]] * (x_g - mean(x_g)) + shift[["bmi_z"]] * (x_b - mean(x_b))
  expect_true(dg$by_covariate$group$common_shift$estimator %in% c("maximal_overlap", "gaussian_mixture") &&
                max(abs(reconstructed - quiet$samp_frac)) <= 1e-12 &&
                is.finite(dg$by_covariate$group$common_shift$adjustment_factor),
              label = "diagnostics: common-shift estimates reproduce samp_frac")
  expect_true(all(c("n_fitting_set", "n_without_usable_fit", "converged", "verification_passes") %in% names(dg$fit)) &&
                identical(dg$first_stage, quiet$settings$first_stage_convergence) && !dg$first_stage$checked,
              label = "diagnostics: fit summary and first-stage check")
})

test_that("print(), summary(), and plot_ancombc3_training()", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  quiet <- verbose_runs()$quiet
  printed <- utils::capture.output(vis <- withVisible(print(quiet)))
  expect_true(!vis$visible && identical(vis$value, quiet) && any(grepl("Discoveries", printed)),
              label = "print() returns the object invisibly")
  sm <- summary(quiet)
  printed <- utils::capture.output(vis <- withVisible(print(sm)))
  expect_true(identical(class(sm), "summary.ancombc3") && !vis$visible && nrow(sm$by_covariate) == 2L &&
                all(c("DA_driven", "DP_driven", "DA_DP", "omnibus_only") %in% names(sm$by_covariate)),
              label = "summary() returns a summary.ancombc3 object, printed invisibly")
  # A result with multigroup tables, assembled from the primary table.
  mg_obj <- quiet
  pair <- data.frame(taxon = rep(c("a", "b"), 3L), comparison = rep(c("B_vs_A", "C_vs_A", "C_vs_B"), each = 2L),
                     da_discovery = c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE), dp_discovery = FALSE,
                     omnibus_discovery = c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE))
  mg_obj$res_global <- data.frame(taxon = c("a", "b"), da_eligible = TRUE, da_discovery = c(TRUE, FALSE),
                                  dp_eligible = TRUE, dp_discovery = FALSE, omnibus_discovery = c(TRUE, FALSE))
  mg_obj$res_pair <- pair
  mg_obj$res_dunn <- pair[1:4, ]
  mg_obj$res_trend <- data.frame(taxon = c("a", "b"), da_discovery = c(TRUE, FALSE), da_pattern = "increasing",
                                 dp_discovery = FALSE, dp_pattern = NA_character_, omnibus_discovery = c(TRUE, FALSE))
  mg_obj$settings$group <- list(name = "grp", levels = c("A", "B", "C"))
  sm_mg <- summary(mg_obj)
  printed <- utils::capture.output(print(mg_obj), print(sm_mg))
  expect_true(identical(sm_mg$multigroup$pairwise$da_discoveries, c(1L, 0L, 2L)) &&
                identical(sm_mg$multigroup$dunnet$comparison, c("B_vs_A", "C_vs_A")) &&
                sm_mg$multigroup$global[["da_discoveries"]] == 1 && any(grepl("Pairwise comparisons", printed)),
              label = "summary() of a result with multigroup tables")
  grDevices::pdf(NULL)
  vis <- withVisible(plot_ancombc3_training(quiet))
  vis_h <- withVisible(plot_ancombc3_training(decreasing, main = "constructed"))
  invisible(grDevices::dev.off())
  expect_true(!vis$visible && identical(vis$value, quiet$training_history) && identical(vis_h$value, decreasing),
              label = "plot_ancombc3_training() on a null device")
  err <- tryCatch(plot_ancombc3_training(list(a = 1)), error = function(e) conditionMessage(e))
  expect_true(is.character(err) && grepl("training_history", err), label = "plot_ancombc3_training() rejects other input")
})
