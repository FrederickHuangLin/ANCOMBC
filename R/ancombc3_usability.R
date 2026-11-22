# =============================================================================
# ANCOM-BC3 reporting: progress display, first-stage convergence check,
# diagnostics, and the methods of the class "ancombc3"
# =============================================================================
#
# None of the functions of this file changes a number of the analysis. The
# progress display (.ac3_progress_reporter()) is a callback of the stage
# functions; it uses no random numbers and leaves the state of the
# random-number stream unchanged. The convergence check
# (ancombc3_vae_convergence()) and the diagnostics (.ac3_main_diagnostics())
# summarize objects that the analysis has already computed.

# -----------------------------------------------------------------------------
# Progress display
# -----------------------------------------------------------------------------

# Progress reporter of one long step of ancombc3(): kind "epochs" (training of
# the ZIPLN-VAE, one update per epoch) or "taxon_fits" (the taxon-wise fits of
# the second stage; each call of fit_ancombc3_da_taxa() is one pass over the
# taxa and starts with an update at zero). The bar is drawn with the package
# cli when it is installed and with utils::txtProgressBar() otherwise. Returns
# a list with update, the callback of the stage function, and done, which
# closes the bar. The state of the random-number stream is restored after
# every call of the display functions.
.ac3_progress_reporter <- function(kind = c("epochs", "taxon_fits")) {
  kind <- match.arg(kind)
  use_cli <- requireNamespace("cli", quietly = TRUE)
  state <- new.env(parent = baseenv())
  state$pass <- 0L
  state$id <- NULL
  state$bar <- NULL
  keep_seed <- function(expr_fun) withr::with_preserve_seed(expr_fun())
  label <- function() {
    if (kind == "epochs") "Training the ZIPLN-VAE" else sprintf("Taxon-wise fits, pass %d", state$pass)
  }
  close_bar <- function() {
    if (use_cli && !is.null(state$id)) {
      cli::cli_progress_done(id = state$id, .envir = state)
    } else if (!is.null(state$bar)) {
      close(state$bar)
    }
    state$id <- NULL
    state$bar <- NULL
  }
  open_bar <- function(total) {
    if (use_cli) {
      state$id <- cli::cli_progress_bar(
        name = label(), total = total, clear = TRUE, .auto_close = FALSE, .envir = state,
        format = "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} {cli::pb_status}"
      )
    } else {
      message(label())
      state$bar <- utils::txtProgressBar(min = 0, max = total, style = 3)
    }
  }
  update <- function(done, total, ...) {
    keep_seed(function() {
      new_pass <- kind == "taxon_fits" && done == 0L
      if (new_pass) {
        close_bar()
        state$pass <- state$pass + 1L
      }
      if (is.null(state$id) && is.null(state$bar)) open_bar(total)
      extra <- list(...)
      status <- if (kind == "epochs" && length(extra) >= 2L) {
        sprintf("loss %.4f / %.4f", extra[[1L]], extra[[2L]])
      } else {
        ""
      }
      if (use_cli) {
        cli::cli_progress_update(set = done, total = total, status = status, id = state$id,
                                 .envir = state)
      } else {
        utils::setTxtProgressBar(state$bar, done)
      }
      if (kind == "epochs" && done >= total) close_bar()
    })
    invisible(NULL)
  }
  done <- function() {
    keep_seed(close_bar)
    invisible(NULL)
  }
  list(update = update, done = done)
}

# Progress reporter of ancombc3(): NULL (no display) unless verbose is TRUE
# and the session is interactive, so that the logs of batch jobs contain the
# step messages only.
.ac3_main_progress <- function(verbose, kind) {
  if (!isTRUE(verbose) || !interactive()) return(NULL)
  .ac3_progress_reporter(kind)
}

# -----------------------------------------------------------------------------
# First-stage convergence check
# -----------------------------------------------------------------------------

# Convergence check of the training of the ZIPLN-VAE from its training history
# (train_zipln_vae(): one row per epoch with epoch, train_loss, val_loss, and
# kl_weight). The validation loss of an epoch is the loss of the withheld
# samples after the epoch, evaluated with the Kullback-Leibler weight of the
# epoch, so that losses are comparable only among the epochs after the
# warm-up (kl_weight = 1). The check compares the mean validation loss of the
# last window epochs, m_last, with that of the window epochs before them,
# m_previous: the relative decrease (m_previous - m_last) / |m_previous|
# exceeding tolerance indicates that the validation loss was still decreasing.
# The check is made only when at least min_after_warmup epochs have
# kl_weight = 1 (the 2 window epochs then lie after the warm-up); otherwise
# checked is FALSE. Without withheld samples (fewer than 11 samples),
# val_loss is the mean training loss of the epoch, and the check uses it.
# Returns checked, still_decreasing, relative_decrease, mean_previous,
# mean_last, n_after_warmup, window, tolerance, and min_after_warmup.
ancombc3_vae_convergence <- function(history, window = 10L, tolerance = 0.01, min_after_warmup = 30L) {
  out <- list(checked = FALSE, still_decreasing = FALSE, relative_decrease = NA_real_,
              mean_previous = NA_real_, mean_last = NA_real_, n_after_warmup = 0L,
              window = as.integer(window), tolerance = tolerance,
              min_after_warmup = as.integer(min_after_warmup))
  if (is.null(history) || nrow(history) == 0L) return(out)
  history <- history[order(history$epoch), , drop = FALSE]
  after <- history$kl_weight >= 1
  out$n_after_warmup <- sum(after)
  if (out$n_after_warmup < min_after_warmup || out$n_after_warmup < 2L * window) return(out)
  v <- history$val_loss[after]
  k <- length(v)
  last <- v[(k - window + 1L):k]
  previous <- v[(k - 2L * window + 1L):(k - window)]
  if (any(!is.finite(c(last, previous)))) return(out)
  out$checked <- TRUE
  out$mean_previous <- mean(previous)
  out$mean_last <- mean(last)
  out$relative_decrease <- (out$mean_previous - out$mean_last) / abs(out$mean_previous)
  out$still_decreasing <- isTRUE(out$relative_decrease > tolerance)
  out
}

# Warning of ancombc3() when the convergence check (ancombc3_vae_convergence())
# finds that the validation loss was still decreasing. Returns TRUE when the
# warning is issued, invisibly.
.ac3_main_convergence_warning <- function(convergence) {
  if (!isTRUE(convergence$still_decreasing)) return(invisible(FALSE))
  warning(sprintf(paste0(
    "The validation loss of the ZIPLN-VAE was still decreasing at the end of training (mean over the ",
    "last %d epochs %.2f%% below that of the %d epochs before); a larger epochs may be appropriate."),
    convergence$window, 100 * convergence$relative_decrease, convergence$window), call. = FALSE)
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Training history plot
# -----------------------------------------------------------------------------

#' @title Training history of the first stage of ANCOM-BC3
#'
#' @description Draws the training and validation losses of the first stage
#' of \code{\link{ancombc3}} against the epoch. A vertical line marks the end
#' of the Kullback-Leibler warm-up.
#'
#' @param x the result of \code{\link{ancombc3}} or its element
#' \code{training_history}.
#' @param ... further arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return the training history, a \code{data.frame}, invisibly.
#'
#' @seealso \code{\link{ancombc3}}
#'
#' @examples
#' history <- data.frame(epoch = 1:60,
#'                       train_loss = 100 + 50 * exp(-(1:60) / 10),
#'                       val_loss = 105 + 50 * exp(-(1:60) / 10),
#'                       kl_weight = pmin(1, (1:60) / 40))
#' plot_ancombc3_training(history)
#'
#' @author Huang Lin
#'
#' @export
plot_ancombc3_training <- function(x, ...) {
  history <- if (inherits(x, "ancombc3") || (is.list(x) && !is.data.frame(x))) x$training_history else x
  if (!is.data.frame(history) || !all(c("epoch", "train_loss", "val_loss", "kl_weight") %in% names(history))) {
    stop("x must be the result of ancombc3() or its element training_history.")
  }
  losses <- c(history$train_loss, history$val_loss)
  args <- list(...)
  defaults <- list(x = history$epoch, y = history$train_loss, type = "l", lwd = 2, col = "#1f77b4",
                   ylim = range(losses[is.finite(losses)]), xlab = "Epoch", ylab = "Loss",
                   main = "ZIPLN-VAE training history")
  do.call(graphics::plot, utils::modifyList(defaults, args))
  graphics::lines(history$epoch, history$val_loss, lwd = 2, lty = 2, col = "#d62728")
  warmup_end <- history$epoch[which(history$kl_weight >= 1)[1L]]
  has_warmup <- length(warmup_end) == 1L && is.finite(warmup_end) && warmup_end > min(history$epoch)
  if (has_warmup) graphics::abline(v = warmup_end, lty = 3, col = "grey40")
  graphics::legend(
    "topright", bty = "n",
    legend = c("Training loss", "Validation loss", if (has_warmup) "End of KL warm-up"),
    lty = c(1, 2, if (has_warmup) 3), lwd = c(2, 2, if (has_warmup) 1),
    col = c("#1f77b4", "#d62728", if (has_warmup) "grey40")
  )
  invisible(history)
}

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

# Common-shift status of one tested column from the scored DA table of
# score_ancombc3_da(): the selected estimator, the selected estimate delta and
# its variance, the maximal-overlap and Gaussian-mixture estimates, the number
# of taxa in the null set of the maximal-overlap estimate, the number of
# eligible coefficients, the selection status, and the empirical-null
# adjustment factor. NULL when the column has no scored row.
.ac3_main_shift_status <- function(scored, column) {
  s <- scored[scored$term == column, , drop = FALSE]
  if (nrow(s) == 0L) return(NULL)
  first <- function(name) if (name %in% names(s)) s[[name]][1L] else NA
  gm <- isTRUE(first("gaussian_mixture_selected"))
  list(
    estimator = if (gm) "gaussian_mixture" else "maximal_overlap",
    delta = first("delta"),
    var_delta = first("var_delta"),
    maximal_overlap_delta = first("maximal_overlap_delta"),
    gaussian_mixture_delta = first("gaussian_mixture_delta"),
    n_null_set = first("maximal_overlap_n_null_set"),
    n_eligible_coefficients = first("n_eligible_coefficients"),
    maximal_overlap_status = first("maximal_overlap_status"),
    shift_selection_status = first("shift_selection_status"),
    adjustment_factor = first("adjustment_factor")
  )
}

# A taxon-wise fit record with a usable fit: a converged fit (ok) with a
# covariance matrix (vcov). A converged fit whose penalized observed
# information is not positive definite has no vcov, and its DA estimate is NA.
.ac3_usable_fit <- function(f) isTRUE(f$ok) && !is.null(f$vcov)

# Diagnostics of ancombc3(), assembled from the objects of the analysis:
# per tested column (by_covariate), the DA eligibility table (from
# ancombc3_da_eligibility(), or from ancombc3_group_eligibility() for a group
# column), the common-shift status, and the DP eligibility counts by reason;
# fit, the convergence summary of the second-stage fit; first_stage, the
# convergence check of the ZIPLN-VAE; and, with cluster, the number of
# subjects and the distribution of the number of samples per subject.
.ac3_main_diagnostics <- function(tested_cols, da_sets, by_term, scored, dp_z, fit, mg, convergence,
                                  subject = NULL, taxa) {
  by_covariate <- lapply(stats::setNames(tested_cols, tested_cols), function(col) {
    da_table <- NULL
    if (!is.null(mg) && col %in% mg$columns) {
      i <- match(col, mg$columns)
      el <- mg$eligibility
      sup <- el$support
      k <- match(taxa, rownames(el$da_pair))
      da_table <- data.frame(
        taxon = taxa, term = col, criterion = "multigroup",
        da_eligible = unname(el$da_pair[k, i]),
        common_criteria_pass = unname(el$da_common[k]),
        n_positive_reference = unname(sup$n_positive[k, 1L]),
        n_positive_level = unname(sup$n_positive[k, i + 1L]),
        at_risk_mass_reference = unname(sup$at_risk_mass[k, 1L]),
        at_risk_mass_level = unname(sup$at_risk_mass[k, i + 1L]),
        row.names = NULL, stringsAsFactors = FALSE
      )
    } else if (!is.null(da_sets$diagnostics)) {
      da_table <- da_sets$diagnostics[da_sets$diagnostics$term == col, , drop = FALSE]
      rownames(da_table) <- NULL
    }
    dp_rows <- dp_z[dp_z$test_term == col, , drop = FALSE]
    dp_reason <- if (nrow(dp_rows) > 0L) table(dp_rows$dp_eligibility_reason) else table(character(0))
    list(
      n_da_eligible = length(by_term[[col]]),
      da_eligibility = da_table,
      common_shift = .ac3_main_shift_status(scored, col),
      n_dp_eligible = sum(as.logical(dp_rows$dp_eligible) %in% TRUE),
      dp_eligibility = stats::setNames(as.integer(dp_reason), names(dp_reason))
    )
  })
  fit_ok <- vapply(fit$fits, .ac3_usable_fit, logical(1))
  out <- list(
    by_covariate = by_covariate,
    fit = list(
      n_taxa = length(fit$fits),
      n_fitting_set = length(fit$fitting_set),
      n_without_usable_fit = sum(!fit_ok),
      taxa_without_usable_fit = names(fit$fits)[!fit_ok],
      converged = fit$converged,
      start_selected = fit$start_selected,
      pass1_iterations = fit$pass1_iterations,
      pass2_iterations = fit$pass2_iterations,
      verification_passes = fit$verification_passes,
      verification_replaced = fit$verification_replaced,
      cross_taxon_penalty = !is.null(fit$cross_taxon_penalty)
    ),
    first_stage = convergence
  )
  if (!is.null(subject)) {
    size <- table(subject)
    out$cluster <- list(n_subjects = length(size), n_samples = length(subject),
                        samples_per_subject = table(as.integer(size), dnn = "samples per subject"))
  }
  out
}

# -----------------------------------------------------------------------------
# Methods of the class "ancombc3"
# -----------------------------------------------------------------------------

# Discoveries per tested covariate (and level of a group column) in res.
.ac3_discovery_counts <- function(res) {
  key <- if ("level" %in% names(res)) {
    ifelse(is.na(res$level), res$covariate, paste0(res$covariate, " (", res$level, ")"))
  } else {
    res$covariate
  }
  keys <- unique(key)
  data.frame(
    covariate = keys,
    da_eligible = vapply(keys, function(k) sum(res$da_eligible[key == k]), integer(1)),
    dp_eligible = vapply(keys, function(k) sum(res$dp_eligible[key == k]), integer(1)),
    da_discoveries = vapply(keys, function(k) sum(res$da_discovery[key == k]), integer(1)),
    dp_discoveries = vapply(keys, function(k) sum(res$dp_discovery[key == k]), integer(1)),
    omnibus_discoveries = vapply(keys, function(k) sum(res$omnibus_discovery[key == k]), integer(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#' @title Print and summary methods for ANCOM-BC3 results
#'
#' @description \code{print} shows the numbers of samples and taxa, the
#' tested covariates, and the DA, DP, and omnibus discoveries of each tested
#' covariate. \code{summary} adds the numbers of eligible taxa, the component
#' labels of the omnibus discoveries, the discoveries of the multigroup
#' tables, the convergence check of the first stage, and the computation time.
#'
#' @param x,object the result of \code{\link{ancombc3}}; for the
#' \code{print} method of the summary, the result of \code{summary}.
#' @param ... not used.
#'
#' @return \code{print} returns its argument invisibly. \code{summary}
#' returns an object of class \code{summary.ancombc3}, a \code{list}.
#'
#' @seealso \code{\link{ancombc3}}
#'
#' @examples
#' # See the examples of ancombc3().
#' methods(class = "ancombc3")
#'
#' @author Huang Lin
#'
#' @name ancombc3-methods
#' @export
print.ancombc3 <- function(x, ...) {
  st <- x$settings
  cat("ANCOM-BC3 analysis\n")
  cat(sprintf("  %d samples, %d taxa; alpha = %.2f; multiplicity adjustment %s\n",
              st$n_samples, st$n_taxa, st$alpha, st$p_adj_method %||% "BH"))
  if (length(st$tested) > 0L) cat("  Tested covariates: ", paste(st$tested, collapse = ", "), "\n", sep = "")
  if (!is.null(st$group)) {
    cat(sprintf("  Group: %s (levels %s; reference %s)\n", st$group$name,
                paste(st$group$levels, collapse = ", "), st$group$levels[1L]))
  }
  if (!is.null(st$cluster)) cat(sprintf("  Cluster: %s (%d subjects)\n", st$cluster, st$n_subjects))
  cat(sprintf("  Engine: %s\n", st$engine %||% "r"))
  cat("Discoveries (DA / DP / omnibus):\n")
  counts <- .ac3_discovery_counts(x$res)
  for (i in seq_len(nrow(counts))) {
    cat(sprintf("  %-30s %4d / %4d / %4d\n", counts$covariate[i], counts$da_discoveries[i],
                counts$dp_discoveries[i], counts$omnibus_discoveries[i]))
  }
  tabs <- intersect(c("res_global", "res_pair", "res_dunn", "res_trend"), names(x))
  cat("Tables: x$res", if (length(tabs) > 0L) paste0(", x$", tabs, collapse = ""), "\n", sep = "")
  cat("Use summary(x) for eligibility, component labels, and diagnostics.\n")
  invisible(x)
}

#' @rdname ancombc3-methods
#' @export
summary.ancombc3 <- function(object, ...) {
  res <- object$res
  counts <- .ac3_discovery_counts(res)
  key <- if ("level" %in% names(res)) {
    ifelse(is.na(res$level), res$covariate, paste0(res$covariate, " (", res$level, ")"))
  } else {
    res$covariate
  }
  label_levels <- c("DA_driven", "DP_driven", "DA_DP", "omnibus_only")
  labels <- do.call(rbind, lapply(counts$covariate, function(k) {
    tab <- table(factor(res$component_label[key == k & res$omnibus_discovery], levels = label_levels))
    as.data.frame(as.list(stats::setNames(as.integer(tab), label_levels)))
  }))
  counts <- cbind(counts, labels)
  multigroup <- list()
  if (!is.null(object$res_global)) {
    g <- object$res_global
    multigroup$global <- c(da_eligible = sum(g$da_eligible), da_discoveries = sum(g$da_discovery),
                           dp_eligible = sum(g$dp_eligible), dp_discoveries = sum(g$dp_discovery),
                           omnibus_discoveries = sum(g$omnibus_discovery))
  }
  comparisons <- function(tab) {
    cmp <- unique(tab$comparison)
    data.frame(
      comparison = cmp,
      da_discoveries = vapply(cmp, function(k) sum(tab$da_discovery[tab$comparison == k]), integer(1)),
      dp_discoveries = vapply(cmp, function(k) sum(tab$dp_discovery[tab$comparison == k]), integer(1)),
      omnibus_discoveries = vapply(cmp, function(k) sum(tab$omnibus_discovery[tab$comparison == k]),
                                   integer(1)),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  if (!is.null(object$res_pair)) multigroup$pairwise <- comparisons(object$res_pair)
  if (!is.null(object$res_dunn)) multigroup$dunnet <- comparisons(object$res_dunn)
  if (!is.null(object$res_trend)) {
    tr <- object$res_trend
    multigroup$trend <- list(
      da = table(factor(tr$da_pattern[tr$da_discovery])),
      dp = table(factor(tr$dp_pattern[tr$dp_discovery])),
      omnibus_discoveries = sum(tr$omnibus_discovery)
    )
  }
  st <- object$settings
  out <- list(
    n_samples = st$n_samples, n_taxa = st$n_taxa, alpha = st$alpha,
    p_adj_method = st$p_adj_method %||% "BH",
    tested = st$tested, group = st$group$name, cluster = st$cluster, n_subjects = st$n_subjects,
    engine = st$engine %||% "r",
    by_covariate = counts,
    multigroup = if (length(multigroup) > 0L) multigroup else NULL,
    first_stage = st$first_stage_convergence,
    elapsed_seconds = st$elapsed_seconds,
    stage_seconds = st$stage_seconds
  )
  class(out) <- "summary.ancombc3"
  out
}

#' @rdname ancombc3-methods
#' @export
print.summary.ancombc3 <- function(x, ...) {
  cat("ANCOM-BC3 analysis summary\n")
  cat(sprintf("  %d samples, %d taxa; alpha = %.2f; multiplicity adjustment %s\n",
              x$n_samples, x$n_taxa, x$alpha, x$p_adj_method))
  if (!is.null(x$group)) cat("  Group: ", x$group, "\n", sep = "")
  if (!is.null(x$cluster)) cat(sprintf("  Cluster: %s (%d subjects)\n", x$cluster, x$n_subjects))
  cat(sprintf("  Engine: %s\n", x$engine %||% "r"))
  cat("\nEligible taxa and discoveries by covariate (component labels of the omnibus discoveries):\n")
  print(x$by_covariate, row.names = FALSE)
  if (!is.null(x$multigroup)) {
    mgp <- x$multigroup
    if (!is.null(mgp$global)) {
      cat("\nGlobal test:\n")
      print(mgp$global)
    }
    if (!is.null(mgp$pairwise)) {
      cat("\nPairwise comparisons (mixed directional FDR):\n")
      print(mgp$pairwise, row.names = FALSE)
    }
    if (!is.null(mgp$dunnet)) {
      cat("\nComparisons with the reference (Dunnett type, mixed directional FDR):\n")
      print(mgp$dunnet, row.names = FALSE)
    }
    if (!is.null(mgp$trend)) {
      cat(sprintf("\nTrend test: %d omnibus discoveries; DA discoveries by pattern: %s; DP: %s\n",
                  mgp$trend$omnibus_discoveries,
                  if (length(mgp$trend$da)) paste(names(mgp$trend$da), mgp$trend$da, collapse = ", ") else "none",
                  if (length(mgp$trend$dp)) paste(names(mgp$trend$dp), mgp$trend$dp, collapse = ", ") else "none"))
    }
  }
  fs <- x$first_stage
  if (!is.null(fs)) {
    cat("\nFirst stage: ")
    if (isTRUE(fs$checked)) {
      cat(sprintf("relative decrease of the mean validation loss over the last %d epochs %.4f (tolerance %.2f)%s\n",
                  fs$window, fs$relative_decrease, fs$tolerance,
                  if (isTRUE(fs$still_decreasing)) "; still decreasing, a larger epochs may be appropriate" else ""))
    } else {
      cat(sprintf("convergence not checked (%d epochs after the KL warm-up; %d required)\n",
                  fs$n_after_warmup, fs$min_after_warmup))
    }
  }
  if (!is.null(x$stage_seconds)) {
    cat("\nComputation time (seconds): ",
        paste(sprintf("%s %.1f", names(x$stage_seconds), x$stage_seconds), collapse = ", "),
        sprintf("; total %.1f\n", x$elapsed_seconds), sep = "")
  } else if (!is.null(x$elapsed_seconds)) {
    cat(sprintf("\nComputation time: %.1f seconds\n", x$elapsed_seconds))
  }
  invisible(x)
}
