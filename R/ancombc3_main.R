# =============================================================================
# ANCOM-BC3 user-level function
# =============================================================================
#
# ancombc3() carries out, in order, the first stage
# (ancombc3_fitted_vae_offsets()), the DA eligibility
# (ancombc3_da_eligibility()), the second-stage DA fit (fit_ancombc3_da()) and
# its scoring (score_ancombc3_da()), the DP eligibility (dp_eligibility()) and
# the DP test (fit_ancombc3_dp_zero_pattern()), and the omnibus test
# (omnibus_pvalue(), ancombc3_component_labels()), with the settings of
# ancombc3_default_settings(). Every estimator input is observed: the counts,
# the covariates, and the library sizes N_i = sum_j y_ij.

#' @title Method settings of ANCOM-BC3
#'
#' @description Returns the method settings of \code{\link{ancombc3}}: the
#' ridge penalties, the convergence criteria of the second-stage fit, the
#' settings of the common-shift estimator, and the thresholds of the DA and DP
#' eligibility criteria. Individual values are replaced through the argument
#' \code{settings} of \code{ancombc3}.
#'
#' @return a named \code{list} of the method settings.
#'
#' @seealso \code{\link{ancombc3}}
#'
#' @examples
#' names(ancombc3_default_settings())
#'
#' @author Huang Lin
#'
#' @export
ancombc3_default_settings <- function() {
  list(
    ridge_direct = 1e-6,
    ridge_dp = 1e-3,
    d_tol = 0.01,
    d_max_iter = 50L,
    n_starts = 10L,
    maxit_taxon = 500L,
    d_bound = 10,
    da_maximal_overlap_min_count = 25L,
    da_maximal_overlap_bandwidths = seq(0.5, 3.0, by = 0.1),
    da_offset_variance = "cluster",
    min_observed_positive = 10,
    min_positive_per_group = 5,
    min_continuous_var_ratio = 0.25,
    min_observed_prevalence = 0.10,
    min_at_risk_eff_n = 10,
    min_at_risk_design_support = 1,
    dp_min_total_absence = 1,
    dp_min_total_at_risk = 5,
    dp_min_side_absence = 0.5,
    dp_min_side_at_risk = 2,
    dp_min_weighted_variance = 0.01
  )
}

# Right-hand side of a formula argument as a string: a one-sided formula, a
# string with or without "~", or a character vector of covariate names.
.ac3_main_rhs <- function(f, name) {
  if (inherits(f, "formula")) {
    if (length(f) != 2L) stop(name, " must be a one-sided formula, such as ~ age + group.")
    return(paste(deparse(f[[2L]], width.cutoff = 500L), collapse = " "))
  }
  if (is.character(f) && length(f) >= 1L && all(nzchar(f))) {
    if (length(f) > 1L) return(paste(f, collapse = " + "))
    return(trimws(sub("^\\s*~", "", f)))
  }
  stop(name, " must be a one-sided formula, a string, or a character vector of covariate names.")
}

.ac3_main_say <- function(verbose, ...) if (isTRUE(verbose)) message(...)

.ac3_main_list <- function(x, k = 5L) {
  out <- paste(utils::head(x, k), collapse = ", ")
  if (length(x) > k) out <- paste0(out, ", ... (", length(x), " in total)")
  out
}

# Multiplicity adjustments supported by p_adj_method.
.ac3_p_adjust_methods <- c("BH", "BY", "holm", "hochberg", "bonferroni", "none")

#' @title Analysis of Compositions of Microbiomes with Bias Correction 3
#' (ANCOM-BC3)
#'
#' @description Tests the association of each taxon with each tested covariate
#' through differential abundance (DA), differential prevalence (DP), and an
#' omnibus test that combines the two. The counts follow a zero-inflated
#' Poisson-log-normal model. DA is a difference in the abundance of a taxon
#' among the samples in which it is structurally present, after correction of
#' the sample-specific bias. DP is a difference in the probability of
#' structural presence, after the zeros expected from sampling are accounted
#' for. Multigroup comparisons (global test, pairwise comparisons, Dunnett's
#' type of test, and trend test) and repeated measures are supported.
#'
#' @details The analysis has two stages. The first stage fits a zero-inflated
#' Poisson-log-normal variational autoencoder to the count table with the
#' package \pkg{torch} and provides sample-by-taxon offsets. The second stage
#' fits a Poisson-log-normal model to each taxon, estimates the sample-specific
#' bias, and provides the DA and DP tests.
#'
#' The package \pkg{torch} and its LibTorch library are required and are not
#' installed with \pkg{ANCOMBC}. Install them with
#' \code{install.packages("torch")} followed by \code{torch::install_torch()}.
#'
#' The input must contain read counts. Relative abundances and base counts are
#' not valid input.
#'
#' Samples with a missing value of a model covariate, and samples or taxa
#' without a positive count, are dropped with a message. Fewer than ten
#' samples give an error, and fewer than 50 samples give a warning.
#'
#' A tested covariate is continuous (numeric with more than two distinct
#' values) or binary (numeric, logical, character, or factor with exactly two
#' distinct values). A binary covariate is coded 0 for its first level and 1
#' for its second level; the coding is returned in
#' \code{settings$covariate_coding}. A factor with three or more levels is
#' analyzed through \code{group}. Tested covariates are not allowed in
#' interaction terms.
#'
#' The DA, DP, and omnibus families are adjusted separately within each tested
#' covariate. No adjustment is made across covariates.
#'
#' Differences from \code{\link{ancombc2}}: the tested covariates are named in
#' \code{tested}; repeated measures are specified by \code{cluster} in place of
#' \code{rand_formula}; structural zeros are analyzed by the DP test in place
#' of \code{struc_zero}; no pseudo-count is used; and the taxonomic level is
#' not an argument, so the counts are aggregated before the analysis. The
#' defaults of \code{alpha} (0.10) and \code{p_adj_method} (\code{"BH"}) differ
#' from those of \code{ancombc2}.
#'
#' @param data the count table. A \code{matrix} or \code{data.frame} of
#' non-negative integer read counts, a \code{phyloseq} object, or a
#' \code{SummarizedExperiment} object (including a
#' \code{TreeSummarizedExperiment}). For a \code{matrix} or \code{data.frame},
#' the row and column names are the taxon and sample names.
#' @param taxa_are_rows logical. Whether the taxa are in the rows of a
#' \code{matrix} or \code{data.frame} \code{data}. Default is \code{TRUE}.
#' @param assay.type alias for \code{assay_name}.
#' @param assay_name character. Name of the count assay of a
#' \code{SummarizedExperiment} object. Default is \code{"counts"}.
#' @param meta_data a \code{data.frame} of sample metadata with one row per
#' sample, identified by its row names or by the column \code{sample_col}.
#' Ignored for \code{phyloseq} and \code{SummarizedExperiment} input, whose
#' sample data are used.
#' @param sample_col character. Optional name of the column of
#' \code{meta_data} that holds the sample names. Default is \code{NULL} (the
#' row names).
#' @param fix_formula the fixed effects of the abundance model: a one-sided
#' formula (\code{~ age + group}), a character string (\code{"age + group"}),
#' or a character vector of covariate names. The intercept is always included.
#' @param prv_formula the fixed effects of the prevalence model, in the same
#' forms. Default is \code{fix_formula}.
#' @param tested character vector. Names of the tested covariates, each a term
#' of both formulas. It may be omitted when \code{group} is given.
#' @param group character. Optional name of a factor or character covariate
#' with three or more levels for the multigroup analysis. Its first level is
#' the reference level. It is a term of both formulas and is not in
#' \code{tested}. Default is \code{NULL}.
#' @param cluster character. For repeated measures, the name of the column of
#' \code{meta_data} that identifies the subject of each sample. The standard
#' errors and tests then use cluster-robust covariances with the subjects as
#' clusters. At least 10 subjects are required. Default is \code{NULL}
#' (independent samples).
#' @param p_adj_method character. Multiplicity adjustment of the DA, DP, and
#' omnibus q-values: \code{"BH"} (default), \code{"BY"}, \code{"holm"},
#' \code{"hochberg"}, \code{"bonferroni"}, or \code{"none"}. The screening step
#' of the mixed directional FDR procedure of the pairwise and Dunnett's type of
#' comparisons uses the Benjamini-Hochberg procedure.
#' @param alpha numeric. Nominal false discovery rate. Default is 0.10.
#' @param n_cl numeric. Number of worker processes for the taxon-wise
#' second-stage fits. The results do not depend on it. Default is 1.
#' @param verbose logical. Whether to report the data checks and the progress
#' of the steps. Default is \code{TRUE}.
#' @param global logical. Whether to perform the global test of equality of the
#' levels of \code{group}. Default is \code{FALSE}.
#' @param pairwise logical. Whether to perform all pairwise comparisons of the
#' levels of \code{group}. Default is \code{FALSE}.
#' @param dunnet logical. Whether to perform the comparisons of every level of
#' \code{group} with the reference level (Dunnett's type of test). Default is
#' \code{FALSE}.
#' @param trend logical. Whether to perform the trend test over the ordered
#' levels of \code{group}. Default is \code{FALSE}.
#' @param mdfdr_control a named list with the element \code{fwer_ctrl_method},
#' the method of \code{\link[stats]{p.adjust}} applied within a taxon in the
#' mixed directional FDR procedure of the pairwise comparisons. Default is
#' \code{list(fwer_ctrl_method = "holm")}.
#' @param trend_control a named list for the trend test. \code{contrast} is a
#' named list of \eqn{(K - 1) \times (K - 1)} constraint matrices \eqn{A} for
#' the alternative \eqn{A b \ge 0}, where \eqn{K} is the number of levels of
#' \code{group}, or \code{NULL} for the monotone increasing and decreasing
#' patterns. \code{node} holds one position per matrix: the last position for a
#' monotone pattern and the position of the peak for an umbrella pattern.
#' \code{B} is the number of draws from the null distribution. Default is
#' \code{list(contrast = NULL, node = NULL, B = 1000)}.
#' @param epochs numeric. Number of training epochs of the first stage.
#' Default is 100.
#' @param seed numeric. Random seed of the first stage and of the null draws of
#' the trend test. The random-number state of the caller is not changed.
#' Default is 1.
#' @param device character. The \pkg{torch} device of the first stage,
#' \code{"cpu"} (default) or \code{"cuda"}.
#' @param engine character. Evaluation of the Poisson-log-normal quadrature of
#' the second stage: \code{"compiled"} (default) uses the compiled routine of
#' the package, and \code{"r"} uses R code that performs the same operations.
#' The results of the two agree up to rounding error.
#' @param settings a named list that replaces elements of
#' \code{\link{ancombc3_default_settings}}.
#'
#' @return an object of class \code{ancombc3}, a \code{list} with components:
#'         \itemize{
#'         \item{ \code{res}, a \code{data.frame} with one row per analyzed
#'         taxon and tested covariate (and per non-reference level of
#'         \code{group}):
#'         \itemize{
#'         \item{ \code{taxon}, \code{covariate}, \code{covariate_type}, and,
#'         with \code{group}, \code{level}.}
#'         \item{ \code{da_eligible}: DA eligibility.}
#'         \item{ \code{da_estimate}: bias-corrected coefficient of the
#'         covariate in the model of the log at-risk mean abundance (natural
#'         logarithm): the log fold change of level 1 versus level 0 for a
#'         binary covariate, and the change per unit for a continuous
#'         covariate.}
#'         \item{ \code{da_se}: standard error of \code{da_estimate}.}
#'         \item{ \code{da_statistic}, \code{da_df}: DA test statistic
#'         \code{da_estimate / da_se} and its degrees of freedom.}
#'         \item{ \code{da_p}, \code{da_q}, \code{da_discovery}: two-sided
#'         p-value from the t distribution, q-value over the DA-eligible taxa
#'         of the covariate, and \code{da_q <= alpha}.}
#'         \item{ \code{dp_eligible}: DP eligibility.}
#'         \item{ \code{dp_estimate}: log odds ratio of structural presence for
#'         the covariate; a positive value is a higher prevalence.}
#'         \item{ \code{dp_statistic}: DP chi-square statistic with one degree
#'         of freedom.}
#'         \item{ \code{dp_p}, \code{dp_q}, \code{dp_discovery}: p-value,
#'         q-value over the DP-eligible taxa of the covariate, and
#'         \code{dp_q <= alpha}.}
#'         \item{ \code{omnibus_p}: \code{min(1, 2 * min(da_p, dp_p))} for a
#'         taxon eligible for both tests; otherwise the p-value of its eligible
#'         test.}
#'         \item{ \code{omnibus_q}, \code{omnibus_discovery}: q-value over the
#'         taxa eligible for at least one test, and
#'         \code{omnibus_q <= alpha}.}
#'         \item{ \code{component_label}: \code{"DA_driven"},
#'         \code{"DP_driven"}, \code{"DA_DP"}, or \code{"omnibus_only"} for an
#'         omnibus discovery; \code{"not_discovered"} or
#'         \code{"not_in_universe"} (eligible for neither test) otherwise.}
#'         }
#'         The DA columns other than \code{da_eligible} are \code{NA} for a
#'         DA-ineligible taxon, and the DP columns other than
#'         \code{dp_eligible} are \code{NA} for a DP-ineligible taxon.}
#'         \item{ \code{res_global}, a \code{data.frame} with one row per
#'         taxon for the global test of \code{group}, when \code{global = TRUE}.}
#'         \item{ \code{res_pair}, a \code{data.frame} with one row per taxon
#'         and pairwise comparison, when \code{pairwise = TRUE}. The comparison
#'         \code{"<l>_vs_<k>"} has the estimate of level \code{l} minus level
#'         \code{k}.}
#'         \item{ \code{res_dunn}, a \code{data.frame} with one row per taxon
#'         and comparison with the reference level, when \code{dunnet = TRUE}.}
#'         \item{ \code{res_trend}, a \code{data.frame} with one row per taxon
#'         for the trend test, when \code{trend = TRUE}.}
#'         \item{ \code{sample_bias}, the estimated sample-specific bias, named
#'         by sample.}
#'         \item{ \code{samp_frac}, the total sample bias: \code{sample_bias}
#'         plus the contribution of the estimated common shifts of the tested
#'         covariates (natural-log scale).}
#'         \item{ \code{bias_correct_log_table}, a \code{matrix} of
#'         bias-corrected log abundances in the orientation of \code{data}:
#'         the log count minus \code{samp_frac} for a positive count, and the
#'         posterior mean log intensity of an at-risk taxon minus
#'         \code{samp_frac} for a zero.}
#'         \item{ \code{struc_absence_prob}, a \code{matrix} in the same
#'         orientation of the posterior probabilities of structural absence:
#'         0 for a positive count.}
#'         \item{ \code{sigma}, the estimated log-scale standard deviation of
#'         each taxon.}
#'         \item{ \code{training_history}, a \code{data.frame} with one row per
#'         epoch of the first stage: \code{epoch}, \code{train_loss},
#'         \code{val_loss}, and \code{kl_weight}.}
#'         \item{ \code{diagnostics}, a \code{list}: \code{by_covariate} (DA
#'         eligibility table, common-shift status, and DP eligibility counts
#'         of each tested covariate), \code{fit} (convergence summary of the
#'         second stage, with the taxa without a usable fit),
#'         \code{first_stage} (convergence check of the first stage), and, with
#'         \code{cluster}, the numbers of subjects and of samples per subject.}
#'         \item{ \code{settings}, a \code{list}: the call, the formulas, the
#'         coding of the tested covariates, the method settings, the engine,
#'         the dropped samples and taxa, the software versions, and the
#'         computation time.}
#'         }
#'
#' @seealso \code{\link{ancombc2}} \code{\link{ancombc3_default_settings}}
#' \code{\link{plot_ancombc3_training}}
#'
#' @examplesIf requireNamespace("torch", quietly = TRUE) && torch::torch_is_installed()
#' counts <- read.csv(system.file("extdata", "ancombc3_toy_counts.csv",
#'                                package = "ANCOMBC"),
#'                    row.names = 1, check.names = FALSE)
#' meta_data <- read.csv(system.file("extdata", "ancombc3_toy_metadata.csv",
#'                                   package = "ANCOMBC"),
#'                       row.names = 1)
#' meta_data$group <- factor(meta_data$group, levels = c("control", "treated"))
#'
#' out <- ancombc3(data = counts, meta_data = meta_data,
#'                 fix_formula = ~ group + bmi_z,
#'                 tested = c("group", "bmi_z"),
#'                 epochs = 50, verbose = FALSE)
#' out
#' head(out$res)
#'
#' @author Huang Lin
#'
#' @export
ancombc3 <- function(data,
                     taxa_are_rows = TRUE,
                     assay.type = assay_name,
                     assay_name = "counts",
                     meta_data = NULL,
                     sample_col = NULL,
                     fix_formula,
                     prv_formula = fix_formula,
                     tested,
                     group = NULL,
                     cluster = NULL,
                     p_adj_method = "BH",
                     alpha = 0.10,
                     n_cl = 1L,
                     verbose = TRUE,
                     global = FALSE,
                     pairwise = FALSE,
                     dunnet = FALSE,
                     trend = FALSE,
                     mdfdr_control = list(fwer_ctrl_method = "holm"),
                     trend_control = list(contrast = NULL, node = NULL, B = 1000L),
                     epochs = 100L,
                     seed = 1L,
                     device = "cpu",
                     engine = c("compiled", "r"),
                     settings = list()) {
  .ac3_require_torch()
  assay_name <- assay.type
  t_start <- proc.time()[["elapsed"]]
  call <- match.call()
  counts <- data
  metadata <- meta_data
  engine <- match.arg(engine)

  # --- Settings ---------------------------------------------------------------
  cfg <- ancombc3_default_settings()
  if (length(settings) > 0L) {
    nm_set <- names(settings)
    if (is.null(nm_set) || anyNA(nm_set) || any(!nzchar(nm_set))) {
      stop("Every element of settings must be named, such as settings = list(d_tol = 0.01). ",
           "See names(ancombc3_default_settings()).")
    }
    if (anyDuplicated(nm_set)) {
      stop("settings has duplicated names: ", paste(unique(nm_set[duplicated(nm_set)]), collapse = ", "), ".")
    }
  }
  unknown <- setdiff(names(settings), names(cfg))
  if (length(unknown) > 0L) {
    stop("Unknown element(s) of settings: ", paste(unknown, collapse = ", "),
         ". See names(ancombc3_default_settings()).")
  }
  cfg[names(settings)] <- settings
  if (!is.character(p_adj_method) || length(p_adj_method) != 1L || !p_adj_method %in% .ac3_p_adjust_methods) {
    stop("p_adj_method must be one of: ", paste(.ac3_p_adjust_methods, collapse = ", "), ".")
  }
  epochs <- suppressWarnings(as.integer(epochs))
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1).")
  }
  if (length(epochs) != 1L || is.na(epochs) || epochs < 1L) stop("epochs must be a positive integer.")
  n_cl_int <- suppressWarnings(as.integer(n_cl))
  if (length(n_cl_int) != 1L || is.na(n_cl_int) || n_cl_int < 1L) stop("n_cl must be a positive integer.")
  # settings$da_offset_variance is one of a fixed set of values, matched
  # exactly.
  offset_variance <- cfg$da_offset_variance
  offset_variance_choices <- c("cluster", "always", "never")
  if (!is.character(offset_variance) || length(offset_variance) != 1L ||
      !offset_variance %in% offset_variance_choices) {
    stop("settings$da_offset_variance must be one of: ",
         paste0("\"", offset_variance_choices, "\"", collapse = ", "), ".", call. = FALSE)
  }
  if (missing(fix_formula)) stop("fix_formula is required.")
  mg <- ancombc3_group_arguments(group, global, pairwise, dunnet, trend, trend_control, mdfdr_control)
  if (!is.null(mg) && (missing(tested) || length(tested) == 0L)) tested <- character(0)
  if (missing(tested) || !is.character(tested) || (length(tested) == 0L && is.null(mg))) {
    stop("tested must name at least one covariate of meta_data.")
  }
  if (anyDuplicated(tested)) stop("tested contains duplicated covariate names.")
  x_rhs <- .ac3_main_rhs(fix_formula, "fix_formula")
  r_rhs <- .ac3_main_rhs(prv_formula, "prv_formula")

  # --- Count table and metadata -----------------------------------------------
  is_container <- inherits(counts, "phyloseq") || inherits(counts, "SummarizedExperiment")
  if (!is_container) {
    if (!is.matrix(counts) && !is.data.frame(counts)) {
      stop("data must be a matrix, a data frame, a phyloseq object, or a SummarizedExperiment.")
    }
    if (is.data.frame(counts)) {
      non_numeric <- names(counts)[!vapply(counts, is.numeric, logical(1))]
      if (length(non_numeric) > 0L) {
        stop("data has non-numeric columns (", .ac3_main_list(non_numeric),
             "). Remove taxonomy or other annotation columns and keep the counts only.")
      }
      counts <- as.matrix(counts)
    }
    if (!is.numeric(counts)) stop("data must be numeric.")
    otu <- if (isTRUE(taxa_are_rows)) t(counts) else counts
    if (is.null(rownames(otu))) {
      stop("data has no sample names (", if (isTRUE(taxa_are_rows)) "column" else "row",
           " names). Sample names are required to match the counts with meta_data.")
    }
    if (is.null(colnames(otu))) colnames(otu) <- paste0("taxon_", seq_len(ncol(otu)))
    if (is.null(metadata)) stop("meta_data is required for a matrix or data frame count table.")
    metadata <- as.data.frame(metadata)
    if (!is.null(sample_col)) {
      if (!sample_col %in% names(metadata)) stop("sample_col '", sample_col, "' is not a column of meta_data.")
      ids <- as.character(metadata[[sample_col]])
      if (anyDuplicated(ids)) stop("meta_data column '", sample_col, "' has duplicated sample names.")
      rownames(metadata) <- ids
    }
    missing_meta <- setdiff(rownames(otu), rownames(metadata))
    if (length(missing_meta) > 0L) {
      stop(length(missing_meta), " sample(s) of data have no row in meta_data: ",
           .ac3_main_list(missing_meta),
           ". The sample names of data must match the row names of meta_data (or sample_col).")
    }
    extra_meta <- setdiff(rownames(metadata), rownames(otu))
    if (length(extra_meta) > 0L) {
      .ac3_main_say(verbose, "Dropped ", length(extra_meta),
                    " sample(s) of meta_data that are not in data.")
    }
    metadata <- metadata[rownames(otu), , drop = FALSE]
    input <- ancombc3_normalize_input(otu, meta_data = metadata, taxa_are_rows = FALSE)
  } else {
    input <- ancombc3_normalize_input(counts, taxa_are_rows = taxa_are_rows, assay_name = assay_name)
  }
  y <- input$y
  metadata <- input$metadata
  if (anyDuplicated(rownames(y))) stop("data has duplicated sample names.")
  if (anyNA(colnames(y)) || any(!nzchar(colnames(y)))) stop("data has missing or empty taxon names.")
  if (anyDuplicated(colnames(y))) stop("data has duplicated taxon names.")
  if (any(abs(y - round(y)) > 1e-8)) {
    stop("data must contain non-negative integers (counts); relative abundances or normalized values are not accepted.")
  }
  y <- round(y)

  # --- Covariates ---------------------------------------------------------------
  x_vars <- all.vars(stats::as.formula(paste("~", x_rhs)))
  r_vars <- all.vars(stats::as.formula(paste("~", r_rhs)))
  model_vars <- unique(c(x_vars, r_vars))
  absent <- setdiff(model_vars, names(metadata))
  if (length(absent) > 0L) stop("Covariate(s) not found in meta_data: ", paste(absent, collapse = ", "))
  x_terms <- attr(stats::terms(stats::as.formula(paste("~", x_rhs))), "term.labels")
  r_terms <- attr(stats::terms(stats::as.formula(paste("~", r_rhs))), "term.labels")
  for (tn in tested) {
    if (!tn %in% names(metadata)) stop("Tested covariate '", tn, "' is not a column of meta_data.")
    if (!tn %in% x_terms) stop("Tested covariate '", tn, "' is not a term of fix_formula.")
    if (!tn %in% r_terms) stop("Tested covariate '", tn, "' is not a term of prv_formula.")
    interactions <- c(x_terms, r_terms)[grepl(":", c(x_terms, r_terms), fixed = TRUE)]
    in_interaction <- vapply(strsplit(interactions, ":", fixed = TRUE), function(v) tn %in% v, logical(1))
    if (any(in_interaction)) {
      stop("Tested covariate '", tn, "' appears in an interaction term; this is not supported.")
    }
  }
  if (!is.null(mg)) .ac3_mg_check_term(group, tested, metadata, x_terms, r_terms)
  incomplete <- !stats::complete.cases(metadata[, model_vars, drop = FALSE])
  if (any(incomplete)) {
    .ac3_main_say(verbose, "Dropped ", sum(incomplete), " sample(s) with missing covariate values: ",
                  .ac3_main_list(rownames(y)[incomplete]), ".")
  }
  dropped_samples <- rownames(y)[incomplete]
  y <- y[!incomplete, , drop = FALSE]
  metadata <- metadata[!incomplete, , drop = FALSE]

  zero_samples <- rowSums(y) == 0
  if (any(zero_samples)) {
    .ac3_main_say(verbose, "Dropped ", sum(zero_samples), " sample(s) with no positive count: ",
                  .ac3_main_list(rownames(y)[zero_samples]), ".")
    dropped_samples <- c(dropped_samples, rownames(y)[zero_samples])
    y <- y[!zero_samples, , drop = FALSE]
    metadata <- metadata[!zero_samples, , drop = FALSE]
  }
  zero_taxa <- colSums(y) == 0
  dropped_taxa <- colnames(y)[zero_taxa]
  if (any(zero_taxa)) {
    .ac3_main_say(verbose, "Dropped ", sum(zero_taxa), " taxon/taxa with no positive count: ",
                  .ac3_main_list(dropped_taxa), ".")
    y <- y[, !zero_taxa, drop = FALSE]
  }
  if (ncol(y) < 2L) stop("Fewer than two taxa have a positive count.")

  # Tested covariates: two distinct values give a 0/1 indicator, more than two
  # numeric values a continuous covariate.
  coding <- list()
  for (tn in tested) {
    v <- metadata[[tn]]
    lev <- if (is.factor(v)) levels(droplevels(v)) else sort(unique(v))
    if (length(lev) < 2L) stop("Tested covariate '", tn, "' is constant over the analyzed samples.")
    if (length(lev) == 2L) {
      lev_chr <- as.character(lev)
      metadata[[tn]] <- as.numeric(as.character(v) == lev_chr[2L])
      coding[[tn]] <- sprintf("binary: 0 = %s, 1 = %s", lev_chr[1L], lev_chr[2L])
    } else if (is.numeric(v)) {
      coding[[tn]] <- "continuous: as supplied"
    } else {
      stop("Tested covariate '", tn, "' has ", length(lev), " levels. Tested covariates must be ",
           "binary or continuous; analyze a factor with three or more levels with the argument group.")
    }
  }

  if (!is.null(mg)) {
    mg <- ancombc3_group_setup(mg, metadata, min_samples_per_level = cfg$min_positive_per_group)
    metadata <- mg$metadata
    mg$metadata <- NULL
    coding[[group]] <- mg$coding
  }

  # --- Sample-size checks -------------------------------------------------------
  n <- nrow(y)
  if (n < max(cfg$min_observed_positive, 2 * cfg$min_positive_per_group)) {
    stop("Too few samples (", n, "): the DA eligibility requires at least ",
         max(cfg$min_observed_positive, 2 * cfg$min_positive_per_group), ".")
  }
  for (tn in tested) {
    if (identical(substr(coding[[tn]], 1L, 6L), "binary")) {
      sizes <- table(factor(metadata[[tn]], levels = c(0, 1)))
      if (min(sizes) < cfg$min_positive_per_group) {
        stop("Tested covariate '", tn, "' has a group with ", min(sizes), " sample(s); the DA ",
             "eligibility requires at least ", cfg$min_positive_per_group, " positive counts per group.")
      }
    }
  }
  if (n < 50L) {
    warning("The analysis has ", n, " samples; the operating characteristics of ANCOM-BC3 were ",
            "evaluated with 50 or more samples.", call. = FALSE)
  }
  # Repeated measures: the subject of each analyzed sample, and the group sizes
  # in subjects of the tested covariates that are constant within subject.
  subject <- NULL
  if (!is.null(cluster)) {
    subject <- .ac3_cluster_subjects(metadata, cluster, verbose)
    .ac3_cluster_check_groups(metadata, tested, coding, subject, cfg$min_positive_per_group)
    if (!is.null(mg)) .ac3_mg_check_subjects(mg, metadata, subject, cfg$min_positive_per_group)
  }

  # --- Designs ------------------------------------------------------------------
  design <- ancombc3_make_design(metadata = metadata, x_formula = x_rhs, test_terms = c(tested, group),
                                 r_formula = r_rhs)
  if (!is.null(mg)) .ac3_mg_check_design(mg, design)
  for (nm in c("x", "r")) {
    mat <- cbind(1, design[[nm]])
    if (qr(mat)$rank < ncol(mat)) {
      stop("The ", if (nm == "x") "abundance" else "prevalence",
           " design matrix is rank deficient; remove collinear covariates.")
    }
    if (n <= ncol(mat) + 1L) stop("Too few samples for the number of covariates.")
  }
  tested_cols <- names(design$tests)

  # The taxa are identified by their names throughout: the first-stage
  # outputs carry the dimnames of y, and every later step matches taxa by name.
  storage.mode(y) <- "double"
  .ac3_main_say(verbose, sprintf(
    "Analyzing %d samples and %d taxa (dropped %d samples and %d taxa); tested: %s; engine: %s.",
    n, ncol(y), length(dropped_samples), length(dropped_taxa), paste(c(tested, group), collapse = ", "),
    engine))
  stage_seconds <- c(first_stage = NA_real_, da_eligibility = NA_real_, second_stage = NA_real_,
                     dp_test = NA_real_, multigroup = NA_real_)
  since <- function(t0) proc.time()[["elapsed"]] - t0

  # --- First stage ------------------------------------------------------------
  .ac3_main_say(verbose, sprintf("First stage: training the ZIPLN-VAE (%d epochs).", epochs))
  t_stage <- proc.time()[["elapsed"]]
  pb <- .ac3_main_progress(verbose, "epochs")
  offsets <- ancombc3_fitted_vae_offsets(y = y, x = design$x, r = design$r, tested_terms = tested_cols,
                                         epochs = epochs, device = device, seed = seed,
                                         cluster = subject, progress = if (is.null(pb)) NULL else pb$update)
  if (!is.null(pb)) pb$done()
  stage_seconds[["first_stage"]] <- since(t_stage)
  history <- offsets$training_history
  convergence <- ancombc3_vae_convergence(history)
  .ac3_main_say(verbose, sprintf(
    "First stage completed in %.1f s: %d epochs; final training loss %.4f, validation loss %.4f (%d withheld samples).",
    stage_seconds[["first_stage"]], nrow(history), history$train_loss[nrow(history)],
    history$val_loss[nrow(history)], offsets$n_validation))
  .ac3_main_convergence_warning(convergence)

  # --- DA -----------------------------------------------------------------------
  t_stage <- proc.time()[["elapsed"]]
  da_sets <- ancombc3_da_eligibility(
    at_risk_prob = offsets$at_risk_prob, y = y, x = design$x,
    test_terms = design$tests[tested_cols],
    min_observed_positive = cfg$min_observed_positive,
    min_positive_per_group = cfg$min_positive_per_group,
    min_continuous_var_ratio = cfg$min_continuous_var_ratio,
    min_observed_prevalence = cfg$min_observed_prevalence,
    min_at_risk_eff_n = cfg$min_at_risk_eff_n,
    min_at_risk_design_support = cfg$min_at_risk_design_support,
    cluster = subject
  )
  if (!is.null(mg)) {
    # Group columns: the eligible taxa of a column are those of the comparison
    # of its level with the reference level (ancombc3_group_eligibility()).
    mg$eligibility <- ancombc3_group_eligibility(
      group = metadata[[group]], y = y, x = design$x, r = design$r, columns = mg$columns,
      at_risk_prob = offsets$at_risk_prob,
      structural_absence_weight = offsets$structural_absence_weight,
      min_positive_per_group = cfg$min_positive_per_group,
      min_observed_prevalence = cfg$min_observed_prevalence,
      min_at_risk_eff_n = cfg$min_at_risk_eff_n,
      min_at_risk_design_support = cfg$min_at_risk_design_support,
      dp_min_total_absence = cfg$dp_min_total_absence, dp_min_total_at_risk = cfg$dp_min_total_at_risk,
      dp_min_side_absence = cfg$dp_min_side_absence, dp_min_side_at_risk = cfg$dp_min_side_at_risk,
      cluster = subject
    )
    da_sets$by_term[mg$columns] <- mg$eligibility$da_by_column
    da_sets$fitting_set <- unique(unlist(da_sets$by_term, use.names = FALSE))
  }
  by_term <- da_sets$by_term
  stage_seconds[["da_eligibility"]] <- since(t_stage)
  .ac3_main_say(verbose, sprintf(
    "DA eligibility: %s of %d taxa.", paste(sprintf("%s %d", tested_cols, lengths(by_term[tested_cols])),
                                            collapse = ", "), ncol(y)))
  if (length(da_sets$fitting_set) == 0L) {
    stop("No taxon is DA-eligible for any tested covariate; the sample-specific bias cannot be estimated.")
  }
  .ac3_main_say(verbose, sprintf(
    "Second stage: taxon-wise fits of %d taxa and the sample-specific bias (%d taxa in the fitting set).",
    ncol(y), length(da_sets$fitting_set)))
  t_stage <- proc.time()[["elapsed"]]
  pb <- .ac3_main_progress(verbose, "taxon_fits")
  fit <- fit_ancombc3_da(
    y = y, x = design$x, r = design$r, sigma_vae = offsets$sigma_vae,
    abundance_offset = offsets$abundance_offset, structural_offset = offsets$structural_offset,
    sampling_zero_prob = offsets$sampling_zero_prob, at_risk_prob = offsets$at_risk_prob,
    test_terms = tested_cols, fitting_set = da_sets$fitting_set, eligible_taxa_by_term = by_term,
    ridge = cfg$ridge_direct, maxit_taxon = cfg$maxit_taxon, d_tol = cfg$d_tol,
    d_max_iter = cfg$d_max_iter, n_starts = cfg$n_starts, d_bound = cfg$d_bound,
    n_cores = n_cl,
    cluster = subject, progress = if (is.null(pb)) NULL else pb$update, engine = engine
  )
  if (!is.null(pb)) pb$done()
  # Variance term of the latent term carried by the abundance offset
  # (ancombc3_add_offset_covariance()). settings$da_offset_variance is
  # "cluster" (the term is included when cluster is given, with the subjects
  # as units), "always" (also without cluster, with the samples as units), or
  # "never".
  if (identical(offset_variance, "always") ||
      (identical(offset_variance, "cluster") && !is.null(subject))) {
    fit <- ancombc3_add_offset_covariance(fit, offsets$abundance_offset, cluster = subject)
  }
  scored <- score_ancombc3_da(
    fit$raw, variance_inputs = fit$da_variance_inputs,
    p_adj_method = p_adj_method, eligible_taxa_by_term = by_term,
    maximal_overlap_min_count = cfg$da_maximal_overlap_min_count,
    maximal_overlap_bandwidths = cfg$da_maximal_overlap_bandwidths
  )
  stage_seconds[["second_stage"]] <- since(t_stage)
  # A usable fit is a converged fit with a covariance matrix
  # (.ac3_usable_fit()).
  n_unusable <- sum(!vapply(fit$fits, .ac3_usable_fit, logical(1)))
  .ac3_main_say(verbose, sprintf(
    paste0("Second stage completed in %.1f s: %d taxa in the fitting set, %d taxa without a usable fit; ",
           "%d first-pass and %d final iterations, %d verification passes (%d fits replaced)%s."),
    stage_seconds[["second_stage"]], length(fit$fitting_set), n_unusable, as.integer(fit$pass1_iterations),
    as.integer(fit$pass2_iterations), as.integer(fit$verification_passes),
    as.integer(fit$verification_replaced), if (isTRUE(fit$converged)) "" else "; the alternation did not converge"))
  if (n_unusable > 0L) {
    warning(n_unusable, " of ", length(fit$fits), " taxa have no usable second-stage fit (no converged fit, ",
            "or a converged fit without a positive definite information matrix); their DA estimates ",
            "are NA. The taxa are listed in diagnostics$fit$taxa_without_usable_fit.", call. = FALSE)
  }
  for (tn in tested_cols) {
    sh <- .ac3_main_shift_status(scored, tn)
    if (is.null(sh)) next
    # The null set of the selected estimator (scored$in_null_set): that of the
    # maximal-overlap estimate, or the taxa with a posterior null-component
    # probability of at least 0.5 under the Gaussian-mixture estimate.
    n_null <- sum(scored$in_null_set[scored$term == tn] %in% TRUE)
    .ac3_main_say(verbose, sprintf("Common shift of %s: %.4f (%s estimator; %d taxa in the null set).",
                                   tn, sh$delta, gsub("_", "-", sh$estimator, fixed = TRUE), n_null))
  }

  # --- DP -----------------------------------------------------------------------
  .ac3_main_say(verbose, "DP test.")
  t_stage <- proc.time()[["elapsed"]]
  dp_elig <- dp_eligibility(
    structural_absence_weight = offsets$structural_absence_weight, r = design$r,
    test_terms = tested_cols,
    min_total_absence = cfg$dp_min_total_absence, min_total_at_risk = cfg$dp_min_total_at_risk,
    min_side_absence = cfg$dp_min_side_absence, min_side_at_risk = cfg$dp_min_side_at_risk,
    min_weighted_variance = cfg$dp_min_weighted_variance,
    cluster = subject
  )
  if (!is.null(mg)) dp_elig <- .ac3_mg_dp_eligibility_table(dp_elig, mg$eligibility, mg$columns)
  dp_z <- fit_ancombc3_dp_zero_pattern(
    y = y, r = design$r, f0 = offsets$sampling_zero_prob, test_terms = tested_cols,
    alpha = alpha, p_adj_method = p_adj_method, ridge = cfg$ridge_dp, eligibility = dp_elig,
    da_fit = fit, x = design$x,
    joint_terms = if (is.null(mg)) NULL else mg$columns,
    cluster = subject, engine = engine
  )
  stage_seconds[["dp_test"]] <- since(t_stage)
  .ac3_main_say(verbose, sprintf(
    "DP test completed in %.1f s: DP-eligible taxa %s.", stage_seconds[["dp_test"]],
    paste(vapply(tested_cols, function(tn) {
      sprintf("%s %d", tn, sum(as.logical(dp_z$dp_eligible[dp_z$test_term == tn]) %in% TRUE))
    }, character(1)), collapse = ", ")))

  # --- Results ------------------------------------------------------------------
  taxa <- colnames(y)
  blank <- function(fill) stats::setNames(rep(fill, length(taxa)), taxa)
  rows <- lapply(tested_cols, function(tn) {
    da_eligible <- stats::setNames(taxa %in% by_term[[tn]], taxa)
    est <- se <- stat <- df_da <- p_da <- q_da <- blank(NA_real_)
    da_t <- scored[scored$term == tn, , drop = FALSE]
    if (nrow(da_t) > 0L) {
      est[da_t$taxon] <- da_t$beta_corrected
      se[da_t$taxon] <- da_t$se_scaled
      stat[da_t$taxon] <- da_t$statistic
      df_da[da_t$taxon] <- da_t$df_t
      p_da[da_t$taxon] <- da_t$p_value
    }
    has_da <- da_eligible & is.finite(p_da)
    if (any(has_da)) q_da[has_da] <- stats::p.adjust(p_da[has_da], method = p_adj_method)

    dp_t <- dp_z[dp_z$test_term == tn, , drop = FALSE]
    dp_est <- dp_stat <- p_dp <- q_dp <- blank(NA_real_)
    dp_eligible <- blank(FALSE)
    if (nrow(dp_t) > 0L) {
      dp_est[dp_t$taxon] <- dp_t$estimate
      dp_stat[dp_t$taxon] <- dp_t$statistic
      p_dp[dp_t$taxon] <- dp_t$p_value
      q_dp[dp_t$taxon] <- dp_t$q_value
      dp_eligible[dp_t$taxon] <- as.logical(dp_t$dp_eligible) %in% TRUE
    }
    has_dp <- dp_eligible & is.finite(p_dp)

    p_omni <- omnibus_pvalue(p_da, p_dp, has_da, has_dp)
    q_omni <- blank(NA_real_)
    ok <- (has_da | has_dp) & is.finite(p_omni)
    if (any(ok)) q_omni[ok] <- stats::p.adjust(p_omni[ok], method = p_adj_method)
    omni_disc <- ok & q_omni <= alpha
    omni_disc[is.na(omni_disc)] <- FALSE
    label <- ancombc3_component_labels(discovered = omni_disc, has_da = has_da, has_dp = has_dp,
                                       p_da = p_da, p_dp = p_dp, alpha = alpha,
                                       p_adj_method = p_adj_method)
    na_unless <- function(v, keep) { v[!keep] <- NA; unname(v) }
    data.frame(
      taxon = taxa,
      covariate = design$tests[[tn]]$parent_term,
      covariate_type = design$tests[[tn]]$type,
      da_eligible = unname(da_eligible),
      da_estimate = na_unless(est, has_da),
      da_se = na_unless(se, has_da),
      da_statistic = na_unless(stat, has_da),
      da_df = na_unless(df_da, has_da),
      da_p = na_unless(p_da, has_da),
      da_q = na_unless(q_da, has_da),
      da_discovery = unname(has_da & q_da <= alpha) %in% TRUE,
      dp_eligible = unname(dp_eligible),
      dp_estimate = na_unless(dp_est, has_dp),
      dp_statistic = na_unless(dp_stat, has_dp),
      dp_p = na_unless(p_dp, has_dp),
      dp_q = na_unless(q_dp, has_dp),
      dp_discovery = unname(has_dp & q_dp <= alpha) %in% TRUE,
      omnibus_p = unname(p_omni),
      omnibus_q = unname(q_omni),
      omnibus_discovery = unname(omni_disc),
      component_label = unname(label),
      stringsAsFactors = FALSE
    )
  })
  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  multigroup <- list()
  if (!is.null(mg)) {
    level <- rep(c(rep(NA_character_, length(tested_cols) - length(mg$columns)), mg$levels[-1L]),
                 each = length(taxa))
    results <- cbind(results[1:2], level = level, results[-(1:2)], stringsAsFactors = FALSE)
    if (any(mg$global, mg$pairwise, mg$dunnet, mg$trend)) {
      .ac3_main_say(verbose, "Multigroup analysis of ", group, ".")
      t_stage <- proc.time()[["elapsed"]]
      multigroup <- ancombc3_multigroup(
        y = y, x = design$x, r = design$r, mg = mg, f0 = offsets$sampling_zero_prob, fit = fit,
        scored = scored, dp_z = dp_z, cfg = cfg, alpha = alpha, p_adj_method = p_adj_method,
        seed = seed, cluster = subject, engine = engine
      )
      stage_seconds[["multigroup"]] <- since(t_stage)
      .ac3_main_say(verbose, sprintf("Multigroup analysis completed in %.1f s: %s.", stage_seconds[["multigroup"]],
                                     paste(names(multigroup), collapse = ", ")))
    }
  }

  sample_bias <- stats::setNames(as.numeric(fit$d_hat), rownames(y))
  sigma <- stats::setNames(rep(NA_real_, length(taxa)), taxa)
  sig <- fit$sigma_hat
  sigma[names(sig)] <- as.numeric(sig)

  for (tn in tested) {
    .ac3_main_say(verbose, sprintf(
      "%s: %d DA-eligible, %d DP-eligible taxa; discoveries at q <= %.2f: DA %d, DP %d, omnibus %d.",
      tn, sum(results$da_eligible[results$covariate == tn]),
      sum(results$dp_eligible[results$covariate == tn]), alpha,
      sum(results$da_discovery[results$covariate == tn]),
      sum(results$dp_discovery[results$covariate == tn]),
      sum(results$omnibus_discovery[results$covariate == tn])))
  }
  if (!is.null(mg)) .ac3_mg_say(verbose, results, mg, alpha)

  taxa_in_rows <- if (inherits(counts, "phyloseq")) {
    isTRUE(phyloseq::taxa_are_rows(counts))
  } else if (inherits(counts, "SummarizedExperiment")) {
    TRUE
  } else {
    isTRUE(taxa_are_rows)
  }
  bc <- .ac3_bias_corrected_abundance(y, fit$d_hat, design$x, tested_cols, scored, taxa_in_rows,
                                      mu = fit$mu, pi = fit$pi, f0 = fit$f0, sigma = fit$sigma_hat,
                                      engine = engine)

  diagnostics <- .ac3_main_diagnostics(tested_cols, da_sets, by_term, scored, dp_z, fit, mg, convergence,
                                       subject = subject, taxa = taxa)
  elapsed <- proc.time()[["elapsed"]] - t_start
  .ac3_main_say(verbose, sprintf("Total time %.1f s.", elapsed))
  out <- list(
    res = results,
    sample_bias = sample_bias,
    samp_frac = bc$total_sample_bias,
    bias_correct_log_table = bc$bias_corrected_log_abundance,
    struc_absence_prob = bc$structural_absence_probability,
    sigma = sigma,
    training_history = history,
    diagnostics = diagnostics,
    settings = list(
      call = call,
      fix_formula = paste("~", x_rhs),
      prv_formula = paste("~", r_rhs),
      tested = tested,
      covariate_coding = coding,
      alpha = alpha,
      p_adj_method = p_adj_method,
      epochs = epochs,
      seed = seed,
      device = device,
      n_cl = as.integer(n_cl),
      engine = engine,
      method = cfg,
      n_samples = n,
      n_taxa = ncol(y),
      dropped_samples = dropped_samples,
      dropped_taxa = dropped_taxa,
      library_size = "N_i = sum_j y_ij (observed library size)",
      versions = list(
        R = R.version.string,
        torch = as.character(utils::packageVersion("torch")),
        platform = R.version$platform
      ),
      first_stage_convergence = convergence,
      elapsed_seconds = elapsed,
      stage_seconds = stage_seconds[is.finite(stage_seconds)]
    )
  )
  if (!is.null(mg)) {
    out$settings$group <- mg[c("name", "levels", "columns", "global", "pairwise", "dunnet", "trend",
                               "trend_control", "mdfdr_control")]
    out <- c(out, multigroup)
  }
  if (!is.null(cluster)) {
    out$settings$cluster <- cluster
    out$settings$n_subjects <- length(unique(subject))
  }
  class(out) <- c("ancombc3", "list")
  out
}
