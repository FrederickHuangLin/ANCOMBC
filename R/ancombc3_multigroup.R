# =============================================================================
# ANCOM-BC3 multigroup analysis
# =============================================================================
#
# Analysis of a group factor with K >= 3 levels (argument group of ancombc3()).
# The K - 1 treatment-coded columns of the factor (level k versus the reference
# level, the first level) are tested columns of both designs, so that the
# first stage, the second-stage fit, and the primary DA and DP tests treat each
# of them as a tested covariate. This file adds, for the DA and the DP channel:
#   the eligibility of each comparison of two levels from the support of the
#     levels (ancombc3_group_level_support(), ancombc3_group_eligibility());
#   the comparisons of two non-reference levels (.ac3_mg_da_contrast() for DA,
#     the one-degree-of-freedom DP test on the design with another reference
#     level for DP);
#   the joint working covariance of the K - 1 group coefficients of a taxon
#     (.ac3_mg_da_joint_covariance() for DA, .ac3_mg_dp_joint_test() for DP);
#   the global test of equality of the K levels;
#   the pairwise and the Dunnett-type comparisons with control of the mixed
#     directional false discovery rate (.ac3_mg_mdfdr_pairwise(),
#     .ac3_mg_mdfdr_dunnett());
#   the trend test (ancombc3_trend.R).
# ancombc3_multigroup() assembles the result tables res_global, res_pair,
# res_dunn, and res_trend.
#
# Notation. Levels are indexed 1, ..., K in factor order, level 1 is the
# reference, and b_k (k = 2, ..., K) is the coefficient of level k relative to
# the reference. A comparison (l, k) with k < l has the estimate "level l minus
# level k"; the comparisons are ordered (2, 1), ..., (K, 1), (3, 2), ...,
# (K, K - 1), so that the first K - 1 are the comparisons with the reference.

# -----------------------------------------------------------------------------
# Arguments and group factor
# -----------------------------------------------------------------------------

# Checks of the multigroup arguments of ancombc3(). Returns NULL when group is
# NULL (the requests global, pairwise, dunnet, and trend then give an error),
# and otherwise a list with the group name, the four requests, and the
# completed lists trend_control (contrast, node, B) and mdfdr_control
# (fwer_ctrl_method).
ancombc3_group_arguments <- function(group, global = FALSE, pairwise = FALSE, dunnet = FALSE,
                                     trend = FALSE, trend_control = list(),
                                     mdfdr_control = list()) {
  requests <- list(global = global, pairwise = pairwise, dunnet = dunnet, trend = trend)
  for (nm in names(requests)) {
    if (!is.logical(requests[[nm]]) || length(requests[[nm]]) != 1L || is.na(requests[[nm]])) {
      stop(nm, " must be TRUE or FALSE.")
    }
  }
  if (is.null(group)) {
    if (any(unlist(requests))) {
      stop("global, pairwise, dunnet, and trend require the argument group.")
    }
    return(NULL)
  }
  if (!is.character(group) || length(group) != 1L || is.na(group) || !nzchar(group)) {
    stop("group must be the name of one column of meta_data.")
  }
  tc <- list(contrast = NULL, node = NULL, B = 1000L)
  if (!is.list(trend_control)) stop("trend_control must be a list.")
  unknown <- setdiff(names(trend_control), names(tc))
  if (length(unknown) > 0L) {
    stop("Unknown element(s) of trend_control: ", paste(unknown, collapse = ", "), ".")
  }
  tc[names(trend_control)] <- trend_control
  B <- tc$B
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B != round(B) || B < 1) {
    stop("trend_control$B must be a positive integer.")
  }
  tc$B <- as.integer(B)
  mc <- list(fwer_ctrl_method = "holm")
  if (!is.list(mdfdr_control)) stop("mdfdr_control must be a list.")
  unknown <- setdiff(names(mdfdr_control), names(mc))
  if (length(unknown) > 0L) {
    stop("Unknown element(s) of mdfdr_control: ", paste(unknown, collapse = ", "), ".")
  }
  mc[names(mdfdr_control)] <- mdfdr_control
  if (!is.character(mc$fwer_ctrl_method) || length(mc$fwer_ctrl_method) != 1L ||
      !mc$fwer_ctrl_method %in% stats::p.adjust.methods) {
    stop("mdfdr_control$fwer_ctrl_method must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), ".")
  }
  c(list(name = group), requests, list(trend_control = tc, mdfdr_control = mc))
}

# Checks of the group term: a column of metadata, a term of both formulas, not
# a tested covariate, and not part of an interaction.
.ac3_mg_check_term <- function(group, tested, metadata, x_terms, r_terms) {
  if (!group %in% names(metadata)) stop("group '", group, "' is not a column of meta_data.")
  if (group %in% tested) {
    stop("group '", group, "' is also in tested; a covariate is analyzed either as a tested ",
         "covariate or as the group factor.")
  }
  if (!group %in% x_terms) stop("group '", group, "' is not a term of fix_formula.")
  if (!group %in% r_terms) stop("group '", group, "' is not a term of prv_formula.")
  interactions <- c(x_terms, r_terms)[grepl(":", c(x_terms, r_terms), fixed = TRUE)]
  in_interaction <- vapply(strsplit(interactions, ":", fixed = TRUE), function(v) group %in% v, logical(1))
  if (any(in_interaction)) {
    stop("group '", group, "' appears in an interaction term; this is not supported.")
  }
  invisible(TRUE)
}

# Group factor of the analyzed samples. The column is a factor or a character
# vector (converted to a factor); unused levels are dropped, an ordered factor
# is treated as unordered with the same level order, and treatment contrasts
# with the first level as reference are attached, so that model.matrix()
# returns the columns <group><level> of the non-reference levels. At least
# three levels are required, each with at least min_samples_per_level samples.
# mg is the list of ancombc3_group_arguments(); the result adds levels, columns
# (the design columns of the non-reference levels), reference_column (the name
# of the indicator column of the reference level), coding (a description),
# metadata with the prepared factor, and, when the trend test is requested,
# patterns (ancombc3_trend_patterns(), which checks trend_control against the
# number of levels).
ancombc3_group_setup <- function(mg, metadata, min_samples_per_level = 5) {
  group <- mg$name
  v <- metadata[[group]]
  if (!is.factor(v) && !is.character(v)) {
    stop("group '", group, "' must be a factor or a character column of meta_data.")
  }
  f <- droplevels(if (is.factor(v)) factor(v, levels = levels(v), ordered = FALSE) else factor(v))
  lev <- levels(f)
  if (length(lev) < 3L) {
    stop("group '", group, "' has ", length(lev), " level(s) among the analyzed samples; the ",
         "multigroup analysis requires at least three. Use tested for a binary covariate.")
  }
  sizes <- table(f)
  if (min(sizes) < min_samples_per_level) {
    stop("group '", group, "' has a level (", names(sizes)[which.min(sizes)], ") with ", min(sizes),
         " sample(s); each level requires at least ", min_samples_per_level, ".")
  }
  if (isTRUE(mg$trend)) {
    mg$patterns <- ancombc3_trend_patterns(lev, mg$trend_control$contrast, mg$trend_control$node)
  }
  stats::contrasts(f) <- stats::contr.treatment(lev)
  metadata[[group]] <- f
  mg$levels <- lev
  mg$columns <- paste0(group, lev[-1L])
  mg$reference_column <- paste0(group, lev[1L])
  mg$coding <- sprintf("group: reference = %s; columns = %s", lev[1L], paste(mg$columns, collapse = ", "))
  mg$metadata <- metadata
  mg
}

# Check that the group columns are tested columns of both designs.
.ac3_mg_check_design <- function(mg, design) {
  absent <- c(setdiff(mg$columns, colnames(design$x)), setdiff(mg$columns, colnames(design$r)),
              setdiff(mg$columns, names(design$tests)))
  if (length(absent) > 0L) {
    stop("The design has no column for the level(s) of group '", mg$name, "': ",
         paste(unique(absent), collapse = ", "), ".")
  }
  invisible(TRUE)
}

# Level sizes of the group factor in subjects (repeated measures): when the
# group factor is constant within every subject (subject, the subject of each
# analyzed sample), each level must hold at least min_subjects_per_level
# subjects.
.ac3_mg_check_subjects <- function(mg, metadata, subject, min_subjects_per_level) {
  ids <- .ac3_cluster_index(subject)
  g <- as.integer(metadata[[mg$name]])
  if (!.ac3_cluster_constant_within(g, ids)) return(invisible(TRUE))
  first <- match(seq_len(max(ids)), ids)
  sizes <- table(factor(mg$levels[g[first]], levels = mg$levels))
  if (min(sizes) < min_subjects_per_level) {
    stop("group '", mg$name, "' has a level (", names(sizes)[which.min(sizes)], ") with ", min(sizes),
         " subject(s); each level requires at least ", min_subjects_per_level, " subjects.")
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Level support and eligibility of comparisons
# -----------------------------------------------------------------------------

# Support of each level of the group factor for each taxon.
#   DA: level k is supported when the number of positive counts among its
#     samples is at least min_positive_per_group and its at-risk mass, the sum
#     of the first-stage at-risk probabilities (restricted to [0, 1]) over its
#     samples, is at least min_at_risk_design_support.
#   DP: level k is supported when its structural-absence mass, the sum of the
#     structural-absence weights s_ij^el (restricted to [0, 1]) over its
#     samples, is at least dp_min_side_absence and its at-risk mass, the sum of
#     1 - s_ij^el, is at least dp_min_side_at_risk.
# These are the per-group quantities of ancombc3_da_eligibility_term() (binary
# covariate) and of dp_eligibility() (sides of a binary covariate), evaluated
# for one level. cluster is an optional vector identifying the subject of each
# sample; when supplied, the number of positive counts of a level is the
# number of distinct subjects with at least one positive count among the
# samples of the level. When, in addition, the group factor is constant within
# every subject, the masses of a level are computed over subjects, as for a
# binary tested covariate (ancombc3_cluster.R, item 5): each mass is the sum,
# over the subjects of the level, of the largest weight over the samples of
# the subject. Returns the levels and the taxa x levels matrices da
# and dp (logical), n_positive, at_risk_mass, absence_mass, and
# dp_at_risk_mass.
ancombc3_group_level_support <- function(group, y, at_risk_prob, structural_absence_weight,
                                         min_positive_per_group = 5,
                                         min_at_risk_design_support = 1,
                                         dp_min_side_absence = 0.5,
                                         dp_min_side_at_risk = 2,
                                         cluster = NULL) {
  group <- droplevels(as.factor(group))
  lev <- levels(group)
  y <- as.matrix(y)
  taxa <- .ac3_taxon_names(y)
  colnames(y) <- taxa
  if (length(group) != nrow(y)) stop("group must have one value per sample.")
  at_risk <- .ac3_match_taxa(at_risk_prob, taxa, "at_risk_prob")
  absence <- .ac3_match_taxa(structural_absence_weight, taxa, "structural_absence_weight")
  if (!identical(dim(at_risk), dim(y)) || !identical(dim(absence), dim(y))) {
    stop("at_risk_prob and structural_absence_weight must have the n x p dimensions of y.")
  }
  at_risk <- pmin(pmax(at_risk, 0), 1)
  absence <- pmin(pmax(absence, 0), 1)
  indicator <- outer(as.integer(group), seq_along(lev), "==") * 1
  positive <- (y > 0) * 1
  if (is.null(cluster)) {
    n_positive <- crossprod(positive, indicator)
  } else {
    if (length(cluster) != nrow(y)) stop("cluster must have one value per sample.")
    n_positive <- vapply(seq_along(lev), function(k) {
      rows <- which(as.integer(group) == k)
      by_subject <- rowsum(positive[rows, , drop = FALSE], as.character(cluster[rows]), reorder = FALSE)
      colSums(by_subject > 0)
    }, numeric(ncol(y)))
  }
  if (!is.null(cluster) && .ac3_cluster_constant_within(as.integer(group), .ac3_cluster_index(cluster))) {
    level_mass <- function(w) {
      vapply(seq_along(lev), function(k) {
        rows <- which(as.integer(group) == k)
        colSums(.ac3_cluster_subject_max(w[rows, , drop = FALSE], .ac3_cluster_index(cluster[rows])))
      }, numeric(ncol(y)))
    }
    at_risk_mass <- level_mass(at_risk)
    absence_mass <- level_mass(absence)
    dp_at_risk_mass <- level_mass(1 - absence)
  } else {
    at_risk_mass <- crossprod(at_risk, indicator)
    absence_mass <- crossprod(absence, indicator)
    dp_at_risk_mass <- crossprod(1 - absence, indicator)
  }
  label <- function(mat) {
    mat <- matrix(mat, nrow = length(taxa), ncol = length(lev))
    dimnames(mat) <- list(taxa, lev)
    mat
  }
  list(
    levels = lev,
    da = label(n_positive >= min_positive_per_group & at_risk_mass >= min_at_risk_design_support),
    dp = label(absence_mass >= dp_min_side_absence & dp_at_risk_mass >= dp_min_side_at_risk),
    n_positive = label(n_positive),
    at_risk_mass = label(at_risk_mass),
    absence_mass = label(absence_mass),
    dp_at_risk_mass = label(dp_at_risk_mass)
  )
}

# Comparisons of two levels of a factor with the given levels, in the order
# (2, 1), ..., (K, 1), (3, 2), ..., (K, K - 1): ref and num are the level
# indices, level_ref and level_num the level names, and comparison is
# "<level_num>_vs_<level_ref>".
.ac3_mg_pairs <- function(levels) {
  K <- length(levels)
  ref <- rep(seq_len(K - 1L), times = (K - 1L):1L)
  num <- unlist(lapply(seq_len(K - 1L), function(k) (k + 1L):K))
  data.frame(ref = ref, num = num, level_ref = levels[ref], level_num = levels[num],
             comparison = paste0(levels[num], "_vs_", levels[ref]), stringsAsFactors = FALSE)
}

# Eligibility of every comparison of two levels and of the global and trend
# tests. A comparison of levels k and l is eligible for a taxon when both
# levels are supported (ancombc3_group_level_support()) and the taxon satisfies
# the criteria of the existing eligibility functions that do not depend on the
# split of the samples into groups:
#   DA (ancombc3_da_eligibility_term()): the at-risk effective sample size is
#     at least min_at_risk_eff_n and the positive-count fraction is at least
#     min_observed_prevalence;
#   DP (dp_eligibility()): the expected structural-zero count is at least
#     dp_min_total_absence and the expected at-risk count is at least
#     dp_min_total_at_risk.
# The global and trend tests require these criteria and the support of all K
# levels. x and r are the abundance and prevalence designs and columns the
# design columns of the non-reference levels, in level order. Returns levels;
# pairs (.ac3_mg_pairs()); the taxa x comparisons logical matrices da_pair and
# dp_pair; the logical vectors da_global, dp_global, da_common, and dp_common
# (the criteria that do not depend on the groups); support; and da_by_column,
# the names of the DA-eligible taxa of each group column, which is the
# comparison of its level with the reference.
ancombc3_group_eligibility <- function(group, y, x, r, columns, at_risk_prob,
                                       structural_absence_weight,
                                       min_positive_per_group = 5,
                                       min_observed_prevalence = 0.10,
                                       min_at_risk_eff_n = 10,
                                       min_at_risk_design_support = 1,
                                       dp_min_total_absence = 1,
                                       dp_min_total_at_risk = 5,
                                       dp_min_side_absence = 0.5,
                                       dp_min_side_at_risk = 2,
                                       cluster = NULL) {
  y <- as.matrix(y)
  taxa <- .ac3_taxon_names(y)
  colnames(y) <- taxa
  absence <- .ac3_match_taxa(structural_absence_weight, taxa, "structural_absence_weight")
  support <- ancombc3_group_level_support(
    group = group, y = y, at_risk_prob = at_risk_prob, structural_absence_weight = absence,
    min_positive_per_group = min_positive_per_group,
    min_at_risk_design_support = min_at_risk_design_support,
    dp_min_side_absence = dp_min_side_absence, dp_min_side_at_risk = dp_min_side_at_risk,
    cluster = cluster
  )
  lev <- support$levels
  K <- length(lev)
  if (length(columns) != K - 1L) stop("columns must name the design columns of the K - 1 non-reference levels.")
  da_term <- ancombc3_da_eligibility_term(
    at_risk_prob = at_risk_prob, y = y, x = x,
    term = list(term = columns[1L], type = "binary", coef_name = columns[1L], columns = columns[1L]),
    min_at_risk_eff_n = min_at_risk_eff_n, min_at_risk_design_support = min_at_risk_design_support,
    min_positive_per_group = min_positive_per_group, min_observed_prevalence = min_observed_prevalence
  )
  da_term <- da_term[match(taxa, da_term$taxon), , drop = FALSE]
  da_common <- da_term$at_risk_eff_n >= min_at_risk_eff_n &
    da_term$observed_prevalence >= min_observed_prevalence
  dp_term <- dp_eligibility(
    structural_absence_weight = absence, r = r, test_terms = columns[1L],
    min_total_absence = dp_min_total_absence, min_total_at_risk = dp_min_total_at_risk,
    min_side_absence = dp_min_side_absence, min_side_at_risk = dp_min_side_at_risk
  )
  dp_term <- dp_term[match(taxa, dp_term$taxon), , drop = FALSE]
  dp_common <- dp_term$absence_mass >= dp_min_total_absence &
    dp_term$at_risk_mass >= dp_min_total_at_risk
  pairs <- .ac3_mg_pairs(lev)
  pair_matrix <- function(common, level_support) {
    mat <- common & level_support[, pairs$ref, drop = FALSE] & level_support[, pairs$num, drop = FALSE]
    dimnames(mat) <- list(taxa, pairs$comparison)
    mat
  }
  da_pair <- pair_matrix(da_common, support$da)
  dp_pair <- pair_matrix(dp_common, support$dp)
  reference <- seq_len(K - 1L)
  list(
    levels = lev,
    pairs = pairs,
    da_pair = da_pair,
    dp_pair = dp_pair,
    da_global = stats::setNames(da_common & rowSums(support$da) == K, taxa),
    dp_global = stats::setNames(dp_common & rowSums(support$dp) == K, taxa),
    da_common = stats::setNames(da_common, taxa),
    dp_common = stats::setNames(dp_common, taxa),
    support = support,
    da_by_column = stats::setNames(lapply(reference, function(i) taxa[da_pair[, i]]), columns)
  )
}

# DP eligibility table (dp_eligibility()) with the rows of the group columns
# set to the eligibility of the comparison of the level with the reference:
# dp_eligible, the reason ("eligible", or the failed criteria: the total
# criteria of dp_eligibility() and "low_level_support"), and the smaller of the
# two level masses in min_side_absence_mass and min_side_at_risk_mass.
.ac3_mg_dp_eligibility_table <- function(dp_elig, eligibility, columns) {
  taxa <- rownames(eligibility$dp_pair)
  sup <- eligibility$support
  for (i in seq_along(columns)) {
    rows <- which(dp_elig$test_term == columns[i])
    k <- match(dp_elig$taxon[rows], taxa)
    eligible <- eligibility$dp_pair[k, i]
    total <- vapply(strsplit(dp_elig$dp_eligibility_reason[rows], ";", fixed = TRUE), function(v) {
      paste(v[v %in% c("low_total_absence", "low_total_at_risk")], collapse = ";")
    }, character(1))
    level_ok <- sup$dp[k, 1L] & sup$dp[k, i + 1L]
    reason <- ifelse(level_ok, total, ifelse(nzchar(total), paste(total, "low_level_support", sep = ";"),
                                             "low_level_support"))
    dp_elig$dp_eligible[rows] <- unname(eligible)
    dp_elig$dp_eligibility_reason[rows] <- ifelse(eligible, "eligible", reason)
    dp_elig$min_side_absence_mass[rows] <- pmin(sup$absence_mass[k, 1L], sup$absence_mass[k, i + 1L])
    dp_elig$min_side_at_risk_mass[rows] <- pmin(sup$dp_at_risk_mass[k, 1L], sup$dp_at_risk_mass[k, i + 1L])
  }
  dp_elig
}

# -----------------------------------------------------------------------------
# DA channel
# -----------------------------------------------------------------------------

# Raw rows of the comparison of two non-reference levels, the contrast
# "coefficient of column col_num minus coefficient of column col_ref" of the
# second-stage fit (fit_ancombc3_da()), for every taxon. For taxon j with
# contrast vector c over the abundance coefficients (the coefficient of column
# t of x is at position 1 + t of omega_j, after the intercept):
#   beta_raw = beta_j,num - beta_j,ref;
#   var_raw_model = c' V_j c, with V_j = fit$fits[[j]]$vcov;
#   sensitivity to the sample-specific bias g_j = g_j,num - g_j,ref, the
#     difference of the rows of fit$da_variance_inputs$G (g_j is linear in c);
#   var_raw_bias = max(g_j' V_d^con g_j, 0) and
#     var_raw_marginal = var_raw_model + var_raw_bias, the marginal working
#     variance of ancombc3_add_bias_covariance() for the contrast;
#   when fit$da_variance_inputs holds the offset term A
#     (ancombc3_add_offset_covariance()), the offset term of the contrast,
#     a_j = a_j,num - a_j,ref, the difference of the columns of A for the two
#     group columns (column j of A is linear in g_j), with a zero column for a
#     taxon whose g_j is not finite; var_raw_offset = ||a_j||^2 is added to
#     var_raw_marginal.
# Every other column is copied from the row of the taxon for col_num, so that
# the degrees of freedom of the t reference are those of the taxon (including
# the column df_count with repeated measures). On a fit with cluster, vcov,
# G, and V are the cluster-robust matrices of fit_ancombc3_da(), so the
# contrast variances are cluster-robust; the model-based columns of the
# repeated-measures analysis (var_raw_model_based, se_raw_model_based,
# var_raw_bias_model_based) are set to NA for the contrast. Returns
# raw (one row per taxon, with term = term), G (taxa x samples), and A (units
# x taxa, NULL without the offset term), the inputs of score_ancombc3_da() for
# the derived term.
.ac3_mg_da_contrast <- function(fit, x_columns, col_num, col_ref, term) {
  taxa <- names(fit$fits)
  raw <- fit$raw
  num <- raw[raw$term == col_num, , drop = FALSE]
  num <- num[match(taxa, num$taxon), , drop = FALSE]
  ref <- raw[raw$term == col_ref, , drop = FALSE]
  ref <- ref[match(taxa, ref$taxon), , drop = FALSE]
  idx <- 1L + match(c(col_num, col_ref), x_columns)
  if (anyNA(idx)) stop("Group columns not found in the abundance design.")
  var_model <- vapply(fit$fits, function(f) {
    v <- f$vcov
    if (!isTRUE(f$ok) || is.null(v) || max(idx) > nrow(v)) return(NA_real_)
    val <- v[idx[1L], idx[1L]] + v[idx[2L], idx[2L]] - 2 * v[idx[1L], idx[2L]]
    if (is.finite(val) && val > 0) val else NA_real_
  }, numeric(1))
  G <- fit$da_variance_inputs$G[[col_num]] - fit$da_variance_inputs$G[[col_ref]]
  bias <- rowSums((G %*% fit$da_variance_inputs$V) * G)
  out <- num
  out$term <- term
  out$beta_raw <- num$beta_raw - ref$beta_raw
  out$var_raw_model <- unname(var_model)
  out$se_raw_model <- sqrt(out$var_raw_model)
  out$var_raw_bias <- NA_real_
  out$var_raw_marginal <- NA_real_
  out$se_raw_marginal <- NA_real_
  model_based <- intersect(c("var_raw_model_based", "se_raw_model_based", "var_raw_bias_model_based"),
                           names(out))
  out[model_based] <- NA_real_
  ok <- is.finite(bias) & is.finite(var_model)
  out$var_raw_bias[ok] <- pmax(bias[ok], 0)
  out$var_raw_marginal[ok] <- out$var_raw_model[ok] + out$var_raw_bias[ok]
  A <- NULL
  A_all <- fit$da_variance_inputs$A
  if (!is.null(A_all) && !is.null(A_all[[col_num]]) && !is.null(A_all[[col_ref]])) {
    A <- A_all[[col_num]] - A_all[[col_ref]]
    A[!is.finite(A)] <- 0
    A[, !(rowSums(!is.finite(G)) == 0)] <- 0
    v_off <- colSums(A^2)[match(taxa, fit$da_variance_inputs$taxa)]
    out$var_raw_offset <- NA_real_
    out$var_raw_offset[ok] <- v_off[ok]
    out$var_raw_marginal[ok] <- out$var_raw_marginal[ok] + v_off[ok]
  }
  out$se_raw_marginal[ok] <- sqrt(out$var_raw_marginal[ok])
  rownames(out) <- NULL
  list(raw = out, G = G, A = A)
}

# score_ancombc3_da() with the scoring settings of the primary analysis
# (cfg, the list of ancombc3_default_settings(), and the multiplicity
# adjustment p_adj_method).
.ac3_mg_score <- function(raw, variance_inputs, eligible_taxa_by_term, cfg, p_adj_method) {
  score_ancombc3_da(
    raw, variance_inputs = variance_inputs,
    p_adj_method = p_adj_method, eligible_taxa_by_term = eligible_taxa_by_term,
    maximal_overlap_min_count = cfg$da_maximal_overlap_min_count,
    maximal_overlap_bandwidths = cfg$da_maximal_overlap_bandwidths
  )
}

# Joint working covariance of the corrected group coefficients of every taxon,
# the matrix form of .ac3_da_corrected_variance(). For the group columns k and
# l, with nu^k the correction weights of column k (scored$correction_weight,
# zero for the taxa that do not contribute to the common-shift estimate of the
# column), V_m[k, l] the model covariance of the two coefficients of taxon m
# (fit$fits[[m]]$vcov), g_mk the sensitivity of the coefficient of column k of
# taxon m to the sample-specific bias (row m of fit$da_variance_inputs$G[[k]]),
# gbar_k = sum_m nu_m^k g_mk, and V_d^con = fit$da_variance_inputs$V,
#   Sigma_j[k, l] = (1 - nu_j^k) (1 - nu_j^l) V_j[k, l]
#                   + sum_{m != j} nu_m^k nu_m^l V_m[k, l]
#                   + (g_jk - gbar_k)' V_d^con (g_jl - gbar_l)
#                   + (a_jk - abar_k)' (a_jl - abar_l),
# where the last term is present when fit$da_variance_inputs holds the offset
# term A (ancombc3_add_offset_covariance()), with a_mk column m of
# fit$da_variance_inputs$A[[k]] and abar_k = sum_m nu_m^k a_mk; the covariance
# of the offset terms of the coefficients of taxon m and column k and of taxon
# m' and column l is a_mk' a_m'l. The diagonal element [k, k], with the bias
# term bounded below by zero, is the working variance V_j,primary of
# score_ancombc3_da() for column k. Returns Sigma
# (taxa x columns x columns, NA where the fit of a taxon has no covariance
# matrix).
.ac3_mg_da_joint_covariance <- function(fit, scored, columns, x_columns) {
  taxa <- names(fit$fits)
  p <- length(taxa)
  m <- length(columns)
  idx <- 1L + match(columns, x_columns)
  if (anyNA(idx)) stop("Group columns not found in the abundance design.")
  V_model <- array(NA_real_, dim = c(p, m, m))
  for (j in seq_len(p)) {
    f <- fit$fits[[j]]
    if (isTRUE(f$ok) && !is.null(f$vcov) && max(idx) <= nrow(f$vcov)) V_model[j, , ] <- f$vcov[idx, idx]
  }
  nu <- matrix(0, p, m)
  for (k in seq_len(m)) {
    s <- scored[scored$term == columns[k], , drop = FALSE]
    if (nrow(s) == 0L || is.null(s$correction_weight)) next
    w <- s$correction_weight
    w[!is.finite(w)] <- 0
    nu[match(s$taxon, taxa), k] <- w
  }
  V_d <- fit$da_variance_inputs$V
  centered <- lapply(seq_len(m), function(k) {
    G <- fit$da_variance_inputs$G[[columns[k]]]
    use <- nu[, k] != 0
    gbar <- if (any(use)) drop(crossprod(G[use, , drop = FALSE], nu[use, k])) else numeric(ncol(G))
    sweep(G, 2L, gbar, "-")
  })
  centered_V <- lapply(centered, function(a) a %*% V_d)
  A_all <- fit$da_variance_inputs$A
  centered_A <- NULL
  if (!is.null(A_all) && all(vapply(columns, function(cn) !is.null(A_all[[cn]]), logical(1)))) {
    centered_A <- lapply(seq_len(m), function(k) {
      A <- A_all[[columns[k]]][, match(taxa, fit$da_variance_inputs$taxa), drop = FALSE]
      A[!is.finite(A)] <- 0
      A - drop(A %*% nu[, k])
    })
  }
  Sigma <- array(NA_real_, dim = c(p, m, m), dimnames = list(taxa, columns, columns))
  for (k in seq_len(m)) {
    for (l in k:m) {
      use <- nu[, k] != 0 & nu[, l] != 0
      shared <- sum(nu[use, k] * nu[use, l] * V_model[use, k, l])
      own <- ifelse(use, nu[, k] * nu[, l] * V_model[, k, l], 0)
      model <- (1 - nu[, k]) * (1 - nu[, l]) * V_model[, k, l] + (shared - own)
      bias <- rowSums(centered_V[[k]] * centered[[l]])
      if (k == l) bias <- pmax(bias, 0)
      if (!is.null(centered_A)) bias <- bias + colSums(centered_A[[k]] * centered_A[[l]])
      Sigma[, k, l] <- model + bias
      Sigma[, l, k] <- Sigma[, k, l]
    }
  }
  list(Sigma = Sigma)
}

# Correlation matrix of a covariance matrix, made positive definite: the
# matrix is symmetrized and converted to a correlation matrix R, and when the
# smallest eigenvalue lambda of R is below min_eigen, R is replaced by
# (1 - t) R + t I with t = (min_eigen - lambda) / (1 - lambda), the smallest
# convex combination with the identity whose smallest eigenvalue is min_eigen.
# NULL when an element is not finite or a variance is not positive.
.ac3_mg_correlation <- function(Sigma, min_eigen = 1e-6) {
  if (is.null(Sigma)) return(NULL)
  Sigma <- as.matrix(Sigma)
  if (any(!is.finite(Sigma)) || any(diag(Sigma) <= 0)) return(NULL)
  R <- stats::cov2cor(0.5 * (Sigma + t(Sigma)))
  R[] <- pmin(pmax(R, -1), 1)
  diag(R) <- 1
  if (nrow(R) > 1L) {
    lambda <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
    if (lambda < min_eigen) {
      shrink <- (min_eigen - lambda) / (1 - lambda)
      R <- (1 - shrink) * R + diag(shrink, nrow(R))
    }
  }
  dimnames(R) <- NULL
  R
}

# Global statistic z' R^-1 z of the standardized coefficients z with
# correlation matrix R; with z_k = b_k / se_k and Sigma = S R S, S = diag(se),
# it equals b' Sigma^-1 b.
.ac3_mg_global_statistic <- function(z, R) {
  sum(z * solve(R, z))
}

# Statistics of the DA channel for the multigroup tables. scored is the scored
# table of the primary analysis (score_ancombc3_da() with the group columns as
# terms), x_columns the column names of the abundance design, mg the list of
# the group (levels, columns, eligibility), and cfg and p_adj_method the
# settings of the primary analysis. With pairs = TRUE, the comparisons of two non-reference
# levels are built by .ac3_mg_da_contrast() and scored by score_ancombc3_da()
# with the eligible taxa of each comparison, so that each receives its own
# common-shift estimate, working variance, t reference, and empirical-null
# scale; no model is refitted. Returns, over taxa x comparisons, est
# (corrected coefficient), se (its standard error, the square root of the
# working variance times the empirical-null scale), stat (the statistic
# est / se), df (per taxon, the degrees of freedom of the t reference), p
# (unadjusted p-value), eligible, and has (eligible with a finite p-value);
# and, for the K - 1 comparisons with the reference,
#   b, se_eff (the columns of est and se), and z = b / se_eff (those of stat);
#   corr: per taxon, the correlation matrix R_j of the joint working
#     covariance (.ac3_mg_da_joint_covariance(), .ac3_mg_correlation()). With
#     Omega = diag(empirical-null scales), the working covariance of b is
#     Omega Sigma_j Omega = S R_j S, S = diag(se_eff), and the working null
#     distribution of z is multivariate t with correlation matrix R_j and
#     df[j] degrees of freedom;
#   global: eligible, stat = z' R_j^-1 z, df2 (df), p (the upper tail
#     probability of stat / (K - 1) in the F distribution with K - 1 and df
#     degrees of freedom), and has.
.ac3_mg_da_channel <- function(fit, scored, x_columns, mg, cfg, p_adj_method, pairs = TRUE) {
  taxa <- names(fit$fits)
  p <- length(taxa)
  el <- mg$eligibility
  pr <- el$pairs
  m <- length(mg$columns)
  blank <- matrix(NA_real_, p, nrow(pr), dimnames = list(taxa, pr$comparison))
  est <- se <- stat <- pv <- dft <- blank
  read_scored <- function(tab, term) {
    s <- tab[tab$term == term, , drop = FALSE]
    s[match(taxa, s$taxon), c("beta_corrected", "se_scaled", "statistic", "p_value", "df_t"), drop = FALSE]
  }
  scored_all <- list(scored)
  terms <- mg$columns
  if (isTRUE(pairs) && nrow(pr) > m) {
    other <- (m + 1L):nrow(pr)
    built <- lapply(other, function(i) {
      .ac3_mg_da_contrast(fit, x_columns, col_num = mg$columns[pr$num[i] - 1L],
                          col_ref = mg$columns[pr$ref[i] - 1L], term = pr$comparison[i])
    })
    scored_all[[2L]] <- .ac3_mg_score(
      do.call(rbind, lapply(built, function(b) b$raw)),
      variance_inputs = c(
        list(taxa = fit$da_variance_inputs$taxa,
             G = stats::setNames(lapply(built, function(b) b$G), pr$comparison[other]),
             V = fit$da_variance_inputs$V),
        if (is.null(fit$da_variance_inputs$A)) NULL else
          list(A = stats::setNames(lapply(built, function(b) b$A), pr$comparison[other]))
      ),
      eligible_taxa_by_term = stats::setNames(lapply(other, function(i) taxa[el$da_pair[, i]]),
                                              pr$comparison[other]),
      cfg = cfg, p_adj_method = p_adj_method
    )
    terms <- c(terms, pr$comparison[other])
  }
  for (i in seq_along(terms)) {
    s <- read_scored(scored_all[[if (i <= m) 1L else 2L]], terms[i])
    est[, i] <- s$beta_corrected
    se[, i] <- s$se_scaled
    stat[, i] <- s$statistic
    pv[, i] <- s$p_value
    dft[, i] <- s$df_t
  }
  eligible <- el$da_pair
  has <- eligible & is.finite(pv)
  reference <- seq_len(m)
  b <- est[, reference, drop = FALSE]
  se_eff <- se[, reference, drop = FALSE]
  z <- stat[, reference, drop = FALSE]
  # Degrees of freedom of the taxon: those of its comparisons with the
  # reference, which are equal (the smallest finite value otherwise).
  df <- apply(dft[, reference, drop = FALSE], 1L, function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0L) NA_real_ else min(v)
  })
  jc <- .ac3_mg_da_joint_covariance(fit, scored, mg$columns, x_columns)
  corr <- lapply(seq_len(p), function(j) {
    .ac3_mg_correlation(jc$Sigma[j, , ])
  })
  g_stat <- vapply(seq_len(p), function(j) {
    if (is.null(corr[[j]]) || any(!is.finite(z[j, ]))) return(NA_real_)
    .ac3_mg_global_statistic(z[j, ], corr[[j]])
  }, numeric(1))
  g_has <- unname(el$da_global) & is.finite(g_stat)
  list(
    est = est, se = se, stat = stat, p = pv, eligible = eligible, has = has,
    b = b, se_eff = se_eff, z = z, corr = corr, df = df,
    global = list(eligible = unname(el$da_global), stat = g_stat, df2 = df,
                  p = stats::pf(g_stat / m, df1 = m, df2 = df, lower.tail = FALSE), has = g_has)
  )
}

# -----------------------------------------------------------------------------
# DP channel
# -----------------------------------------------------------------------------

# Joint DP quantities of the coefficients at the positions idx of the
# zero-pattern model of one taxon, computed inside
# fit_ancombc3_dp_zero_pattern() (argument joint_terms) from its unrestricted
# fit (fit, with ok its convergence indicator), the zero indicators Z, the
# sampling-zero probabilities f0, the design D = [1, r], and the second-stage
# quantities ss (.dp_second_stage_taxon(); NULL for a taxon with the supplied
# f0).
#   Global test of the m = length(idx) coefficients: the restricted fit holds
#     them at zero (BFGS from the unrestricted estimate of the other
#     coefficients). The penalized likelihood-ratio statistic T^LR is
#     multiplied by
#       rho = m / tr{(V^DP)^-1 (V^DP + V^prop)},
#     where V^DP and V^prop = G V^theta G' are the blocks of the m coefficients
#     of the inverse information and of the propagated covariance
#     (.dp_two_step_variance()) at the restricted fit. rho is the first-order
#     (Rao-Scott type) generalization of the scalar variance ratio of the
#     one-degree-of-freedom test, to which it reduces for m = 1; rho = 1 for a
#     taxon with the supplied f0 or without the two-step variance. The P value
#     is 1 - F_chi2_m(rho T^LR). When either fit fails, the statistic is 0 and
#     the P value is 1.
#   Covariance of the m estimates: the block of V^DP + V^prop at the
#     unrestricted fit (V^DP alone for a taxon with the supplied f0), or at the
#     restricted fit when it is not available at the unrestricted fit.
# With repeated measures, ids is the subject index of each sample
# (.ac3_cluster_index()), and V^DP is replaced by its cluster-robust version
# V^DP,CR (.ac3_cluster_dp_covariance()) in the same way as in the
# one-degree-of-freedom test (.ac3_cluster_dp_variance()):
#   rho = m / tr{(V^DP)^-1 (V^DP,CR + V^prop)} at the restricted fit, for every
#     taxon (V^prop = 0 for a taxon with the supplied f0), which reduces to the
#     variance ratio of the one-degree-of-freedom test for m = 1; the ratio
#     used is the smaller of this value and the ratio without ids;
#   covariance: the block of V^DP,CR + V^prop.
# When the cluster-robust components are not available at a fit, the
# quantities without cluster are used at that fit.
# Returns estimate, statistic, statistic_lr, p_value, df, variance_ratio,
# covariance (NULL when not available), covariance_source ("unrestricted" or
# "restricted"), and robust (TRUE when the covariance is cluster-robust).
.ac3_mg_dp_joint_test <- function(Z, f0, D, ridge, fit, ok, ss, idx, ids = NULL) {
  m <- length(idx)
  pcol <- ncol(D)
  out <- list(estimate = if (ok) fit$par[idx] else rep(NA_real_, m), statistic = 0,
              statistic_lr = 0, p_value = 1, df = m, variance_ratio = 1, covariance = NULL,
              covariance_source = NA_character_, robust = FALSE)
  null_nll <- .dp_zero_pattern_nll(Z, f0, D, ridge = ridge, fix_idx = idx, fix_vals = rep(0, m))
  init <- if (ok) fit$par[-idx] else rep(0, pcol - m)
  null_fit <- tryCatch(
    stats::optim(init, null_nll, method = "BFGS", control = list(maxit = 500)),
    error = function(e) NULL
  )
  null_ok <- !is.null(null_fit) && null_fit$convergence == 0 && all(is.finite(null_fit$par))
  variance_blocks <- function(coef) {
    if (!is.null(ss)) {
      tv <- .dp_two_step_variance(Z, D, coef, ridge, ss)
      if (is.null(tv)) return(NULL)
      blocks <- list(V_dp = tv$V_dp[idx, idx, drop = FALSE], V_prop = tv$V_prop[idx, idx, drop = FALSE])
    } else {
      info <- .dp_zero_pattern_information(Z, f0, D, coef, ridge)
      if (is.null(info$H)) return(NULL)
      blocks <- list(V_dp = chol2inv(chol(info$H))[idx, idx, drop = FALSE], V_prop = matrix(0, m, m))
    }
    total <- blocks$V_dp + blocks$V_prop
    if (any(!is.finite(total)) || any(diag(total) <= 0)) return(NULL)
    blocks
  }
  robust_blocks <- function(coef) {
    cr <- .ac3_cluster_dp_covariance(Z, f0, D, coef, ridge, ids)
    if (is.null(cr)) return(NULL)
    V_prop <- matrix(0, m, m)
    if (!is.null(ss)) {
      G <- cr$V_dp %*% crossprod(cr$dscore_df0, ss$jacobian)
      V_full <- G %*% ss$vcov_theta %*% t(G)
      V_prop <- (0.5 * (V_full + t(V_full)))[idx, idx, drop = FALSE]
    }
    blocks <- list(V_dp = cr$V_dp[idx, idx, drop = FALSE], V_cr = cr$V_cr[idx, idx, drop = FALSE],
                   V_prop = V_prop)
    total <- blocks$V_cr + blocks$V_prop
    if (any(!is.finite(total)) || any(diag(total) <= 0) || any(diag(blocks$V_dp) <= 0)) return(NULL)
    blocks
  }
  null_blocks <- NULL
  null_robust <- NULL
  if (null_ok) {
    coef0 <- numeric(pcol)
    coef0[-idx] <- null_fit$par
    null_blocks <- variance_blocks(coef0)
    if (!is.null(ids)) null_robust <- robust_blocks(coef0)
  }
  rho <- 1
  if (!is.null(ss) && !is.null(null_blocks)) {
    excess <- tryCatch(sum(diag(solve(null_blocks$V_dp, null_blocks$V_prop))), error = function(e) NA_real_)
    if (is.finite(excess)) rho <- m / (m + max(excess, 0))
  }
  if (!is.null(null_robust)) {
    spread <- tryCatch(sum(diag(solve(null_robust$V_dp, null_robust$V_cr + null_robust$V_prop))),
                       error = function(e) NA_real_)
    if (is.finite(spread) && spread > 0) rho <- min(rho, m / spread)
  }
  if (ok && null_ok && is.finite(fit$value) && is.finite(null_fit$value)) {
    out$statistic_lr <- max(2 * ((-fit$value) - (-null_fit$value)), 0)
    out$statistic <- rho * out$statistic_lr
    out$p_value <- stats::pchisq(out$statistic, df = m, lower.tail = FALSE)
  }
  out$variance_ratio <- rho
  # Covariance of the m estimates at the restricted fit.
  null_total <- if (!is.null(null_robust)) {
    null_robust$V_cr + null_robust$V_prop
  } else if (!is.null(null_blocks)) {
    null_blocks$V_dp + null_blocks$V_prop
  } else NULL
  out$covariance_null <- if (is.null(null_total)) NULL else 0.5 * (null_total + t(null_total))
  if (!is.null(ids)) {
    robust <- if (ok) robust_blocks(fit$par) else NULL
    out$covariance_source <- "unrestricted"
    if (is.null(robust)) {
      robust <- null_robust
      out$covariance_source <- "restricted"
    }
    if (!is.null(robust)) {
      total <- robust$V_cr + robust$V_prop
      out$covariance <- 0.5 * (total + t(total))
      out$robust <- TRUE
      return(out)
    }
  }
  blocks <- if (ok) variance_blocks(fit$par) else NULL
  out$covariance_source <- "unrestricted"
  if (is.null(blocks)) {
    blocks <- null_blocks
    out$covariance_source <- "restricted"
  }
  if (is.null(blocks)) {
    out$covariance_source <- NA_character_
  } else {
    total <- blocks$V_dp + blocks$V_prop
    out$covariance <- 0.5 * (total + t(total))
  }
  out
}

# Statistics of the DP channel for the multigroup tables. dp_z is the result of
# fit_ancombc3_dp_zero_pattern() for the primary analysis with
# joint_terms = mg$columns: its rows give the one-degree-of-freedom tests of
# the comparisons with the reference, and its attribute joint the global test
# and the covariance of the group coefficients (.ac3_mg_dp_joint_test()). With
# pairs = TRUE, the comparison of two non-reference levels k < l is the
# existing one-degree-of-freedom test of the column of level l in the
# prevalence design in which level k is the reference: the column of level k is
# replaced by the indicator of the original reference level, and
# fit_ancombc3_dp_zero_pattern() is applied to that design with the same
# sampling-zero probabilities and second-stage fit (and, with repeated
# measures, the subject of each sample in cluster, so that these tests use the
# cluster-robust variance ratio of the primary DP test). Returns, over taxa x
# comparisons, est, stat (chi-square statistic T^DP), p, eligible, and has;
# and, for the K - 1 comparisons with the reference,
#   z: signed roots sign(est) sqrt(T^DP);
#   se_eff: effective standard errors |est| / |z| (the square root of the
#     diagonal of the covariance when z is zero or not finite);
#   cov_null: per taxon, the covariance of the group coefficients at the
#     restricted fit (NULL where not available);
#   se_trend, corr_trend: the standard errors and correlation matrices of the
#     trend test (below);
#   corr: per taxon, the correlation matrix of the covariance of the group
#     coefficients (.ac3_mg_correlation());
#   global: eligible, stat, p, and has.
.ac3_mg_dp_channel <- function(y, r, x, f0, da_fit, dp_z, mg, cfg, alpha, pairs = TRUE, cluster = NULL,
                               engine = "r") {
  y <- as.matrix(y)
  taxa <- .ac3_taxon_names(y)
  p <- length(taxa)
  el <- mg$eligibility
  pr <- el$pairs
  m <- length(mg$columns)
  blank <- matrix(NA_real_, p, nrow(pr), dimnames = list(taxa, pr$comparison))
  est <- stat <- pv <- blank
  read_dp <- function(tab, term) {
    s <- tab[tab$test_term == term, , drop = FALSE]
    s[match(taxa, s$taxon), c("estimate", "statistic", "p_value"), drop = FALSE]
  }
  for (i in seq_len(m)) {
    s <- read_dp(dp_z, mg$columns[i])
    est[, i] <- s$estimate
    stat[, i] <- s$statistic
    pv[, i] <- s$p_value
  }
  if (isTRUE(pairs) && m >= 2L) {
    r <- as.matrix(r)
    reference_indicator <- 1 - rowSums(r[, mg$columns, drop = FALSE])
    for (k in 2:m) {
      ref_col <- mg$columns[k - 1L]
      num_cols <- mg$columns[k:m]
      r_k <- r
      r_k[, ref_col] <- reference_indicator
      colnames(r_k)[colnames(r_k) == ref_col] <-
        utils::tail(make.unique(c(setdiff(colnames(r), ref_col), mg$reference_column)), 1L)
      tab <- fit_ancombc3_dp_zero_pattern(
        y = y, r = r_k, f0 = f0, test_terms = num_cols, alpha = alpha, p_adj_method = "BH",
        ridge = cfg$ridge_dp, eligibility = NULL,
        da_fit = da_fit, x = x, cluster = cluster, engine = engine
      )
      for (l in (k + 1L):(m + 1L)) {
        i <- which(pr$ref == k & pr$num == l)
        s <- read_dp(tab, mg$columns[l - 1L])
        est[, i] <- s$estimate
        stat[, i] <- s$statistic
        pv[, i] <- s$p_value
      }
    }
  }
  eligible <- el$dp_pair
  has <- eligible & is.finite(pv)
  reference <- seq_len(m)
  b <- est[, reference, drop = FALSE]
  z <- ifelse(is.finite(b), sign(b), 0) * sqrt(pmax(stat[, reference, drop = FALSE], 0))
  se_eff <- abs(b) / abs(z)
  joint <- attr(dp_z, "joint")
  corr <- vector("list", p)
  cov_null <- vector("list", p)
  g_stat <- g_p <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    jt <- joint[[taxa[j]]]
    if (is.null(jt)) next
    if (!is.null(jt$covariance_null)) cov_null[[j]] <- jt$covariance_null
    g_stat[j] <- jt$statistic
    g_p[j] <- jt$p_value
    if (is.null(jt$covariance)) next
    corr[[j]] <- .ac3_mg_correlation(jt$covariance)
    replace <- !(is.finite(se_eff[j, ]) & se_eff[j, ] > 0)
    se_eff[j, replace] <- sqrt(diag(jt$covariance))[replace]
  }
  g_has <- unname(el$dp_global) & is.finite(g_p)
  # Covariance of the trend test: the covariance of the group coefficients at
  # the restricted fit (se_trend its standard errors, corr_trend its
  # correlation matrix); a taxon without it uses se_eff and corr.
  se_trend <- se_eff
  corr_trend <- corr
  for (j in seq_len(p)) {
    if (is.null(cov_null[[j]])) next
    R_null <- .ac3_mg_correlation(cov_null[[j]])
    if (is.null(R_null)) next
    se_trend[j, ] <- sqrt(diag(cov_null[[j]]))
    corr_trend[[j]] <- R_null
  }
  list(
    est = est, stat = stat, p = pv, eligible = eligible, has = has,
    b = b, se_eff = se_eff, z = z, corr = corr, cov_null = cov_null,
    se_trend = se_trend, corr_trend = corr_trend,
    global = list(eligible = unname(el$dp_global), stat = g_stat, p = g_p, has = g_has)
  )
}

# -----------------------------------------------------------------------------
# Mixed directional false discovery rate
# -----------------------------------------------------------------------------
# Two-stage procedures for directional decisions on several comparisons per
# taxon (Guo, Sarkar, and Peddada, 2010, Biometrics 66:485-492). Stage 1: each
# taxon with at least one eligible comparison has a screening p-value for the
# hypothesis that all its eligible comparisons are null; the Benjamini-Hochberg
# procedure at level alpha over the n_tax screened taxa selects R taxa.
# Stage 2: within a selected taxon, the p-values of its eligible comparisons
# are adjusted for the family-wise error rate at level R alpha / n_tax;
# equivalently, the adjusted p-values are multiplied by n_tax / R and compared
# with alpha. The reported q is this product, capped at one, for a selected
# taxon and one for a screened taxon that is not selected; a comparison is a
# discovery when q <= alpha, and its direction is the sign of the estimate.

# Stage 1 for a vector of screening p-values (NA for the taxa that are not
# screened). Returns q_screen, selected, n_tax, and n_selected.
.ac3_mg_mdfdr_screen <- function(p_screen, alpha) {
  screened <- is.finite(p_screen)
  q_screen <- rep(NA_real_, length(p_screen))
  if (any(screened)) q_screen[screened] <- stats::p.adjust(p_screen[screened], method = "BH")
  selected <- screened & q_screen <= alpha
  selected[is.na(selected)] <- FALSE
  list(q_screen = q_screen, selected = selected, n_tax = sum(screened), n_selected = sum(selected))
}

# Pairwise comparisons. p is the taxa x comparisons matrix of unadjusted
# p-values and has flags the eligible comparisons with a p-value. The
# screening p-value of a taxon is the p-value of its global test (p_global)
# when the taxon has one (has_global: all levels supported), and otherwise the
# Bonferroni-adjusted minimum min(1, m_j min p) over its m_j eligible
# comparisons. Within a selected taxon, the p-values of its eligible
# comparisons are adjusted by stats::p.adjust() with fwer_ctrl_method (Holm by
# default). Returns q (NA for the comparisons without a p-value), p_screen,
# q_screen, selected, n_tax, and n_selected.
.ac3_mg_mdfdr_pairwise <- function(p, has, p_global, has_global, alpha, fwer_ctrl_method = "holm") {
  p <- as.matrix(p)
  n_comparisons <- rowSums(has)
  p_min <- suppressWarnings(apply(ifelse(has, p, NA_real_), 1L, min, na.rm = TRUE))
  p_screen <- ifelse(has_global %in% TRUE & is.finite(p_global), p_global, pmin(1, n_comparisons * p_min))
  p_screen[n_comparisons == 0L] <- NA_real_
  screen <- .ac3_mg_mdfdr_screen(p_screen, alpha)
  q <- matrix(NA_real_, nrow(p), ncol(p), dimnames = dimnames(p))
  q[has] <- 1
  for (j in which(screen$selected)) {
    e <- which(has[j, ])
    q[j, e] <- pmin(1, stats::p.adjust(p[j, e], method = fwer_ctrl_method) * screen$n_tax / screen$n_selected)
  }
  c(list(q = q, p_screen = p_screen), screen)
}

# Rectangle probability P(|Z_l| <= bound for all l) in dimension m = 2 or 3,
# for Z ~ N(0, R) (df = Inf) or Z multivariate t with correlation matrix R and
# df degrees of freedom (a positive integer), by inclusion and exclusion over
# the 2^m orthant probabilities P(Z_l <= s_l bound for all l), s_l = +1 or -1,
# each computed by the deterministic bivariate and trivariate algorithm of
# Genz (2004) (mvtnorm::pmvnorm() and mvtnorm::pmvt() with algorithm TVPACK).
.ac3_mg_rectangle_exact <- function(bound, R, df = Inf) {
  m <- nrow(R)
  total <- 0
  for (pattern in 0:(2^m - 1L)) {
    negative <- bitwAnd(pattern, 2^(0:(m - 1L))) > 0
    upper <- ifelse(negative, -bound, bound)
    orthant <- if (is.finite(df)) {
      mvtnorm::pmvt(lower = rep(-Inf, m), upper = upper, df = as.integer(df), corr = R,
                    algorithm = mvtnorm::TVPACK(abseps = 1e-12))[[1L]]
    } else {
      mvtnorm::pmvnorm(lower = rep(-Inf, m), upper = upper, corr = R,
                       algorithm = mvtnorm::TVPACK(abseps = 1e-12))[[1L]]
    }
    total <- total + (if (sum(negative) %% 2L == 0L) orthant else -orthant)
  }
  min(max(total, 0), 1)
}

# Single-step adjusted p-values of the many-to-one (Dunnett-type) comparisons
# of one taxon: for standardized statistics z with correlation matrix R and
# the events A_l = {|Z_l| > |z_k|}, where Z ~ N(0, R) (df = Inf) or Z is
# multivariate t with correlation matrix R and df degrees of freedom,
#   p_k = P(A_1 or ... or A_m) = 1 - P(|Z_l| <= |z_k| for all l).
# The probability is computed
#   in closed form for one comparison, the two-sided tail probability of z_k;
#   by .ac3_mg_rectangle_exact() for two and three comparisons;
#   for four or more comparisons, from the decomposition
#       p_k = P(A_1 or A_2 or A_3) + sum_{l >= 4} P(A_l and not A_1, ..., A_{l-1}),
#     whose first term is computed by .ac3_mg_rectangle_exact() and whose
#     further terms equal 2 P(Z_l > |z_k|, |Z_i| <= |z_k| for i < l) by the
#     symmetry of the distribution. Each of these is computed by the
#     quasi-Monte Carlo algorithm of Genz and Bretz (mvtnorm::pmvnorm() or
#     mvtnorm::pmvt(), at
#     most 25,000 points, absolute tolerance 0.0005 times the marginal
#     p-value), so that the numerical error is small relative to p_k for
#     small and large p_k alike. The algorithm uses random shifts, which are
#     generated from a fixed seed so that the result is a deterministic
#     function of z and R; the state of the random-number stream of the caller
#     is restored. The reported error bounds of the algorithm are added to
#     p_k.
# The result is restricted to its exact bounds, the marginal two-sided p-value
# from below and its Bonferroni multiple from above. With R = NULL
# (correlation not available), the Bonferroni bound is returned.
.ac3_mg_dunnett_p <- function(z, R = NULL, df = Inf) {
  m <- length(z)
  a <- abs(z)
  if (length(df) != 1L || is.na(df) || df <= 0) df <- Inf
  marginal <- 2 * stats::pt(-a, df = df)
  upper <- pmin(1, m * marginal)
  if (m == 1L) return(marginal)
  if (is.null(R)) return(upper)
  if (m <= 3L) {
    p <- vapply(a, function(bound) {
      tryCatch(1 - .ac3_mg_rectangle_exact(bound, R, df), error = function(e) NA_real_)
    }, numeric(1))
  } else {
    withr::local_seed(20261110L, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                       .rng_sample_kind = "Rejection")
    first <- seq_len(3L)
    p <- vapply(seq_len(m), function(k) {
      bound <- a[k]
      tryCatch({
        total <- 1 - .ac3_mg_rectangle_exact(bound, R[first, first], df)
        for (l in 4:m) {
          region <- list(
            lower = c(rep(-bound, l - 1L), bound), upper = c(rep(bound, l - 1L), Inf),
            corr = R[seq_len(l), seq_len(l)],
            algorithm = mvtnorm::GenzBretz(maxpts = 25000L, abseps = max(1e-12, 5e-4 * marginal[k]))
          )
          term <- if (is.finite(df)) {
            do.call(mvtnorm::pmvt, c(region, list(df = as.integer(df))))
          } else {
            do.call(mvtnorm::pmvnorm, region)
          }
          error <- attr(term, "error")
          total <- total + 2 * (term[[1L]] + (if (is.null(error) || !is.finite(error)) 0 else error))
        }
        total
      }, error = function(e) NA_real_)
    }, numeric(1))
  }
  p[!is.finite(p)] <- upper[!is.finite(p)]
  pmin(pmax(p, marginal), upper)
}

# Dunnett-type comparisons with the reference. z is the taxa x (K - 1) matrix
# of standardized statistics, has flags the supported comparisons with a
# statistic, and corr holds the correlation matrix of each taxon (NULL where
# not available). For each taxon, the adjusted p-values of its supported
# comparisons are computed by .ac3_mg_dunnett_p() under N(0, R_j) restricted to
# those comparisons, or under the multivariate t distribution with df[j]
# degrees of freedom when df is given; the screening p-value is their minimum, the p-value of
# the maximum absolute statistic. Returns p_adjusted, q (min(1, p_adjusted
# n_tax / R) for a selected taxon, one otherwise, NA without a statistic),
# p_screen, q_screen, selected, n_tax, and n_selected.
.ac3_mg_mdfdr_dunnett <- function(z, has, corr, alpha, df = NULL) {
  z <- as.matrix(z)
  has <- has & is.finite(z)
  p_adjusted <- matrix(NA_real_, nrow(z), ncol(z), dimnames = dimnames(z))
  p_screen <- rep(NA_real_, nrow(z))
  for (j in which(rowSums(has) > 0L)) {
    e <- which(has[j, ])
    R <- corr[[j]]
    if (!is.null(R)) R <- R[e, e, drop = FALSE]
    p_adjusted[j, e] <- .ac3_mg_dunnett_p(z[j, e], R, if (is.null(df)) Inf else df[j])
    p_screen[j] <- min(p_adjusted[j, e])
  }
  screen <- .ac3_mg_mdfdr_screen(p_screen, alpha)
  q <- matrix(NA_real_, nrow(z), ncol(z), dimnames = dimnames(z))
  q[has] <- 1
  sel <- which(screen$selected)
  if (length(sel) > 0L) {
    scaled <- pmin(1, p_adjusted[sel, , drop = FALSE] * screen$n_tax / screen$n_selected)
    q[sel, ] <- ifelse(has[sel, , drop = FALSE], scaled, NA_real_)
  }
  c(list(p_adjusted = p_adjusted, q = q, p_screen = p_screen), screen)
}

# -----------------------------------------------------------------------------
# Result tables
# -----------------------------------------------------------------------------

.ac3_mg_na_unless <- function(v, keep) {
  v[!keep] <- NA
  unname(v)
}

# Omnibus columns of one comparison or test: omnibus_pvalue() applied to the
# DA and DP p-values, the adjusted p-values over
# the taxa with at least one of them (p_adj_method), and the discoveries at
# level alpha.
.ac3_mg_omnibus <- function(p_da, p_dp, has_da, has_dp, p_adj_method, alpha) {
  p_omni <- omnibus_pvalue(unname(p_da), unname(p_dp), unname(has_da), unname(has_dp))
  q_omni <- rep(NA_real_, length(p_omni))
  ok <- (has_da | has_dp) & is.finite(p_omni)
  if (any(ok)) q_omni[ok] <- stats::p.adjust(p_omni[ok], method = p_adj_method)
  discovery <- ok & q_omni <= alpha
  discovery[is.na(discovery)] <- FALSE
  data.frame(omnibus_p = p_omni, omnibus_q = q_omni, omnibus_discovery = unname(discovery))
}

# Table of the comparisons which (indices into mg$eligibility$pairs), one row
# per taxon and comparison. da and dp are the channel statistics; p_da, q_da,
# p_dp, and q_dp are taxa x length(which) matrices of the reported p-values and
# of the mixed-directional-FDR adjusted p-values. The omnibus columns use the
# unadjusted p-values of the comparison. The columns of a channel other than
# its eligibility indicator are NA for a taxon without a p-value in that
# channel. p_adj_method is the multiplicity adjustment of the omnibus columns.
.ac3_mg_comparison_table <- function(taxa, mg, which, da, dp, p_da, q_da, p_dp, q_dp, p_adj_method,
                                     alpha) {
  pr <- mg$eligibility$pairs
  direction <- function(discovery, estimate) {
    ifelse(discovery, ifelse(estimate > 0, "positive", "negative"), NA_character_)
  }
  rows <- lapply(seq_along(which), function(ii) {
    i <- which[ii]
    has_da <- unname(da$has[, i])
    has_dp <- unname(dp$has[, i])
    disc_da <- (has_da & q_da[, ii] <= alpha) %in% TRUE
    disc_dp <- (has_dp & q_dp[, ii] <= alpha) %in% TRUE
    data.frame(
      taxon = taxa,
      group = mg$name,
      comparison = pr$comparison[i],
      level_num = pr$level_num[i],
      level_ref = pr$level_ref[i],
      da_eligible = unname(da$eligible[, i]),
      da_estimate = .ac3_mg_na_unless(da$est[, i], has_da),
      da_se = .ac3_mg_na_unless(da$se[, i], has_da),
      da_statistic = .ac3_mg_na_unless(da$stat[, i], has_da),
      da_p = .ac3_mg_na_unless(p_da[, ii], has_da),
      da_q = .ac3_mg_na_unless(q_da[, ii], has_da),
      da_discovery = disc_da,
      da_direction = unname(direction(disc_da, da$est[, i])),
      dp_eligible = unname(dp$eligible[, i]),
      dp_estimate = .ac3_mg_na_unless(dp$est[, i], has_dp),
      dp_statistic = .ac3_mg_na_unless(dp$stat[, i], has_dp),
      dp_p = .ac3_mg_na_unless(p_dp[, ii], has_dp),
      dp_q = .ac3_mg_na_unless(q_dp[, ii], has_dp),
      dp_discovery = disc_dp,
      dp_direction = unname(direction(disc_dp, dp$est[, i])),
      .ac3_mg_omnibus(da$p[, i], dp$p[, i], has_da, has_dp, p_adj_method, alpha),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Adjusted p-values over the taxa with a p-value (has) by the given method.
.ac3_mg_adjust <- function(p, has, method) {
  q <- rep(NA_real_, length(p))
  if (any(has)) q[has] <- stats::p.adjust(p[has], method = method)
  q
}

# Table of the global tests, one row per taxon. The DA, DP, and omnibus
# q-values use p_adj_method (the p_adj_method of ancombc3()), as in the
# primary analysis, each over the taxa eligible for the test.
.ac3_mg_global_table <- function(taxa, mg, da, dp, p_adj_method, alpha) {
  m <- length(mg$columns)
  has_da <- da$global$has
  has_dp <- dp$global$has
  q_da <- .ac3_mg_adjust(da$global$p, has_da, p_adj_method)
  q_dp <- .ac3_mg_adjust(dp$global$p, has_dp, p_adj_method)
  data.frame(
    taxon = taxa,
    group = mg$name,
    da_eligible = da$global$eligible,
    da_statistic = .ac3_mg_na_unless(da$global$stat, has_da),
    da_df = ifelse(has_da, m, NA_integer_),
    da_df2 = .ac3_mg_na_unless(da$global$df2, has_da),
    da_p = .ac3_mg_na_unless(da$global$p, has_da),
    da_q = q_da,
    da_discovery = (has_da & q_da <= alpha) %in% TRUE,
    dp_eligible = dp$global$eligible,
    dp_statistic = .ac3_mg_na_unless(dp$global$stat, has_dp),
    dp_df = ifelse(has_dp, m, NA_integer_),
    dp_p = .ac3_mg_na_unless(dp$global$p, has_dp),
    dp_q = q_dp,
    dp_discovery = (has_dp & q_dp <= alpha) %in% TRUE,
    .ac3_mg_omnibus(da$global$p, dp$global$p, has_da, has_dp, p_adj_method, alpha),
    stringsAsFactors = FALSE
  )
}

# Table of the trend tests, one row per taxon. trend is the result of
# ancombc3_trend_test() for the channels da and dp. The fit columns hold the
# constrained estimates of the K - 1 non-reference levels under the reported
# pattern. The q-values are adjusted as in .ac3_mg_global_table().
.ac3_mg_trend_table <- function(taxa, mg, da, dp, trend, p_adj_method, alpha) {
  channel <- function(res, eligible, method, prefix) {
    has <- eligible & is.finite(res$p_value)
    q <- .ac3_mg_adjust(res$p_value, has, method)
    fit <- res$fit
    fit[!has, ] <- NA_real_
    colnames(fit) <- paste0(prefix, "_fit_", mg$levels[-1L])
    out <- data.frame(
      eligible = eligible,
      statistic = .ac3_mg_na_unless(res$statistic, has),
      p = .ac3_mg_na_unless(res$p_value, has),
      q = q,
      discovery = (has & q <= alpha) %in% TRUE,
      pattern = .ac3_mg_na_unless(res$pattern, has),
      stringsAsFactors = FALSE
    )
    names(out) <- paste0(prefix, "_", names(out))
    list(table = cbind(out, as.data.frame(fit, check.names = FALSE, stringsAsFactors = FALSE)), has = has)
  }
  t_da <- channel(trend$da, da$global$eligible, p_adj_method, "da")
  t_dp <- channel(trend$dp, dp$global$eligible, p_adj_method, "dp")
  out <- cbind(
    data.frame(taxon = taxa, group = mg$name, stringsAsFactors = FALSE),
    t_da$table, t_dp$table,
    .ac3_mg_omnibus(trend$da$p_value, trend$dp$p_value, t_da$has, t_dp$has, p_adj_method, alpha)
  )
  rownames(out) <- NULL
  out
}

# Multigroup result tables of ancombc3(). y, x, and r are the counts and the
# two designs; mg is the list of the group (ancombc3_group_arguments(),
# ancombc3_group_setup(), with eligibility = ancombc3_group_eligibility());
# f0 is the supplied sampling-zero probability of the DP test; fit and scored
# are the second-stage fit, also used by the DP tests of the comparisons, and
# its scored table; dp_z is the DP result with joint_terms = mg$columns; cfg,
# alpha, p_adj_method, and seed are the settings, the level, the multiplicity
# adjustment, and the seed of the analysis; cluster is the subject of each
# sample with repeated measures (NULL otherwise), used by the DP tests of the
# comparisons of two non-reference levels. With cluster, the DA quantities are
# cluster-robust through fit (fit_ancombc3_da()) and the joint DP quantities
# through dp_z. Returns the requested tables among
#   res_global: global tests of equality of the K levels;
#   res_pair: all K (K - 1) / 2 comparisons, with q from
#     .ac3_mg_mdfdr_pairwise() (the p columns are the unadjusted p-values);
#   res_dunn: the K - 1 comparisons with the reference, with q from
#     .ac3_mg_mdfdr_dunnett() (the p columns are the single-step adjusted
#     p-values);
#   res_trend: trend tests (ancombc3_trend_test()) with the estimates and the
#     effective covariance of the comparisons with the reference.
# engine selects the evaluation of the quadrature of the second-stage f0 of
# the DP comparisons (fit_ancombc3_dp_zero_pattern()).
ancombc3_multigroup <- function(y, x, r, mg, f0, fit, scored, dp_z, cfg, alpha, p_adj_method, seed,
                                cluster = NULL, engine = "r") {
  taxa <- .ac3_taxon_names(as.matrix(y))
  m <- length(mg$columns)
  reference <- seq_len(m)
  da <- .ac3_mg_da_channel(fit, scored, colnames(x), mg, cfg, p_adj_method, pairs = isTRUE(mg$pairwise))
  dp <- .ac3_mg_dp_channel(y, r, x, f0, fit, dp_z, mg, cfg, alpha, pairs = isTRUE(mg$pairwise),
                           cluster = cluster, engine = engine)
  out <- list()
  if (isTRUE(mg$global)) out$res_global <- .ac3_mg_global_table(taxa, mg, da, dp, p_adj_method, alpha)
  if (isTRUE(mg$pairwise)) {
    method <- mg$mdfdr_control$fwer_ctrl_method
    md_da <- .ac3_mg_mdfdr_pairwise(da$p, da$has, da$global$p, da$global$has, alpha, method)
    md_dp <- .ac3_mg_mdfdr_pairwise(dp$p, dp$has, dp$global$p, dp$global$has, alpha, method)
    out$res_pair <- .ac3_mg_comparison_table(taxa, mg, seq_len(nrow(mg$eligibility$pairs)), da, dp,
                                             da$p, md_da$q, dp$p, md_dp$q, p_adj_method, alpha)
  }
  if (isTRUE(mg$dunnet)) {
    dn_da <- .ac3_mg_mdfdr_dunnett(da$z, da$has[, reference, drop = FALSE], da$corr, alpha, df = da$df)
    dn_dp <- .ac3_mg_mdfdr_dunnett(dp$z, dp$has[, reference, drop = FALSE], dp$corr, alpha)
    out$res_dunn <- .ac3_mg_comparison_table(taxa, mg, reference, da, dp, dn_da$p_adjusted, dn_da$q,
                                             dn_dp$p_adjusted, dn_dp$q, p_adj_method, alpha)
  }
  if (isTRUE(mg$trend)) {
    trend <- ancombc3_trend_test(
      channels = list(
        da = list(b = da$b, se = da$se_eff, corr = da$corr, eligible = da$global$has, df = da$df),
        dp = list(b = dp$b, se = dp$se_trend, corr = dp$corr_trend, eligible = dp$global$has)
      ),
      patterns = mg$patterns, B = mg$trend_control$B, seed = seed
    )
    out$res_trend <- .ac3_mg_trend_table(taxa, mg, da, dp, trend, p_adj_method, alpha)
  }
  out
}

# Progress message of ancombc3() for the group columns of the primary table.
.ac3_mg_say <- function(verbose, results, mg, alpha) {
  for (i in seq_along(mg$columns)) {
    rows <- results$covariate == mg$name & results$level %in% mg$levels[i + 1L]
    .ac3_main_say(verbose, sprintf(
      "%s, %s versus %s: %d DA-eligible, %d DP-eligible taxa; discoveries at q <= %.2f: DA %d, DP %d, omnibus %d.",
      mg$name, mg$levels[i + 1L], mg$levels[1L], sum(results$da_eligible[rows]),
      sum(results$dp_eligible[rows]), alpha, sum(results$da_discovery[rows]),
      sum(results$dp_discovery[rows]), sum(results$omnibus_discovery[rows])))
  }
  invisible(NULL)
}
