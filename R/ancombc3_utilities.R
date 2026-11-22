# =============================================================================
# ANCOM-BC3 data import, designs, and DA eligibility
# =============================================================================
#
# ancombc3_normalize_input() reads a count table and its sample metadata, and
# ancombc3_make_design() builds the abundance design x, the prevalence design
# r, and the tested terms. reference_test_columns() sets the tested covariates
# to their reference values, the sample means, for the decoding of the
# first-stage offsets. ancombc3_da_eligibility() applies the DA eligibility
# criteria to every tested covariate and returns the eligible taxa of each
# covariate and their union, the fitting set used to estimate the
# sample-specific bias.
#
# Taxon matching. The taxa of an analysis are the column names of the count
# matrix y, or T1, ..., Tp when y has no column names (.ac3_taxon_names()).
# Every other sample-by-taxon matrix and every per-taxon vector is matched to
# these taxa by name when it carries taxon names, and by position otherwise
# (.ac3_match_taxa(), .ac3_match_taxa_vector()).

`%||%` <- function(a, b) if (is.null(a)) b else a

# Taxon names of a count matrix: its column names, or T1, ..., Tp.
.ac3_taxon_names <- function(y) {
  colnames(y) %||% paste0("T", seq_len(ncol(y)))
}

# Columns of the sample-by-taxon matrix mat for the taxa taxa, in that order,
# with column names taxa. A matrix with column names is matched by name and
# must contain every taxon of taxa; further columns are dropped. A matrix
# without column names is matched by position and must have one column per
# taxon. name is the argument name used in the error messages.
.ac3_match_taxa <- function(mat, taxa, name) {
  mat <- as.matrix(mat)
  if (is.null(colnames(mat))) {
    if (ncol(mat) != length(taxa)) stop(name, " must have one column per taxon.")
    colnames(mat) <- taxa
    return(mat)
  }
  if (anyDuplicated(colnames(mat))) stop(name, " has duplicated taxon names.")
  absent <- setdiff(taxa, colnames(mat))
  if (length(absent) > 0L) {
    stop(name, " has no column for ", length(absent), " taxon/taxa of the counts: ",
         paste(utils::head(absent, 5L), collapse = ", "), ".")
  }
  mat[, taxa, drop = FALSE]
}

# Elements of the per-taxon vector v for the taxa taxa, in that order, with
# names taxa: by name when v is named, by position otherwise.
.ac3_match_taxa_vector <- function(v, taxa, name) {
  if (is.null(names(v))) {
    if (length(v) != length(taxa)) stop(name, " must have one value per taxon.")
    return(stats::setNames(v, taxa))
  }
  if (anyDuplicated(names(v))) stop(name, " has duplicated taxon names.")
  absent <- setdiff(taxa, names(v))
  if (length(absent) > 0L) {
    stop(name, " has no value for ", length(absent), " taxon/taxa of the counts: ",
         paste(utils::head(absent, 5L), collapse = ", "), ".")
  }
  v[taxa]
}

# Count table and sample metadata of ancombc3(), with the samples in rows. data
# is a phyloseq object (its OTU table and sample data), a SummarizedExperiment
# (the assay assay_name and colData), or a matrix or data frame of counts, with
# taxa in rows when taxa_are_rows is TRUE, whose sample metadata is meta_data.
# Samples and taxa without names are named S1, ..., Sn and T1, ..., Tp, and
# metadata without row names is matched to the samples by position. Returns
# y (n x p) and metadata, with the rows of metadata in the order of the rows
# of y.
ancombc3_normalize_input <- function(data,
                                     meta_data = NULL,
                                     taxa_are_rows = TRUE,
                                     assay_name = "counts") {
  if (inherits(data, "phyloseq")) {
    otu <- methods::as(phyloseq::otu_table(data), "matrix")
    if (phyloseq::taxa_are_rows(data)) otu <- t(otu)
    metadata <- methods::as(phyloseq::sample_data(data), "data.frame")
  } else if (inherits(data, "SummarizedExperiment")) {
    otu <- as.matrix(SummarizedExperiment::assay(data, assay_name))
    otu <- t(otu)
    metadata <- as.data.frame(SummarizedExperiment::colData(data))
  } else {
    otu <- as.matrix(data)
    metadata <- meta_data
    if (taxa_are_rows) otu <- t(otu)
  }

  storage.mode(otu) <- "double"
  if (anyNA(otu)) stop("Count table contains missing values.")
  if (any(otu < 0)) stop("Count table contains negative values.")

  if (is.null(rownames(otu))) rownames(otu) <- paste0("S", seq_len(nrow(otu)))
  if (is.null(colnames(otu))) colnames(otu) <- paste0("T", seq_len(ncol(otu)))

  metadata <- as.data.frame(metadata)
  if (nrow(metadata) != nrow(otu)) {
    stop("meta_data must have one row per sample after orientation normalization.")
  }
  # Metadata with automatic row names (for example the colData of a
  # SummarizedExperiment without column names) is matched by position.
  if (.row_names_info(metadata) < 0L) rownames(metadata) <- rownames(otu)
  metadata <- metadata[rownames(otu), , drop = FALSE]

  list(
    y = otu,
    metadata = metadata
  )
}

# Designs of one analysis. x_formula and r_formula are the right-hand sides
# (strings) of the abundance and prevalence formulas; test_terms are the
# tested terms (default: the terms of x_formula). Returns metadata, the
# abundance design x and the prevalence design r without intercept columns,
# and tests, one element per tested column of r (term, parent_term, type,
# coef_name, columns).
ancombc3_make_design <- function(metadata,
                                 x_formula,
                                 test_terms = NULL,
                                 r_formula = x_formula) {
  metadata <- as.data.frame(metadata)
  x <- model.matrix(stats::as.formula(paste("~", x_formula)), data = metadata)
  r_formula_obj <- stats::as.formula(paste("~", r_formula))
  r <- model.matrix(r_formula_obj, data = metadata)
  # Term of each column of r: the index into the term labels of r_formula
  # (0 for the intercept), recorded before the intercept column is removed.
  r_assign <- attr(r, "assign")
  r_term_labels <- attr(stats::terms(r_formula_obj), "term.labels")
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  r_keep <- colnames(r) != "(Intercept)"
  r <- r[, r_keep, drop = FALSE]
  r_assign <- r_assign[r_keep]

  if (is.null(test_terms)) {
    test_terms <- attr(stats::terms(stats::as.formula(paste("~", x_formula))), "term.labels")
  }

  tests <- list()
  for (term in test_terms) {
    value <- metadata[[term]]
    # Columns of the prevalence design generated by the term; a test term that
    # is not a term of r_formula is matched to a column of the same name.
    term_index <- match(term, r_term_labels)
    columns <- if (is.na(term_index)) character(0) else colnames(r)[r_assign == term_index]
    if (length(columns) == 0L && term %in% colnames(r)) columns <- term
    if (length(columns) == 0L) {
      stop("Could not match test term '", term, "' to columns of the prevalence design.")
    }
    term_type <- if (is.numeric(value) && length(unique(value[is.finite(value)])) == 2L) {
      "binary"
    } else if (is.numeric(value)) {
      "continuous"
    } else {
      "categorical"
    }
    for (column in columns) {
      tests[[column]] <- list(
        term = column,
        parent_term = term,
        type = if (term_type == "categorical" && length(columns) > 1L) "binary" else term_type,
        coef_name = column,
        columns = column
      )
    }
  }

  list(
    metadata = metadata,
    x = x,
    r = r,
    tests = tests
  )
}

# Copy of a design matrix with each named column set to its reference value,
# the sample mean.
reference_test_columns <- function(design, columns) {
  reference <- as.matrix(design)
  for (column in columns) reference[, column] <- mean(reference[, column])
  dimnames(reference) <- dimnames(design)
  reference
}

# At-risk-weighted sum of squares of a continuous covariate,
# sum_i w_i (x_i - xbar_w)^2 with xbar_w the w-weighted mean, over the samples
# with a finite covariate value and a positive, finite weight.
weighted_continuous_support <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2L) return(0)
  x <- x[ok]
  w <- w[ok]
  w_mean <- sum(w * x) / sum(w)
  sum(w * (x - w_mean)^2)
}

# DA eligibility of every taxon for one tested term, a list with the elements
# term, type ("continuous", "binary", or a vector contrast otherwise), and
# columns, as in ancombc3_make_design()$tests. at_risk_prob is the n x p matrix
# of first-stage at-risk probabilities pi^full, y the counts, and x the
# abundance design. The taxa are those of y (.ac3_taxon_names(): its column
# names, or T1, ..., Tp when it has none); at_risk_prob is matched to them by
# name when it has column names, and by position otherwise.
# A taxon is eligible when it satisfies three groups of criteria.
#   Fitted at-risk criteria: the at-risk effective sample size
#     sum_i pi_ij^full is at least min_at_risk_eff_n, and the at-risk design
#     support is at least min_at_risk_design_support. The design support is
#     the at-risk-weighted sum of squares of a continuous covariate, the
#     smaller of the two group sums of the at-risk probabilities for a binary
#     covariate, and the smallest eigenvalue of the at-risk-weighted
#     cross-product matrix of the design with an intercept for a vector
#     contrast.
#   Observed-data criteria: for a continuous covariate, at least
#     min_observed_positive positive counts and a covariate variance among the
#     positive samples of at least min_continuous_var_ratio times its
#     full-sample variance; for a binary covariate, at least
#     min_positive_per_group positive counts in each group; for a vector
#     contrast, at least min_observed_positive positive counts and an observed
#     design support (the smallest eigenvalue of the cross-product matrix over
#     the positive samples) of at least min_at_risk_design_support.
#   Positive-count fraction: n^-1 sum_i 1{y_ij > 0} is at least
#     min_observed_prevalence.
# The default thresholds are those of the manuscript (Methods, table of
# eligibility criteria). Returns one row per taxon with the eligibility
# indicator da_eligible and the quantities compared with the thresholds.
# With repeated measures, cluster is the subject of each sample. For a binary
# covariate that is constant within every subject, the group quantities are
# then computed over subjects (.ac3_cluster_binary_support()): the number of
# positive counts of a group is the number of its subjects with at least one
# positive count, and the group sum of the at-risk probabilities is the sum
# over its subjects of the largest at-risk probability of the subject. The
# other criteria are unchanged.
ancombc3_da_eligibility_term <- function(at_risk_prob,
                                         y,
                                         x,
                                         term,
                                         min_at_risk_eff_n = 10,
                                         min_at_risk_design_support = 1,
                                         min_observed_positive = 10,
                                         min_positive_per_group = 5,
                                         min_continuous_var_ratio = 0.25,
                                         min_observed_prevalence = 0.10,
                                         cluster = NULL) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  taxa <- .ac3_taxon_names(y)
  colnames(y) <- taxa
  at_risk_prob <- .ac3_match_taxa(at_risk_prob, taxa, "at_risk_prob")
  if (!identical(dim(at_risk_prob), dim(y))) stop("at_risk_prob and y must have the same n x p dimensions.")

  columns <- term$columns
  if (!all(columns %in% colnames(x))) {
    stop("Term columns not found in x: ", paste(setdiff(columns, colnames(x)), collapse = ", "))
  }

  # At-risk probabilities, restricted to [0, 1].
  q_at_risk <- pmin(pmax(at_risk_prob, 0), 1)
  observed_positive <- y > 0
  n_eff <- colSums(q_at_risk)
  observed_positive_n <- colSums(observed_positive)
  at_risk_support <- numeric(ncol(y))
  observed_support <- numeric(ncol(y))
  min_positive_by_group <- rep(NA_real_, ncol(y))

  if (term$type == "continuous") {
    x_term <- as.numeric(x[, columns[1L]])
    full_var <- stats::var(x_term)
    if (!is.finite(full_var) || full_var <= 0) full_var <- NA_real_
    for (j in seq_len(ncol(y))) {
      at_risk_support[j] <- weighted_continuous_support(x_term, q_at_risk[, j])
      pos <- observed_positive[, j]
      observed_support[j] <- if (sum(pos) >= 3L && is.finite(full_var)) stats::var(x_term[pos]) / full_var else 0
    }
    observed_pass <- observed_positive_n >= min_observed_positive &
      observed_support >= min_continuous_var_ratio
  } else if (term$type == "binary") {
    x_term <- as.numeric(x[, columns[1L]])
    group0 <- x_term <= 0.5
    group1 <- x_term > 0.5
    for (j in seq_len(ncol(y))) {
      at_risk_support[j] <- min(sum(q_at_risk[group0, j]), sum(q_at_risk[group1, j]))
      min_positive_by_group[j] <- min(sum(observed_positive[group0, j]), sum(observed_positive[group1, j]))
      observed_support[j] <- min_positive_by_group[j]
    }
    if (!is.null(cluster) && .ac3_cluster_constant_within(x_term, .ac3_cluster_index(cluster))) {
      by_subject <- .ac3_cluster_binary_support(q_at_risk, observed_positive, group1,
                                                .ac3_cluster_index(cluster))
      at_risk_support <- by_subject$at_risk_support
      min_positive_by_group <- by_subject$min_positive_by_group
      observed_support <- min_positive_by_group
    }
    observed_pass <- min_positive_by_group >= min_positive_per_group
  } else {
    x_term <- cbind("(Intercept)" = 1, x[, columns, drop = FALSE])
    for (j in seq_len(ncol(y))) {
      xw <- x_term * sqrt(q_at_risk[, j])
      at_risk_support[j] <- min(eigen(crossprod(xw), symmetric = TRUE, only.values = TRUE)$values)
      pos <- observed_positive[, j]
      observed_support[j] <- if (sum(pos) > ncol(x_term)) {
        min(eigen(crossprod(x_term[pos, , drop = FALSE]), symmetric = TRUE, only.values = TRUE)$values)
      } else {
        0
      }
    }
    observed_pass <- observed_positive_n >= min_observed_positive &
      observed_support >= min_at_risk_design_support
  }

  observed_prevalence <- colSums(observed_positive) / nrow(y)
  at_risk_pass <- n_eff >= min_at_risk_eff_n &
    at_risk_support >= min_at_risk_design_support
  eligible <- observed_pass & at_risk_pass &
    observed_prevalence >= min_observed_prevalence

  data.frame(
    taxon = taxa,
    term = term$term,
    term_type = term$type,
    da_eligible = as.logical(eligible),
    at_risk_pass = as.logical(at_risk_pass),
    observed_support_pass = as.logical(observed_pass),
    at_risk_eff_n = as.numeric(n_eff),
    at_risk_design_support = as.numeric(at_risk_support),
    observed_positive_n = as.numeric(observed_positive_n),
    observed_design_support = as.numeric(observed_support),
    min_positive_by_group = as.numeric(min_positive_by_group),
    observed_prevalence = as.numeric(observed_prevalence),
    row.names = NULL
  )
}

# DA eligibility of every tested covariate. test_terms is either a character
# vector of columns of x or a list of term specifications (elements term, type,
# and columns, as in ancombc3_make_design()$tests). A column given by name is a
# binary term when it takes exactly the two values 0 and 1, and a continuous
# term otherwise. The further arguments (...) are the thresholds of
# ancombc3_da_eligibility_term() and its argument cluster. Returns
#   by_term: the names of the eligible taxa of each tested covariate (the
#     taxa of y, matched as in ancombc3_da_eligibility_term()), in the column
#     order of y, named by the names of test_terms (by term$term for an
#     unnamed specification);
#   fitting_set: the union of by_term, the taxa used to estimate the
#     sample-specific bias (the argument fitting_set of fit_ancombc3_da());
#   diagnostics: the tables of ancombc3_da_eligibility_term(), stacked by row.
ancombc3_da_eligibility <- function(at_risk_prob, y, x, test_terms, ...) {
  x <- as.matrix(x)
  if (is.character(test_terms)) {
    test_terms <- lapply(stats::setNames(test_terms, test_terms), function(column) {
      if (!column %in% colnames(x)) stop("Tested column not found in x: ", column)
      values <- unique(x[is.finite(x[, column]), column])
      type <- if (length(values) == 2L && all(values %in% c(0, 1))) "binary" else "continuous"
      list(term = column, type = type, coef_name = column, columns = column)
    })
  }
  if (!is.list(test_terms) || length(test_terms) == 0L) {
    stop("test_terms must identify at least one tested term.")
  }
  term_names <- names(test_terms) %||% rep("", length(test_terms))
  unnamed <- !nzchar(term_names)
  term_names[unnamed] <- vapply(test_terms[unnamed], function(term) as.character(term$term), character(1))
  by_term <- stats::setNames(vector("list", length(test_terms)), term_names)
  diagnostics <- vector("list", length(test_terms))
  for (k in seq_along(test_terms)) {
    elig <- ancombc3_da_eligibility_term(at_risk_prob = at_risk_prob, y = y, x = x,
                                         term = test_terms[[k]], ...)
    by_term[[k]] <- elig$taxon[elig$da_eligible]
    diagnostics[[k]] <- elig
  }
  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL
  list(
    by_term = by_term,
    fitting_set = unique(unlist(by_term, use.names = FALSE)),
    diagnostics = diagnostics
  )
}
