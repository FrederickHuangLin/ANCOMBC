# Comma-separated list of names for a message. A list longer than max_n is
# shown by its first max_n names and the number of the remaining names.
.name_list = function(x, max_n = 10) {
    x = as.character(x)
    n = length(x)
    if (n <= max_n) return(paste(x, collapse = ", "))
    paste0(paste(x[seq_len(max_n)], collapse = ", "),
           ", ... and ", n - max_n, " more")
}

# The function to extract off-diagonol elements
.odiag = function(x) x[col(x) != row(x)]

# The function to calculate the variance of difference
.var_diff = function(x) {
    sum(diag(x)) - sum(.odiag(x))
}

# The function of combination
.combn_fun = function(x, fun, sep) {
    y = c(x, utils::combn(x, 2, FUN = fun))
    combn_mat = utils::combn(names(x), 2)
    combn_name = paste(combn_mat[2, ], combn_mat[1, ], sep = sep)
    names(y) = c(names(x), combn_name)
    return(y)
}

.combn_fun2 = function(x, fun, sep) {
    combn_mat = utils::combn(colnames(x), 2)
    y = numeric(ncol(combn_mat))
    for (i in seq(ncol(combn_mat))) {
        idx = c(combn_mat[2, i], combn_mat[1, i])
        y[i] = fun(x[idx, idx])
    }
    y = c(diag(x), y)
    combn_name = paste(combn_mat[2, ], combn_mat[1, ], sep = sep)
    names(y) = c(colnames(x), combn_name)
    return(y)
}

# The mdFDR correction
.mdfdr = function(global_test = c("pairwise", "dunnet"),
                  W, dof, fwer_ctrl_method, ...) {

    input_list = list(...)

    # The total number of null hypotheses rejected in the global test
    if (global_test == "pairwise") {
      if (is.null(input_list$rand_formula)) {
        res_screen = .ancombc_global_F(x = input_list$x,
                                       group = input_list$group,
                                       beta_hat = input_list$beta_hat,
                                       vcov_hat = input_list$vcov_hat,
                                       dof = dof,
                                       p_adj_method = "BH",
                                       alpha = input_list$alpha,
                                       global_p = input_list$global_p)
      } else {
        res_screen = .ancombc_global_LRT(full_model = input_list$full_model,
                                         fix_formula = input_list$fix_formula,
                                         rand_formula = input_list$rand_formula,
                                         control = input_list$control,
                                         x = input_list$x,
                                         group = input_list$group,
                                         y = input_list$y,
                                         meta_data = input_list$meta_data,
                                         p_adj_method = "BH",
                                         alpha = input_list$alpha)
      }
    } else {
        res_screen = .dunn_global(x = input_list$x, group = input_list$group,
                                  W = W,
                                  B = input_list$B,
                                  dof = dof,
                                  p_adj_method = "BH",
                                  alpha = input_list$alpha,
                                  dunnet_null = input_list$dunnet_null)
    }
    R = sum(res_screen$diff_abn)

    # P-values for pairwise tests
    p_val = 2 * (pt(abs(W), df = dof, lower.tail = FALSE))

    # Only consider R significant taxa with regards to the global test
    screen_ind = res_screen$diff_abn
    p_val = p_val * screen_ind
    p_val[p_val == 0] = 1
    p_val[is.na(p_val)] = 1

    # Adjust pairwise p-values at level of R * alpha / d
    n_tax = nrow(W)
    q_val = t(apply(p_val, 1, function(x)
        p.adjust(x, method = fwer_ctrl_method, n = length(x) * n_tax / R)))

    output = list(p_val = p_val, q_val = q_val)
    return(output)
}

# Estimate coefficients under constraints
.safe_inverse_spd = function(A, ridge = 1e-8) {
    A <- (A + t(A)) / 2

    tryCatch(
        chol2inv(chol(A)),
        error = function(e) chol2inv(chol(A + ridge * diag(nrow(A))))
    )
}

.constrain_est = function(beta_hat, vcov_hat, contrast) {
    # A taxon with a coefficient that is not estimable has no constrained
    # estimate
    if (anyNA(beta_hat) || anyNA(vcov_hat)) return(rep(NA_real_, length(beta_hat)))

    P = .safe_inverse_spd(vcov_hat)
    
    Dmat = 2 * P
    Dmat = (Dmat + t(Dmat)) / 2
    diag(Dmat) = diag(Dmat) + 1e-10

    dvec = as.numeric(2 * P %*% beta_hat)
    Amat = t(contrast)
    bvec = rep(0, nrow(contrast))

    fit = try(
        quadprog::solve.QP(
            Dmat = Dmat,
            dvec = dvec,
            Amat = Amat,
            bvec = bvec,
            meq = 0
        ),
        silent = TRUE
    )

    if (inherits(fit, "try-error")) {
        rep(0, length(beta_hat))
    } else {
        as.numeric(fit$solution)
    }
}

# Compute the l_infty norm for a pattern
.l_infty = function(beta_opt, node) {
    l = max(abs(beta_opt[node]),
            abs(beta_opt[node] - beta_opt[length(beta_opt)]),
            na.rm = TRUE)
    return(l)
}

# Generate random variables from the poisson log-normal distribution
.rplnm = function(mu, sigma, n, N) {
    d = length(mu)
    y = MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
    x = N * exp(y)
    otu_table = matrix(rpois(n = n * d, lambda = x), nrow = n)
    return(otu_table)
}

# Regularize eigenvalues
.regularize_eigenvalues = function(mat){

  # Decomposition
  decomp = eigen(mat)

  # Replace zero and negative eigenvalues with the minimum positive eigenvalue
  decomp$values[decomp$values <= 0] = min(decomp$values[decomp$values > 0])

  # Ratio of eigenvalues to the largest eigenvalue
  ratios = decomp$values[1]/decomp$values

  # Difference between consecutive ratios
  diff_ratio = vapply(2:length(ratios),
                      FUN = function(x){ ratios[x]-ratios[x-1] },
                      FUN.VALUE = numeric(1))

  # Find the eigenvalue that corresponds to the maximum difference (the largest gap)
  delta = decomp$values[which(diff_ratio == max(diff_ratio,na.rm = TRUE))]

  # Delta must be between 0.01 -1
  i = 1
  while(delta <0.01){
    delta = decomp$values[which(diff_ratio == max(diff_ratio, na.rm = TRUE)) - i]
    i= i + 1
  }
  if(delta > 1)
    delta = 1

  # Winsorize eigenvalues to the delta
  decomp$values[decomp$values < delta] = delta

  # Reconstruct the covariance matrix
  cov_mat_pos = decomp$vectors %*% diag(decomp$values) %*% solve(decomp$vectors)
  rownames(cov_mat_pos) = rownames(mat)
  colnames(cov_mat_pos) = colnames(mat)

  # Make the matrix symmetric
  cov_mat_pos[lower.tri(cov_mat_pos)] = t(cov_mat_pos)[lower.tri(cov_mat_pos)]

  return(cov_mat_pos)
}

# Check if a matrix is positive semi-definite
.is_psd = function(matrix) {
    # Get eigenvalues
    eigenvals <- eigen(matrix, symmetric = TRUE)$values
    # Check if all eigenvalues are non-negative (within numerical precision)
    all(eigenvals > -1e-10)
}


# Validate a method_control list. Entries that are not supplied take the values
# in defaults; an entry that is not in defaults is an error.
.check_method_control = function(method_control, defaults) {
    if (is.null(method_control)) method_control = list()
    if (!is.list(method_control) ||
        (length(method_control) > 0 && is.null(names(method_control)))) {
        stop("`method_control` should be a named list", call. = FALSE)
    }
    unknown = setdiff(names(method_control), names(defaults))
    if (length(unknown) > 0) {
        stop("Unknown entries in `method_control`: ",
             paste(unknown, collapse = ", "),
             "\nThe available entries are: ",
             paste(names(defaults), collapse = ", "), call. = FALSE)
    }
    out = defaults
    out[names(method_control)] = method_control
    choices = list(vcov = c("pooled", "observed"),
                   global_p = c("two_sided", "upper"),
                   failed_taxa = c("include", "exclude"),
                   dunnet_null = c("all", "estimable"))
    for (nm in intersect(names(out), names(choices))) {
        v = out[[nm]]
        if (!is.character(v) || length(v) != 1 || !v %in% choices[[nm]]) {
            stop(sprintf("`method_control$%s` should be one of: %s", nm,
                         paste0("\"", choices[[nm]], "\"", collapse = ", ")),
                 call. = FALSE)
        }
    }
    if ("sens_parallel" %in% names(out)) {
        v = out$sens_parallel
        if (!is.logical(v) || length(v) != 1 || is.na(v)) {
            stop("`method_control$sens_parallel` should be TRUE or FALSE",
                 call. = FALSE)
        }
    }
    out
}

# Logical index of model-matrix columns belonging to the main effect of group.
# Uses the formula terms and assign attribute on model.matrix() output when
# available so unrelated covariates whose names contain `group` and interaction
# terms are excluded. Falls back to substring matching for plain matrices.
.group_coef_ind = function(x, group) {
    assign_idx = attr(x, "assign")
    tt = attr(x, "terms")
    if (!is.null(assign_idx) && !is.null(tt)) {
        term_labels = attr(tt, "term.labels")
        term_order = attr(tt, "order")
        is_group_term = term_order == 1L & vapply(term_labels, function(v) {
            vars = all.vars(str2lang(v))
            length(vars) == 1L && identical(vars, group)
        }, FUN.VALUE = logical(1))
        return(assign_idx %in% which(is_group_term))
    }
    covariates = colnames(x)
    grepl(group, covariates, fixed = TRUE) & !grepl(":", covariates, fixed = TRUE)
}

# Columns of the pairwise comparison results selected by pairs, a list of
# pairs of levels of the group variable. Each comparison is the later level in
# levels(group) minus the earlier level, as in the full set of pairwise
# comparisons, whose columns are the group coefficients (level versus the
# reference level) followed by their differences "<later>_<earlier>".
.pair_columns = function(x, group, group_levels, pairs) {
    if (is.null(pairs)) return(NULL)
    covariates = colnames(x)
    group_ind = .group_coef_ind(x, group)
    coef_name = covariates[group_ind]
    names(coef_name) = group_levels[-1]
    unique(vapply(pairs, function(pr) {
        idx = sort(match(pr, group_levels))
        if (idx[1] == 1L) {
            coef_name[[group_levels[idx[2]]]]
        } else {
            paste(coef_name[[group_levels[idx[2]]]],
                  coef_name[[group_levels[idx[1]]]], sep = "_")
        }
    }, FUN.VALUE = character(1)))
}
