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
    y = vector(mode = "numeric")
    for (i in seq(ncol(combn_mat))) {
        idx = c(combn_mat[2, i], combn_mat[1, i])
        y = c(y, fun(x[idx, idx]))
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
                                       alpha = input_list$alpha)
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
                                  alpha = input_list$alpha)
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
.constrain_est = function(beta_hat, vcov_hat, contrast, solver) {
    beta_opt = CVXR::Variable(rows = length(beta_hat), cols = 1, name = "beta")
    obj = CVXR::Minimize(CVXR::matrix_frac(beta_opt - beta_hat, vcov_hat))
    cons = suppressMessages(contrast %*% beta_opt >= 0)
    problem = CVXR::Problem(objective = obj, constraints = list(cons))

    suppressMessages(result <- try(CVXR::solve(problem, solver = solver),
                                   silent = TRUE))

    if (inherits(result, "try-error")) {
        beta_opt = rep(0, length(beta_hat))
    } else {
        beta_opt = as.numeric(result$getValue(beta_opt))
    }
    return(beta_opt)
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

# Get the p-values for the sensitivity analysis
.get_p = function(y, data, formula, group, n_levels, pairwise, global, trend) {
  tformula = paste0("y ~ ", formula)
  df = data.frame(y = y, data)
  lm_fit = stats::lm(formula(tformula), data = df)
  summ = summary(lm_fit)
  p_val = summ$coefficients[, "Pr(>|t|)"]
  p_val[p_val == 0] = 2e-16
  names(p_val) = rownames(summ$coefficients)

  if (pairwise) {
    mcp_arg = paste0(group, ' = "Tukey"')
    comparison = multcomp::glht(lm_fit, linfct = eval(parse(text = paste0("multcomp::mcp(", mcp_arg, ")"))))
    summ = summary(comparison, test = multcomp::adjusted("none"))
    pair_p_val = summ$test$pvalues
    pair_p_val[pair_p_val == 0] = 2e-16
    names(pair_p_val) = paste0(summ$focus, names(pair_p_val))
    pair_p_val = pair_p_val[-(seq_len(n_levels - 1))]
    p_val = c(p_val, pair_p_val)
  }

  if (global | trend) {
    anova_fit = anova(lm_fit)
    group_p_val = anova_fit$`Pr(>F)`[grepl(group, rownames(anova_fit))]
    if (group_p_val == 0) group_p_val = 2e-16
    if (global) p_val = c(p_val, global = group_p_val)
    if (trend) p_val = c(p_val, trend = group_p_val)
  }

  return(p_val)
}

# Internal wrappers for mia::agglomerateByRank/mergeRows
.merge_features = function(x, merge.by, ...) {
    # Check if merge.by parameter belongs to taxonomyRanks
    if (is.character(merge.by) && length(merge.by) == 1 && merge.by %in% mia::taxonomyRanks(x)) {
        # Merge using agglomerateByRank
        x = mia::agglomerateByRank(x, rank = merge.by, ...)
    } else {
        # Merge using mia::mergeRows
        f = factor(SummarizedExperiment::rowData(x)[, merge.by])
        x = mia::mergeRows(x, f = f, ...)
    }
    return(x)
}


