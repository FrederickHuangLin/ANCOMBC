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
.safe_inverse_spd = function(A, ridge = 1e-8) {
    A <- (A + t(A)) / 2

    tryCatch(
        chol2inv(chol(A)),
        error = function(e) chol2inv(chol(A + ridge * diag(nrow(A))))
    )
}

.constrain_est = function(beta_hat, vcov_hat, contrast, solver) {
    beta_hat_mat = matrix(beta_hat, ncol = 1)
    vcov_hat_inv = .safe_inverse_spd(vcov_hat)

    beta_opt = CVXR::Variable(shape = c(length(beta_hat), 1),
                              name = "beta")
    obj = CVXR::Minimize(CVXR::quad_form(beta_opt - beta_hat_mat, vcov_hat_inv))
    cons = suppressMessages(contrast %*% beta_opt >= 0)
    problem = CVXR::Problem(objective = obj, constraints = list(cons))

    suppressMessages(opt_val <- try(CVXR::psolve(problem, solver = solver),
                                    silent = TRUE))

    if (inherits(opt_val, "try-error")) {
        beta_opt_val = rep(0, length(beta_hat))
    } else {
        beta_opt_val = as.numeric(CVXR::value(beta_opt))
    }
    return(beta_opt_val)
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

