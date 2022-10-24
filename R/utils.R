# The function to extract off-diagonol elements
odiag = function(x) x[col(x) != row(x)]

# The function to calculate the variance of difference
var_diff = function(x) {
    sum(diag(x)) - sum(odiag(x))
}

# The function of combination
combn_fun = function(x, fun, sep) {
    y = c(x, utils::combn(x, 2, FUN = fun))
    combn_mat = utils::combn(names(x), 2)
    combn_name = paste(combn_mat[2, ], combn_mat[1, ], sep = sep)
    names(y) = c(names(x), combn_name)
    return(y)
}

combn_fun2 = function(x, fun, sep) {
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
mdfdr = function(global_test = c("pairwise", "dunnet"),
                 W, fwer_ctrl_method, ...) {

    input_list = list(...)

    # The total number of null hypotheses rejected in the global test
    if (global_test == "pairwise") {
        res_screen = ancombc_global(x = input_list$x, group = input_list$group,
                                    beta_hat = input_list$beta_hat,
                                    vcov_hat = input_list$vcov_hat,
                                    p_adj_method = "BH",
                                    alpha = input_list$alpha)
    } else {
        res_screen = dunn_global(x = input_list$x, group = input_list$group,
                                 W = W,
                                 B = input_list$B, p_adj_method = "BH",
                                 alpha = input_list$alpha)
    }
    R = sum(res_screen$diff_abn)

    # P-values for pairwise tests
    p_val = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)

    # Only consider R significant taxa with regards to the global test
    screen_ind = res_screen$diff_abn
    p_val = p_val * screen_ind
    p_val[p_val == 0] = 1

    # Adjust pairwise p-values at level of R * alpha / d
    n_tax = nrow(W)
    q_val = t(apply(p_val, 1, function(x)
        p.adjust(x, method = fwer_ctrl_method, n = length(x) * n_tax / R)))

    output = list(p_val = p_val, q_val = q_val)
    return(output)
}

# Estimate coefficients under constraints
constrain_est = function(beta_hat, vcov_hat, contrast, solver) {
    beta_opt = CVXR::Variable(rows = length(beta_hat), cols = 1, name = "beta")
    obj = CVXR::Minimize(CVXR::matrix_frac(beta_opt - beta_hat, vcov_hat))
    cons = contrast %*% beta_opt >= 0
    problem = CVXR::Problem(objective = obj, constraints = list(cons))

    suppressWarnings(result <- try(CVXR::solve(problem, solver = solver),
                                   silent = TRUE))

    if (inherits(result, "try-error")) {
        beta_opt = rep(0, length(beta_hat))
    } else {
        beta_opt = as.numeric(result$getValue(beta_opt))
    }
    return(beta_opt)
}

# Compute the l_infty norm for a pattern
l_infty = function(beta_opt, node) {
    l = max(abs(beta_opt[node]),
            abs(beta_opt[node] - beta_opt[length(beta_opt)]),
            na.rm = TRUE)
    return(l)
}

# Generate random variables from the poisson log-normal distribution
rplnm = function(mu, sigma, n, N) {
    d = length(mu)
    y = MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
    x = N * exp(y)
    otu_table = matrix(rpois(n = n * d, lambda = x), nrow = n)
    return(otu_table)
}






