fit_summary = function(y, x, beta, var_hat, delta_em, var_delta, conserve) {
    n_taxa = nrow(y)

    beta_hat = beta
    beta_hat[, -1] = t(t(beta_hat[, -1]) - delta_em)

    if (conserve) {
        # Account for the variance of delta_hat
        se_hat = sqrt(sweep(var_hat, 2, c(0, var_delta), "+") +
                          2 * sqrt(sweep(var_hat, 2, c(0, var_delta), "*")))
    }else{ se_hat = sqrt(var_hat) }

    d_hat = matrix(NA, nrow = nrow(y), ncol = ncol(y))
    for (i in seq_len(n_taxa)) {
        d_hat[i, ] = y[i, ] - x %*% beta_hat[i, ]
    }
    d_hat = colMeans(d_hat, na.rm = TRUE)

    # Remove uninformative intercept column
    beta_hat = beta_hat[, setdiff(colnames(beta_hat), "(Intercept)"),
                        drop = FALSE]
    se_hat = se_hat[, setdiff(colnames(se_hat), "(Intercept)"),
                    drop = FALSE]

    fiuo_fit = list(beta_hat = beta_hat, se_hat = se_hat, d_hat = d_hat)
    return(fiuo_fit)
}
