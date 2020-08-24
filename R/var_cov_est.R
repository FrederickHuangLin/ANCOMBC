var_cov_est = function(x, e, n_taxa) {
    covariates = colnames(x)
    n_covariates = length(covariates)
    n_samp = nrow(x)
    XTX_inv = MASS::ginv(t(x[complete.cases(x), ]) %*% x[complete.cases(x), ])
    var_cov_hat = vector(mode = "list", length = n_taxa) # Covariances
    var_hat = matrix(NA, nrow = n_taxa, ncol = n_covariates) # Variances
    for (i in seq_len(n_taxa)) {
        sigma2_xxT = matrix(0, ncol = n_covariates, nrow = n_covariates)
        for (j in seq_len(n_samp)) {
            sigma2_xxT_j = e[i, j]^2 * x[j, ] %*% t(x[j, ])
            sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0
            sigma2_xxT = sigma2_xxT + sigma2_xxT_j
        }
        var_cov_hat[[i]] = XTX_inv %*% sigma2_xxT %*% XTX_inv
        rownames(var_cov_hat[[i]]) = covariates
        colnames(var_cov_hat[[i]]) = covariates
        var_hat[i, ] = diag(var_cov_hat[[i]])
    }
    fiuo_var_cov = list(var_cov_hat = var_cov_hat, var_hat = var_hat)
    return(fiuo_var_cov)
}
