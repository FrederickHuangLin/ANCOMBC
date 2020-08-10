para_est = function(y, meta_data, formula, tol, max_iter) {
    options(na.action = "na.pass") # Keep NA's in rows of x
    x = model.matrix(formula(paste0("~", formula)), data = meta_data)
    options(na.action = "na.omit") # Switch it back
    taxa_id = rownames(y); n_taxa = nrow(y)
    samp_id = colnames(y); n_samp = ncol(y)
    covariates = colnames(x); n_covariates = length(covariates)

    # Sampling fractions
    d = rep(0, n_samp)
    tformula = formula(paste0("y ~ ", formula))
    fits = lapply(seq_len(n_taxa), function(i) {
        df = data.frame(y = unlist(y[i, ]) - d, meta_data)
        return(lm(tformula, data = df))
    })
    # Regression coefficients
    beta = lapply(fits, function(i) {
        beta_i = rep(NA, length(covariates)) # prevent errors of missing values
        coef_i = coef(i)
        beta_i[match(names(coef_i), covariates)] = coef_i
        return(beta_i)
    })
    beta = Reduce('rbind', beta)

    # Iterative least square
    iterNum = 0; epsilon = 100
    while (epsilon > tol & iterNum < max_iter) {
        # Updating beta
        fits = lapply(seq_len(n_taxa), function(i) {
            df = data.frame(y = unlist(y[i, ]) - d, meta_data)
            return(lm(tformula, data = df))
        })
        beta_new = lapply(fits, function(i) {
            beta_i = rep(NA, length(covariates))
            coef_i = coef(i)
            beta_i[match(names(coef_i), covariates)] = coef_i
            return(beta_i)
        })
        beta_new = Reduce('rbind', beta_new)

        # Updating d
        y_hat = lapply(fits, function(i) {
            y_hat_i = rep(NA, n_samp)
            fit_i = fitted(i)
            y_hat_i[match(names(fit_i), samp_id)] = fit_i
            return(y_hat_i)

        })
        y_hat = Reduce('rbind', y_hat)
        d_new = colMeans(y - y_hat, na.rm = TRUE)

        # Iteration
        epsilon = sqrt(sum((beta_new - beta)^2, na.rm = TRUE) +
                           sum((d_new - d)^2, na.rm = TRUE))
        iterNum = iterNum + 1
        beta = beta_new; d = d_new
    }

    # Regression residuals
    y_hat = lapply(fits, function(i) {
        y_hat_i = rep(NA, n_samp)
        fit_i = fitted(i)
        y_hat_i[match(names(fit_i), samp_id)] = fit_i
        return(y_hat_i)

    })
    y_hat = Reduce('rbind', y_hat)
    e = t(t(y - y_hat) - d)

    # Variance-covariance matrices of coefficients
    fiuo_var_cov = var_cov_est(x, e, n_taxa)
    var_cov_hat = fiuo_var_cov$var_cov_hat; var_hat = fiuo_var_cov$var_hat

    colnames(beta) = covariates; rownames(beta) = taxa_id
    names(d) = samp_id
    names(var_cov_hat) = taxa_id
    colnames(var_hat) = covariates; rownames(var_hat) = taxa_id

    fiuo_para = list(beta = beta, d = d, e = e,
                     var_cov_hat = var_cov_hat, var_hat = var_hat)
    return(fiuo_para)
}
