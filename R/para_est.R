para_est = function(x, y, meta_data, formula, tol, max_iter) {
    taxa_id = rownames(y)
    n_taxa = nrow(y)
    samp_id = colnames(y)
    n_samp = ncol(y)
    covariates = colnames(x)
    n_covariates = length(covariates)
    tformula = formula(paste0("y ~ ", formula))

    # Test for over-parameterization
    lm_smoke = lm(tformula, data = data.frame(y = unlist(y[1, ]), meta_data))

    if (any(is.na(lm_smoke$coefficients))) {
        stop_txt = sprintf(paste("Estimation failed for the following covariates:",
                                 paste(names(which(is.na(lm_smoke$coefficients))), collapse = ", "),
                                 "Consider removing these covariates",
                                 sep = "\n"))
        stop(stop_txt, call. = FALSE)
    }

    if (lm_smoke$df.residual == 0) {
        stop_txt = sprintf(paste("No residual degrees of freedom! The model is over-parameterized",
                                 "A more parsimonious model is needed",
                                 sep = "\n"))
        stop(stop_txt, call. = FALSE)
    }

    # Sampling fractions
    d = rep(0, n_samp)
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
    iterNum = 0
    epsilon = 100
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
        beta = beta_new
        d = d_new
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

    colnames(beta) = covariates
    rownames(beta) = taxa_id
    names(d) = samp_id
    names(var_cov_hat) = taxa_id
    colnames(var_hat) = covariates
    rownames(var_hat) = taxa_id

    fiuo_para = list(beta = beta, d = d, e = e,
                     var_cov_hat = var_cov_hat, var_hat = var_hat)
    return(fiuo_para)
}
