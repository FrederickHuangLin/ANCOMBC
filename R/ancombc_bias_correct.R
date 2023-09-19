# Iterative MLE
.iter_mle = function(x, y, meta_data, formula, theta = NULL,
                     tol, max_iter, verbose = FALSE) {
    tax_id = rownames(y)
    n_tax = nrow(y)
    samp_id = colnames(y)
    n_samp = ncol(y)
    fix_eff = colnames(x)
    n_fix_eff = length(fix_eff)
    tformula = formula(paste0("y_crt ~ ", formula))

    # Test for over-parameterization
    lm_smoke = stats::lm(formula = tformula,
                         data = data.frame(y_crt = rnorm(n = n_samp), meta_data))

    if (any(is.na(lm_smoke$coefficients))) {
        stop_txt = sprintf(paste("Estimation failed for the following covariates:",
                                 paste(names(which(is.na(lm_smoke$coefficients))), collapse = ", "),
                                 "Please ensure that these covariates do not have missing values and check for multicollinearity before re-estimating the model",
                                 sep = "\n"))
        stop(stop_txt, call. = FALSE)
    }

    if (lm_smoke$df.residual == 0) {
        stop_txt = sprintf(paste("No residual degrees of freedom! The model is over-parameterized",
                                 "Please consider a more parsimonious model",
                                 sep = "\n"))
        stop(stop_txt, call. = FALSE)
    }

    if (is.null(theta)) {
        # Sampling fractions
        theta = rep(0, n_samp)

        # ML fits
        fits = lapply(seq_len(n_tax), function(i) {
            df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
            suppressWarnings(fit <- try(stats::lm(tformula, data = df),
                                        silent = TRUE))
            if (inherits(fit, "try-error")) {fit = NA}
            return(fit)
        })

        # Degree of freedom
        dof = NULL

        # Degree of freedom
        dof = NULL

        # Degree of freedom
        dof = NULL

        # Coefficients
        empty_coef = rep(NA, n_fix_eff)
        names(empty_coef) = fix_eff
        beta = lapply(fits, function(i) {
            beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
            coef_i = if (inherits(i, "lm")) {
              stats::coef(i)
            } else {
              empty_coef
            }
            beta_i[match(names(coef_i), fix_eff)] = coef_i
            return(beta_i)
        })
        beta = do.call("rbind", beta)

        # Iterative least square
        iterNum = 0
        epsilon = 100
        empty_fitted = rep(NA, n_samp)
        names(empty_fitted) = samp_id
        while (epsilon > tol & iterNum < max_iter) {
            # Updating beta
            fits = lapply(seq_len(n_tax), function(i) {
                df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
                suppressWarnings(fit <- try(stats::lm(tformula, data = df),
                                            silent = TRUE))
                if (inherits(fit, "try-error")) {fit = NA}
                return(fit)
            })

            beta_new = lapply(fits, function(i) {
                beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
                coef_i = if (inherits(i, "lm")) {
                  stats::coef(i)
                } else {
                  empty_coef
                }
                beta_i[match(names(coef_i), fix_eff)] = coef_i
                return(beta_i)
            })
            beta_new = do.call("rbind", beta_new)

            # Updating theta
            y_crt_hat = lapply(fits, function(i) {
                y_crt_hat_i = rep(0, n_samp)
                fitted_i = if (inherits(i, "lm")) {
                  stats::fitted(i)
                } else {
                  empty_fitted
                }
                y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
                return(y_crt_hat_i)
            })
            y_crt_hat = do.call("rbind", y_crt_hat)
            theta_new = colMeans(y - y_crt_hat, na.rm = TRUE)

            # Iteration
            epsilon = sqrt(sum((beta_new - beta)^2, na.rm = TRUE) +
                               sum((theta_new - theta)^2, na.rm = TRUE))
            iterNum = iterNum + 1
            beta = beta_new
            theta = theta_new

            if (verbose) {
                txt = sprintf(paste0("ML iteration = ", iterNum,
                                     ", epsilon = ", signif(epsilon, 2)))
                message(txt)
            }
        }

        # Variance-covariance matrices
        y_crt_hat = lapply(fits, function(i) {
          y_crt_hat_i = rep(0, n_samp)
          fitted_i = if (inherits(i, "lm")) {
            stats::fitted(i)
          } else {
            empty_fitted
          }
          y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
          return(y_crt_hat_i)
        })
        y_crt_hat = do.call("rbind", y_crt_hat)
        eps = t(t(y - y_crt_hat) - theta)

        XTX_inv = MASS::ginv(t(x[complete.cases(x), ]) %*% x[complete.cases(x), ])
        vcov_hat = vector(mode = "list", length = n_tax)
        var_hat = matrix(NA, nrow = n_tax, ncol = n_fix_eff)
        for (i in seq_len(n_tax)) {
          sigma2_xxT = matrix(0, ncol = n_fix_eff, nrow = n_fix_eff)
          for (j in seq_len(n_samp)) {
            sigma2_xxT_j = eps[i, j]^2 * x[j, ] %*% t(x[j, ])
            sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0.1
            sigma2_xxT = sigma2_xxT + sigma2_xxT_j
          }
          vcov_hat[[i]] = XTX_inv %*% sigma2_xxT %*% XTX_inv
          rownames(vcov_hat[[i]]) = fix_eff
          colnames(vcov_hat[[i]]) = fix_eff
          var_hat[i, ] = diag(vcov_hat[[i]])
        }
    } else {
        # ML fits
        fits = lapply(seq_len(n_tax), function(i) {
            df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
            suppressWarnings(fit <- try(stats::lm(tformula, data = df),
                                        silent = TRUE))
            if (inherits(fit, "try-error")) {fit = NA}
            return(fit)
        })

        # Degree of freedom
        dof = vapply(fits, function(i) {
          if (inherits(i, "lm")) {
            summary(i)$df[2]
          } else {
            999L
          }
        }, FUN.VALUE = integer(1))
        dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)

        # Degree of freedom
        dof = vapply(fits, function(i) {
          if (inherits(i, "lm")) {
            summary(i)$df[2]
          } else {
            999L
          }
        }, FUN.VALUE = integer(1))
        dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)

        # Degree of freedom
        dof = vapply(fits, function(i) {
          if (inherits(i, "lm")) {
            summary(i)$df[2]
          } else {
            999L
          }
        }, FUN.VALUE = integer(1))
        dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)

        # Coefficients
        empty_coef = rep(NA, n_fix_eff)
        names(empty_coef) = fix_eff

        beta = lapply(fits, function(i) {
            beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
            coef_i = if (inherits(i, "lm")) {
              stats::coef(i)
            } else {
              empty_coef
            }
            beta_i[match(names(coef_i), fix_eff)] = coef_i
            return(beta_i)
        })
        beta = do.call("rbind", beta)

        # Variance-covariance matrices
        empty_fitted = rep(NA, n_samp)
        names(empty_fitted) = samp_id

        y_crt_hat = lapply(fits, function(i) {
          y_crt_hat_i = rep(0, n_samp)
          fitted_i = if (inherits(i, "lm")) {
            stats::fitted(i)
          } else {
            empty_fitted
          }
          y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
          return(y_crt_hat_i)
        })
        y_crt_hat = do.call("rbind", y_crt_hat)
        eps = t(t(y - y_crt_hat) - theta)

        XTX_inv = MASS::ginv(t(x[complete.cases(x), ]) %*% x[complete.cases(x), ])
        vcov_hat = vector(mode = "list", length = n_tax)
        var_hat = matrix(NA, nrow = n_tax, ncol = n_fix_eff)
        for (i in seq_len(n_tax)) {
          sigma2_xxT = matrix(0, ncol = n_fix_eff, nrow = n_fix_eff)
          for (j in seq_len(n_samp)) {
            sigma2_xxT_j = eps[i, j]^2 * x[j, ] %*% t(x[j, ])
            sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0.1
            sigma2_xxT = sigma2_xxT + sigma2_xxT_j
          }
          vcov_hat[[i]] = XTX_inv %*% sigma2_xxT %*% XTX_inv
          rownames(vcov_hat[[i]]) = fix_eff
          colnames(vcov_hat[[i]]) = fix_eff
          var_hat[i, ] = diag(vcov_hat[[i]])
        }
    }

    if (!is.null(dof)) {
      colnames(dof) = fix_eff
      rownames(dof) = tax_id
    }
    colnames(beta) = fix_eff
    rownames(beta) = tax_id
    names(theta) = samp_id
    names(vcov_hat) = tax_id
    colnames(var_hat) = fix_eff
    rownames(var_hat) = tax_id

    output = list(dof = dof, beta = beta, theta = theta,
                  vcov_hat = vcov_hat, var_hat = var_hat)
    return(output)
}

# Iterative REML
.iter_remle = function(x, y, meta_data, fix_formula, rand_formula,
                       lme_control = lme_control, theta = NULL,
                       tol, max_iter, verbose = FALSE) {
    tax_id = rownames(y)
    n_tax = nrow(y)
    samp_id = colnames(y)
    n_samp = ncol(y)
    fix_eff = colnames(x)
    n_fix_eff = length(fix_eff)
    tformula = formula(paste0("y_crt ~ ", fix_formula, "+ ", rand_formula))

    # Test for over-parameterization
    lm_smoke = stats::lm(formula = formula(paste0("y ~ ", fix_formula)),
                         data = data.frame(y = rnorm(n = n_samp), meta_data))

    if (any(is.na(lm_smoke$coefficients))) {
      stop_txt = sprintf(paste("Estimation failed for the following covariates:",
                               paste(names(which(is.na(lm_smoke$coefficients))), collapse = ", "),
                               "Please ensure that these covariates do not have missing values and check for multicollinearity before re-estimating the model",
                               sep = "\n"))
      stop(stop_txt, call. = FALSE)
    }

    if (lm_smoke$df.residual == 0) {
      stop_txt = sprintf(paste("No residual degrees of freedom! The model is over-parameterized",
                               "Please consider a more parsimonious model",
                               sep = "\n"))
      stop(stop_txt, call. = FALSE)
    }

    # Test for the fitting of linear mixed-effects model
    tryCatch({
        # Try to run the lmerTest model
        result <- lmerTest::lmer(formula = tformula,
                                 data = data.frame(y_crt = y[1, ], meta_data),
                                 control = lme_control)
    },
    error = function(e) {
        # This block will be executed if there's an error in the above code
        message <- sprintf(paste("Encountering the error for `lmerTest` package.",
                                 "Please try to select one of your taxa and use its raw counts to fix the same linear mixed-effects model using `lmerTest` without the `ANCOMBC` package.",
                                 "Load all necessary packages EXCEPT `ANCOMBC`, and see if the error arises due to package incompatibility or other issues.",
                                 "The error message from `lmerTest` is as follows:",
                                 e$message, sep = "\n"))
        stop(message, call. = FALSE)
    })

    # Estimate sample-specific biases
    if (is.null(theta)) {
        # Initial values
        theta = rep(0, n_samp)

        # REML fits
        fits = lapply(seq_len(n_tax), function(i) {
            df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
            fit = tryCatch(
                {
                    suppressWarnings(suppressMessages(
                        lmerTest::lmer(tformula, data = df, control = lme_control)
                    ))
                },
                error = function(e) {
                    NA
                }
            )
            return(fit)
        })

        # Degree of freedom
        dof = NULL

        # Coefficients
        empty_coef = rep(NA, n_fix_eff)
        names(empty_coef) = fix_eff
        beta = lapply(fits, function(i) {
          beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
          if (inherits(i, "lmerModLmerTest")) {
            summ_i = summary(i)
            coef_i = summ_i$coefficients[, "Estimate"]
          } else {
            coef_i = empty_coef
          }
          beta_i[match(names(coef_i), fix_eff)] = coef_i
          return(beta_i)
        })
        beta = do.call("rbind", beta)

        # Iterative REML
        iterNum = 0
        epsilon = 100
        empty_fitted = rep(NA, n_samp)
        names(empty_fitted) = samp_id
        while (epsilon > tol & iterNum < max_iter) {
            # Updating beta
            fits = lapply(seq_len(n_tax), function(i) {
              df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
              fit = tryCatch(
                  {
                      suppressWarnings(suppressMessages(
                          lmerTest::lmer(tformula,
                                         data = df,
                                         control = lme_control)
                      ))
                  },
                  error = function(e) {
                      NA
                  }
              )
              return(fit)
            })

            beta_new = lapply(fits, function(i) {
                beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
                if (inherits(i, "lmerModLmerTest")) {
                  summ_i = summary(i)
                  coef_i = summ_i$coefficients[, "Estimate"]
                } else {
                  coef_i = empty_coef
                }
                beta_i[match(names(coef_i), fix_eff)] = coef_i
                return(beta_i)
            })
            beta_new = do.call("rbind", beta_new)

            # Updating theta
            y_crt_hat = lapply(fits, function(i) {
              y_crt_hat_i = rep(0, n_samp)
              fitted_i = if (inherits(i, "lmerModLmerTest")) {
                stats::fitted(i)
              } else {
                empty_fitted
              }
              y_crt_hat_i[match(names(fitted_i), samp_id)] = fitted_i
              return(y_crt_hat_i)
            })
            y_crt_hat = do.call("rbind", y_crt_hat)
            theta_new = colMeans(y - y_crt_hat, na.rm = TRUE)

            # Iteration
            epsilon = sqrt(sum((beta_new - beta)^2, na.rm = TRUE) +
                               sum((theta_new - theta)^2, na.rm = TRUE))
            iterNum = iterNum + 1
            beta = beta_new
            theta = theta_new

            if (verbose) {
                txt = sprintf(paste0("REML iteration = ", iterNum,
                                     ", epsilon = ", signif(epsilon, 2)))
                message(txt)
            }
        }

        # Residuals
        empty_resid = rep(NA, n_samp)
        names(empty_resid) = samp_id
        eps = lapply(fits, function(i) {
            eps_i = rep(0, n_samp)
            if (inherits(i, "lmerModLmerTest")) {
              summ_i = summary(i)
              resid_i = summ_i$residuals
            } else {
              resid_i = empty_resid
            }
            eps_i[match(names(resid_i), samp_id)] = resid_i
            return(eps_i)
        })
        eps = do.call("rbind", eps)

        # Variance-covariance matrices
        empty_vcov = matrix(NA, nrow = n_fix_eff, ncol = n_fix_eff)
        colnames(empty_vcov) = fix_eff
        rownames(empty_vcov) = fix_eff
        vcov_hat = lapply(fits, function(i) {
            Sigma_hat_i = diag(0.1, nrow = n_fix_eff)
            colnames(Sigma_hat_i) = fix_eff
            rownames(Sigma_hat_i) = fix_eff
            if (inherits(i, "lmerModLmerTest")) {
              summ_i = summary(i)
              vcov_hat_i = as.matrix(summ_i$vcov)
            } else {
              vcov_hat_i = empty_vcov
            }
            Sigma_hat_i[match(rownames(vcov_hat_i), fix_eff),
                        match(colnames(vcov_hat_i), fix_eff)] = vcov_hat_i
            return(Sigma_hat_i)
        })
        var_hat = lapply(vcov_hat, function(i) {
            return(diag(i))
        })
        var_hat = do.call("rbind", var_hat)
    } else {
        # REML fits
        fits = lapply(seq_len(n_tax), function(i) {
          df = data.frame(y_crt = unlist(y[i, ]) - theta, meta_data)
          fit = tryCatch(
              {
                  suppressWarnings(suppressMessages(
                      lmerTest::lmer(tformula, data = df, control = lme_control)
                  ))
              },
              error = function(e) {
                  NA
              }
          )
          return(fit)
        })

        # Degree of freedom
        dof = lapply(fits, function(i) {
          if (inherits(i, "lmerModLmerTest")) {
            summary(i)$coefficients[, "df"]
          } else {
            rep(999, n_fix_eff)
          }
        })
        dof = do.call("rbind", dof)

        # Coefficients
        empty_coef = rep(NA, n_fix_eff)
        names(empty_coef) = fix_eff
        beta = lapply(fits, function(i) {
          beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
          if (inherits(i, "lmerModLmerTest")) {
            summ_i = summary(i)
            coef_i = summ_i$coefficients[, "Estimate"]
          } else {
            coef_i = empty_coef
          }
          beta_i[match(names(coef_i), fix_eff)] = coef_i
          return(beta_i)
        })
        beta = do.call("rbind", beta)

        # Residuals
        empty_resid = rep(NA, n_samp)
        names(empty_resid) = samp_id
        eps = lapply(fits, function(i) {
          eps_i = rep(0, n_samp)
          if (inherits(i, "lmerModLmerTest")) {
            summ_i = summary(i)
            resid_i = summ_i$residuals
          } else {
            resid_i = empty_resid
          }
          eps_i[match(names(resid_i), samp_id)] = resid_i
          return(eps_i)
        })
        eps = do.call("rbind", eps)

        # Variance-covariance matrices
        empty_vcov = matrix(NA, nrow = n_fix_eff, ncol = n_fix_eff)
        colnames(empty_vcov) = fix_eff
        rownames(empty_vcov) = fix_eff
        vcov_hat = lapply(fits, function(i) {
          Sigma_hat_i = diag(0.1, nrow = n_fix_eff)
          colnames(Sigma_hat_i) = fix_eff
          rownames(Sigma_hat_i) = fix_eff
          if (inherits(i, "lmerModLmerTest")) {
            summ_i = summary(i)
            vcov_hat_i = as.matrix(summ_i$vcov)
          } else {
            vcov_hat_i = empty_vcov
          }
          Sigma_hat_i[match(rownames(vcov_hat_i), fix_eff),
                      match(colnames(vcov_hat_i), fix_eff)] = vcov_hat_i
          return(Sigma_hat_i)
        })
        var_hat = lapply(vcov_hat, function(i) {
          return(diag(i))
        })
        var_hat = do.call("rbind", var_hat)
    }

    if (!is.null(dof)) {
      colnames(dof) = fix_eff
      rownames(dof) = tax_id
    }
    colnames(beta) = fix_eff
    rownames(beta) = tax_id
    names(theta) = samp_id
    names(vcov_hat) = tax_id
    colnames(var_hat) = fix_eff
    rownames(var_hat) = tax_id
    rownames(eps) = tax_id

    output = list(fits = fits, beta = beta, theta = theta, eps = eps,
                  dof = dof, vcov_hat = vcov_hat, var_hat = var_hat)
    return(output)
}

# E-M algorithm
.bias_em = function(beta, var_hat, tol, max_iter) {
    beta = beta[!is.na(beta)]
    nu0 = var_hat
    nu0 = nu0[!is.na(nu0)]

    if (any(nu0 == 0)) {
        stop_txt = sprintf(paste("Zero variances have been detected for the following taxa:",
                                 paste(names(which(nu0 == 0)), collapse = ", "),
                                 "Please remove these taxa or select a more parsimonious model",
                                 sep = "\n"))
        stop(stop_txt, call. = FALSE)
    }

    # Initials
    pi0_0 = 0.75
    pi1_0 = 0.125
    pi2_0 = 0.125
    delta_0 = mean(beta[beta >= quantile(beta, 0.25, na.rm = TRUE)&
                            beta <= quantile(beta, 0.75, na.rm = TRUE)],
                   na.rm = TRUE)
    if(is.na(delta_0)) delta_0 = mean(beta, na.rm = TRUE)
    l1_0 = mean(beta[beta < quantile(beta, 0.125, na.rm = TRUE)],
                na.rm = TRUE)
    if(is.na(l1_0)) l1_0 = min(beta, na.rm = TRUE)
    l2_0 = mean(beta[beta > quantile(beta, 0.875, na.rm = TRUE)],
                na.rm = TRUE)
    if(is.na(l2_0)) l2_0 = max(beta, na.rm = TRUE)
    kappa1_0 = var(beta[beta < quantile(beta, 0.125, na.rm = TRUE)],
                   na.rm = TRUE)
    if(is.na(kappa1_0)|kappa1_0 == 0) kappa1_0 = 1
    kappa2_0 = var(beta[beta > quantile(beta, 0.875, na.rm = TRUE)],
                   na.rm = TRUE)
    if(is.na(kappa2_0)|kappa2_0 == 0) kappa2_0 = 1

    # Apply E-M algorithm
    # Store all paras in vectors/matrices
    pi0_vec = pi0_0
    pi1_vec = pi1_0
    pi2_vec = pi2_0
    delta_vec = delta_0
    l1_vec = l1_0
    l2_vec = l2_0
    kappa1_vec = kappa1_0
    kappa2_vec = kappa2_0
    n_tax = length(beta)

    # E-M iteration
    iterNum = 0
    epsilon = 100
    while (epsilon > tol & iterNum < max_iter) {
        # Current value of paras
        pi0 = pi0_vec[length(pi0_vec)]
        pi1 = pi1_vec[length(pi1_vec)]
        pi2 = pi2_vec[length(pi2_vec)]
        delta = delta_vec[length(delta_vec)]
        l1 = l1_vec[length(l1_vec)]
        l2 = l2_vec[length(l2_vec)]
        kappa1 = kappa1_vec[length(kappa1_vec)]
        kappa2 = kappa2_vec[length(kappa2_vec)]

        # E-step
        pdf0 = vapply(seq(n_tax), function(i)
            dnorm(beta[i], delta, sqrt(nu0[i])),
            FUN.VALUE = double(1))
        pdf1 = vapply(seq(n_tax), function(i)
            dnorm(beta[i], delta + l1, sqrt(nu0[i] + kappa1)),
            FUN.VALUE = double(1))
        pdf2 = vapply(seq(n_tax), function(i)
            dnorm(beta[i], delta + l2, sqrt(nu0[i] + kappa2)),
            FUN.VALUE = double(1))
        r0i = pi0*pdf0/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
        r0i[is.na(r0i)] = 0
        r1i = pi1*pdf1/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
        r1i[is.na(r1i)] = 0
        r2i = pi2*pdf2/(pi0*pdf0 + pi1*pdf1 + pi2*pdf2)
        r2i[is.na(r2i)] = 0

        # M-step
        pi0_new = mean(r0i, na.rm = TRUE)
        pi1_new = mean(r1i, na.rm = TRUE)
        pi2_new = mean(r2i, na.rm = TRUE)
        delta_new = sum(r0i*beta/nu0 + r1i*(beta-l1)/(nu0+kappa1) +
                            r2i*(beta-l2)/(nu0+kappa2), na.rm = TRUE)/
            sum(r0i/nu0 + r1i/(nu0+kappa1) + r2i/(nu0+kappa2), na.rm = TRUE)
        l1_new = min(sum(r1i*(beta-delta)/(nu0+kappa1), na.rm = TRUE)/
                         sum(r1i/(nu0+kappa1), na.rm = TRUE), 0)
        if (is.na(l1_new)) l1_new = 0
        l2_new = max(sum(r2i*(beta-delta)/(nu0+kappa2), na.rm = TRUE)/
                         sum(r2i/(nu0+kappa2), na.rm = TRUE), 0)
        if (is.na(l2_new)) l2_new = 0

        # Nelder-Mead simplex algorithm for kappa1 and kappa2
        obj_kappa1 = function(x){
            log_pdf = log(vapply(seq(n_tax), function(i)
                dnorm(beta[i], delta+l1, sqrt(nu0[i]+x)),
                FUN.VALUE = double(1)))
            log_pdf[is.infinite(log_pdf)] = 0
            -sum(r1i*log_pdf, na.rm = TRUE)
        }
        kappa1_new = nloptr::neldermead(x0 = kappa1,
                                        fn = obj_kappa1, lower = 0)$par

        obj_kappa2 = function(x){
            log_pdf = log(vapply(seq(n_tax), function(i)
                dnorm(beta[i], delta+l2, sqrt(nu0[i]+x)),
                FUN.VALUE = double(1)))
            log_pdf[is.infinite(log_pdf)] = 0
            -sum(r2i*log_pdf, na.rm = TRUE)
        }
        kappa2_new = nloptr::neldermead(x0 = kappa2,
                                        fn = obj_kappa2, lower = 0)$par

        # Merge to the paras vectors/matrices
        pi0_vec = c(pi0_vec, pi0_new)
        pi1_vec = c(pi1_vec, pi1_new)
        pi2_vec = c(pi2_vec, pi2_new)
        delta_vec = c(delta_vec, delta_new)
        l1_vec = c(l1_vec, l1_new)
        l2_vec = c(l2_vec, l2_new)
        kappa1_vec = c(kappa1_vec, kappa1_new)
        kappa2_vec = c(kappa2_vec, kappa2_new)

        # Calculate the new epsilon
        epsilon = sqrt((pi0_new-pi0)^2 + (pi1_new-pi1)^2 + (pi2_new-pi2)^2 +
                           (delta_new-delta)^2 + (l1_new-l1)^2 + (l2_new-l2)^2 +
                           (kappa1_new-kappa1)^2 + (kappa2_new-kappa2)^2)
        iterNum = iterNum + 1
    }

    # The EM estimator of bias
    delta_em = delta_new

    # The WLS estimator of bias
    pi1 = pi1_new
    pi2 = pi2_new
    l1 = l1_new
    l2 = l2_new
    kappa1 = kappa1_new
    kappa2 = kappa2_new
    # Cluster 0
    C0 = which(beta >= quantile(beta, pi1, na.rm = TRUE) &
                   beta < quantile(beta, 1 - pi2, na.rm = TRUE))
    # Cluster 1
    C1 = which(beta < quantile(beta, pi1, na.rm = TRUE))
    # Cluster 2
    C2 = which(beta >= quantile(beta, 1 - pi2, na.rm = TRUE))
    # Numerator of the WLS estimator
    nu = nu0
    nu[C1] = nu[C1] + kappa1
    nu[C2] = nu[C2] + kappa2
    wls_deno = sum(1 / nu)
    # Denominator of the WLS estimator
    wls_nume = 1 / nu
    wls_nume[C0] = (wls_nume * beta)[C0]
    wls_nume[C1] = (wls_nume * (beta - l1))[C1]
    wls_nume[C2] = (wls_nume * (beta - l2))[C2]
    wls_nume = sum(wls_nume)

    delta_wls = wls_nume / wls_deno

    # Estimate the variance of bias
    var_delta = 1 / wls_deno
    if (is.na(var_delta)) var_delta = 0

    output = c(delta_em = delta_em,
               delta_wls = delta_wls,
               var_delta = var_delta)
}




