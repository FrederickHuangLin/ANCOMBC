# Fit the fixed-effects design x to every taxon's theta-adjusted response. Taxa
# are grouped by missing-value pattern and each pattern is solved with one
# multi-response QR. A rank-deficient group is refitted taxon by taxon with
# stats::lm(), which drops an absent factor level and fails on a single-level
# factor; tformula and meta_data define that per-taxon model.
#
# Returns, aligned to the taxa of Ymat:
#   beta    taxa x covariate coefficients, NA for an unfitted taxon
#   fitted  taxa x sample fitted values, 0 for a sample dropped from the fit
#   dof     residual degrees of freedom per taxon, 999 when not estimable
.lm_fit_all = function(x, Ymat, meta_data, tformula) {
    n_tax = nrow(Ymat)
    n_samp = ncol(Ymat)
    p = ncol(x)
    fix_eff = colnames(x)
    tax_id = rownames(Ymat)
    samp_id = colnames(Ymat)

    x_ok = stats::complete.cases(x)
    beta = matrix(NA_real_, nrow = n_tax, ncol = p,
                  dimnames = list(tax_id, fix_eff))
    fitted = matrix(NA_real_, nrow = n_tax, ncol = n_samp,
                    dimnames = list(tax_id, samp_id))
    dof = rep(999L, n_tax)
    names(dof) = tax_id

    # Single-taxon fit. Ymat[i, ] is theta-adjusted, so it is the response.
    fit_one = function(i) {
        df = data.frame(y_crt = Ymat[i, ], meta_data)
        fit = suppressWarnings(try(stats::lm(tformula, data = df), silent = TRUE))
        if (inherits(fit, "lm")) {
            bi = rep(0, p)
            ci = stats::coef(fit)
            bi[match(names(ci), fix_eff)] = ci
            beta[i, ] <<- bi
            fi = rep(0, n_samp)
            fv = stats::fitted(fit)
            fi[match(names(fv), samp_id)] = fv
            fitted[i, ] <<- fi
            dof[i] <<- fit$df.residual
        }
        # An unfitted taxon keeps beta and fitted at NA and dof at 999
    }

    # use[i, j] indicates that sample j is usable for taxon i: the response is
    # finite and the design row is complete. Taxa sharing a usable-sample
    # pattern are solved together.
    if (all(x_ok) && all(is.finite(Ymat))) {
        use = matrix(TRUE, nrow = n_tax, ncol = n_samp)
        groups = list(seq_len(n_tax))
    } else {
        use = is.finite(Ymat) & matrix(x_ok, nrow = n_tax, ncol = n_samp,
                                       byrow = TRUE)
        keys = do.call(paste0, asplit(use * 1L, 2L))
        groups = split(seq_len(n_tax), factor(keys, levels = unique(keys)))
    }

    for (idx in groups) {
        rows = use[idx[1L], ]
        if (!any(rows)) {
            # No usable samples
            for (i in idx) fit_one(i)
            next
        }
        xr = x[rows, , drop = FALSE]
        Yr = t(Ymat[idx, rows, drop = FALSE])  # n_used x length(idx)
        fit = stats::lm.fit(xr, Yr)

        if (fit$rank < ncol(xr)) {
            # Rank-deficient design
            for (i in idx) fit_one(i)
            next
        }

        co = fit$coefficients
        if (is.null(dim(co))) co = matrix(co, ncol = 1L)
        beta[idx, ] = t(co)

        fit_vals = fit$fitted.values
        if (is.null(dim(fit_vals))) fit_vals = matrix(fit_vals, ncol = 1L)
        fitted[idx, rows] = t(fit_vals)
        fitted[idx, !rows] = 0
        dof[idx] = fit$df.residual
    }

    list(beta = beta, fitted = fitted, dof = dof)
}

# Sandwich (Huber-White) variance estimator. For taxon i,
# V_i = (X'X)^- ( sum_j eps_ij^2 x_j x_j' ) (X'X)^-. Any entry of a term
# involving a missing value is set to 0.1. The per-sample outer products
# x_j x_j' do not depend on the taxon, so they are computed once, stored as the
# rows of an n_samp x p^2 matrix, and accumulated over samples for all taxa
# together.
.sandwich_vcov = function(x, eps, fix_eff) {
    n_tax = nrow(eps)
    n_samp = ncol(eps)
    p = ncol(x)

    x_cc = x[stats::complete.cases(x), ]
    XTX_inv = MASS::ginv(t(x_cc) %*% x_cc)

    # Row j of XX is the outer product x_j x_j' in column-major order.
    XX = matrix(NA_real_, nrow = n_samp, ncol = p * p)
    for (j in seq_len(n_samp)) XX[j, ] = as.vector(x[j, ] %*% t(x[j, ]))

    eps2 = eps^2
    vcov_hat = vector(mode = "list", length = n_tax)
    var_hat = matrix(NA, nrow = n_tax, ncol = p)
    dn = list(fix_eff, fix_eff)

    # Accumulate taxa in blocks to bound the size of the temporaries
    block = max(1L, as.integer(2^17 / (p * p)))
    for (start in seq.int(1L, n_tax, by = block)) {
        idx = seq.int(start, min(start + block - 1L, n_tax))
        sigma2_xxT = matrix(0, nrow = length(idx), ncol = p * p)
        for (j in seq_len(n_samp)) {
            term_j = outer(eps2[idx, j], XX[j, ])
            term_j[is.na(term_j)] = 0.1
            sigma2_xxT = sigma2_xxT + term_j
        }
        for (k in seq_along(idx)) {
            v_i = XTX_inv %*% matrix(sigma2_xxT[k, ], nrow = p, ncol = p) %*%
                XTX_inv
            dimnames(v_i) = dn
            vcov_hat[[idx[k]]] = v_i
            var_hat[idx[k], ] = diag(v_i)
        }
    }
    list(vcov_hat = vcov_hat, var_hat = var_hat)
}

# Iterative MLE
.iter_mle = function(x, y, meta_data, formula, theta = NULL,
                     tol, max_iter, verbose = FALSE) {
    y = as.matrix(y)
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
        beta = .lm_fit_all(x, sweep(y, 2, theta, "-"), meta_data, tformula)$beta

        # Degree of freedom
        dof = NULL

        # Iterative least squares
        iterNum = 0
        epsilon = 100
        y_crt_hat = NULL
        while (epsilon > tol & iterNum < max_iter) {
            # Updating beta
            fit_res = .lm_fit_all(x, sweep(y, 2, theta, "-"), meta_data, tformula)
            beta_new = fit_res$beta
            y_crt_hat = fit_res$fitted

            # Updating theta
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

        # Residuals
        if (is.null(y_crt_hat)) {
            y_crt_hat = .lm_fit_all(x, sweep(y, 2, theta, "-"), meta_data, tformula)$fitted
        }
        eps = t(t(y - y_crt_hat) - theta)

        # Variance-covariance matrices
        sw = .sandwich_vcov(x, eps, fix_eff)
        vcov_hat = sw$vcov_hat
        var_hat = sw$var_hat
    } else {
        # ML fits
        fit_res = .lm_fit_all(x, sweep(y, 2, theta, "-"), meta_data, tformula)
        beta = fit_res$beta
        y_crt_hat = fit_res$fitted

        # Degree of freedom
        dof = fit_res$dof
        dof = matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)

        # Residuals
        eps = t(t(y - y_crt_hat) - theta)

        # Variance-covariance matrices
        sw = .sandwich_vcov(x, eps, fix_eff)
        vcov_hat = sw$vcov_hat
        var_hat = sw$var_hat
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

# Build the lme4 model structure for one pattern of usable samples. The
# structure does not depend on the response values, so it is shared by every
# taxon and every iteration with that pattern. theta0 and lambdat0 hold the
# starting values of the covariance parameters, which lme4 overwrites in place
# while optimizing. keep records the samples retained by na.action.
.lmer_struct = function(tformula, data, lme_control) {
    lmod = lme4::lFormula(formula = tformula, data = data,
                          control = lme_control)
    lmod$formula = NULL
    na_act = attr(lmod$fr, "na.action")
    keep = seq_len(nrow(data))
    if (!is.null(na_act)) keep = keep[-unclass(na_act)]
    lmod$keep = keep
    lmod$theta0 = lmod$reTrms$theta
    lmod$lambdat0 = lmod$reTrms$Lambdat@x
    return(lmod)
}

# Fit one response with a pre-built model structure. The covariance parameters
# are restored to their starting values in a private copy, so the fit does not
# depend on the fits that precede it.
.lmer_refit = function(lmod, y_crt, lme_control) {
    theta0 = lmod$theta0
    lambdat0 = lmod$lambdat0
    lmod$keep = NULL
    lmod$theta0 = NULL
    lmod$lambdat0 = NULL
    lmod$fr[[1L]] = y_crt
    lambdat = lmod$reTrms$Lambdat
    lambdat@x = lambdat0[seq_along(lambdat0)]
    lmod$reTrms$Lambdat = lambdat
    lmod$reTrms$theta = theta0[seq_along(theta0)]

    devfun = do.call(lme4::mkLmerDevfun,
                     c(lmod, list(start = NULL, verbose = 0L,
                                  control = lme_control)))
    rho = environment(devfun)
    n_obs = nrow(lmod$fr)
    n_par = length(rho$lower)
    calc_derivs = lme_control$calc.derivs
    if (is.null(calc_derivs)) {
        calc_derivs = n_obs < lme_control$checkConv$check.conv.nobsmax &&
            n_par < lme_control$checkConv$check.conv.nparmax
    }
    opt = lme4::optimizeLmer(devfun, optimizer = lme_control$optimizer,
                             restart_edge = lme_control$restart_edge,
                             boundary.tol = lme_control$boundary.tol,
                             control = lme_control$optCtrl, verbose = 0L,
                             start = NULL, calc.derivs = calc_derivs,
                             force.calc.derivs = isTRUE(lme_control$calc.derivs),
                             use.last.params = lme_control$use.last.params)
    conv = lme4::checkConv(attr(opt, "derivs"), opt$par,
                           ctrl = lme_control$checkConv,
                           lbound = rho$lower, ubound = rho$upper,
                           nobs = n_obs, ndim = n_par)
    fit = lme4::mkMerMod(rho, opt, lmod$reTrms, fr = lmod$fr,
                         mc = quote(lme4::lmer()), lme4conv = conv)
    return(fit)
}

# Group the taxa of y by their pattern of non-missing samples and build one
# model structure per group. A group whose structure cannot be built is recorded
# as NULL and its taxa are reported as unfitted.
.lmer_struct_all = function(y, meta_data, tformula, lme_control) {
    n_tax = nrow(y)
    obs = !is.na(y)
    if (all(obs)) {
        grp = rep(1L, n_tax)
    } else {
        keys = do.call(paste0, asplit(obs * 1L, 2L))
        grp = match(keys, unique(keys))
    }
    n_grp = max(grp)
    structs = vector(mode = "list", length = n_grp)
    for (g in seq_len(n_grp)) {
        df = data.frame(y_crt = y[match(g, grp), ], meta_data)
        structs[[g]] = tryCatch(
            suppressWarnings(suppressMessages(
                .lmer_struct(tformula, df, lme_control)
            )),
            error = function(e) NULL)
    }
    return(list(grp = grp, structs = structs))
}

# Quantities taken from one fitted model. A model whose variance-covariance
# matrix is not positive definite raises an error and the taxon is reported as
# unfitted.
.remle_extract = function(fit) {
    output = list(coef = lme4::fixef(fit),
                  fitted = stats::fitted(fit),
                  eps = stats::residuals(fit, "pearson", scaled = TRUE),
                  vcov = as.matrix(stats::vcov(fit)))
    return(output)
}

# Fit every taxon at the current sampling fractions. A taxon that cannot be
# fitted contributes NA. A fixed effect or sample absent from a taxon's fit
# contributes 0.
.remle_fit_all = function(struct_list, y, theta, tformula, meta_data,
                          fix_eff, samp_id, lme_control) {
    n_tax = nrow(y)
    n_samp = ncol(y)
    n_fix_eff = length(fix_eff)
    grp = struct_list$grp
    structs = struct_list$structs

    beta = matrix(NA_real_, nrow = n_tax, ncol = n_fix_eff)
    fitted = matrix(NA_real_, nrow = n_tax, ncol = n_samp)
    eps = matrix(NA_real_, nrow = n_tax, ncol = n_samp)
    vcov_hat = vector(mode = "list", length = n_tax)
    var_hat = matrix(NA_real_, nrow = n_tax, ncol = n_fix_eff)
    empty_vcov = diag(0.1, nrow = n_fix_eff)
    colnames(empty_vcov) = fix_eff
    rownames(empty_vcov) = fix_eff
    theta_ok = !anyNA(theta)

    for (i in seq_len(n_tax)) {
        lmod = structs[[grp[i]]]
        y_crt = y[i, ] - theta
        fit_i = NULL
        if (theta_ok) {
            if (!is.null(lmod)) {
                fit_i = tryCatch(
                    suppressWarnings(suppressMessages(
                        .remle_extract(.lmer_refit(lmod, y_crt[lmod$keep],
                                                   lme_control))
                    )),
                    error = function(e) NULL)
            }
        } else {
            # A sampling fraction that could not be estimated adds missing
            # responses and changes the model frame, so the taxon is fitted
            # separately
            fit_i = tryCatch(
                suppressWarnings(suppressMessages(
                    .remle_extract(lme4::lmer(tformula,
                                              data = data.frame(y_crt = y_crt,
                                                                meta_data),
                                              control = lme_control))
                )),
                error = function(e) NULL)
        }

        vcov_i = empty_vcov
        if (is.null(fit_i)) {
            vcov_i[] = NA_real_
        } else {
            beta_i = rep(0, n_fix_eff)
            beta_i[match(names(fit_i$coef), fix_eff)] = fit_i$coef
            beta[i, ] = beta_i

            fitted_i = rep(0, n_samp)
            fitted_i[match(names(fit_i$fitted), samp_id)] = fit_i$fitted
            fitted[i, ] = fitted_i

            eps_i = rep(0, n_samp)
            eps_i[match(names(fit_i$eps), samp_id)] = fit_i$eps
            eps[i, ] = eps_i

            vcov_i[match(rownames(fit_i$vcov), fix_eff),
                   match(colnames(fit_i$vcov), fix_eff)] = fit_i$vcov
        }
        vcov_hat[[i]] = vcov_i
        var_hat[i, ] = diag(vcov_i)
    }

    output = list(beta = beta, fitted = fitted, eps = eps,
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

        # Model structures, one per pattern of usable samples. Subtracting
        # theta does not change which responses are missing, so the structures
        # are reused by every iteration
        struct_list = .lmer_struct_all(y = y, meta_data = meta_data,
                                       tformula = tformula,
                                       lme_control = lme_control)

        # REML fits
        para = .remle_fit_all(struct_list = struct_list, y = y, theta = theta,
                              tformula = tformula, meta_data = meta_data,
                              fix_eff = fix_eff, samp_id = samp_id,
                              lme_control = lme_control)

        # Degree of freedom
        dof = NULL

        # Coefficients
        beta = para$beta

        # Iterative REML
        iterNum = 0
        epsilon = 100
        while (epsilon > tol & iterNum < max_iter) {
            # Updating beta. The first iteration uses the fits at theta = 0
            if (iterNum > 0) {
                para = .remle_fit_all(struct_list = struct_list, y = y,
                                      theta = theta, tformula = tformula,
                                      meta_data = meta_data, fix_eff = fix_eff,
                                      samp_id = samp_id,
                                      lme_control = lme_control)
            }
            beta_new = para$beta

            # Updating theta
            y_crt_hat = para$fitted
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

        # Residuals and variance-covariance matrices
        fits = NULL
        eps = para$eps
        vcov_hat = para$vcov_hat
        var_hat = para$var_hat
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

        # Model summaries
        summ_list = lapply(fits, function(i) {
          if (inherits(i, "lmerModLmerTest")) summary(i) else NULL
        })

        # Degree of freedom
        dof = lapply(summ_list, function(i) {
          if (is.null(i)) rep(999, n_fix_eff) else i$coefficients[, "df"]
        })
        dof = do.call("rbind", dof)

        # Coefficients
        empty_coef = rep(NA, n_fix_eff)
        names(empty_coef) = fix_eff
        beta = lapply(summ_list, function(i) {
          beta_i = rep(0, length(fix_eff)) # prevent errors of missing values
          coef_i = if (is.null(i)) empty_coef else i$coefficients[, "Estimate"]
          beta_i[match(names(coef_i), fix_eff)] = coef_i
          return(beta_i)
        })
        beta = do.call("rbind", beta)

        # Residuals
        empty_resid = rep(NA, n_samp)
        names(empty_resid) = samp_id
        eps = lapply(summ_list, function(i) {
          eps_i = rep(0, n_samp)
          resid_i = if (is.null(i)) empty_resid else i$residuals
          eps_i[match(names(resid_i), samp_id)] = resid_i
          return(eps_i)
        })
        eps = do.call("rbind", eps)

        # Variance-covariance matrices
        empty_vcov = matrix(NA, nrow = n_fix_eff, ncol = n_fix_eff)
        colnames(empty_vcov) = fix_eff
        rownames(empty_vcov) = fix_eff
        vcov_hat = lapply(summ_list, function(i) {
          Sigma_hat_i = diag(0.1, nrow = n_fix_eff)
          colnames(Sigma_hat_i) = fix_eff
          rownames(Sigma_hat_i) = fix_eff
          vcov_hat_i = if (is.null(i)) empty_vcov else as.matrix(i$vcov)
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
    nu0 = var_hat
    neither_na = !(is.na(beta) | is.na(nu0))
    beta = beta[neither_na]
    nu0 = nu0[neither_na]

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
    # Nelder-Mead options
    nm_opts = nloptr::nl.opts(list())
    nm_opts["algorithm"] = "NLOPT_LN_NELDERMEAD"

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
        pdf0 = dnorm(beta, delta, sqrt(nu0))
        pdf1 = dnorm(beta, delta + l1, sqrt(nu0 + kappa1))
        pdf2 = dnorm(beta, delta + l2, sqrt(nu0 + kappa2))
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
            log_pdf = log(dnorm(beta, delta + l1, sqrt(nu0 + x)))
            log_pdf[is.infinite(log_pdf)] = 0
            -sum(r1i*log_pdf, na.rm = TRUE)
        }
        kappa1_new = nloptr::nloptr(x0 = kappa1, eval_f = obj_kappa1,
                                    lb = 0, ub = NULL,
                                    opts = nm_opts)$solution

        obj_kappa2 = function(x){
            log_pdf = log(dnorm(beta, delta + l2, sqrt(nu0 + x)))
            log_pdf[is.infinite(log_pdf)] = 0
            -sum(r2i*log_pdf, na.rm = TRUE)
        }
        kappa2_new = nloptr::nloptr(x0 = kappa2, eval_f = obj_kappa2,
                                    lb = 0, ub = NULL,
                                    opts = nm_opts)$solution

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




