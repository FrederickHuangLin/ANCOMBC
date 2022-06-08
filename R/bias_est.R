# E-M algorithm for estimating the bias term
bias_est = function(beta, var_hat, tol, max_iter) {
    delta_em = rep(NA, ncol(beta) - 1)
    delta_wls = rep(NA, ncol(beta) - 1)
    var_delta = rep(NA, ncol(beta) - 1)
    for (i in seq_along(delta_em)) {
        # Ignore the intercept
        Delta = beta[, i + 1]
        Delta = Delta[!is.na(Delta)]
        nu0 = var_hat[, i + 1]
        nu0 = nu0[!is.na(nu0)]

        if (any(nu0 == 0)) {
            stop_txt = sprintf(paste("Zero variances detected for the following taxa:",
                                     paste(names(which(nu0 == 0)), collapse = ", "),
                                     "Consider removing these taxa or choosing a more parsimonious model",
                                     sep = "\n"))
            stop(stop_txt, call. = FALSE)
        }

        # Initials
        pi0_0 = 0.75
        pi1_0 = 0.125
        pi2_0 = 0.125
        delta_0 = mean(Delta[Delta >= quantile(Delta, 0.25, na.rm = TRUE)&
                                 Delta <= quantile(Delta, 0.75, na.rm = TRUE)],
                       na.rm = TRUE)
        if(is.na(delta_0)) delta_0 = mean(Delta, na.rm = TRUE)
        l1_0 = mean(Delta[Delta < quantile(Delta, 0.125, na.rm = TRUE)],
                    na.rm = TRUE)
        if(is.na(l1_0)) l1_0 = min(Delta, na.rm = TRUE)
        l2_0 = mean(Delta[Delta > quantile(Delta, 0.875, na.rm = TRUE)],
                    na.rm = TRUE)
        if(is.na(l2_0)) l2_0 = max(Delta, na.rm = TRUE)
        kappa1_0 = var(Delta[Delta < quantile(Delta, 0.125, na.rm = TRUE)],
                       na.rm = TRUE)
        if(is.na(kappa1_0)|kappa1_0 == 0) kappa1_0 = 1
        kappa2_0 = var(Delta[Delta > quantile(Delta, 0.875, na.rm = TRUE)],
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
        n_taxa = length(Delta)

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
            pdf0 = vapply(seq(n_taxa), function(i)
                dnorm(Delta[i], delta, sqrt(nu0[i])), FUN.VALUE = double(1))
            pdf1 = vapply(seq(n_taxa), function(i)
                dnorm(Delta[i], delta + l1, sqrt(nu0[i] + kappa1)),
                FUN.VALUE = double(1))
            pdf2 = vapply(seq(n_taxa), function(i)
                dnorm(Delta[i], delta + l2, sqrt(nu0[i] + kappa2)),
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
            delta_new = sum(r0i*Delta/nu0 + r1i*(Delta-l1)/(nu0+kappa1) +
                                r2i*(Delta-l2)/(nu0+kappa2), na.rm = TRUE)/
                sum(r0i/nu0 + r1i/(nu0+kappa1) + r2i/(nu0+kappa2), na.rm = TRUE)
            l1_new = min(sum(r1i*(Delta-delta)/(nu0+kappa1), na.rm = TRUE)/
                             sum(r1i/(nu0+kappa1), na.rm = TRUE), 0)
            l2_new = max(sum(r2i*(Delta-delta)/(nu0+kappa2), na.rm = TRUE)/
                             sum(r2i/(nu0+kappa2), na.rm = TRUE), 0)

            # Nelder-Mead simplex algorithm for kappa1 and kappa2
            obj_kappa1 = function(x){
                log_pdf = log(vapply(seq(n_taxa), function(i)
                    dnorm(Delta[i], delta+l1, sqrt(nu0[i]+x)),
                    FUN.VALUE = double(1)))
                log_pdf[is.infinite(log_pdf)] = 0
                -sum(r1i*log_pdf, na.rm = TRUE)
            }
            kappa1_new = nloptr::neldermead(x0 = kappa1,
                                            fn = obj_kappa1, lower = 0)$par

            obj_kappa2 = function(x){
                log_pdf = log(vapply(seq(n_taxa), function(i)
                    dnorm(Delta[i], delta+l2, sqrt(nu0[i]+x)),
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
        delta_em[i] = delta_new

        # The WLS estimator of bias
        pi1 = pi1_new
        pi2 = pi2_new
        l1 = l1_new
        l2 = l2_new
        kappa1 = kappa1_new
        kappa2 = kappa2_new
        # Cluster 0
        C0 = which(Delta >= quantile(Delta, pi1, na.rm = TRUE) &
                       Delta < quantile(Delta, 1 - pi2, na.rm = TRUE))
        # Cluster 1
        C1 = which(Delta < quantile(Delta, pi1, na.rm = TRUE))
        # Cluster 2
        C2 = which(Delta >= quantile(Delta, 1 - pi2, na.rm = TRUE))
        # Numerator of the WLS estimator
        nu = nu0
        nu[C1] = nu[C1] + kappa1
        nu[C2] = nu[C2] + kappa2
        wls_deno = sum(1 / nu)
        # Denominator of the WLS estimator
        wls_nume = 1 / nu
        wls_nume[C0] = (wls_nume * Delta)[C0]
        wls_nume[C1] = (wls_nume * (Delta - l1))[C1]
        wls_nume[C2] = (wls_nume * (Delta - l2))[C2]
        wls_nume = sum(wls_nume)

        delta_wls[i] = wls_nume / wls_deno

        # Estimate the variance of bias
        var_delta[i] = 1 / wls_deno
        if (is.na(var_delta[i])) var_delta[i] = 0
    }

    fiuo_bias = list(delta_em = delta_em, delta_wls = delta_wls,
                     var_delta = var_delta)
}

