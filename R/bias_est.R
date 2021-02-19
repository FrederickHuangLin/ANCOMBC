# E-M algorithm for estimating the bias term
bias_est = function(beta, var_hat, tol, max_iter, n_taxa) {
    delta_em = rep(NA, ncol(beta) - 1)
    delta_wls = rep(NA, ncol(beta) - 1)
    var_delta = rep(NA, ncol(beta) - 1)
    for (i in seq_along(delta_em)) {
        # Ignore the intercept
        Delta = beta[, i + 1]
        Delta = Delta[!is.na(Delta)]
        nu0 = var_hat[, i + 1]
        nu0 = nu0[!is.na(nu0)]

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
        fiuo_em = em_iter(Delta, nu0, pi0_0, pi1_0, pi2_0, delta_0,
                          l1_0, l2_0, kappa1_0, kappa2_0, tol, max_iter)

        # The EM estimator of bias
        delta_em[i] = fiuo_em$delta

        # The WLS estimator of bias
        pi1 = fiuo_em$pi1
        pi2 = fiuo_em$pi2
        l1 = fiuo_em$l1
        l2 = fiuo_em$l2
        kappa1 = fiuo_em$kappa1
        kappa2 = fiuo_em$kappa2
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

