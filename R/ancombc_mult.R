# ANCOM-BC global test
ancombc_global = function(x, group, beta_hat, vcov_hat, p_adj_method, alpha){
    tax_id = rownames(beta_hat)
    n_tax = nrow(beta_hat)
    covariates = colnames(x)
    output = data.frame(matrix(NA, nrow = n_tax, ncol = 4))
    rownames(output) = tax_id
    colnames(output) = c("W", "p_val", "q_val", "diff_abn")

    # Loop over the parameters of interest
    group_ind = grepl(group, covariates)
    beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
    vcov_hat_sub = lapply(vcov_hat, function(x) {
        x = x[group_ind, group_ind, drop = FALSE]
    })

    for (i in seq_len(n_tax)) {
        # Loop over taxa
        beta_hat_sub_i = beta_hat_sub[i, ]
        vcov_hat_sub_i = vcov_hat_sub[[i]]
        A = diag(x = 1, nrow = length(beta_hat_sub_i))
        W_global = t(A %*% beta_hat_sub_i) %*%
            MASS::ginv(A %*% vcov_hat_sub_i %*% t(A)) %*%
            (A %*% beta_hat_sub_i)
        p_global = 2 * min(pchisq(W_global, df = length(beta_hat_sub_i),
                                  lower.tail = TRUE),
                           pchisq(W_global, df = length(beta_hat_sub_i),
                                  lower.tail = FALSE))
        output[i, "W"] = W_global
        output[i, "p_val"] = p_global
    }
    # Model summary
    q_global = p.adjust(output[, "p_val"], method = p_adj_method)
    q_global[is.na(q_global)] = 1
    diff_global = q_global < alpha & !is.na(q_global)

    output$q_val = q_global
    output$diff_abn = diff_global
    return(output)
}

# ANCOM-BC pairwise test
ancombc_pair = function(x, group, beta_hat, var_hat, vcov_hat,
                        fwer_ctrl_method, alpha) {
    covariates = colnames(x)

    # Subset the parameters of interest
    group_ind = grepl(group, covariates)
    beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
    vcov_hat_sub = lapply(vcov_hat, function(x) {
        x[group_ind, group_ind, drop = FALSE]
    })

    # Run the combination function to obtain pairwise comparison results
    beta_hat_pair = t(apply(beta_hat_sub, 1, function(x)
        combn_fun(x, fun = base::diff, sep = "_")))
    var_hat_pair = t(vapply(vcov_hat_sub, function(x)
        combn_fun2(x, fun = var_diff, sep = "_"),
        FUN.VALUE = double(ncol(beta_hat_pair))))
    rownames(var_hat_pair) = rownames(beta_hat_pair)
    se_hat_pair = sqrt(var_hat_pair)
    W_pair = beta_hat_pair/se_hat_pair

    # Obtain p-values and mdFDR adjusted p-values
    p_q_pair = mdfdr(global_test = "pairwise", W = W_pair,
                     fwer_ctrl_method = fwer_ctrl_method,
                     x = x, group = group,
                     beta_hat = beta_hat,
                     vcov_hat = vcov_hat,
                     alpha = alpha)
    p_hat_pair = p_q_pair$p_val
    q_hat_pair = p_q_pair$q_val
    diff_pair = ifelse(q_hat_pair < alpha, TRUE, FALSE)

    output = list(beta = beta_hat_pair, se = se_hat_pair,
                  W = W_pair, p_val = p_hat_pair,
                  q_val = q_hat_pair, diff_abn = diff_pair)
    return(output)
}

# ANCOM-BC Dunnet's type of test
dunn_global = function(x, group, W, B, p_adj_method, alpha) {
    covariates = colnames(x)
    group_ind = grepl(group, covariates)
    n_group = sum(group_ind)
    n_tax = nrow(W)
    tax_id = rownames(W)
    output = data.frame(matrix(NA, nrow = n_tax, ncol = 4))
    rownames(output) = tax_id
    colnames(output) = c("W", "p_val", "q_val", "diff_abn")

    W_global = apply(W, 1, function(x) max(abs(x), na.rm = TRUE))

    W_global_null = matrix(NA, nrow = n_tax, ncol = B)
    for (b in seq_len(B)) {
        W_null_b = matrix(rnorm(n_tax * n_group), nrow = n_tax)
        W_global_null_b = apply(W_null_b, 1, function(x)
            max(abs(x), na.rm = TRUE))
        W_global_null[, b] = W_global_null_b
    }
    p_global = 1/B * apply(W_global_null > W_global, 1, function(x)
        sum(x, na.rm = TRUE))

    q_global = p.adjust(p_global, method = p_adj_method)
    q_global[is.na(q_global)] = 1
    diff_global = q_global < alpha & !is.na(q_global)

    output$W = W_global
    output$p_val = p_global
    output$q_val = q_global
    output$diff_abn = diff_global
    return(output)
}

ancombc_dunn = function(x, group, beta_hat, var_hat,
                        B, fwer_ctrl_method, alpha) {
    covariates = colnames(x)

    # Subset the parameters of interest
    group_ind = grepl(group, covariates)
    beta_hat_dunn = beta_hat[, group_ind]
    var_hat_dunn = var_hat[, group_ind]
    se_hat_dunn = sqrt(var_hat_dunn)
    W_dunn = beta_hat_dunn/se_hat_dunn

    # Obtain p-values and mdFDR adjusted p-values
    p_q_dunn = mdfdr(global_test = "dunnet", W = W_dunn,
                     fwer_ctrl_method = fwer_ctrl_method,
                     x = x, group = group, B = B, alpha = alpha)
    p_hat_dunn = p_q_dunn$p_val
    q_hat_dunn = p_q_dunn$q_val
    diff_dunn = ifelse(q_hat_dunn < alpha, TRUE, FALSE)

    output = list(beta = beta_hat_dunn, se = se_hat_dunn,
                  W = W_dunn, p_val = p_hat_dunn,
                  q_val = q_hat_dunn, diff_abn = diff_dunn)
    return(output)
}

# ANCOM-BC trend test
ancombc_trend = function(x, group, beta_hat, var_hat, vcov_hat,
                         p_adj_method, alpha,
                         trend_control = list(contrast = NULL,
                                              node = NULL,
                                              solver = "ECOS",
                                              B = 100)){
    tax_id = rownames(beta_hat)
    n_tax = nrow(beta_hat)
    covariates = colnames(x)

    group_ind = grepl(group, covariates)
    n_group = sum(group_ind)
    beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
    var_hat_sub = var_hat[, group_ind, drop = FALSE]
    vcov_hat_sub = lapply(vcov_hat, function(x) {
        x = x[group_ind, group_ind, drop = FALSE]
    })

    contrast = trend_control$contrast
    node = trend_control$node
    solver = trend_control$solver
    B = trend_control$B

    n_trend = length(contrast)
    trend_name = names(contrast)

    fun_list = list(constrain_est, l_infty)

    beta_hat_opt_all = foreach(i = seq_len(n_tax), .combine = rbind) %dorng% {
        beta_hat_opt = unlist(lapply(X = contrast,
                                     FUN = fun_list[[1]],
                                     beta_hat = beta_hat_sub[i, ],
                                     vcov_hat = vcov_hat_sub[[i]],
                                     solver = solver))
    }

    l = matrix(NA, nrow = n_tax, ncol = n_trend)
    for (i in seq_len(n_trend)) {
        beta_hat_opt_i = beta_hat_opt_all[, grepl(trend_name[i], colnames(beta_hat_opt_all))]
        node_i = node[[i]]
        l[, i] = apply(beta_hat_opt_i, 1, function(x) fun_list[[2]](x, node_i))
    }
    W_trend = apply(l, 1, function(x) max(x, na.rm = TRUE))
    names(W_trend) = tax_id
    opt_trend = apply(l, 1, function(x) trend_name[which.max(x)])
    beta_hat_trend = matrix(NA, nrow = n_tax, ncol = n_group)
    for (i in seq_len(n_tax)) {
        beta_hat_trend[i, ] = beta_hat_opt_all[i, grepl(opt_trend[i], colnames(beta_hat_opt_all))]
    }
    rownames(beta_hat_trend) = tax_id
    colnames(beta_hat_trend) = colnames(beta_hat_sub)

    # Generate the null distribution of W
    rng = rngtools::RNGseq(B * n_tax, 1234)
    b = r = NULL
    W_trend_null = foreach(b = seq_len(B), .combine = 'cbind') %:%
        foreach(i = seq_len(n_tax), r = rng[(b - 1) * n_tax + seq_len(n_tax)], .combine = 'c') %dopar% {
            rngtools::setRNG(r)
            beta_null = rnorm(n_group) * sqrt(var_hat_sub[i, ])
            beta_null_opt = unlist(lapply(X = contrast,
                                          FUN = fun_list[[1]],
                                          beta_hat = beta_null,
                                          vcov_hat = vcov_hat_sub[[i]],
                                          solver = solver))
            l_null = unlist(lapply(X = seq_len(n_trend),
                                   FUN = function(j) {
                                       l = fun_list[[2]](beta_opt = beta_null_opt[grepl(trend_name[j], names(beta_null_opt))],
                                                         node = node[[j]])
                                       return(l)
                                   }))
            W_null = max(l_null, na.rm = TRUE)
        }
    W_trend_null = as.matrix(W_trend_null)

    p_trend = 1/B * apply(W_trend_null > W_trend, 1, function(x)
        sum(x, na.rm = TRUE))

    q_trend = p.adjust(p_trend, method = p_adj_method)
    q_trend[is.na(q_trend)] = 1
    diff_trend = q_trend < alpha & !is.na(q_trend)

    output = list(beta = beta_hat_trend,
                  se = sqrt(var_hat_sub),
                  W = W_trend, p_val = p_trend,
                  q_val = q_trend, diff_abn = diff_trend)
    return(output)
}
