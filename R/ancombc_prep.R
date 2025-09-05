# Filter data by prevalence and library size
.data_core = function(data, meta_data, prv_cut, lib_cut,
                      tax_keep = NULL, samp_keep = NULL) {
    feature_table = data

    # Discard taxa with prevalences < prv_cut
    if (is.null(tax_keep)) {
        prevalence = apply(feature_table, 1, function(x)
            sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
        tax_keep = which(prevalence >= prv_cut)
    }else if (length(tax_keep) == 0) {
        stop("All taxa contain structural zeros", call. = FALSE)
    } else {
        # Discard taxa with structural zeros
        feature_table = feature_table[tax_keep, , drop = FALSE]
        prevalence = apply(feature_table, 1, function(x)
            sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
        tax_keep = which(prevalence >= prv_cut)
    }

    if (length(tax_keep) > 0) {
        feature_table = feature_table[tax_keep, , drop = FALSE]
    } else {
        stop("No taxa remain under the current cutoff", call. = FALSE)
    }

    # Discard samples with library sizes < lib_cut
    if (is.null(samp_keep)) {
        lib_size = colSums(feature_table, na.rm = TRUE)
        samp_keep = which(lib_size >= lib_cut)
    }
    if (length(samp_keep) > 0){
        feature_table = feature_table[, samp_keep, drop = FALSE]
        meta_data = meta_data[samp_keep, , drop = FALSE]
    } else {
        stop("No samples remain under the current cutoff", call. = FALSE)
    }

    output = list(feature_table = feature_table,
                  meta_data = meta_data,
                  tax_keep = tax_keep,
                  samp_keep = samp_keep)
    return(output)
}

# Identify structural zeros
.get_struc_zero = function(data, meta_data, group, neg_lb) {
    feature_table = data
    tax_name = rownames(data)
    group_data = factor(meta_data[, group])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1
    n_tax = nrow(feature_table)
    n_group = nlevels(group_data)

    p_hat = matrix(NA, nrow = n_tax, ncol = n_group)
    rownames(p_hat) = rownames(feature_table)
    colnames(p_hat) = levels(group_data)
    for (i in seq_len(n_tax)) {
        p_hat[i, ] = tapply(present_table[i, ], group_data,
                            function(x) mean(x, na.rm = TRUE))
    }

    samp_size = matrix(NA, nrow = n_tax, ncol = n_group)
    rownames(samp_size) = rownames(feature_table)
    colnames(samp_size) = levels(group_data)
    for (i in seq_len(n_tax)) {
        samp_size[i, ] = tapply(as.matrix(feature_table)[i, ], group_data,
                                function(x) length(x[!is.na(x)]))
    }

    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    output = (p_hat == 0)
    # Shall we classify a taxon as a structural zero by its negative lower bound?
    if (neg_lb) output[p_hat_lo <= 0] = TRUE

    output = cbind(tax_name, output)
    colnames(output) = c("taxon",
                         paste0("structural_zero (", group,
                                " = ", colnames(output)[-1], ")"))
    output = data.frame(output, check.names = FALSE, row.names = NULL)
    output[, -1] = apply(output[, -1], 2, as.logical)
    return(output)
}

# ANCOM-BC2 core function
.ancombc2_core = function(data, aggregate_data = NULL,
                          meta_data = NULL,
                          fix_formula, rand_formula = NULL,
                          p_adj_method = "holm", pseudo = 0,
                          s0_perc = 0.05, group = NULL,
                          alpha = 0.05, verbose = TRUE,
                          global = FALSE, pairwise = FALSE,
                          dunnet = FALSE, trend = FALSE,
                          iter_control = list(tol = 1e-2,
                                              max_iter = 20,
                                              verbose = FALSE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm",
                                               B = 100),
                          trend_control = list(contrast = NULL,
                                               node = NULL,
                                               solver = "ECOS",
                                               B = 100)){

    # 1. Data pre-processing
    O1 = data + pseudo
    O2 = aggregate_data + pseudo

    n_tax = nrow(O2)
    tax_name = rownames(O2)

    # 2. Estimation of the sample-specific biases
    options(na.action = "na.pass") # Keep NA's in rows of x
    x = stats::model.matrix(formula(paste0("~", fix_formula)), data = meta_data)
    options(na.action = "na.omit") # Switch it back
    fix_eff = colnames(x)
    n_fix_eff = length(fix_eff)

    if (nrow(O1) < 50 & pseudo == 0) {
        warn_txt = sprintf(paste0("The number of taxa used for estimating ",
                                  "sample-specific biases is: ",
                                  nrow(O1),
                                  "\nA large number of taxa (>50) is required ",
                                  "for the consistent estimation of biases"))
        warning(warn_txt, call. = FALSE)
    }
    o1 = log(O1)
    o1[is.infinite(o1)] = NA
    y1 = o1 - rowMeans(o1, na.rm = TRUE)

    # Obtain initial estimates
    if (verbose) {
        message("Obtaining initial estimates ...")
    }

    if (is.null(rand_formula)) {
        para1 = .iter_mle(x = x, y = y1, meta_data = meta_data,
                          formula = fix_formula,
                          theta = NULL, tol = iter_control$tol,
                          max_iter = iter_control$max_iter,
                          verbose = iter_control$verbose)
    } else {
        para1 = .iter_remle(x = x, y = y1, meta_data = meta_data,
                            fix_formula = fix_formula,
                            rand_formula = rand_formula,
                            lme_control = lme_control,
                            theta = NULL, tol = iter_control$tol,
                            max_iter = iter_control$max_iter,
                            verbose = iter_control$verbose)
    }
    beta1 = para1$beta
    var_hat1 = para1$var_hat

    # Apply E-M algorithm
    if (verbose) {
        message("Estimating sample-specific biases ...")
    }
    fun_list = list(.bias_em)

    bias1 = foreach(i = seq_len(ncol(beta1)), .combine = rbind) %dorng% {
        output = fun_list[[1]](beta = beta1[, i],
                               var_hat = var_hat1[, i],
                               tol = em_control$tol,
                               max_iter = em_control$max_iter)
    }
    bias1 = data.frame(bias1, row.names = fix_eff, check.names = FALSE)
    colnames(bias1) = c("delta_em", "delta_wls", "var_delta")

    delta_em = bias1$delta_em
    delta_wls = bias1$delta_wls
    var_delta = bias1$var_delta

    # Obtain the final estimates for sample-specific biases
    beta1 = t(t(beta1) - delta_em)
    theta_hat = matrix(NA, nrow = nrow(y1), ncol = ncol(y1))
    for (i in seq_len(nrow(y1))) {
        theta_hat[i, ] = y1[i, ] - base::rowSums(t(t(x) * beta1[i, ]), na.rm = TRUE)
    }
    theta_hat = colMeans(theta_hat, na.rm = TRUE)
    names(theta_hat) = colnames(y1)

    if (any(is.na(theta_hat))) {
        warn_txt = sprintf(paste("Estimation of sampling fractions failed for the following samples:",
                                 paste(names(which(is.na(theta_hat))), collapse = ", "),
                                 "These samples may have an excessive number of zero values",
                                 sep = "\n"))
        warning(warn_txt, call. = FALSE)
    }

    # 3. Obtain unbiased estimates
    o2 = log(O2)
    o2[is.infinite(o2)] = NA
    y2 = o2 - rowMeans(o2, na.rm = TRUE)
    y_bias_crt = data.frame(t(t(y2) - theta_hat), check.names = FALSE)
    if (is.null(rand_formula)) {
        para2 = .iter_mle(x = x, y = y2, meta_data = meta_data,
                          formula = fix_formula,
                          theta = theta_hat, tol = iter_control$tol,
                          max_iter = iter_control$max_iter,
                          verbose = iter_control$verbose)
    } else {
        para2 = .iter_remle(x = x, y = y2, meta_data = meta_data,
                            fix_formula = fix_formula,
                            rand_formula = rand_formula,
                            lme_control = lme_control,
                            theta = theta_hat, tol = iter_control$tol,
                            max_iter = iter_control$max_iter,
                            verbose = iter_control$verbose)
    }
    beta_hat = para2$beta
    var_hat = para2$var_hat
    dof = para2$dof

    # Account for the variance of delta
    var_hat = sweep(var_hat, 2, var_delta, "+") +
        2 * sqrt(sweep(var_hat, 2, var_delta, "*"))

    # Add a small positive constant to stabilize the variance
    if (is.null(s0_perc)) {
        s02 = 0
    } else {
        s02 = apply(var_hat, 2, function(x)
            stats::quantile(x, s0_perc, na.rm = TRUE))
    }
    var_hat = t(t(var_hat) + s02)
    var_hat[is.na(beta_hat)] = NA
    se_hat = sqrt(var_hat)
    vcov_hat = lapply(seq_len(n_tax), function(i) {
        diag(para2$vcov_hat[[i]]) = var_hat[i, ]
        return(para2$vcov_hat[[i]])
    })

    # 4. Primary results
    if (verbose) {
        message("ANCOM-BC2 primary results ...")
    }
    W = beta_hat/se_hat
    p_hat = 2 * pt(abs(W), df = dof, lower.tail = FALSE)
    p_hat[is.na(p_hat)] = 1
    q_hat = apply(p_hat, 2, function(x) p.adjust(x, method = p_adj_method))
    diff_abn = q_hat <= alpha & !is.na(q_hat)

    beta_prim = data.frame(beta_hat, check.names = FALSE)
    se_prim = data.frame(se_hat, check.names = FALSE)
    W_prim = data.frame(W, check.names = FALSE)
    p_prim = data.frame(p_hat, check.names = FALSE)
    q_prim = data.frame(q_hat, check.names = FALSE)
    diff_prim = data.frame(diff_abn, check.names = FALSE)
    colnames(beta_prim) = paste0("lfc_", colnames(beta_hat))
    colnames(se_prim) = paste0("se_", colnames(se_hat))
    colnames(W_prim) = paste0("W_", colnames(W))
    colnames(p_prim) = paste0("p_", colnames(p_hat))
    colnames(q_prim) = paste0("q_", colnames(q_hat))
    colnames(diff_prim) = paste0("diff_", colnames(diff_abn))
    res = do.call("cbind", list(data.frame(taxon = tax_name),
                                beta_prim, se_prim, W_prim,
                                p_prim, q_prim, diff_prim))
    rownames(res) = NULL

    # 5. Results of global test
    if (global) {
        if (verbose) {
            message("ANCOM-BC2 global test ...")
        }
        if (is.null(rand_formula)) {
            res_global = .ancombc_global_F(x = x, group = group,
                                           beta_hat = beta_hat,
                                           vcov_hat = vcov_hat,
                                           dof = dof,
                                           p_adj_method = p_adj_method,
                                           alpha = alpha)
        } else {
            res_global = .ancombc_global_LRT(full_model = para2$fits,
                                             fix_formula = fix_formula,
                                             rand_formula = rand_formula,
                                             control = lme_control,
                                             x = x, group = group,
                                             y = y_bias_crt,
                                             meta_data = meta_data,
                                             p_adj_method = p_adj_method,
                                             alpha = alpha)
        }
        rownames(res_global) = NULL
    } else { res_global = NULL }

    # 6. Results of multiple pairwise comparisons
    if (pairwise) {
        if (verbose) {
            message("ANCOM-BC2 multiple pairwise comparisons ...")
        }
        res_pair = .ancombc_pair(x = x, group = group,
                                 beta_hat = beta_hat,
                                 var_hat = var_hat,
                                 vcov_hat = vcov_hat,
                                 dof = dof,
                                 fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                                 alpha = alpha,
                                 full_model = para2$fits,
                                 fix_formula = fix_formula,
                                 rand_formula = rand_formula,
                                 control = lme_control,
                                 y = y_bias_crt,
                                 meta_data = meta_data)
        beta_pair = data.frame(res_pair$beta, check.names = FALSE)
        se_pair = data.frame(res_pair$se, check.names = FALSE)
        W_pair = data.frame(res_pair$W, check.names = FALSE)
        p_pair = data.frame(res_pair$p_val, check.names = FALSE)
        q_pair = data.frame(res_pair$q_val, check.names = FALSE)
        diff_pair = data.frame(res_pair$diff_abn, check.names = FALSE)

        # Directional test summary
        colnames(beta_pair) = paste0("lfc_", colnames(beta_pair))
        colnames(se_pair) = paste0("se_", colnames(se_pair))
        colnames(W_pair) = paste0("W_", colnames(W_pair))
        colnames(p_pair) = paste0("p_", colnames(p_pair))
        colnames(q_pair) = paste0("q_", colnames(q_pair))
        colnames(diff_pair) = paste0("diff_", colnames(diff_pair))
        res_pair = do.call("cbind", list(data.frame(taxon = tax_name),
                                         beta_pair, se_pair, W_pair,
                                         p_pair, q_pair, diff_pair))
        pair_col_name = gsub("lfc_", "", colnames(beta_pair))
        rownames(res_pair) = NULL
    } else {
        res_pair = NULL
    }

    # 7. Results of Dunnet's type of test
    if (dunnet) {
        if (verbose) {
            message("ANCOM-BC2 multiple pairwise comparisons against the reference group ...")
        }
        res_dunn = .ancombc_dunn(x = x, group = group, beta_hat = beta_hat,
                                 var_hat = var_hat, dof = dof,
                                 fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                                 B = mdfdr_control$B, alpha = alpha)
        beta_dunn = data.frame(res_dunn$beta, check.names = FALSE)
        se_dunn = data.frame(res_dunn$se, check.names = FALSE)
        W_dunn = data.frame(res_dunn$W, check.names = FALSE)
        p_dunn = data.frame(res_dunn$p_val, check.names = FALSE)
        q_dunn = data.frame(res_dunn$q_val, check.names = FALSE)
        diff_dunn = data.frame(res_dunn$diff_abn, check.names = FALSE)

        # Directional test summary
        colnames(beta_dunn) = paste0("lfc_", colnames(beta_dunn))
        colnames(se_dunn) = paste0("se_", colnames(se_dunn))
        colnames(W_dunn) = paste0("W_", colnames(W_dunn))
        colnames(p_dunn) = paste0("p_", colnames(p_dunn))
        colnames(q_dunn) = paste0("q_", colnames(q_dunn))
        colnames(diff_dunn) = paste0("diff_", colnames(diff_dunn))
        res_dunn = do.call("cbind", list(data.frame(taxon = tax_name),
                                         beta_dunn, se_dunn, W_dunn,
                                         p_dunn, q_dunn, diff_dunn))
        dunn_col_name = gsub("lfc_", "", colnames(beta_dunn))
        rownames(res_dunn) = NULL
    } else {
        res_dunn = NULL
    }

    # 8. Results of pattern analysis
    if (trend) {
        if (verbose) {
            message("ANCOM-BC2 pattern analysis ...")
        }
        res_trend = .ancombc_trend(
            x = x, group = group, beta_hat = beta_hat,
            var_hat = var_hat, vcov_hat = vcov_hat,
            p_adj_method = p_adj_method, alpha = alpha,
            trend_control = trend_control)
        beta_trend = res_trend$beta
        se_trend = res_trend$se
        W_trend = res_trend$W
        p_trend = res_trend$p_val
        q_trend = res_trend$q_val
        diff_trend = res_trend$diff_abn

        # Directional test summary
        colnames(beta_trend) = paste0("lfc_", colnames(beta_trend))
        colnames(se_trend) = paste0("se_", colnames(se_trend))
        res_trend = cbind(data.frame(taxon = tax_name),
                          beta_trend, se_trend,
                          data.frame(W = W_trend, p_val = p_trend,
                                     q_val = q_trend, diff_abn = diff_trend))
        rownames(res_trend) = NULL
    } else {
        res_trend = NULL
    }

    # 9. Outputs
    out = list(feature_table = O2,
               bias_correct_log_table = y_bias_crt,
               samp_frac = theta_hat,
               delta_em = delta_em,
               delta_wls = delta_wls,
               res = res,
               res_global = res_global,
               res_pair = res_pair,
               res_dunn = res_dunn,
               res_trend = res_trend)

    return(out)
}

# ANCOM-BC3 core function
.ancombc3_core = function(data, aggregate_data = NULL,
                          meta_data = NULL,
                          fix_formula, rand_formula = NULL,
                          p_adj_method = "holm", pseudo = 0,
                          s0_perc = 0.05, group = NULL,
                          alpha = 0.05, verbose = TRUE,
                          global = FALSE, pairwise = FALSE,
                          dunnet = FALSE, trend = FALSE,
                          iter_control = list(tol = 1e-2,
                                              max_iter = 20,
                                              verbose = FALSE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm",
                                               B = 100),
                          trend_control = list(contrast = NULL,
                                               node = NULL,
                                               solver = "ECOS",
                                               B = 100)){

    # 1. Data pre-processing
    O1 = data + pseudo
    O2 = aggregate_data + pseudo

    n_tax = nrow(O2)
    tax_name = rownames(O2)

    # 2. Estimation of the sample-specific biases
    options(na.action = "na.pass") # Keep NA's in rows of x
    x = stats::model.matrix(formula(paste0("~", fix_formula)), data = meta_data)
    options(na.action = "na.omit") # Switch it back
    fix_eff = colnames(x)
    n_fix_eff = length(fix_eff)

    if (nrow(O1) < 50 & pseudo == 0) {
        warn_txt = sprintf(paste0("Taxa with prevalence >= 0.50 across samples",
                                  " will be used for estimating sample-specific",
                                  " biases.",
                                  "\nThere are ",
                                  nrow(O1),
                                  " taxa with prevalence >= 0.50",
                                  "\nA large number of taxa (> 50) is required ",
                                  "for the consistent estimation of biases"))
        warning(warn_txt, call. = FALSE)
    }
    o1 = log(O1)
    o1[is.infinite(o1)] = NA
    y1 = o1 - rowMeans(o1, na.rm = TRUE)

    # Obtain initial estimates
    if (verbose) {
        message("Obtaining initial estimates ...")
    }

    if (is.null(rand_formula)) {
        para1 = .iter_mle(x = x, y = y1, meta_data = meta_data,
                          formula = fix_formula,
                          theta = NULL, tol = iter_control$tol,
                          max_iter = iter_control$max_iter,
                          verbose = iter_control$verbose)
    } else {
        para1 = .iter_remle(x = x, y = y1, meta_data = meta_data,
                            fix_formula = fix_formula,
                            rand_formula = rand_formula,
                            lme_control = lme_control,
                            theta = NULL, tol = iter_control$tol,
                            max_iter = iter_control$max_iter,
                            verbose = iter_control$verbose)
    }
    beta1 = para1$beta
    var_hat1 = para1$var_hat

    # Apply E-M algorithm
    if (verbose) {
        message("Estimating sample-specific biases ...")
    }
    fun_list = list(.bias_em)

    bias1 = foreach(i = seq_len(ncol(beta1)), .combine = rbind) %dorng% {
        output = fun_list[[1]](beta = beta1[, i],
                               var_hat = var_hat1[, i],
                               tol = em_control$tol,
                               max_iter = em_control$max_iter)
    }
    bias1 = data.frame(bias1, row.names = fix_eff, check.names = FALSE)
    colnames(bias1) = c("delta_em", "delta_wls", "var_delta")

    delta_em = bias1$delta_em
    delta_wls = bias1$delta_wls
    var_delta = bias1$var_delta

    # Obtain the final estimates for sample-specific biases
    beta1 = t(t(beta1) - delta_em)
    theta_hat = matrix(NA, nrow = nrow(y1), ncol = ncol(y1))
    for (i in seq_len(nrow(y1))) {
        theta_hat[i, ] = y1[i, ] - base::rowSums(t(t(x) * beta1[i, ]), na.rm = TRUE)
    }
    theta_hat = colMeans(theta_hat, na.rm = TRUE)
    names(theta_hat) = colnames(y1)

    if (any(is.na(theta_hat))) {
        warn_txt = sprintf(paste("Estimation of sampling fractions failed for the following samples:",
                                 paste(names(which(is.na(theta_hat))), collapse = ", "),
                                 "These samples may have an excessive number of zero values",
                                 sep = "\n"))
        warning(warn_txt, call. = FALSE)
    }

    # 3. Obtain unbiased estimates
    o2 = log(O2)
    o2[is.infinite(o2)] = NA
    y2 = o2 - rowMeans(o2, na.rm = TRUE)
    y_bias_crt = data.frame(t(t(y2) - theta_hat), check.names = FALSE)
    if (is.null(rand_formula)) {
        para2 = .iter_mle(x = x, y = y2, meta_data = meta_data,
                          formula = fix_formula,
                          theta = theta_hat, tol = iter_control$tol,
                          max_iter = iter_control$max_iter,
                          verbose = iter_control$verbose)
    } else {
        para2 = .iter_remle(x = x, y = y2, meta_data = meta_data,
                            fix_formula = fix_formula,
                            rand_formula = rand_formula,
                            lme_control = lme_control,
                            theta = theta_hat, tol = iter_control$tol,
                            max_iter = iter_control$max_iter,
                            verbose = iter_control$verbose)
    }
    beta_hat = para2$beta
    var_hat = para2$var_hat
    dof = para2$dof

    # Account for the variance of delta
    var_hat = sweep(var_hat, 2, var_delta, "+") +
        2 * sqrt(sweep(var_hat, 2, var_delta, "*"))

    # Add a small positive constant to stabilize the variance
    if (is.null(s0_perc)) {
        s02 = 0
    } else {
        s02 = apply(var_hat, 2, function(x)
            stats::quantile(x, s0_perc, na.rm = TRUE))
    }
    var_hat = t(t(var_hat) + s02)
    var_hat[is.na(beta_hat)] = NA
    se_hat = sqrt(var_hat)
    vcov_hat = lapply(seq_len(n_tax), function(i) {
        diag(para2$vcov_hat[[i]]) = var_hat[i, ]
        return(para2$vcov_hat[[i]])
    })

    # 4. Primary results
    if (verbose) {
        message("ANCOM-BC3 DAA primary results ...")
    }
    W = beta_hat/se_hat
    p_hat = 2 * pt(abs(W), df = dof, lower.tail = FALSE)
    p_hat[is.na(p_hat)] = 1
    q_hat = apply(p_hat, 2, function(x) p.adjust(x, method = p_adj_method))
    diff_abn = q_hat <= alpha & !is.na(q_hat)

    beta_prim = data.frame(beta_hat, check.names = FALSE)
    se_prim = data.frame(se_hat, check.names = FALSE)
    W_prim = data.frame(W, check.names = FALSE)
    p_prim = data.frame(p_hat, check.names = FALSE)
    q_prim = data.frame(q_hat, check.names = FALSE)
    diff_prim = data.frame(diff_abn, check.names = FALSE)
    colnames(beta_prim) = paste0("lfc_", colnames(beta_hat))
    colnames(se_prim) = paste0("se_", colnames(se_hat))
    colnames(W_prim) = paste0("W_", colnames(W))
    colnames(p_prim) = paste0("p_", colnames(p_hat))
    colnames(q_prim) = paste0("q_", colnames(q_hat))
    colnames(diff_prim) = paste0("diff_", colnames(diff_abn))
    res = do.call("cbind", list(data.frame(taxon = tax_name),
                                beta_prim, se_prim, W_prim,
                                p_prim, q_prim, diff_prim))
    rownames(res) = NULL

    # 5. Results of global test
    if (global) {
        if (verbose) {
            message("ANCOM-BC3 DAA global test ...")
        }
        if (is.null(rand_formula)) {
            res_global = .ancombc_global_F(x = x, group = group,
                                           beta_hat = beta_hat,
                                           vcov_hat = vcov_hat,
                                           dof = dof,
                                           p_adj_method = p_adj_method,
                                           alpha = alpha)
        } else {
            res_global = .ancombc_global_LRT(full_model = para2$fits,
                                             fix_formula = fix_formula,
                                             rand_formula = rand_formula,
                                             control = lme_control,
                                             x = x, group = group,
                                             y = y_bias_crt,
                                             meta_data = meta_data,
                                             p_adj_method = p_adj_method,
                                             alpha = alpha)
        }
        rownames(res_global) = NULL
    } else { res_global = NULL }

    # 6. Results of multiple pairwise comparisons
    if (pairwise) {
        if (verbose) {
            message("ANCOM-BC3 DAA multiple pairwise comparisons ...")
        }
        res_pair = .ancombc_pair(x = x, group = group,
                                 beta_hat = beta_hat,
                                 var_hat = var_hat,
                                 vcov_hat = vcov_hat,
                                 dof = dof,
                                 fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                                 alpha = alpha,
                                 full_model = para2$fits,
                                 fix_formula = fix_formula,
                                 rand_formula = rand_formula,
                                 control = lme_control,
                                 y = y_bias_crt,
                                 meta_data = meta_data)
        beta_pair = data.frame(res_pair$beta, check.names = FALSE)
        se_pair = data.frame(res_pair$se, check.names = FALSE)
        W_pair = data.frame(res_pair$W, check.names = FALSE)
        p_pair = data.frame(res_pair$p_val, check.names = FALSE)
        q_pair = data.frame(res_pair$q_val, check.names = FALSE)
        diff_pair = data.frame(res_pair$diff_abn, check.names = FALSE)

        # Directional test summary
        colnames(beta_pair) = paste0("lfc_", colnames(beta_pair))
        colnames(se_pair) = paste0("se_", colnames(se_pair))
        colnames(W_pair) = paste0("W_", colnames(W_pair))
        colnames(p_pair) = paste0("p_", colnames(p_pair))
        colnames(q_pair) = paste0("q_", colnames(q_pair))
        colnames(diff_pair) = paste0("diff_", colnames(diff_pair))
        res_pair = do.call("cbind", list(data.frame(taxon = tax_name),
                                         beta_pair, se_pair, W_pair,
                                         p_pair, q_pair, diff_pair))
        pair_col_name = gsub("lfc_", "", colnames(beta_pair))
        rownames(res_pair) = NULL
    } else {
        res_pair = NULL
    }

    # 7. Results of Dunnet's type of test
    if (dunnet) {
        if (verbose) {
            message("ANCOM-BC3 DAA multiple pairwise comparisons against the reference group ...")
        }
        res_dunn = .ancombc_dunn(x = x, group = group, beta_hat = beta_hat,
                                 var_hat = var_hat, dof = dof,
                                 fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                                 B = mdfdr_control$B, alpha = alpha)
        beta_dunn = data.frame(res_dunn$beta, check.names = FALSE)
        se_dunn = data.frame(res_dunn$se, check.names = FALSE)
        W_dunn = data.frame(res_dunn$W, check.names = FALSE)
        p_dunn = data.frame(res_dunn$p_val, check.names = FALSE)
        q_dunn = data.frame(res_dunn$q_val, check.names = FALSE)
        diff_dunn = data.frame(res_dunn$diff_abn, check.names = FALSE)

        # Directional test summary
        colnames(beta_dunn) = paste0("lfc_", colnames(beta_dunn))
        colnames(se_dunn) = paste0("se_", colnames(se_dunn))
        colnames(W_dunn) = paste0("W_", colnames(W_dunn))
        colnames(p_dunn) = paste0("p_", colnames(p_dunn))
        colnames(q_dunn) = paste0("q_", colnames(q_dunn))
        colnames(diff_dunn) = paste0("diff_", colnames(diff_dunn))
        res_dunn = do.call("cbind", list(data.frame(taxon = tax_name),
                                         beta_dunn, se_dunn, W_dunn,
                                         p_dunn, q_dunn, diff_dunn))
        dunn_col_name = gsub("lfc_", "", colnames(beta_dunn))
        rownames(res_dunn) = NULL
    } else {
        res_dunn = NULL
    }

    # 8. Results of pattern analysis
    if (trend) {
        if (verbose) {
            message("ANCOM-BC3 DAA pattern analysis ...")
        }
        res_trend = .ancombc_trend(
            x = x, group = group, beta_hat = beta_hat,
            var_hat = var_hat, vcov_hat = vcov_hat,
            p_adj_method = p_adj_method, alpha = alpha,
            trend_control = trend_control)
        beta_trend = res_trend$beta
        se_trend = res_trend$se
        W_trend = res_trend$W
        p_trend = res_trend$p_val
        q_trend = res_trend$q_val
        diff_trend = res_trend$diff_abn

        # Directional test summary
        colnames(beta_trend) = paste0("lfc_", colnames(beta_trend))
        colnames(se_trend) = paste0("se_", colnames(se_trend))
        res_trend = cbind(data.frame(taxon = tax_name),
                          beta_trend, se_trend,
                          data.frame(W = W_trend, p_val = p_trend,
                                     q_val = q_trend, diff_abn = diff_trend))
        rownames(res_trend) = NULL
    } else {
        res_trend = NULL
    }

    # 9. Outputs
    out = list(feature_table = O2,
               bias_correct_log_table = y_bias_crt,
               samp_frac = theta_hat,
               delta_em = delta_em,
               delta_wls = delta_wls,
               res = res,
               res_global = res_global,
               res_pair = res_pair,
               res_dunn = res_dunn,
               res_trend = res_trend)

    return(out)
}
