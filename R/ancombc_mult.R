# ANCOM-BC global test
.ancombc_global_F = function(x, group, beta_hat, vcov_hat,
                             dof = NULL, p_adj_method, alpha){
    tax_id = rownames(beta_hat)
    n_tax = nrow(beta_hat)
    covariates = colnames(x)
    output = data.frame(matrix(NA, nrow = n_tax, ncol = 5))
    colnames(output) = c("taxon", "W", "p_val", "q_val", "diff_abn")
    output$taxon = tax_id

    # Loop over the parameters of interest
    group_ind = grepl(group, covariates)
    beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
    vcov_hat_sub = lapply(vcov_hat, function(x) {
        x = x[group_ind, group_ind, drop = FALSE]
    })

    if (is.null(dof)) {
      for (i in seq_len(n_tax)) {
        # Loop over taxa
        beta_hat_sub_i = beta_hat_sub[i, ]
        vcov_hat_sub_i = vcov_hat_sub[[i]]
        A = diag(x = 1, nrow = length(beta_hat_sub_i))

        suppressWarnings(W_global <- try(t(A %*% beta_hat_sub_i) %*%
                                           MASS::ginv(A %*% vcov_hat_sub_i %*% t(A)) %*%
                                           (A %*% beta_hat_sub_i),
                                         silent = TRUE))

        if (inherits(W_global, "try-error")) {
          output[i, "W"] = NA
          output[i, "p_val"] = 1
        } else {
          p_global = 2 * min(pchisq(W_global, df = length(beta_hat_sub_i),
                                    lower.tail = TRUE),
                             pchisq(W_global, df = length(beta_hat_sub_i),
                                    lower.tail = FALSE))
          output[i, "W"] = W_global
          output[i, "p_val"] = p_global
        }
      }
    } else {
      for (i in seq_len(n_tax)) {
        # Loop over taxa
        beta_hat_sub_i = beta_hat_sub[i, ]
        vcov_hat_sub_i = vcov_hat_sub[[i]]
        dof_i = unique(dof[i, ])
        A = diag(x = 1, nrow = length(beta_hat_sub_i))

        suppressWarnings(W_global <- try(t(A %*% beta_hat_sub_i) %*%
                                           MASS::ginv(A %*% vcov_hat_sub_i %*% t(A)) %*%
                                           (A %*% beta_hat_sub_i),
                                         silent = TRUE))

        if (inherits(W_global, "try-error")) {
          output[i, "W"] = NA
          output[i, "p_val"] = 1
        } else {
          p_global = 2 * min(pf(W_global,
                                df1 = length(beta_hat_sub_i),
                                df2 = dof_i,
                                lower.tail = TRUE),
                             pf(W_global,
                                df1 = length(beta_hat_sub_i),
                                df2 = dof_i,
                                lower.tail = FALSE))
          output[i, "W"] = W_global
          output[i, "p_val"] = p_global
        }
      }
    }

    # Model summary
    q_global = p.adjust(output[, "p_val"], method = p_adj_method)
    q_global[is.na(q_global)] = 1
    diff_global = q_global <= alpha & !is.na(q_global)

    output$q_val = q_global
    output$diff_abn = diff_global
    return(output)
}

.ancombc_global_LRT = function(full_model, fix_formula, rand_formula,
                               control, x, group,
                               y, meta_data, p_adj_method, alpha){
  tax_id = rownames(y)
  n_tax = nrow(y)
  covariates = colnames(x)
  output = data.frame(matrix(NA, nrow = n_tax, ncol = 5))
  colnames(output) = c("taxon", "W", "p_val", "q_val", "diff_abn")
  output$taxon = tax_id

  # Perform LRT
  reduce_fix_formula = gsub(pattern = paste0(" \\+ ", group),
                            replacement = "", x = fix_formula)
  reduce_formula = formula(paste0("y ~ ",
                                  reduce_fix_formula,
                                  "+ ", rand_formula))

  reduced_model = lapply(seq_len(n_tax), function(i) {
    df = data.frame(y = unlist(y[i, ]), meta_data)
    fit = try(suppressMessages(lmerTest::lmer(reduce_formula,
                                              data = df,
                                              control = control)),
              silent = TRUE)
    if (inherits(fit, "try-error")) {fit = NA}
    return(fit)
  })

  W_p_global = lapply(seq_len(n_tax), function(i) {
    model_comparison = try(suppressMessages(anova(full_model[[i]], reduced_model[[i]])),
              silent = TRUE)
    if (inherits(model_comparison, "try-error")) {
      output = c(W = NA, p = 1)
    } else {
      output = c(W = model_comparison$Chisq[2],
                 p = model_comparison$`Pr(>Chisq)`[2])
    }
    return(output)
  })
  W_p_global = do.call("rbind", W_p_global)
  W_global = W_p_global[, "W"]
  p_global = W_p_global[, "p"]

  # Model summary
  q_global = p.adjust(p_global, method = p_adj_method)
  q_global[is.na(q_global)] = 1
  diff_global = q_global <= alpha & !is.na(q_global)

  output$W = W_global
  output$p_val = p_global
  output$q_val = q_global
  output$diff_abn = diff_global
  return(output)
}

# ANCOM-BC multiple pairwise comparisons
.ancombc_pair = function(x, group, beta_hat, var_hat, vcov_hat, dof,
                         fwer_ctrl_method, alpha, full_model,
                         fix_formula, rand_formula, control, y, meta_data) {
    covariates = colnames(x)

    # Subset the parameters of interest
    group_ind = grepl(group, covariates)
    beta_hat_sub = beta_hat[, group_ind, drop = FALSE]
    vcov_hat_sub = lapply(vcov_hat, function(x) {
        x[group_ind, group_ind, drop = FALSE]
    })
    dof_group = dof[, group_ind]

    # Run the combination function to obtain pairwise comparison results
    beta_hat_pair = t(apply(beta_hat_sub, 1, function(x)
        .combn_fun(x, fun = base::diff, sep = "_")))
    var_hat_pair = t(vapply(vcov_hat_sub, function(x)
        .combn_fun2(x, fun = .var_diff, sep = "_"),
        FUN.VALUE = double(ncol(beta_hat_pair))))
    rownames(var_hat_pair) = rownames(beta_hat_pair)
    se_hat_pair = sqrt(var_hat_pair)
    W_pair = beta_hat_pair/se_hat_pair

    # Obtain p-values and mdFDR adjusted p-values
    p_q_pair = .mdfdr(global_test = "pairwise",
                      W = W_pair,
                      dof = dof_group,
                      fwer_ctrl_method = fwer_ctrl_method,
                      x = x, group = group,
                      beta_hat = beta_hat,
                      vcov_hat = vcov_hat,
                      alpha = alpha,
                      full_model = full_model,
                      fix_formula = fix_formula,
                      rand_formula = rand_formula,
                      control = control,
                      y = y,
                      meta_data = meta_data)
    p_hat_pair = p_q_pair$p_val
    q_hat_pair = p_q_pair$q_val
    diff_pair = ifelse(q_hat_pair <= alpha, TRUE, FALSE)

    output = list(beta = beta_hat_pair, se = se_hat_pair,
                  W = W_pair, p_val = p_hat_pair,
                  q_val = q_hat_pair, diff_abn = diff_pair)
    return(output)
}

# ANCOM-BC Dunnet's type of test
.dunn_global = function(x, group, W, B, dof, p_adj_method, alpha) {
    covariates = colnames(x)
    group_ind = grepl(group, covariates)
    n_group = sum(group_ind)
    n_tax = nrow(W)
    tax_id = rownames(W)
    output = data.frame(matrix(NA, nrow = n_tax, ncol = 5))
    colnames(output) = c("taxon", "W", "p_val", "q_val", "diff_abn")
    output$taxon = tax_id

    suppressWarnings(W_global <- apply(W, 1, function(x) max(abs(x), na.rm = TRUE)))

    W_global_null = matrix(NA, nrow = n_tax, ncol = B)
    for (b in seq_len(B)) {
        W_null_b = matrix(unlist(apply(dof, seq_len(2), function(df) rt(1, df = df))),
                          nrow = nrow(dof), ncol = ncol(dof))
        W_global_null_b = apply(W_null_b, 1, function(x)
            max(abs(x), na.rm = TRUE))
        W_global_null[, b] = W_global_null_b
    }
    p_global = 1/B * apply(W_global_null > W_global, 1, function(x)
        sum(x, na.rm = TRUE))

    q_global = p.adjust(p_global, method = p_adj_method)
    q_global[is.na(q_global)] = 1
    diff_global = q_global <= alpha & !is.na(q_global)

    output$W = W_global
    output$p_val = p_global
    output$q_val = q_global
    output$diff_abn = diff_global
    return(output)
}

.ancombc_dunn = function(x, group, beta_hat, var_hat, dof,
                         B, fwer_ctrl_method, alpha) {
    covariates = colnames(x)

    # Subset the parameters of interest
    group_ind = grepl(group, covariates)
    beta_hat_dunn = beta_hat[, group_ind]
    var_hat_dunn = var_hat[, group_ind]
    se_hat_dunn = sqrt(var_hat_dunn)
    W_dunn = beta_hat_dunn/se_hat_dunn
    dof_dunn = dof[, group_ind]

    # Obtain p-values and mdFDR adjusted p-values
    p_q_dunn = .mdfdr(global_test = "dunnet", W = W_dunn, dof = dof_dunn,
                      fwer_ctrl_method = fwer_ctrl_method,
                      x = x, group = group, B = B, alpha = alpha)
    p_hat_dunn = p_q_dunn$p_val
    q_hat_dunn = p_q_dunn$q_val
    diff_dunn = ifelse(q_hat_dunn <= alpha, TRUE, FALSE)

    output = list(beta = beta_hat_dunn, se = se_hat_dunn,
                  W = W_dunn, p_val = p_hat_dunn,
                  q_val = q_hat_dunn, diff_abn = diff_dunn)
    return(output)
}

# ANCOM-BC pattern analysis
.ancombc_trend = function(x, group, beta_hat, var_hat, vcov_hat,
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

    fun_list = list(.constrain_est, .l_infty)

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

    ident_mat = diag(1, nrow = n_group)
    var_hat_sub_dup = var_hat_sub[, rep(seq_len(ncol(var_hat_sub)), n_trend)]
    var_hat_sub_dup[is.na(var_hat_sub_dup)] = 1
    b = NULL
    W_trend_null = foreach(b = seq_len(B), .combine = 'cbind') %dorng%
      {
        set.seed(b)
        beta_null = matrix(rnorm(n_group * n_tax), nrow = n_tax)
        beta_null_opt = t(apply(beta_null, 1, function(x) {
          beta_null_opt_x = unlist(lapply(X = contrast,
                                          FUN = fun_list[[1]],
                                          beta_hat = x,
                                          vcov_hat = ident_mat,
                                          solver = solver))
          return(beta_null_opt_x)
        }))
        beta_null_opt = beta_null_opt * sqrt(var_hat_sub_dup)

        beta_null_opt_list = split(beta_null_opt, row(beta_null_opt))
        beta_null_opt_names = colnames(beta_null_opt)
        l_null = lapply(beta_null_opt_list, function(x) {
          l_null_x = unlist(lapply(X = seq_len(n_trend),
                                   FUN = function(j) {
                                     l = fun_list[[2]](beta_opt = x[grepl(trend_name[j], beta_null_opt_names)],
                                                       node = node[[j]])
                                     return(l)
                                   }))
          return(l_null_x)
        })
        l_null = do.call(rbind, l_null)
        W_null = apply(l_null, 1, function(x) max(x, na.rm = TRUE))
      }
    W_trend_null = as.matrix(W_trend_null)

    p_trend = 1/B * apply(W_trend_null > W_trend, 1, function(x)
        sum(x, na.rm = TRUE))

    q_trend = p.adjust(p_trend, method = p_adj_method)
    q_trend[is.na(q_trend)] = 1
    diff_trend = q_trend <= alpha & !is.na(q_trend)

    output = list(beta = beta_hat_trend,
                  se = sqrt(var_hat_sub),
                  W = W_trend, p_val = p_trend,
                  q_val = q_trend, diff_abn = diff_trend)
    return(output)
}
