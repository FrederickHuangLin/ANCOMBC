res_combine_zero = function(x, group, struc_zero, zero_ind, alpha,
                            global, res, res_global) {
    covariates = colnames(x); n_covariates = length(covariates)

    # Set p/q-values of structural zeros to be 0s.
    if (struc_zero) {
        group_ind = grepl(group, covariates)
        zero_mask = 1 - apply(zero_ind, 1, function(x) any(x == 1))
        res$p_val[, group_ind] = res$p_val[, group_ind] * zero_mask
        res$q_val[, group_ind] = res$q_val[, group_ind] * zero_mask
        res$diff_abn = ifelse(res$q_val < alpha, TRUE, FALSE)

        # Global test
        if (global) {
            res_global[, "p_val"] = res_global[, "p_val"] * zero_mask
            res_global[, "q_val"] = res_global[, "q_val"] * zero_mask
            res_global[, "diff_abn"] = ifelse(res_global[, "q_val"] < alpha,
                                              TRUE, FALSE)
        }
    }
    fiuo_out = list(res = res, res_global = res_global)
    return(fiuo_out)
}
