res_combine_zero = function(x, group, struc_zero, zero_ind, alpha,
                            global, res, res_global) {
    covariates = setdiff(colnames(x), "(Intercept)")

    # Set p/q-values of structural zeros to be 0s.
    if (struc_zero) {
        group_ind = grepl(group, covariates)
        zero_mask = 1 - apply(zero_ind, 1, function(x)
            sum(x) > 0 & sum(x) < ncol(zero_ind))
        res$p_val[, group_ind] = res$p_val[, group_ind] * zero_mask
        res$q_val[, group_ind] = res$q_val[, group_ind] * zero_mask
        res$diff_abn = data.frame(res$q_val < alpha & !is.na(res$q_val),
                                  check.names = FALSE)

        # Global test
        if (global) {
            zero_mask = 1 - apply(zero_ind, 1, function(x)
                sum(x) > 0 & sum(x) < ncol(zero_ind))
            res_global[, "p_val"] = res_global[, "p_val"] * zero_mask
            res_global[, "q_val"] = res_global[, "q_val"] * zero_mask
            res_global[, "diff_abn"] = res_global[, "q_val"] < alpha &
                !is.na(res_global[, "q_val"])
        }
    }
    fiuo_out = list(res = res, res_global = res_global)
    return(fiuo_out)
}
