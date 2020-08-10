# Identify structural zeros
get_struc_zero = function(feature_table, meta_data, group, neg_lb) {
    group_data = factor(meta_data[, group])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1

    p_hat = t(apply(present_table, 1, function(x)
        unlist(tapply(x, group_data, function(y) mean(y, na.rm = TRUE)))))
    samp_size = t(apply(feature_table, 1, function(x)
        unlist(tapply(x, group_data, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    zero_ind = (p_hat == 0)
    # Do we classify a taxon as a structural zero by its negative lower bound?
    if (neg_lb) zero_ind[p_hat_lo <= 0] = TRUE

    colnames(zero_ind) = paste0("structural_zero (", group,
                                " = ", colnames(zero_ind), ")")
    return(zero_ind)
}
