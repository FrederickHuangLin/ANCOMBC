# Identify structural zeros
get_struc_zero = function(feature_table, meta_data, group, neg_lb) {
    group_data = factor(meta_data[, group])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1
    n_taxa = nrow(feature_table)
    n_group = nlevels(group_data)

    p_hat = matrix(NA, nrow = n_taxa, ncol = n_group)
    rownames(p_hat) = rownames(feature_table)
    colnames(p_hat) = levels(group_data)
    for (i in seq_len(n_taxa)) {
        p_hat[i, ] = tapply(present_table[i, ], group_data,
                            function(x) mean(x, na.rm = TRUE))
    }

    samp_size = matrix(NA, nrow = n_taxa, ncol = n_group)
    rownames(samp_size) = rownames(feature_table)
    colnames(samp_size) = levels(group_data)
    for (i in seq_len(n_taxa)) {
        samp_size[i, ] = tapply(as.matrix(feature_table)[i, ], group_data,
                                function(x) length(x[!is.na(x)]))
    }

    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    zero_ind = (p_hat == 0)
    # Do we classify a taxon as a structural zero by its negative lower bound?
    if (neg_lb) zero_ind[p_hat_lo <= 0] = TRUE

    colnames(zero_ind) = paste0("structural_zero (", group,
                                " = ", colnames(zero_ind), ")")
    return(zero_ind)
}
