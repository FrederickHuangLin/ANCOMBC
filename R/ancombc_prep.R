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
