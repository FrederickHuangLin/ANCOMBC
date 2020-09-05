# Data pre-processing
data_prep = function(phyloseq, group, zero_cut, lib_cut, global = global) {
    feature_table = as(otu_table(phyloseq), "matrix")
    feature_table = data.frame(feature_table, check.names = FALSE)
    meta_data = as(sample_data(phyloseq), "data.frame")
    # Drop unused levels
    meta_data[] = lapply(meta_data, function(x)
        if(is.factor(x)) factor(x) else x)
    # Check the group variable
    if (is.null(group)) {
        if (global) {
            stop("Please specify the group variable for the global test.")
        }
    } else {
        n_level = length(unique(meta_data[, group]))
        if (n_level < 2) {
            stop("The group variable should have >= 2 categories.")
        } else if (n_level < 3) {
            global = FALSE
            warning("The multi-group comparison will be deactivated as the group variable has < 3 categories.")
        }
    }

    # Discard taxa with zeros >= zero_cut
    zero_prop = apply(feature_table, 1, function(x)
        sum(x == 0, na.rm = TRUE)/length(x[!is.na(x)]))
    tax_del = which(zero_prop >= zero_cut)
    if (length(tax_del) > 0) {
        feature_table = feature_table[- tax_del, ]
    }

    # Discard samples with library size < lib_cut
    lib_size = colSums(feature_table, na.rm = TRUE)
    if(any(lib_size < lib_cut)){
        subj_del = which(lib_size < lib_cut)
        feature_table = feature_table[, - subj_del]
        meta_data = meta_data[- subj_del, ]
    }
    fiuo_prep = list(feature_table = feature_table,
                     meta_data = meta_data,
                     global = global)
    return(fiuo_prep)
}

