# Data pre-processing
data_prep = function(phyloseq, group, zero_cut, lib_cut, global = global) {
    feature_table = abundances(phyloseq)
    meta_data = meta(phyloseq)
    # Drop unused levels
    meta_data[] = lapply(meta_data, function(x)
        if(is.factor(x)) factor(x) else x)
    # Check the group variable
    if (is.null(group)) {
        if (global) {
            stop("Please specify the group variable for the global test",
                 call. = FALSE)
        }
    } else {
        # Check the number of groups
        n_level = length(unique(meta_data[, group]))
        if (n_level < 2) {
            stop("The group variable should have >= 2 categories",
                 call. = FALSE)
        } else if (n_level < 3) {
            global = FALSE
            warning("The multi-group comparison will be deactivated as the group variable has < 3 categories",
                    call. = FALSE)
        }

        # Check the sample size per group
        size_per_group = tapply(meta_data[, group], meta_data[, group], length)
        if (any(size_per_group < 2)) {
            stop_txt = sprintf(paste("Sample size per group should be >= 2",
                                     "Small sample size detected for the following group(s): ",
                                     paste(names(size_per_group)[which(size_per_group < 2)], collapse = " "),
                                     sep = "\n"))
            stop(stop_txt, call. = FALSE)
        } else if (any(size_per_group < 5)) {
            warning_txt = sprintf(paste("Small sample size detected for the following group(s): ",
                                        paste(names(size_per_group)[which(size_per_group < 5)], collapse = " "),
                                        "ANCOM-BC results would be unstable when the sample size is < 5 per group",
                                        sep = "\n"))
            warning(warning_txt, call. = FALSE)
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
        feature_table = feature_table[, - subj_del, drop = FALSE]
        meta_data = meta_data[- subj_del, , drop = FALSE]
    }
    fiuo_prep = list(feature_table = feature_table,
                     meta_data = meta_data,
                     global = global)
    return(fiuo_prep)
}

