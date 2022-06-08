# Data pre-processing
data_qc = function(meta_data, group, global) {
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
            warn_txt = sprintf(paste("The group variable has < 3 categories",
                                     "The multi-group comparison will be deactivated (global = FALSE)",
                                     sep = "\n"))
            warning(warn_txt, call. = FALSE)
        }

        # Check the sample size per group
        size_per_group = tapply(meta_data[, group], meta_data[, group], length)
        if (any(size_per_group < 2)) {
            stop_txt = sprintf(paste("Sample size per group should be >= 2",
                                     "Small sample size detected for the following group(s): ",
                                     paste(names(size_per_group)[which(size_per_group < 2)], collapse = ", "),
                                     sep = "\n"))
            stop(stop_txt, call. = FALSE)
        } else if (any(size_per_group < 5)) {
            warn_txt = sprintf(paste("Small sample size detected for the following group(s): ",
                                     paste(names(size_per_group)[which(size_per_group < 5)], collapse = ", "),
                                     "ANCOM-BC results would be unstable when the sample size is < 5 per group",
                                     sep = "\n"))
            warning(warn_txt, call. = FALSE)
        }
    }

    fiuo_qc = list(meta_data = meta_data,
                   global = global)
    return(fiuo_qc)
}

