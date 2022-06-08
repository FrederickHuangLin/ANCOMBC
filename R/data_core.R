# Data pre-processing
data_core = function(pseq, prv_cut, lib_cut,
                     tax_keep = NULL, samp_keep = NULL) {
    feature_table = abundances(pseq)
    meta_data = meta(pseq)

    # Discard taxa with prevalences < prv_cut
    if (is.null(tax_keep)) {
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

    fiuo_core = list(feature_table = feature_table,
                     meta_data = meta_data)
    return(fiuo_core)
}

