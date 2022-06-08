abn_est = function(pseqs, pseudo, prv_cut, lib_cut) {
    pseq1 = pseqs[[1]]
    pseq2 = pseqs[[2]]

    if (! all(sample_names(pseq1) == sample_names(pseq2))) {
        stop("Sample names of two phyloseq objects does not match.")
    }

    # Sampling fraction difference estimation
    O1 = data_core(pseq1, prv_cut, lib_cut,
                   tax_keep = NULL, samp_keep = NULL)$feature_table
    s_diff_hat = s_diff_est(O1)

    # Data pre-processing
    samp_keep = names(s_diff_hat)
    O2 = data_core(pseq2, prv_cut, lib_cut,
                   tax_keep = NULL, samp_keep = samp_keep)$feature_table
    O2 = O2 + pseudo
    o = log(O2)
    o[is.infinite(o)] = NA
    n = ncol(o)
    d = nrow(o)

    # Absolute abundance estimation
    y_hat = o - rowMeans(o, na.rm = TRUE) - t(replicate(d, s_diff_hat))

    abs_list = list(s_diff_hat = s_diff_hat, y_hat = y_hat)
    return(abs_list)
}
