# Sampling fraction difference estimation
.s_diff_est = function(feature_table) {
    if (nrow(feature_table) < 50) {
        warn_txt = sprintf(paste("The number of taxa used for estimating sample-specific biases is: ",
                                 nrow(feature_table),
                                 "A large number of taxa (> 50) is required for the statistical consistency",
                                 sep = "\n"))
        warning(warn_txt, call. = FALSE)
    }
    o = log(feature_table)

    o[is.infinite(o)] = NA
    o_center = o - rowMeans(o, na.rm = TRUE)

    # Estimate weights
    wt = apply(o_center, 1, function(x) 1/var(x, na.rm = TRUE))
    o_center = o_center[is.finite(wt), ]
    wt = wt[is.finite(wt)]

    # Estimate sampling fraction difference
    s_diff_hat = apply(o_center, 2, function(x) {
        weighted.mean(x, wt, na.rm = TRUE)}
    )

    return(s_diff_hat)
}

# Bias-corrected abundance estimation
.abn_est = function(tse, tax_level, assay_name, pseudo, prv_cut, lib_cut) {
    # Sampling fraction difference estimation
    core1 = .data_core(tse = tse, tax_level = tax_level,
                       assay_name = assay_name, alt = FALSE,
                       prv_cut = prv_cut, lib_cut = lib_cut,
                       tax_keep = NULL, samp_keep = NULL)
    O1 = core1$feature_table
    s_diff_hat = .s_diff_est(O1)

    # Data pre-processing
    samp_keep = names(s_diff_hat)
    core2 = .data_core(tse = tse, tax_level = tax_level,
                       assay_name = assay_name, alt = TRUE,
                       prv_cut = prv_cut, lib_cut = lib_cut,
                       tax_keep = NULL, samp_keep = samp_keep)
    O2 = core2$feature_table
    O2 = O2 + pseudo
    o = log(O2)
    o[is.infinite(o)] = NA
    n = ncol(o)
    d = nrow(o)

    # Bias-corrected abundance estimation
    y_hat = o - rowMeans(o, na.rm = TRUE) - t(replicate(d, s_diff_hat))

    abs_list = list(s_diff_hat = s_diff_hat, y_hat = y_hat)
    return(abs_list)
}

