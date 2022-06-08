s_diff_est = function(feature_table) {
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
