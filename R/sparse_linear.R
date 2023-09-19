.sparse_linear = function(mat, wins_quant, method, soft, thresh_len,
                          n_cv, thresh_hard, max_p) {
    # Thresholding
    mat_thresh = function(mat, th, soft){
        mat_sign = sign(mat)
        mat_th = mat
        mat_th[abs(mat) <= th] = 0
        if (soft) {
            mat_th[abs(mat) > th] = abs(mat_th[abs(mat) > th]) - th
            mat_th = mat_th * mat_sign
        }
        return(mat_th)
    }

    # Threshold loss function
    thresh_loss = function(mat1, mat2, method, th, soft) {
        corr1 = cor(mat1, method = method, use = "pairwise.complete.obs")
        corr2 = cor(mat2, method = method, use = "pairwise.complete.obs")
        corr_diff = mat_thresh(corr1, th, soft) - corr2
        corr_diff[is.na(corr_diff)] = 0
        loss = norm(corr_diff, type = "F")
        return(loss)
    }

    # Filtering based on p-values
    p_filter = function(mat, mat_p, max_p){
        ind_p = mat_p
        ind_p[mat_p > max_p] = 0
        ind_p[mat_p <= max_p] = 1

        mat_filter = mat * ind_p
        return(mat_filter)
    }

    # Sort taxa
    sort_taxa = sort(colnames(mat))
    mat = mat[, sort_taxa]

    # Winsorization
    mat = apply(mat, 2, function(x)
        DescTools::Winsorize(x, probs = wins_quant, na.rm = TRUE))

    # Co-occurrence matrix
    mat_occur = mat
    mat_occur[mat_occur != 0] = 1
    mat_occur[mat_occur == 0] = 0
    mat_occur[is.na(mat_occur)] = 0

    df_occur = as.data.frame(mat_occur)
    df_occur$sample_id = rownames(df_occur)
    df_occur_long = stats::reshape(df_occur,
                                   direction = "long",
                                   varying = list(colnames(mat_occur)),
                                   v.names = "occur",
                                   idvar = "sample_id",
                                   times = colnames(mat_occur),
                                   new.row.names = seq_len(nrow(df_occur)*ncol(df_occur)))
    names(df_occur_long)[names(df_occur_long) == "time"] = "taxon"
    df_occur_long = df_occur_long[df_occur_long$occur == 1, ]

    mat_cooccur = matrix(0, nrow = ncol(mat_occur), ncol = ncol(mat_occur))
    rownames(mat_cooccur) = colnames(mat_occur)
    colnames(mat_cooccur) = colnames(mat_occur)

    mat_cooccur_comp = crossprod(table(df_occur_long[, seq_len(2)]))
    idx = base::match(colnames(mat_cooccur_comp), colnames(mat_cooccur))
    mat_cooccur[idx, idx] = mat_cooccur_comp
    diag(mat_cooccur) = colSums(mat_occur)

    if (any(mat_cooccur < 10)) {
        warn_txt = sprintf(paste("There are some pairs of taxa that have insufficient (< 10) overlapping samples",
                                 "Proceed with caution since the point estimates for these pairs are unstable",
                                 "For pairs of taxa with no overlapping samples, the point estimates will be replaced with 0s,",
                                 "and the corresponding p-values will be replaced with 1s",
                                 "Please check `mat_cooccur` for details about the co-occurrence pattern",
                                 sep = "\n"))
        warning(warn_txt)
    }

    # Sample size for training and test sets
    n = dim(mat)[1]
    n1 = n - floor(n/log(n))
    n2 = n - n1
    d = dim(mat)[2]

    # Correlation matrix
    corr_list = suppressWarnings(Hmisc::rcorr(x = mat, type = method))
    corr = corr_list$r
    corr[mat_cooccur < 2] = 0
    corr[is.infinite(corr)] = 0

    # Cross-Validation
    max_thresh = max(abs(corr[corr != 1]), na.rm = TRUE)
    thresh_grid = seq(from = 0, to = max_thresh, length.out = thresh_len)

    loss_mat = foreach(i = seq_len(n_cv), .combine = rbind) %dorng% {
        index = sample(seq_len(n), size = n1, replace = FALSE)
        mat1 = mat[index,]
        mat2 = mat[-index,]
        loss = vapply(thresh_grid, FUN = thresh_loss,
                      mat1 = mat1, mat2 = mat2,
                      method = method, soft = soft,
                      FUN.VALUE = double(1))
    }

    # Correlation matrix after thresholding
    loss_vec = colMeans(loss_mat)
    thresh_opt = thresh_grid[which.min(loss_vec)]
    corr_th = mat_thresh(mat = corr, th = thresh_opt, soft = soft)
    corr_th = mat_thresh(mat = corr_th, th = thresh_hard, soft = FALSE)

    # Correlation matrix after filtering
    corr_p = corr_list$P
    diag(corr_p) = 0
    corr_p[mat_cooccur < 2] = 1
    corr_p[is.na(corr_p)] = 1
    corr_p[is.infinite(corr_p)] = 1
    corr_fl = p_filter(mat = corr, mat_p = corr_p, max_p = max_p)
    corr_fl = mat_thresh(mat = corr_fl, th = thresh_hard, soft = FALSE)

    # Output
    result = list(cv_error = loss_vec,
                  thresh_grid = thresh_grid,
                  thresh_opt = thresh_opt,
                  mat_cooccur = mat_cooccur,
                  corr = corr,
                  corr_p = corr_p,
                  corr_th = corr_th,
                  corr_fl = corr_fl)
    return(result)
}
