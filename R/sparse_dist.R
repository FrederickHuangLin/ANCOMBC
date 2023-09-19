.sparse_dist = function(mat, wins_quant, R, thresh_hard, max_p) {
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

    # Calculation
    d = dim(mat)[2]
    taxanames = colnames(mat)
    comb = function(...) {
        mapply('rbind', ..., SIMPLIFY = FALSE)
    }

    idx = NULL
    dcorr_list = foreach(idx = seq_len(d - 1), .combine = 'comb', .multicombine = TRUE, .packages = "energy") %dorng% {
        dcorr_idx = rep(NA, d)
        p_val_idx = rep(NA, d)

        mat_x = mat[!is.na(mat[, idx]), , drop = FALSE]
        x = mat_x[, idx]

        if (length(x) == 1) {
          dcorr_idx = rep(0, d)
          p_val_idx = rep(1, d)
        } else {
          # Distance correlation
          dcorr_idx[(idx + 1):d] = apply(mat_x[, (idx + 1):d, drop = FALSE], 2,
                                         function(y) {
                                           z = x[!is.na(y)]
                                           y = y[!is.na(y)]
                                           dcor(z, y, index = 1.0)
                                         })

          # P-values
          p_val_idx[(idx + 1):d] = apply(mat_x[, (idx + 1):d, drop = FALSE], 2,
                                         function(y) {
                                           z = x[!is.na(y)]
                                           y = y[!is.na(y)]
                                           dcor.test(z, y, index = 1.0, R = R)$p.value
                                         })
        }

        list(dcorr_idx, p_val_idx)
        }

    dcorr = rbind(dcorr_list[[1]], rep(NA, d))
    dcorr_p = rbind(dcorr_list[[2]], rep(NA, d))
    # Symmetrize the matrix
    dcorr[lower.tri(dcorr)] = t(dcorr)[lower.tri(dcorr)]
    diag(dcorr) = 1
    dcorr[mat_cooccur < 2] = 0
    dcorr[is.infinite(dcorr)] = 0
    dcorr_p[lower.tri(dcorr_p)] = t(dcorr_p)[lower.tri(dcorr_p)]
    diag(dcorr_p) = 0
    dcorr_p[mat_cooccur < 2] = 1
    dcorr_p[is.na(dcorr_p)] = 1
    dcorr_p[is.infinite(dcorr_p)] = 1
    dimnames(dcorr) = list(taxanames, taxanames)
    dimnames(dcorr_p) = list(taxanames, taxanames)

    dcorr_fl = p_filter(mat = dcorr, mat_p = dcorr_p, max_p = max_p)
    dcorr_fl = mat_thresh(mat = dcorr_fl, th = thresh_hard, soft = FALSE)

    # Output
    result = list(mat_cooccur = mat_cooccur,
                  dcorr = dcorr,
                  dcorr_p = dcorr_p,
                  dcorr_fl = dcorr_fl)
    return(result)
}
