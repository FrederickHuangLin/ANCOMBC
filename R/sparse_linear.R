# Threshold loss over a grid of thresholding and penalty parameters, for one
# cross-validation split with training correlation matrix corr1 and validation
# correlation matrix corr2. The entries of corr1 are sorted by absolute value,
# so the sums over the entries retained at a given threshold are suffix sums.
.thresh_loss_grid = function(corr1, corr2, soft, thresh, alpha) {
    a = as.vector(corr1)
    b = as.vector(corr2)

    # Entries with a missing value contribute 0 to the squared Frobenius norm
    # at every threshold
    ok = !is.na(a) & !is.na(b)
    a = a[ok]
    b = b[ok]
    o = order(abs(a))
    a = a[o]
    b = b[o]
    a_abs = abs(a)
    a_b = a - b
    m = length(a)

    # Number of entries set to 0 at each threshold
    k = findInterval(thresh, a_abs)

    # Entries set to 0 contribute b^2, entries retained contribute (a - b)^2
    pre_b = c(0, cumsum(b * b))
    suf_ab = c(rev(cumsum(rev(a_b * a_b))), 0)
    fro_sq = pre_b[k + 1] + suf_ab[k + 1]
    if (soft) {
        # Soft thresholding shifts a retained entry by the threshold
        suf_g = c(rev(cumsum(rev(sign(a) * a_b))), 0)
        fro_sq = fro_sq - 2 * thresh * suf_g[k + 1] + thresh^2 * (m - k)
    }
    fro_sq[fro_sq < 0] = 0
    loss = 0.5 * sqrt(fro_sq)

    pos = which(alpha > 0)
    if (length(pos) > 0) {
        if (anyNA(corr1) || anyNA(corr2)) {
            loss[pos] = NA
        } else {
            epsilon = 1e-10  # Small constant to prevent division by zero
            weight = 1/abs(corr2 + epsilon)
            weight[corr2 == 0] = 0
            diag(weight) = 0
            w_abs = abs(as.vector(weight)[o])
            if (soft) {
                suf_wa = c(rev(cumsum(rev(w_abs * a_abs))), 0)
                suf_w = c(rev(cumsum(rev(w_abs))), 0)
                penalty = suf_wa[k[pos] + 1] - thresh[pos] * suf_w[k[pos] + 1]
            } else {
                suf_wa = c(rev(cumsum(rev(w_abs * a_abs))), 0)
                penalty = suf_wa[k[pos] + 1]
            }
            loss[pos] = loss[pos] + alpha[pos] * penalty
        }
    }
    return(loss)
}

.sparse_linear = function(mat, wins_quant, method, soft, alpha_grid,
                          thresh_len, n_cv, thresh_hard, max_p) {
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
        DescTools::Winsorize(x, val = stats::quantile(x, probs = wins_quant, na.rm = TRUE)))

    # Co-occurrence matrix
    mat_occur = mat
    mat_occur[mat_occur != 0] = 1
    mat_occur[mat_occur == 0] = 0
    mat_occur[is.na(mat_occur)] = 0

    # The number of samples in which both taxa are present is the cross-product
    # of the presence indicators
    mat_cooccur = crossprod(mat_occur)
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

    # Regularization of the covariance matrix
    if(method == "spearman"){
        # Convert to rank
        mat = apply(mat,2,function(x) {
            r = rank(x, na.last = NA)
            x[!is.na(x)] = r
            return(x)
        }
        )
    }
    # Covariance matrix
    cov_mat = stats::cov(mat, use = "pairwise.complete.obs")
    cov_mat[is.na(cov_mat)] = 0

    # Regularize the covariance matrix
    cov_mat_pos = .regularize_eigenvalues(cov_mat)
    cov_mat_pos[mat_cooccur < 2] = 0
    cov_mat_pos[is.infinite(cov_mat_pos)] = 0

    # Check if it is positive semi-definite, if not repeat the regularization process
    while(!.is_psd(cov_mat_pos)) {
        cov_mat_pos = .regularize_eigenvalues(cov_mat_pos)
        cov_mat_pos[mat_cooccur < 2] = 0
        cov_mat_pos[is.infinite(cov_mat_pos)] = 0
    }
    # Convert to correlation coefficient
    corr_reg = cov2cor(cov_mat_pos)

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
    if (is.null(alpha_grid)) alpha_grid = 0
    param_grid = expand.grid(thresh = thresh_grid, alpha = alpha_grid)

    # The splits are evaluated in an environment detached from the package
    # namespace. The loss over the whole parameter grid is obtained in one pass
    # per split.
    loop_env = .detached_env(mat = mat, n = n, n1 = n1, method = method,
                             soft = soft, n_cv = n_cv,
                             grid_thresh = param_grid$thresh,
                             grid_alpha = param_grid$alpha,
                             thresh_loss_grid = .thresh_loss_grid)
    environment(loop_env$thresh_loss_grid) = loop_env

    loss_mat = eval(quote(
        foreach(i = seq_len(n_cv), .combine = rbind) %dorng% {
        # Create training and validation splits
        index = sample(seq_len(n), size = n1, replace = FALSE)
        mat1 = mat[index,]
        mat2 = mat[-index,]

        # Correlations depend on the split only, not on the thresholding or
        # penalty grid
        corr1 = stats::cor(mat1, method = method, use = "pairwise.complete.obs")
        corr2 = stats::cor(mat2, method = method, use = "pairwise.complete.obs")

        thresh_loss_grid(corr1 = corr1, corr2 = corr2, soft = soft,
                         thresh = grid_thresh, alpha = grid_alpha)
        }), loop_env)

    # Calculate mean loss across CV folds
    mean_losses = colMeans(loss_mat)

    # Find optimal parameters
    opt_index = which.min(mean_losses)
    thresh_opt = param_grid$thresh[opt_index]
    alpha_opt = param_grid$alpha[opt_index]

    # Apply optimal thresholding
    corr = stats::cor(mat, method = method, use = "pairwise.complete.obs")
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
    result = list(cv_error = mean_losses,
                  thresh_grid = thresh_grid,
                  thresh_opt = thresh_opt,
                  alpha_opt = alpha_opt,
                  mat_cooccur = mat_cooccur,
                  corr = corr,
                  corr_p = corr_p,
                  corr_th = corr_th,
                  corr_fl = corr_fl,
                  corr_reg = corr_reg)
    return(result)
}
