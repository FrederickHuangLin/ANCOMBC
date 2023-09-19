# Construct TSE object
.tse_construct = function(data, assay_name, tax_level, phyloseq) {
    if (!is.null(data)) {
        # Check data types
        if (!inherits(data, c("phyloseq", "SummarizedExperiment",
                              "TreeSummarizedExperiment"))) {
            stop_txt = paste0("The input data should be in one of the ",
                              "following format: ",
                              "`phyloseq`, `SummarizedExperiment`, ",
                              "`TreeSummarizedExperiment`")
            stop(stop_txt, call. = FALSE)
        }

        if (inherits(data, "phyloseq")) {
            tse = mia::makeTreeSummarizedExperimentFromPhyloseq(data)
            assay_name = "counts"
        } else {
            tse = data
            assay_name = assay_name
        }

        # Check feature metadata
        if (ncol(SummarizedExperiment::rowData(tse)) == 0) {
            tax_tab = matrix(rownames(tse), ncol = 1)
            rownames(tax_tab) = rownames(tse)
            colnames(tax_tab) = c("Species")
            tax_tab = S4Vectors::DataFrame(tax_tab)
            SummarizedExperiment::rowData(tse) = tax_tab
        }

        # Check if agglomeration should be performed
        if (is.null(tax_level)) {
            tax_level = "ASV"
            tse_alt = tse
        } else {
            tse_alt = .merge_features(tse, tax_level)
        }
        SingleCellExperiment::altExp(tse, tax_level) = tse_alt
    } else if (!is.null(phyloseq)) {
        tse = mia::makeTreeSummarizedExperimentFromPhyloseq(phyloseq)
        assay_name = "counts"

        if (ncol(SummarizedExperiment::rowData(tse)) == 0) {
            tax_tab = matrix(rownames(tse), ncol = 1)
            rownames(tax_tab) = rownames(tse)
            colnames(tax_tab) = c("Species")
            tax_tab = S4Vectors::DataFrame(tax_tab)
            SummarizedExperiment::rowData(tse) = tax_tab
        }

        if (is.null(tax_level)) {
            tax_level = "ASV"
            tse_alt = tse
        } else {
            tse_alt = .merge_features(tse, tax_level)
        }
        SingleCellExperiment::altExp(tse, tax_level) = tse_alt
    } else {
        stop_txt = paste0("The input data are missing. ",
                          "Please specify either `data` or `phyloseq`")
        stop(stop_txt, call. = FALSE)
    }

    tse_obj = list(tse = tse, assay_name = assay_name,
                   tax_level = tax_level, tse_alt = tse_alt)

    return(tse_obj)
}

# Filter data by prevalence and library size
.data_core = function(tse = tse, tax_level, assay_name = assay_name,
                      alt = FALSE, prv_cut, lib_cut,
                      tax_keep = NULL, samp_keep = NULL) {
    if (alt) {
        tse_alt = SingleCellExperiment::altExp(tse, tax_level)
        feature_table = SummarizedExperiment::assay(tse_alt, assay_name)
        meta_data = SummarizedExperiment::colData(tse_alt)
    } else {
        feature_table = SummarizedExperiment::assay(tse, assay_name)
        meta_data = SummarizedExperiment::colData(tse)
    }

    # Discard taxa with prevalences < prv_cut
    if (is.null(tax_keep)) {
        prevalence = apply(feature_table, 1, function(x)
            sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
        tax_keep = which(prevalence >= prv_cut)
    }else if (length(tax_keep) == 0) {
        stop("All taxa contain structural zeros", call. = FALSE)
    } else {
        # Discard taxa with structural zeros
        feature_table = feature_table[tax_keep, , drop = FALSE]
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

    output = list(feature_table = feature_table,
                  meta_data = meta_data,
                  tax_keep = tax_keep,
                  samp_keep = samp_keep)
    return(output)
}

# Metadata and arguments check
.data_qc = function(meta_data, formula, group, struc_zero,
                    global, pairwise, dunnet,
                    mdfdr_control, trend, trend_control) {
    # Drop unused levels
    meta_data[] = lapply(meta_data, function(x)
        if(is.factor(x)) factor(x) else x)

    # Check if all covariates specified in the formula are columns in meta_data
    vars = unlist(strsplit(formula, split = "\\s*\\+\\s*"))
    missing_vars = vars[!vars %in% colnames(meta_data)]
    if(length(missing_vars) > 0) {
        stop("The following variables specified are not in the meta data: ",
             paste(missing_vars, collapse = ", "))
    }

    # Check the group variable
    if (is.null(group)) {
        if (any(c(global, pairwise, dunnet, trend))) {
            stop_txt = paste0("Group variable is required for the multi-group comparison \n",
                             "`group` is `NULL` while some of the arguments ",
                             "(`global`, `pairwise`, `dunnet`, `trend`) are `TRUE`")
            stop(stop_txt, call. = FALSE)
        }
        if (struc_zero) {
            stop_txt = paste0("Please specify the group variable for detecting structural zeros \n",
                              "Otherwise, set `struc_zero = FALSE` to proceed")
            stop(stop_txt, call. = FALSE)
        }
    } else {
        meta_data[, group] = as.factor(meta_data[, group])
        # Check the number of groups
        n_level = nlevels(meta_data[, group])
        if (n_level < 2) {
            stop("The group variable should have >= 2 categories",
                 call. = FALSE)
        } else if (n_level < 3) {
            global = FALSE
            pairwise = FALSE
            dunnet = FALSE
            trend = FALSE
            warn_txt = paste0("The group variable has < 3 categories \n",
                              "The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated")
            warning(warn_txt, call. = FALSE)
        }

        # Check the mdfdr setting for pairwise and dunnet's tests
        if (pairwise | dunnet) {
            if (is.null(mdfdr_control)) {
                stop("Please specify `mdfdr_control` for pairwise or dunnet's test",
                     call. = FALSE)
            }
        }

        # Check contrast matrices and nodes for trend test
        if (trend) {
            if (is.null(trend_control)) {
                stop("Please specify the `trend_control` parameter for the trend test.",
                     call. = FALSE)
            }
            if (is.null(trend_control$contrast)) {
                stop("Please specify the contrast matrices for the trend test.",
                     call. = FALSE)
            }
            if (is.null(trend_control$node)) {
                stop("Please specify the nodes for the trend test",
                     call. = FALSE)
            }
            if (length(trend_control$contrast) != length(trend_control$node)) {
                stop("The number of nodes should match the number of contrast matrices",
                     call. = FALSE)
            }
            sq_mat_check = vapply(trend_control$contrast, function(x)
                nrow(x) == ncol(x), FUN.VALUE = logical(1))
            if (any(sq_mat_check == FALSE)) {
                stop("The contrast matrices for the trend test should be square matrices",
                     call. = FALSE)
            }
            dim_mat_check = vapply(trend_control$contrast, function(x)
                nrow(x), FUN.VALUE = integer(1))
            if (any(dim_mat_check != (n_level - 1))) {
                stop_txt = paste0("The contrast matrices for the trend test should be square matrices ",
                                  "with dimension #group - 1 \n",
                                  "The number of groups in current data is: ",
                                  n_level)

                stop(stop_txt, call. = FALSE)
            }

            n_trend = length(trend_control$contrast)
            if (is.null(names(trend_control$contrast))) {
                names(trend_control$contrast) = paste0("trend", seq_len(n_trend))
                names(trend_control$node) = paste0("trend", seq_len(n_trend))
            }
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
                                     "Variance estimation would be unstable when the sample size is < 5 per group",
                                     sep = "\n"))
            warning(warn_txt, call. = FALSE)
        }
    }

    output = list(meta_data = meta_data,
                  global = global,
                  pairwise = pairwise,
                  dunnet = dunnet,
                  trend = trend,
                  trend_control = trend_control)
    return(output)
}

# Identify structural zeros
.get_struc_zero = function(tse, tax_level, assay_name,
                           alt = FALSE, group, neg_lb) {
    if (alt) {
        tse_alt = SingleCellExperiment::altExp(tse, tax_level)
        feature_table = SummarizedExperiment::assay(tse_alt, assay_name)
        meta_data = SummarizedExperiment::colData(tse_alt)
        tax_name = rownames(tse_alt)
    } else {
        feature_table = SummarizedExperiment::assay(tse, assay_name)
        meta_data = SummarizedExperiment::colData(tse)
        tax_name = rownames(tse)
    }
    group_data = factor(meta_data[, group])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1
    n_tax = nrow(feature_table)
    n_group = nlevels(group_data)

    p_hat = matrix(NA, nrow = n_tax, ncol = n_group)
    rownames(p_hat) = rownames(feature_table)
    colnames(p_hat) = levels(group_data)
    for (i in seq_len(n_tax)) {
        p_hat[i, ] = tapply(present_table[i, ], group_data,
                            function(x) mean(x, na.rm = TRUE))
    }

    samp_size = matrix(NA, nrow = n_tax, ncol = n_group)
    rownames(samp_size) = rownames(feature_table)
    colnames(samp_size) = levels(group_data)
    for (i in seq_len(n_tax)) {
        samp_size[i, ] = tapply(as.matrix(feature_table)[i, ], group_data,
                                function(x) length(x[!is.na(x)]))
    }

    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

    output = (p_hat == 0)
    # Shall we classify a taxon as a structural zero by its negative lower bound?
    if (neg_lb) output[p_hat_lo <= 0] = TRUE

    output = cbind(tax_name, output)
    colnames(output) = c("taxon",
                         paste0("structural_zero (", group,
                                " = ", colnames(output)[-1], ")"))
    output = data.frame(output, check.names = FALSE, row.names = NULL)
    output[, -1] = apply(output[, -1], 2, as.logical)
    return(output)
}
