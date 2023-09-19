#' @title Sparse estimation of linear correlations among microbiomes
#'
#' @description Obtain the sparse correlation matrix for linear correlations
#' between taxa. The current version of \code{secom_linear} function supports
#' either of the three correlation coefficients: Pearson, Spearman, and
#' Kendall's \eqn{\tau}.
#'
#' @param data a list of the input data. The \code{data} parameter should be a
#' list containing input data objects, which can be either \code{phyloseq}
#' or \code{TreeSummarizedExperiment} objects. Each object within the list
#' consists of a feature table (microbial count table), a sample metadata table,
#' a taxonomy table (optional), and a phylogenetic tree (optional).
#' Ensure that the row names of the metadata table match the sample names in the
#' feature table, and the row names of the taxonomy table match the taxon
#' (feature) names in the feature table. For detailed information, refer to
#' \code{?phyloseq::phyloseq} or
#' \code{?TreeSummarizedExperiment::TreeSummarizedExperiment}.
#' It is recommended to use low taxonomic levels, such as OTU or species level,
#' as the estimation of sampling fractions requires a large number of taxa.
#' If working with multiple ecosystems, such as gut and tongue, stack the data
#' by specifying the list of input data as
#' \code{data = list(gut = tse1, tongue = tse2)}.
#' @param assay_name character. Name of the feature table within the data object
#' (only applicable if the data object is a \code{(Tree)SummarizedExperiment}).
#' Default is "counts".
#' See \code{?SummarizedExperiment::assay} for more details.
#' @param assay.type alias for \code{assay_name}.
#' @param tax_level character. The taxonomic level of interest. The input data
#' can be agglomerated at different taxonomic levels based on your research
#' interest. Default is NULL, i.e., do not perform agglomeration, and the
#' SECOM anlysis will be performed at the lowest taxonomic level of the
#' input \code{data}.
#' @param rank alias for \code{tax_level}.
#' @param pseudo numeric. Add pseudo-counts to the data.
#' Default is 0 (no pseudo-counts).
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' (the proportion of samples in which the taxon is present)
#' less than \code{prv_cut} will be excluded in the analysis. For example,
#' if there are 100 samples, and a taxon has nonzero counts present in less than
#' 100*prv_cut samples, it will not be considered in the analysis.
#' Default is 0.50.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis. Default is 1000.
#' @param corr_cut numeric. To avoid false positives caused by taxa with small
#' variances, taxa with Pearson correlation coefficients greater than
#' \code{corr_cut} with the estimated sample-specific bias will be flagged.
#' When taxa are flagged, the pairwise correlation coefficient between them will
#' be set to 0s. Default is 0.5.
#' @param wins_quant a numeric vector of probabilities with values between
#' 0 and 1. Replace extreme values in the abundance data with less
#' extreme values. Default is \code{c(0.05, 0.95)}. For details,
#' see \code{?DescTools::Winsorize}.
#' @param method character. It indicates which correlation coefficient is to be
#' computed. It can be either "pearson" or "spearman".
#' @param soft logical. \code{TRUE} indicates that soft thresholding is applied
#' to achieve the sparsity of the correlation matrix. \code{FALSE} indicates
#' that hard thresholding is applied to achieve the sparsity of the correlation
#' matrix. Default is \code{FALSE}.
#' @param thresh_len numeric. Grid-search is implemented to find the optimal
#' values over \code{thresh_len} thresholds for the thresholding operator.
#' Default is 100.
#' @param n_cv numeric. The fold number in cross validation.
#' Default is 10 (10-fold cross validation).
#' @param thresh_hard Numeric. Pairwise correlation coefficients
#' (in their absolute value) that are less than or equal to \code{thresh_hard}
#' will be set to 0. Default is 0.3.
#' @param max_p numeric. Obtain the sparse correlation matrix by
#' p-value filtering. Pairwise correlation coefficients with p-value greater
#' than \code{max_p} will be set to 0s. Default is 0.005.
#' @param n_cl numeric. The number of nodes to be forked. For details, see
#' \code{?parallel::makeCluster}. Default is 1 (no parallel computing).
#'
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{s_diff_hat}, a numeric vector of estimated
#'         sample-specific biases.}
#'         \item{ \code{y_hat}, a matrix of bias-corrected abundances}
#'         \item{ \code{cv_error}, a numeric vector of cross-validation error
#'         estimates, which are the Frobenius norm differences between
#'         correlation matrices using training set and validation set,
#'         respectively.}
#'         \item{ \code{thresh_grid}, a numeric vector of thresholds
#'         in the cross-validation.}
#'         \item{ \code{thresh_opt}, numeric. The optimal threshold through
#'         cross-validation.}
#'         \item{ \code{mat_cooccur}, a matrix of taxon-taxon co-occurrence
#'         pattern. The number in each cell represents the number of complete
#'         (nonzero) samples for the corresponding pair of taxa.}
#'         \item{ \code{corr}, the sample correlation matrix (using the measure
#'         specified in \code{method}) computed using the bias-corrected
#'         abundances \code{y_hat}.}
#'         \item{ \code{corr_p}, the p-value matrix corresponding to the sample
#'         correlation matrix \code{corr}.}
#'         \item{ \code{corr_th}, the sparse correlation matrix obtained by
#'         thresholding based on the method specified in \code{soft}.}
#'         \item{ \code{corr_fl}, the sparse correlation matrix obtained by
#'         p-value filtering based on the cutoff specified in \code{max_p}.}
#'         }
#'
#' @seealso \code{\link{secom_dist}}
#'
#' @examples
#' library(ANCOMBC)
#' data(dietswap, package = "microbiome")
#' tse = mia::makeTreeSummarizedExperimentFromPhyloseq(dietswap)
#'
#' # subset to baseline
#' tse = tse[, tse$timepoint == 1]
#'
#' set.seed(123)
#' res_linear = secom_linear(data = list(tse), assay_name = "counts",
#'                           tax_level = "Phylum", pseudo = 0,
#'                           prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
#'                           wins_quant = c(0.05, 0.95), method = "pearson",
#'                           soft = FALSE, thresh_len = 20, n_cv = 10,
#'                           thresh_hard = 0.3, max_p = 0.005, n_cl = 2)
#'
#' corr_th = res_linear$corr_th
#' corr_fl = res_linear$corr_fl
#'
#' @author Huang Lin
#'
#' @importFrom mia makeTreeSummarizedExperimentFromPhyloseq taxonomyRanks agglomerateByRank
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom energy dcor dcor.test
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng%
#' @importFrom gtools smartbind
#' @importFrom Hmisc rcorr
#' @importFrom DescTools Winsorize
#' @importFrom Rdpack reprompt
#'
#' @export
secom_linear = function(data, assay.type = assay_name, assay_name = "counts",
                        rank = tax_level, tax_level = NULL,
                        pseudo = 0, prv_cut = 0.5, lib_cut = 1000,
                        corr_cut = 0.5, wins_quant = c(0.05, 0.95),
                        method = c("pearson", "spearman"),
                        soft = FALSE, thresh_len = 100, n_cv = 10,
                        thresh_hard = 0, max_p = 0.005, n_cl = 1) {
    # Check for aliases
    if (!is.null(assay.type)) {
        assay_name = assay.type
    }

    if (!is.null(rank)) {
        tax_level = rank
    }

    # ===========Sampling fraction and absolute abundance estimation============
    if (length(data) == 1) {
        tse_obj = .tse_construct(data = data[[1]], assay_name = assay_name[1],
                                 tax_level = tax_level[1], phyloseq = NULL)
        abn_list = .abn_est(tse = tse_obj$tse, tax_level = tse_obj$tax_level,
                            assay_name = tse_obj$assay_name, pseudo = pseudo,
                            prv_cut = prv_cut, lib_cut = lib_cut)
        s_diff_hat = abn_list$s_diff_hat
        y_hat = abn_list$y_hat
    } else {
        if (is.null(names(data))) names(data) = paste0("data", seq_along(data))

        # Check common samples
        samp_names = lapply(data, function(x) colnames(x))
        samp_common = Reduce(intersect, samp_names)
        samp_txt = sprintf(paste0("Number of common samples ",
                                  "across datasets: ",
                                  length(samp_common)))
        message(samp_txt)
        if (length(samp_common) < 10) {
            stop_txt = paste0("Insufficient common samples: ",
                              "Multi-dataset computation not recommended")
            stop(stop_txt)
        }

        # Rename taxa
        tse_list = lapply(seq_along(data), function(i) {
            tse_obj = .tse_construct(data = data[[i]], assay_name = assay_name[i],
                                     tax_level = tax_level[i], phyloseq = NULL)
            return(tse_obj)
        })

        for (i in seq_along(tse_list)) {
            rownames(SingleCellExperiment::altExp(tse_list[[i]]$tse,
                                                  tse_list[[i]]$tax_level)) =
                paste(names(data)[[i]],
                      rownames(SingleCellExperiment::altExp(tse_list[[i]]$tse,
                                                            tse_list[[i]]$tax_level)),
                      sep = " - ")
        }

        abn_list = lapply(seq_along(tse_list), function(i) {
            .abn_est(tse = tse_list[[i]]$tse,
                     tax_level = tse_list[[i]]$tax_level,
                     assay_name = assay_name[i], pseudo = pseudo,
                     prv_cut = prv_cut, lib_cut = lib_cut)
        })
        s_diff_hat = lapply(abn_list, function(x) x$s_diff_hat)
        y_hat = do.call(gtools::smartbind, lapply(abn_list, function(x) as.data.frame(x$y_hat)))
        y_hat_rownames = do.call(c, lapply(abn_list, function(x) rownames(x$y_hat)))
        y_hat = as.matrix(y_hat)
        rownames(y_hat) = y_hat_rownames
    }

    # =================Sparse estimation on linear correlations=================
    if (n_cl > 1) {
      cl = parallel::makeCluster(n_cl)
      doParallel::registerDoParallel(cl)
    } else {
      foreach::registerDoSEQ()
    }

    if (method %in% c("pearson", "spearman")) {
        res_corr = .sparse_linear(mat = t(y_hat), wins_quant, method, soft,
                                  thresh_len, n_cv, thresh_hard, max_p)
    } else {
        stop_txt = paste0("The specified correlation coefficient type is not valid \n",
                          "Please choose either 'pearson' or 'spearman' as the type of correlation coefficient")
        stop(stop_txt, call. = FALSE)
    }

    if (n_cl > 1) {
      parallel::stopCluster(cl)
    }

    # To prevent FP from taxa with extremely small variances
    if (length(data) == 1) {
        corr_s = cor(cbind(s_diff_hat, t(y_hat)),
                     use = "pairwise.complete.obs")[1, -1]
        fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
        fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
        fp_ind = (fp_ind1 * fp_ind2 == 1)
        diag(fp_ind) = FALSE
        res_corr$corr[fp_ind] = 0
        res_corr$corr_th[fp_ind] = 0
        res_corr$corr_fl[fp_ind] = 0
        res_corr$corr_p[fp_ind] = 1
    } else {
        for (i in seq_along(data)) {
            df_s = data.frame(sample_id = names(s_diff_hat[[i]]),
                              s = s_diff_hat[[i]])
            rownames(df_s) = NULL
            df_y = data.frame(sample_id = rownames(t(y_hat)), t(y_hat),
                              check.names = FALSE)
            rownames(df_y) = NULL

            df_merge = df_y
            df_merge$s = df_s$s[match(df_y$sample_id, df_s$sample_id)]
            df_merge$sample_id = NULL
            df_merge = df_merge[c('s', setdiff(names(df_merge), 's'))]

            corr_s = cor(df_merge, use = "pairwise.complete.obs")[1, -1]
            fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
            fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
            fp_ind = (fp_ind1 * fp_ind2 == 1)
            diag(fp_ind) = FALSE
            res_corr$corr[fp_ind] = 0
            res_corr$corr_th[fp_ind] = 0
            res_corr$corr_fl[fp_ind] = 0
            res_corr$corr_p[fp_ind] = 1
        }
    }

    # ==================================Outputs=================================
    res = c(list(s_diff_hat = s_diff_hat, y_hat = y_hat), res_corr)
    return(res)
}
