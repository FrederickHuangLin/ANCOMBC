#' @title Analysis of Compositions of Microbiomes with Bias Correction
#' (ANCOM-BC)
#'
#' @description Determine taxa whose absolute abundances, per unit volume, of
#' the ecosystem (e.g., gut) are significantly different with changes in the
#' covariate of interest (e.g., group). The current version of
#' \code{ancombc} function implements Analysis of Compositions of Microbiomes
#' with Bias Correction (ANCOM-BC) in cross-sectional data while allowing
#' for covariate adjustment.
#'
#' @details A taxon is considered to have structural zeros in some (>=1)
#' groups if it is completely (or nearly completely) missing in these groups.
#' For instance, suppose there are three groups: g1, g2, and g3.
#' If the counts of taxon A in g1 are 0 but nonzero in g2 and g3,
#' then taxon A will be considered to contain structural zeros in g1.
#' In this example, taxon A is declared to be differentially abundant between
#' g1 and g2, g1 and g3, and consequently, it is globally differentially
#' abundant with respect to this group variable.
#' Such taxa are not further analyzed using ANCOM-BC, but the results are
#' summarized in the overall summary. For more details about the structural
#' zeros, please go to the
#' \href{https://doi.org/10.3389/fmicb.2017.02114}{ANCOM-II} paper.
#' Setting \code{neg_lb = TRUE} indicates that you are using both criteria
#' stated in section 3.2 of
#' \href{https://doi.org/10.3389/fmicb.2017.02114}{ANCOM-II}
#' to detect structural zeros; otherwise, the algorithm will only use the
#' equation 1 in section 3.2 for declaring structural zeros. Generally, it is
#' recommended to set \code{neg_lb = TRUE} when the sample size per group is
#' relatively large (e.g. > 30).
#'
#' @param data the input data. The \code{data} parameter should be either a
#' \code{phyloseq} or a \code{TreeSummarizedExperiment} object, which
#' consists of a feature table (microbial count table), a sample metadata table,
#' a taxonomy table (optional), and a phylogenetic tree (optional).
#' Ensure that the row names of the metadata table match the sample names in the
#' feature table, and the row names of the taxonomy table match the taxon
#' (feature) names in the feature table. For detailed information, refer to
#' \code{?phyloseq::phyloseq} or
#' \code{?TreeSummarizedExperiment::TreeSummarizedExperiment}.
#' @param assay_name character. Name of the count table in the data object
#' (only applicable if data object is a \code{(Tree)SummarizedExperiment}).
#' Default is "counts".
#' See \code{?SummarizedExperiment::assay} for more details.
#' @param assay.type alias for \code{assay_name}.
#' @param tax_level character. The taxonomic level of interest. The input data
#' can be agglomerated at different taxonomic levels based on your research
#' interest. Default is NULL, i.e., do not perform agglomeration, and the
#' ANCOM-BC anlysis will be performed at the lowest taxonomic level of the
#' input \code{data}.
#' @param rank alias for \code{tax_level}.
#' @param phyloseq a \code{phyloseq} object. Will be deprecated.
#' @param formula the character string expresses how microbial absolute
#' abundances for each taxon depend on the variables in metadata. When
#' specifying the \code{formula}, make sure to include the \code{group} variable
#' in the formula if it is not NULL.
#' @param p_adj_method character. method to adjust p-values. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{?stats::p.adjust} for more details.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' (the proportion of samples in which the taxon is present)
#' less than \code{prv_cut} will be excluded in the analysis. For example,
#' if there are 100 samples, and a taxon has nonzero counts present in less than
#' 100*prv_cut samples, it will not be considered in the analysis.
#' Default is 0.10.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis. Default is 0, i.e. do not discard any sample.
#' @param group character. the name of the group variable in metadata.
#' The \code{group} parameter should be a character string representing the name
#' of the group variable in the metadata. The \code{group} variable should be
#' discrete, meaning it consists of categorical values. Specifying the
#' \code{group} variable is required if you are interested in detecting
#' structural zeros and performing global tests. However, if these analyses are
#' not of interest to you, you can leave the \code{group} parameter as NULL.
#' If the \code{group} variable of interest contains only two categories, you
#' can also leave the \code{group} parameter as NULL. Default is NULL.
#' @param struc_zero logical. whether to detect structural zeros based on
#' \code{group}. Default is FALSE.
#' @param neg_lb logical. whether to classify a taxon as a structural zero using
#' its asymptotic lower bound. Default is FALSE.
#' @param tol numeric. the iteration convergence tolerance for the E-M
#' algorithm. Default is 1e-05.
#' @param max_iter numeric. the maximum number of iterations for the E-M
#' algorithm. Default is 100.
#' @param conserve logical. whether to use a conservative variance estimator for
#' the test statistic. It is recommended if the sample size is small and/or
#' the number of differentially abundant taxa is believed to be large.
#' Default is FALSE.
#' @param alpha numeric. level of significance. Default is 0.05.
#' @param global logical. whether to perform the global test. Default is FALSE.
#' @param n_cl numeric. The number of nodes to be forked. For details, see
#' \code{?parallel::makeCluster}. Default is 1 (no parallel computing).
#' @param verbose logical. Whether to generate verbose output during the
#' ANCOM-BC fitting process. Default is FALSE.
#'
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{feature_table}, a \code{data.frame} of pre-processed
#'         (based on \code{prv_cut} and \code{lib_cut}) microbial count table.}
#'         \item{ \code{zero_ind}, a logical \code{data.frame} with TRUE
#'         indicating the taxon is detected to contain structural zeros in
#'         some specific groups.}
#'         \item{ \code{samp_frac}, a numeric vector of estimated sampling
#'         fractions in log scale (natural log).}
#'         \item{ \code{delta_em}, estimated sample-specific biases
#'         through E-M algorithm.}
#'         \item{ \code{delta_wls}, estimated sample-specific biases through
#'         weighted least squares (WLS) algorithm.}
#'         \item{ \code{res},  a \code{list} containing ANCOM-BC primary result,
#'         which consists of:}
#'         \itemize{
#'         \item{ \code{lfc}, a \code{data.frame} of log fold changes
#'         obtained from the ANCOM-BC log-linear (natural log) model.}
#'         \item{ \code{se}, a \code{data.frame} of standard errors (SEs) of
#'         \code{lfc}.}
#'         \item{ \code{W}, a \code{data.frame} of test statistics.
#'         \code{W = lfc/se}.}
#'         \item{ \code{p_val}, a \code{data.frame} of p-values. P-values are
#'         obtained from two-sided Z-test using the test statistic \code{W}.}
#'         \item{ \code{q_val}, a \code{data.frame} of adjusted p-values.
#'         Adjusted p-values are obtained by applying \code{p_adj_method}
#'         to \code{p_val}.}
#'         \item{ \code{diff_abn}, a logical \code{data.frame}. TRUE if the
#'         taxon has \code{q_val} less than \code{alpha}.}
#'         }
#'         \item{ \code{res_global},  a \code{data.frame} containing ANCOM-BC
#'         global test result for the variable specified in \code{group},
#'         each column is:}
#'         \itemize{
#'         \item{ \code{W}, test statistics.}
#'         \item{ \code{p_val}, p-values, which are obtained from two-sided
#'         Chi-square test using \code{W}.}
#'         \item{ \code{q_val}, adjusted p-values. Adjusted p-values are
#'         obtained by applying \code{p_adj_method} to \code{p_val}.}
#'         \item{ \code{diff_abn}, A logical vector. TRUE if the taxon has
#'         \code{q_val} less than \code{alpha}.}
#'         }
#'         }
#'
#' @seealso \code{\link{ancom}} \code{\link{ancombc2}}
#'
#' @examples
#' #===========Build a TreeSummarizedExperiment Object from Scratch=============
#' library(mia)
#'
#' # microbial count table
#' otu_mat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
#' rownames(otu_mat) = paste0("taxon", 1:nrow(otu_mat))
#' colnames(otu_mat) = paste0("sample", 1:ncol(otu_mat))
#' assays = SimpleList(counts = otu_mat)
#'
#' # sample metadata
#' smd = data.frame(group = sample(LETTERS[1:4], size = 10, replace = TRUE),
#'                  row.names = paste0("sample", 1:ncol(otu_mat)),
#'                  stringsAsFactors = FALSE)
#' smd = DataFrame(smd)
#'
#' # taxonomy table
#' tax_tab = matrix(sample(letters, 70, replace = TRUE),
#'                  nrow = nrow(otu_mat), ncol = 7)
#' rownames(tax_tab) = rownames(otu_mat)
#' colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
#'                       "Family", "Genus", "Species")
#' tax_tab = DataFrame(tax_tab)
#'
#' # create TSE
#' tse = TreeSummarizedExperiment(assays = assays,
#'                                colData = smd,
#'                                rowData = tax_tab)
#'
#' # convert TSE to phyloseq
#' pseq = makePhyloseqFromTreeSummarizedExperiment(tse)
#'
#' #========================Run ANCOMBC Using a Real Data=======================
#' library(ANCOMBC)
#' data(atlas1006, package = "microbiome")
#' tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)
#'
#' # subset to baseline
#' tse = tse[, tse$time == 0]
#'
#' # run ancombc function
#' set.seed(123)
#' out = ancombc(data = tse, assay_name = "counts",
#'               tax_level = "Family", phyloseq = NULL,
#'               formula = "age + nationality + bmi_group",
#'               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#'               group = "bmi_group", struc_zero = TRUE, neg_lb = FALSE,
#'               tol = 1e-5, max_iter = 100, conserve = TRUE,
#'               alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
#'
#' res_prim = out$res
#' res_global = out$res_global
#'
#' # to run ancombc using the phyloseq object
#' tse_alt = agglomerateByRank(tse, "Family")
#' pseq = makePhyloseqFromTreeSummarizedExperiment(tse_alt)
#' set.seed(123)
#' out = ancombc(data = NULL, assay_name = NULL,
#'               tax_level = "Family", phyloseq = pseq,
#'               formula = "age + nationality + bmi_group",
#'               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#'               group = "bmi_group", struc_zero = TRUE, neg_lb = FALSE,
#'               tol = 1e-5, max_iter = 100, conserve = TRUE,
#'               alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
#'
#' @author Huang Lin
#'
#' @references
#' \insertRef{kaul2017analysis}{ANCOMBC}
#'
#' \insertRef{lin2020analysis}{ANCOMBC}
#'
#' @rawNamespace import(stats, except = filter)
#' @importFrom mia makeTreeSummarizedExperimentFromPhyloseq taxonomyRanks agglomerateByRank
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng%
#' @importFrom MASS ginv
#' @importFrom nloptr neldermead
#' @importFrom Rdpack reprompt
#'
#' @export
ancombc = function(data = NULL, assay.type = NULL, assay_name = "counts",
                   rank = NULL, tax_level = NULL, phyloseq = NULL,
                   formula, p_adj_method = "holm", prv_cut = 0.10,
                   lib_cut = 0, group = NULL, struc_zero = FALSE,
                   neg_lb = FALSE, tol = 1e-05, max_iter = 100,
                   conserve = FALSE, alpha = 0.05, global = FALSE,
                   n_cl = 1, verbose = FALSE){
    message("'ancombc' has been fully evolved to 'ancombc2'. \n",
            "Explore the enhanced capabilities of our refined method!")

    if (n_cl > 1) {
      cl = parallel::makeCluster(n_cl)
      doParallel::registerDoParallel(cl)
    } else {
      foreach::registerDoSEQ()
    }

    # 1. Data pre-processing
    # Check for aliases
    if (!is.null(assay.type)) {
        assay_name = assay.type
    }

    if (!is.null(rank)) {
        tax_level = rank
    }

    # TSE data construction
    tse_obj = .tse_construct(data = data, assay_name = assay_name,
                             tax_level = tax_level, phyloseq = phyloseq)
    tse = tse_obj$tse
    assay_name = tse_obj$assay_name
    tax_level = tse_obj$tax_level

    # Filter data by prevalence and library size
    core = .data_core(tse = tse, tax_level = tax_level, assay_name = assay_name,
                      alt = TRUE, prv_cut = prv_cut, lib_cut = lib_cut,
                      tax_keep = NULL, samp_keep = NULL)
    feature_table = core$feature_table
    meta_data = core$meta_data
    tax_keep = core$tax_keep
    n_tax = nrow(feature_table)
    tax_name = rownames(feature_table)
    if (n_tax < 10) {
        warn_txt = sprintf(paste("ANCOM-BC results would be unreliable when the number of taxa is too small (e.g. < 10)",
                                 "The number of taxa in the current dataset is: ",
                                 n_tax, sep = "\n"))
        warning(warn_txt, call. = FALSE)
    }

    # Metadata and arguments check
    qc = .data_qc(meta_data = meta_data,
                  formula = formula, group = group,
                  struc_zero = struc_zero, global = global,
                  pairwise = FALSE, dunnet = FALSE,
                  mdfdr_control = NULL, trend = FALSE, trend_control = NULL)
    meta_data = qc$meta_data
    global = qc$global

    # Add pseudocount (1) and take logarithm.
    y = log(feature_table + 1)
    options(na.action = "na.pass") # Keep NA's in rows of x
    x = model.matrix(formula(paste0("~", formula)), data = meta_data)
    options(na.action = "na.omit") # Switch it back
    covariates = colnames(x)
    n_covariates = length(covariates)

    # 2. Identify taxa with structural zeros
    if (struc_zero) {
        if (is.null(group)) {
            stop_txt = paste("Please specify the group variable for",
                             "detecting structural zeros.",
                             "Otherwise, set struc_zero = FALSE to proceed")
            stop(stop_txt, call. = FALSE)
        }
        zero_ind = .get_struc_zero(tse = tse, tax_level = tax_level,
                                   assay_name = assay_name,
                                   alt = TRUE, group = group, neg_lb = neg_lb)
        zero_ind = zero_ind[tax_keep, ]
        rownames(zero_ind) = NULL
    }else{ zero_ind = NULL }

    # 3. Estimation of parameters
    if (verbose) {
        message("Obtaining initial estimates ...")
    }
    para = .iter_mle(x = x, y = y, meta_data = meta_data,
                     formula = formula, theta = NULL, tol = tol,
                     max_iter = max_iter, verbose = FALSE)
    beta = para$beta
    vcov_hat = para$vcov_hat
    var_hat = para$var_hat

    # 4. Estimation of the sample-specific bias
    if (verbose) {
        message("Estimating sample-specific biases ...")
    }
    fun_list = list(.bias_em)
    bias = foreach(i = seq_len(ncol(beta)), .combine = rbind) %dorng% {
        output = fun_list[[1]](beta = beta[, i],
                               var_hat = var_hat[, i],
                               tol = tol,
                               max_iter = max_iter)
    }
    bias = data.frame(bias, row.names = covariates, check.names = FALSE)
    delta_em = bias$delta_em
    delta_wls = bias$delta_wls
    var_delta = bias$var_delta

    # 5. Obtain coefficients, standard errors, and sampling fractions
    beta_hat = beta
    beta_hat = t(t(beta_hat) - delta_em)

    # Account for the variance of delta_hat
    if (conserve) {
        var_hat = sweep(var_hat, 2, var_delta, "+") +
            2 * sqrt(sweep(var_hat, 2, var_delta, "*"))

        vcov_hat = lapply(seq_len(n_tax), function(i) {
            diag(vcov_hat[[i]]) = var_hat[i, ]
            return(vcov_hat[[i]])
        })
        se_hat = sqrt(var_hat)
    }else{ se_hat = sqrt(var_hat) }

    theta_hat = matrix(NA, nrow = n_tax, ncol = ncol(y))
    for (i in seq_len(n_tax)) {
        theta_hat[i, ] = y[i, ] - x %*% beta_hat[i, ]
    }
    theta_hat = colMeans(theta_hat, na.rm = TRUE)

    # 6. Primary results
    if (verbose) {
        message("ANCOM-BC primary results ...")
    }
    W = beta_hat/se_hat
    p = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)
    q = apply(p, 2, function(x) p.adjust(x, method = p_adj_method))
    diff_abn = q <= alpha & !is.na(q)

    beta_prim = cbind(taxon = data.frame(taxon = tax_name),
                      data.frame(beta_hat, check.names = FALSE,
                                 row.names = NULL))
    se_prim = cbind(taxon = data.frame(taxon = tax_name),
                    data.frame(se_hat, check.names = FALSE,
                               row.names = NULL))
    W_prim = cbind(taxon = data.frame(taxon = tax_name),
                   data.frame(W, check.names = FALSE,
                              row.names = NULL))
    p_prim = cbind(taxon = data.frame(taxon = tax_name),
                   data.frame(p, check.names = FALSE,
                              row.names = NULL))
    q_prim = cbind(taxon = data.frame(taxon = tax_name),
                   data.frame(q, check.names = FALSE,
                              row.names = NULL))
    diff_prim = cbind(taxon = data.frame(taxon = tax_name),
                      data.frame(diff_abn, check.names = FALSE,
                                 row.names = NULL))

    res = list(lfc = beta_prim,
               se = se_prim,
               W = W_prim,
               p_val = p_prim,
               q_val = q_prim,
               diff_abn = diff_prim)

    # 7. Global test results
    if (global) {
        if (verbose) {
            message("ANCOM-BC global test ...")
        }
        res_global = .ancombc_global_F(x = x, group = group,
                                       beta_hat = beta_hat,
                                       vcov_hat = vcov_hat,
                                       p_adj_method = p_adj_method,
                                       alpha = alpha)
    } else { res_global = NULL }

    # 8. Combine the information of structural zeros
    # Set p/q-values and SEs of structural zeros to be 0s.
    if (struc_zero) {
        if (verbose) {
            txt = paste0("Merge the information of structural zeros ... \n",
                         "Note that taxa with structural zeros will have ",
                         "0 p/q-values and SEs")
            message(txt)
        }
        zero_idx = as.matrix(zero_ind[, -1])
        group_ind = grepl(group, c("taxon", covariates))
        zero_mask = 1 - abs((zero_idx - zero_idx[, 1]))
        zero_mask = zero_mask[, -1, drop = FALSE]
        res$se[, group_ind] = res$se[, group_ind] * zero_mask
        res$p_val[, group_ind] = res$p_val[, group_ind] * zero_mask
        res$q_val[, group_ind] = res$q_val[, group_ind] * zero_mask
        res$diff_abn[, group_ind] = data.frame(res$q_val[, group_ind] <= alpha &
                                                   !is.na(res$q_val[, group_ind]),
                                  check.names = FALSE)

        # Global test
        if (global) {
            zero_mask = 1 - apply(zero_idx, 1, function(x)
                sum(x) > 0 & sum(x) < ncol(zero_idx))
            res_global[, "p_val"] = res_global[, "p_val"] * zero_mask
            res_global[, "q_val"] = res_global[, "q_val"] * zero_mask
            res_global[, "diff_abn"] = res_global[, "q_val"] <= alpha &
                !is.na(res_global[, "q_val"])
        }
    }

    # 9. Outputs
    out = list(feature_table = feature_table, zero_ind = zero_ind,
               samp_frac = theta_hat, delta_em = delta_em,
               delta_wls = delta_wls, res = res, res_global = res_global)

    if (n_cl > 1) {
      parallel::stopCluster(cl)
    }

    return(out)
}



