#' @title Analysis of Compositions of Microbiomes with Bias Correction 2
#' (ANCOM-BC2)
#'
#' @description Determine taxa whose absolute abundances, per unit volume, of
#' the ecosystem (e.g., gut) are significantly different with changes in the
#' covariate of interest (e.g., group). The current version of
#' \code{ancombc2} function implements Analysis of Compositions of Microbiomes
#' with Bias Correction (ANCOM-BC2) in cross-sectional and repeated measurements
#' data. In addition to the two-group comparison, ANCOM-BC2 also supports
#' testing for continuous covariates and multi-group comparisons,
#' including the global test, pairwise directional test, Dunnett's type of
#' test, and trend test.
#'
#' @details A taxon is considered to have structural zeros in some (>=1)
#' groups if it is completely (or nearly completely) missing in these groups.
#' For instance, suppose there are three groups: g1, g2, and g3.
#' If the counts of taxon A in g1 are 0 but nonzero in g2 and g3,
#' then taxon A will be considered to contain structural zeros in g1.
#' In this example, taxon A is declared to be differentially abundant between
#' g1 and g2, g1 and g3, and consequently, it is globally differentially
#' abundant with respect to this group variable.
#' Such taxa are not further analyzed using ANCOM-BC2, but the results are
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
#' Like other differential abundance analysis methods, ANCOM-BC2 log transforms
#' the observed counts. However, to deal with zero counts, a pseudo-count is
#' added before the log transformation. Several studies have shown that
#' differential abundance results could be sensitive to the choice of
#' pseudo-count
#' (\href{https://doi.org/10.1038/nmeth.2897}{Costea et al. (2014)};
#' \href{https://doi.org/10.1038/nmeth.2898}{Paulson, Bravo, and Pop (2014)}),
#' resulting in an inflated false positive rate. To avoid such false positives,
#' we conduct a sensitivity analysis and provide a sensitivity score for
#' each taxon to determine if a particular taxon is sensitive to the choice of
#' pseudo-count. The larger the score, the more likely the significant
#' result is a false positive.
#'
#' When performning pairwise directional (or Dunnett's type of) test, the mixed
#' directional false discover rate (mdFDR) should be taken into account.
#' The mdFDR is the combination of false discovery rate due to multiple testing,
#' multiple pairwise comparisons, and directional tests within each pairwise
#' comparison. For example, suppose we have five taxa and three experimental
#' groups: g1, g2, and g3. Thus, we are performing five tests corresponding to
#' five taxa. For each taxon, we are also conducting three pairwise comparisons
#' (g1 vs. g2, g2 vs. g3, and g1 vs. g3). Within each pairwise comparison,
#' we wish to determine if the abundance has increased or decreased or did not
#' change (direction of the effect size). Errors could occur in each step.
#' The overall false discovery rate is controlled by the mdFDR methodology we
#' adopted from
#' \href{https://doi.org/10.1111/j.1541-0420.2009.01292.x}{Guo, Sarkar, and Peddada (2010)} and
#' \href{https://doi.org/10.1186/s12859-016-0937-5}{Grandhi, Guo, and Peddada (2016)}.
#'
#' @param data the input data. A
#' \code{phyloseq}, \code{SummarizedExperiment}, or
#' \code{TreeSummarizedExperiment} object, which consists of
#' a feature table (microbial count table), a sample metadata, a
#' taxonomy table (optional), and a phylogenetic tree (optional). The row names
#' of the metadata must match the sample names of the feature table, and the
#' row names of the taxonomy table must match the taxon (feature) names of the
#' feature table. See \code{?phyloseq::phyloseq},
#' \code{?SummarizedExperiment::SummarizedExperiment}, or
#' \code{?TreeSummarizedExperiment::TreeSummarizedExperiment} for more details.
#' It is highly recommended that the input data
#' are in low taxonomic levels, such as OTU or species level, as the estimation
#' of sampling fractions requires a large number of taxa.
#' @param assay_name character. Name of the count table in the data object
#' (only applicable if data object is a \code{(Tree)SummarizedExperiment}).
#' Default is "counts".
#' See \code{?SummarizedExperiment::assay} for more details.
#' @param tax_level character. The taxonomic level of interest. The input data
#' can be agglomerated at different taxonomic levels based on your research
#' interest. Default is NULL, i.e., do not perform agglomeration, and the
#' ANCOM-BC2 anlysis will be performed at the lowest taxonomic level of the
#' input \code{data}.
#' @param fix_formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the fixed effects in metadata.
#' @param rand_formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the random effects in metadata. ANCOM-BC2
#' follows the \code{lmerTest} package in formulating the random effects. See
#' \code{?lmerTest::lmer} for more details. Default is \code{NULL}.
#' @param p_adj_method character. method to adjust p-values. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{?stats::p.adjust} for more details.
#' @param pseudo numeric. Add pseudo-counts to the data.
#' Default is 0 (no pseudo-count addition).
#' @param pseudo_sens logical. Whether to perform the sensitivity analysis to
#' the pseudo-count addition. Default is \code{TRUE}. See \code{Details} for
#' a more comprehensive discussion on this sensitivity analysis.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' less than \code{prv_cut} will be excluded in the analysis. For instance,
#' suppose there are 100 samples, if a taxon has nonzero counts presented in
#' less than 10 samples, it will not be further analyzed. Default is 0.10.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis. Default is 0, i.e. do not discard any sample.
#' @param s0_perc a numerical fraction between 0 and 1. Inspired by
#' \href{https://doi.org/10.1073/pnas.091062498}{Significance
#' Analysis of Microarrays (SAM)} methodology, a small positive constant is
#' added to the denominator of ANCOM-BC2 test statistic corresponding to
#' each taxon to avoid the significance due to extremely small standard errors,
#' especially for rare taxa. This small positive constant is chosen as
#' \code{s0_perc}-th percentile of standard error values for each fixed effect.
#' Default is 0.05 (5th percentile).
#' @param group character. The name of the group variable in metadata.
#' \code{group} should be discrete. Specifying \code{group} is required for
#' detecting structural zeros and performing multi-group comparisons (global
#' test, pairwise directional test, Dunnett's type of test, and trend test).
#' Default is NULL. If the \code{group} of interest contains only two
#' categories, leave it as NULL.
#' @param struc_zero logical. Whether to detect structural zeros based on
#' \code{group}. Default is FALSE. See \code{Details} for
#' a more comprehensive discussion on structural zeros.
#' @param neg_lb logical. Whether to classify a taxon as a structural zero using
#' its asymptotic lower bound. Default is FALSE.
#' @param alpha numeric. Level of significance. Default is 0.05.
#' @param n_cl numeric. The number of nodes to be forked. For details, see
#' \code{?parallel::makeCluster}. Default is 1 (no parallel computing).
#' @param verbose logical. Whether to generate verbose output during the
#' ANCOM-BC2 fitting process. Default is FALSE.
#' @param global logical. Whether to perform the global test. Default is FALSE.
#' @param pairwise logical. Whether to perform the pairwise directional test.
#' Default is FALSE.
#' @param dunnet logical. Whether to perform the Dunnett's type of test.
#' Default is FALSE.
#' @param trend logical. Whether to perform trend test. Default is FALSE.
#' @param iter_control a named list of control parameters for the iterative
#' MLE or RMEL algorithm, including 1) \code{tol}: the iteration convergence
#' tolerance (default is 1e-02), 2) \code{max_iter}: the maximum number of
#' iterations (default is 20), and 3)\code{verbose}: whether to show the verbose
#' output (default is FALSE).
#' @param em_control a named list of control parameters for the E-M algorithm,
#' including 1) \code{tol}: the iteration convergence tolerance
#' (default is 1e-05) and 2) \code{max_iter}: the maximum number of iterations
#' (default is 100).
#' @param lme_control a list of control parameters for mixed model fitting.
#' See \code{?lme4::lmerControl} for details.
#' @param mdfdr_control a named list of control parameters for mixed directional
#' false discover rate (mdFDR), including 1) \code{fwer_ctrl_method}: family
#' wise error (FWER) controlling procedure, such as "holm", "hochberg",
#' "bonferroni", etc (default is "holm") and 2) \code{B}: the number of
#' bootstrap samples (default is 100). Increase \code{B} will lead to a more
#' accurate p-values. See \code{Details} for a more comprehensive discussion on
#' mdFDR.
#' @param trend_control a named list of control parameters for the trend test,
#' including 1) \code{contrast}: the list of contrast matrices for
#' constructing inequalities, 2) \code{node}: the list of positions for the
#' nodal parameter, 3) \code{solver}: a string indicating the solver to use
#' (default is "ECOS"), and 4) \code{B}: the number of bootstrap samples
#' (default is 100). Increase \code{B} will lead to a more accurate p-values.
#' See \code{vignette} for the corresponding trend test examples.
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
#'         \item{ \code{pseudo_sens_tab}, the results of sensitivity analysis
#'         for the pseudo-count addition. }
#'         \item{ \code{res},  a \code{data.frame} containing ANCOM-BC2 primary
#'         result:}
#'         \itemize{
#'         \item{ columns started with \code{lfc}: log fold changes
#'         obtained from the ANCOM-BC2 log-linear (natural log) model.}
#'         \item{ columns started with \code{se}: standard errors (SEs) of
#'         \code{lfc}.}
#'         \item{ columns started with \code{W}: test statistics.
#'         \code{W = lfc/se}.}
#'         \item{ columns started with \code{p}: p-values. P-values are
#'         obtained from two-sided Z-test using the test statistic \code{W}.}
#'         \item{ columns started with \code{q}: adjusted p-values.
#'         Adjusted p-values are obtained by applying \code{p_adj_method}
#'         to \code{p}.}
#'         \item{ columns started with \code{diff}: TRUE if the
#'         taxon is significant (has \code{q} less than \code{alpha}).}
#'         }
#'         \item{ \code{res_global},  a \code{data.frame} containing ANCOM-BC2
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
#'         \item{ \code{res_pair},  a \code{data.frame} containing ANCOM-BC2
#'         pairwise directional test result for the variable specified in
#'         \code{group}:}
#'         \itemize{
#'         \item{ columns started with \code{lfc}: log fold changes.}
#'         \item{ columns started with \code{se}: standard errors (SEs).}
#'         \item{ columns started with \code{W}: test statistics.}
#'         \item{ columns started with \code{p}: p-values.}
#'         \item{ columns started with \code{q}: adjusted p-values.}
#'         \item{ columns started with \code{diff}: TRUE if the
#'         taxon is significant (has \code{q} less than \code{alpha}).}
#'         }
#'         \item{ \code{res_dunn},  a \code{data.frame} containing ANCOM-BC2
#'         Dunnett's type of test result for the variable specified in
#'         \code{group}:}
#'         \itemize{
#'         \item{ columns started with \code{lfc}: log fold changes.}
#'         \item{ columns started with \code{se}: standard errors (SEs).}
#'         \item{ columns started with \code{W}: test statistics.}
#'         \item{ columns started with \code{p}: p-values.}
#'         \item{ columns started with \code{q}: adjusted p-values.}
#'         \item{ columns started with \code{diff}: TRUE if the
#'         taxon is significant (has \code{q} less than \code{alpha}).}
#'         }
#'         \item{ \code{res_trend},  a \code{data.frame} containing ANCOM-BC2
#'         trend test result for the variable specified in
#'         \code{group}:}
#'         \itemize{
#'         \item{ columns started with \code{lfc}: log fold changes.}
#'         \item{ columns started with \code{se}: standard errors (SEs).}
#'         \item{ \code{W}: test statistics.}
#'         \item{ \code{p_val}: p-values.}
#'         \item{ \code{q_val}: adjusted p-values.}
#'         \item{ \code{diff_abn}: TRUE if the
#'         taxon is significant (has \code{q} less than \code{alpha}).}
#'         }
#'         }
#'
#' @seealso \code{\link{ancom}} \code{\link{ancombc}}
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
#' #=======================Run ANCOMBC2 Using a Real Data=======================
#' library(ANCOMBC)
#' data(dietswap)
#'
#' colData(dietswap)$bmi_group = factor(colData(dietswap)$bmi_group,
#'                                      levels = c("obese",
#'                                                 "overweight",
#'                                                 "lean"))
#'
#' set.seed(123)
#' # Note that setting pseudo_sens = FALSE, max_iter = 1, and B = 1 is
#' # only for the sake of speed
#' # Set pseudo_sens = TRUE, and use default or larger values for max_iter and B
#' # for better performance
#' out = ancombc2(data = dietswap, assay_name = "counts", tax_level = "Phylum",
#'                fix_formula = "nationality + timepoint + bmi_group",
#'                rand_formula = "(timepoint | subject)",
#'                p_adj_method = "holm", pseudo = 0, pseudo_sens = FALSE,
#'                prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#'                group = "bmi_group", struc_zero = TRUE, neg_lb = TRUE,
#'                alpha = 0.05, n_cl = 1, verbose = TRUE,
#'                global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
#'                iter_control = list(tol = 1e-2, max_iter = 1, verbose = TRUE),
#'                em_control = list(tol = 1e-5, max_iter = 1),
#'                lme_control = lme4::lmerControl(),
#'                mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
#'                trend_control = list(contrast =
#'                                           list(matrix(c(1, 0, -1, 1),
#'                                                       nrow = 2,
#'                                                       byrow = TRUE)),
#'                                       node = list(2),
#'                                       solver = "ECOS",
#'                                       B = 1))
#' res_prim = out$res
#' res_global = out$res_global
#' res_pair = out$res_pair
#' res_dunn = out$res_dunn
#' res_trend = out$res_trend
#'
#' @author Huang Lin
#'
#' @references
#' \insertRef{kaul2017analysis}{ANCOMBC}
#'
#' \insertRef{lin2020analysis}{ANCOMBC}
#'
#' \insertRef{tusher2001significance}{ANCOMBC}
#'
#' \insertRef{costea2014fair}{ANCOMBC}
#'
#' \insertRef{paulson2014reply}{ANCOMBC}
#'
#' \insertRef{guo2010controlling}{ANCOMBC}
#'
#' \insertRef{grandhi2016multiple}{ANCOMBC}
#'
#' @rawNamespace import(stats, except = filter)
#' @importFrom utils combn
#' @importFrom mia makeTreeSummarizedExperimentFromPhyloseq taxonomyRanks agglomerateByRank
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
#' @importFrom emmeans emmeans contrast
#' @importFrom CVXR Variable Minimize Problem solve matrix_frac
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% %:%
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng%
#' @importFrom rngtools setRNG
#' @importFrom MASS ginv
#' @importFrom nloptr neldermead
#' @importFrom Rdpack reprompt
#'
#' @export
ancombc2 = function(data, assay_name = "counts", tax_level = NULL,
                    fix_formula , rand_formula = NULL,
                    p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                    prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                    group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                    alpha = 0.05, n_cl = 1, verbose = FALSE,
                    global = FALSE, pairwise = FALSE,
                    dunnet = FALSE, trend = FALSE,
                    iter_control = list(tol = 1e-2,
                                        max_iter = 20,
                                        verbose = FALSE),
                    em_control = list(tol = 1e-5, max_iter = 100),
                    lme_control = lme4::lmerControl(),
                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                    trend_control = list(contrast = NULL,
                                         node = NULL,
                                         solver = "ECOS",
                                         B = 100)){
    cl = parallel::makeCluster(n_cl)
    doParallel::registerDoParallel(cl)

    # 1. Data pre-processing
    # TSE data construction
    tse_obj = tse_construct(data = data, assay_name = assay_name,
                            tax_level = tax_level, phyloseq = NULL)
    tse = tse_obj$tse
    assay_name = tse_obj$assay_name
    tax_level = tse_obj$tax_level
    tse_alt = tse_obj$tse_alt

    # Identify taxa with structural zeros
    if (struc_zero) {
        zero_ind = get_struc_zero(tse = tse, tax_level = tax_level,
                                  assay_name = assay_name,
                                  alt = TRUE, group = group, neg_lb = neg_lb)
        # Taxa with structural zeros will be removed from ANCOM-BC2 analyses
        tax_idx = apply(zero_ind[, -1], 1, function(x) all(x == FALSE))
        tax_keep = which(tax_idx)
    }else{
        zero_ind = NULL
        # Taxa with structural zeros will be removed from ANCOM-BC2 analyses
        tax_keep = seq(nrow(tse_alt))}

    # Filter data by prevalence and library size
    core1 = data_core(tse = tse, tax_level = tax_level,
                      assay_name = assay_name, alt = FALSE,
                      prv_cut = prv_cut, lib_cut = lib_cut,
                      tax_keep = NULL, samp_keep = NULL)
    O1 = core1$feature_table
    samp_keep = colnames(O1)

    core2 = data_core(tse = tse, tax_level = tax_level,
                      assay_name = assay_name, alt = TRUE,
                      prv_cut = prv_cut, lib_cut = lib_cut,
                      tax_keep = tax_keep, samp_keep = samp_keep)
    O2 = core2$feature_table
    n_tax = nrow(O2)
    tax_name = rownames(O2)

    # Add pseudo-counts
    O1 = O1 + pseudo
    O2 = O2 + pseudo

    # Metadata and arguments check
    qc = data_qc(meta_data = core2$meta_data, group = group,
                 struc_zero = struc_zero, global = global,
                 pairwise = pairwise, dunnet = dunnet,
                 mdfdr_control = mdfdr_control,
                 trend = trend, trend_control = trend_control)
    meta_data = qc$meta_data
    global = qc$global
    pairwise = qc$pairwise
    dunnet = qc$dunnet
    trend = qc$trend
    trend_control = qc$trend_control

    # 2. Estimation of the sample-specific biases
    options(na.action = "na.pass") # Keep NA's in rows of x
    x = model.matrix(formula(paste0("~", fix_formula)), data = meta_data)
    options(na.action = "na.omit") # Switch it back
    fix_eff = colnames(x)
    n_fix_eff = length(fix_eff)

    if (nrow(O1) < 50) {
        warn_txt = sprintf(paste0("The number of taxa used for estimating ",
                                  "sample-specific biases is: ",
                                  nrow(O1),
                                  "\nA large number of taxa (>50) is required ",
                                  "for the consistent estimation of biases"))
        warning(warn_txt, call. = FALSE)
    }
    o1 = log(O1)
    o1[is.infinite(o1)] = 0
    y1 = o1 - rowMeans(o1, na.rm = TRUE)

    # Obtain initial estimates
    if (verbose) {
        message("Obtaining initial estimates ...")
    }

    if (is.null(rand_formula)) {
        para1 = iter_mle(x = x, y = y1, meta_data = meta_data,
                         formula = fix_formula,
                         theta = NULL, tol = iter_control$tol,
                         max_iter = iter_control$max_iter,
                         verbose = iter_control$verbose)
    } else {
        para1 = iter_remle(x = x, y = y1, meta_data = meta_data,
                           fix_formula = fix_formula,
                           rand_formula = rand_formula,
                           lme_control = lme_control,
                           theta = NULL, tol = iter_control$tol,
                           max_iter = iter_control$max_iter,
                           verbose = iter_control$verbose)
    }
    beta1 = para1$beta
    var_hat1 = para1$var_hat

    # Apply E-M algorithm
    if (verbose) {
        message("Estimating sample-specific biases ...")
    }
    fun_list = list(bias_em)
    bias1 = foreach(i = seq_len(ncol(beta1)), .combine = rbind) %dorng% {
        output = fun_list[[1]](beta = beta1[, i],
                               var_hat = var_hat1[, i],
                               tol = em_control$tol,
                               max_iter = em_control$max_iter)
    }
    bias1 = data.frame(bias1, row.names = fix_eff, check.names = FALSE)
    delta_em = bias1$delta_em
    delta_wls = bias1$delta_wls
    var_delta = bias1$var_delta

    # Obtain the final estimates for sample-specific biases
    beta1 = t(t(beta1) - delta_em)
    theta_hat = matrix(NA, nrow = nrow(y1), ncol = ncol(y1))
    x_impute = x
    x_impute[is.na(x_impute)] = 0
    for (i in seq_len(nrow(y1))) {
        theta_hat[i, ] = y1[i, ] - x_impute %*% beta1[i, ]
    }
    theta_hat = colMeans(theta_hat, na.rm = TRUE)

    # 3. Obtain unbiased estimates
    o2 = log(O2)
    o2[is.infinite(o2)] = 0
    y2 = o2 - rowMeans(o2, na.rm = TRUE)
    if (is.null(rand_formula)) {
        para2 = iter_mle(x = x, y = y2, meta_data = meta_data,
                         formula = fix_formula,
                         theta = theta_hat, tol = iter_control$tol,
                         max_iter = iter_control$max_iter,
                         verbose = iter_control$verbose)
    } else {
        para2 = iter_remle(x = x, y = y2, meta_data = meta_data,
                           fix_formula = fix_formula,
                           rand_formula = rand_formula,
                           lme_control = lme_control,
                           theta = theta_hat, tol = iter_control$tol,
                           max_iter = iter_control$max_iter,
                           verbose = iter_control$verbose)
    }
    beta_hat = para2$beta
    var_hat = para2$var_hat

    # Account for the variance of delta
    var_hat = sweep(var_hat, 2, var_delta, "+") +
        2 * sqrt(sweep(var_hat, 2, var_delta, "*"))

    # Add a small positive constant to stabilize the variance
    s02 = apply(var_hat, 2, function(x)
        stats::quantile(x, s0_perc, na.rm = TRUE))
    var_hat = t(t(var_hat) + s02)
    se_hat = sqrt(var_hat)
    vcov_hat = lapply(seq_len(n_tax), function(i) {
        diag(para2$vcov_hat[[i]]) = var_hat[i, ]
        return(para2$vcov_hat[[i]])
    })

    # 4. Primary results
    if (verbose) {
        message("ANCOM-BC2 primary results ...")
    }
    W = beta_hat/se_hat
    p_hat = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)
    q_hat = apply(p_hat, 2, function(x) p.adjust(x, method = p_adj_method))
    diff_abn = q_hat < alpha & !is.na(q_hat)

    beta_prim = data.frame(beta_hat, check.names = FALSE)
    se_prim = data.frame(se_hat, check.names = FALSE)
    W_prim = data.frame(W, check.names = FALSE)
    p_prim = data.frame(p_hat, check.names = FALSE)
    q_prim = data.frame(q_hat, check.names = FALSE)
    diff_prim = data.frame(diff_abn, check.names = FALSE)
    colnames(beta_prim) = paste0("lfc_", colnames(beta_hat))
    colnames(se_prim) = paste0("se_", colnames(se_hat))
    colnames(W_prim) = paste0("W_", colnames(W))
    colnames(p_prim) = paste0("p_", colnames(p_hat))
    colnames(q_prim) = paste0("q_", colnames(q_hat))
    colnames(diff_prim) = paste0("diff_", colnames(diff_abn))
    res = do.call("cbind", list(data.frame(taxon = tax_name),
                                beta_prim, se_prim, W_prim,
                                p_prim, q_prim, diff_prim))
    rownames(res) = NULL

    # Sensitivity analysis for the pseudo-count addition
    if (pseudo_sens) {
        if (verbose) {
            message("The sensitivity analysis for the pseudo-count addition ...")
        }
        pseudo_list = seq(0, 1, 0.01)
        tformula = paste0("~ ", fix_formula)
        pseudo_sens_tab = foreach(i = seq_len(n_tax), .combine = rbind) %dorng% {
            count1 = core2$feature_table[i, ]
            nlp_tab = vector(mode = "numeric")
            for (j in seq_along(pseudo_list)) {
                pseudo = pseudo_list[j]
                count2 = count1 + pseudo
                log_count = log(count2)
                log_count[is.infinite(log_count)] = 0
                log_resid = log_count - mean(log_count)
                log_resid_crt = log_resid - theta_hat
                df = data.frame(y = log_resid_crt, meta_data)
                lm_fit = stats::lm(formula(paste("y", tformula)), data = df)
                # Negative log p-value
                nlp = -log(summary(lm_fit)$coef[, "Pr(>|t|)"])
                if (global) {
                    anova_fit = anova(lm_fit)
                    nlp_global = -log(anova_fit$`Pr(>F)`[grep(group, rownames(anova_fit))])
                    nlp = c(nlp, global = nlp_global)
                }
                if (pairwise) {
                    emm_fit = emmeans::emmeans(lm_fit, specs = group)
                    emm_pair = emmeans::contrast(emm_fit, "pairwise",
                                                 adjust = "none")
                    nlp_pair = -log(summary(emm_pair)[, "p.value"])
                    names(nlp_pair) = summary(emm_pair)[, "contrast"]
                    nlp = c(nlp, nlp_pair)
                }
                nlp_tab = cbind(nlp_tab, nlp)
            }
            apply(nlp_tab, 1, function(x) {
                ifelse(all(x > -log(alpha)), 0, sd(x))
            })
        }
        pseudo_sens_tab = as.data.frame(pseudo_sens_tab)
        pseudo_sens_tab = cbind(data.frame(taxon = tax_name),
                                pseudo_sens_tab)
        rownames(pseudo_sens_tab) = NULL
    } else {
        pseudo_sens_tab = NULL
    }

    # 5. Results of global test
    if (global) {
        if (verbose) {
            message("ANCOM-BC2 global test ...")
        }
        res_global = ancombc_global(x = x, group = group,
                                    beta_hat = beta_hat,
                                    vcov_hat = vcov_hat,
                                    p_adj_method = p_adj_method,
                                    alpha = alpha)
    } else { res_global = NULL }

    # 6. Results of pairwise test
    if (pairwise) {
        if (verbose) {
            message("ANCOM-BC2 pairwise directional test ...")
        }
        res_pair = ancombc_pair(x = x, group = group,
                                beta_hat = beta_hat,
                                var_hat = var_hat, vcov_hat = vcov_hat,
                                fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                                alpha = alpha)
        beta_pair = data.frame(res_pair$beta, check.names = FALSE)
        se_pair = data.frame(res_pair$se, check.names = FALSE)
        W_pair = data.frame(res_pair$W, check.names = FALSE)
        p_pair = data.frame(res_pair$p_val, check.names = FALSE)
        q_pair = data.frame(res_pair$q_val, check.names = FALSE)
        diff_pair = data.frame(res_pair$diff_abn, check.names = FALSE)

        # Directional test summary
        colnames(beta_pair) = paste0("lfc_", colnames(beta_pair))
        colnames(se_pair) = paste0("se_", colnames(se_pair))
        colnames(W_pair) = paste0("W_", colnames(W_pair))
        colnames(p_pair) = paste0("p_", colnames(p_pair))
        colnames(q_pair) = paste0("q_", colnames(q_pair))
        colnames(diff_pair) = paste0("diff_", colnames(diff_pair))
        res_pair = do.call("cbind", list(data.frame(taxon = tax_name),
                                         beta_pair, se_pair, W_pair,
                                         p_pair, q_pair, diff_pair))
        rownames(res_pair) = NULL
    } else {
        res_pair = NULL
    }

    # 7. Results of Dunnet's type of test
    if (dunnet) {
        if (verbose) {
            message("ANCOM-BC2 Dunnet's type of test ...")
        }
        res_dunn = ancombc_dunn(x = x, group = group,
                                beta_hat = beta_hat, var_hat = var_hat,
                                fwer_ctrl_method = mdfdr_control$fwer_ctrl_method,
                                B = mdfdr_control$B, alpha = alpha)
        beta_dunn = data.frame(res_dunn$beta, check.names = FALSE)
        se_dunn = data.frame(res_dunn$se, check.names = FALSE)
        W_dunn = data.frame(res_dunn$W, check.names = FALSE)
        p_dunn = data.frame(res_dunn$p_val, check.names = FALSE)
        q_dunn = data.frame(res_dunn$q_val, check.names = FALSE)
        diff_dunn = data.frame(res_dunn$diff_abn, check.names = FALSE)

        # Directional test summary
        colnames(beta_dunn) = paste0("lfc_", colnames(beta_dunn))
        colnames(se_dunn) = paste0("se_", colnames(se_dunn))
        colnames(W_dunn) = paste0("W_", colnames(W_dunn))
        colnames(p_dunn) = paste0("p_", colnames(p_dunn))
        colnames(q_dunn) = paste0("q_", colnames(q_dunn))
        colnames(diff_dunn) = paste0("diff_", colnames(diff_dunn))
        res_dunn = do.call("cbind", list(data.frame(taxon = tax_name),
                                         beta_dunn, se_dunn, W_dunn,
                                         p_dunn, q_dunn, diff_dunn))
        rownames(res_dunn) = NULL
    } else {
        res_dunn = NULL
    }

    # 8. Results of trend test
    if (trend) {
        if (verbose) {
            message("ANCOM-BC2 trend test ...")
        }
        res_trend = ancombc_trend(
            x = x, group = group, beta_hat = beta_hat,
            var_hat = var_hat, vcov_hat = vcov_hat,
            p_adj_method = p_adj_method, alpha = alpha,
            trend_control = trend_control)
        beta_trend = res_trend$beta
        se_trend = res_trend$se
        W_trend = res_trend$W
        p_trend = res_trend$p_val
        q_trend = res_trend$q_val
        diff_trend = res_trend$diff_abn

        # Directional test summary
        colnames(beta_trend) = paste0("lfc_", colnames(beta_trend))
        colnames(se_trend) = paste0("se_", colnames(se_trend))
        res_trend = cbind(data.frame(taxon = tax_name),
                          beta_trend, se_trend,
                          data.frame(W = W_trend, p_val = p_trend,
                                     q_val = q_trend, diff_abn = diff_trend))
        rownames(res_trend) = NULL
    } else {
        res_trend = NULL
    }

  # 9. Outputs
  out = list(feature_table = O2,
             zero_ind = zero_ind,
             samp_frac = theta_hat,
             delta_em = delta_em,
             delta_wls = delta_wls,
             pseudo_sens_tab = pseudo_sens_tab,
             res = res,
             res_global = res_global,
             res_pair = res_pair,
             res_dunn = res_dunn,
             res_trend = res_trend)

  parallel::stopCluster(cl)
  return(out)
}



