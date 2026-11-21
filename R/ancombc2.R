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
#' @details Taxa with structural zeros are excluded from the primary result and
#' are reported in \code{zero_ind}. The detection of structural zeros is
#' controlled by \code{struc_zero} and \code{neg_lb}, which are defined in
#' \code{\link{ANCOMBC-concepts}}.
#'
#'Like other differential abundance analysis methods, ANCOM-BC2 applies a log
#'transformation to the observed counts. However, the presence of zero counts
#'poses a challenge, and researchers often consider adding a pseudo-count before
#'the log transformation. However, it has been shown that the choice of
#'pseudo-count can impact the results and lead to an inflated false positive
#'rate (\href{https://doi.org/10.1038/nmeth.2897}{Costea et al. (2014)};
#' \href{https://doi.org/10.1038/nmeth.2898}{Paulson, Bravo, and Pop (2014)}).
#' To address this issue, we conduct a sensitivity analysis to assess the impact
#' of different pseudo-counts on zero counts for each taxon. The
#' \code{conservative} argument specifies the pseudo-counts that are added and
#' the quantity on which the agreement between runs is evaluated. Let
#' \eqn{\alpha} denote the significance level, specified by \code{alpha}, and
#' let \eqn{I(\cdot)}{I(.)} denote the indicator function.
#'
#' Under \code{conservative = TRUE}, the pseudo-counts
#' \eqn{c \in \{0, 0.1, 0.5, 1\}}{c in {0, 0.1, 0.5, 1}} are added to all counts before the estimation
#' of sampling fractions, and the full ANCOM-BC2 algorithm is re-run for each
#' \eqn{c}. Let \eqn{q_c} denote the adjusted p-value of a taxon in run \eqn{c}.
#' The sensitivity score of the taxon is
#' \deqn{s = \frac{1}{4} \sum_{c \in \{0, 0.1, 0.5, 1\}} I(q_c > \alpha).}{s = (1/4) * sum_c I(q_c > alpha).}
#' The taxon passes the sensitivity analysis, \code{passed_ss = TRUE}, when
#' \deqn{s = 0 \quad \mbox{or} \quad s = 1.}{s = 0 or s = 1.}
#'
#' Under \code{conservative = FALSE}, the sampling fractions are estimated once
#' on the complete data, the pseudo-counts
#' \eqn{c \in \{0.01, 0.02, \ldots, 0.5\}}{c in {0.01, 0.02, ..., 0.5}} are added to the zero counts of the
#' bias-corrected data, and the fixed effects are re-estimated by ordinary least
#' squares for each \eqn{c}. Let \eqn{p_c} denote the unadjusted p-value of a
#' taxon in refit \eqn{c}, and let \eqn{p_0} denote the unadjusted p-value of
#' the taxon in the main run. The sensitivity score of the taxon is
#' \deqn{s = \frac{1}{50} \sum_{c \in \{0.01, 0.02, \ldots, 0.5\}} I(p_c > \alpha).}{s = (1/50) * sum_c I(p_c > alpha).}
#' The taxon passes the sensitivity analysis when
#' \deqn{(s = 0 \mbox{ and } p_0 \le \alpha) \quad \mbox{or} \quad (s = 1 \mbox{ and } p_0 > \alpha).}{(s = 0 and p_0 <= alpha) or (s = 1 and p_0 > alpha).}
#'
#' Under both settings, \code{diff_robust} is TRUE when \code{diff_abn} is TRUE
#' and \code{passed_ss} is TRUE. Use the default unless power is the primary
#' concern, since \code{conservative = FALSE} has higher power and a higher
#' false positive rate.
#'
#' The pairwise directional test and the Dunnett's type of test control the
#' mixed directional false discovery rate (mdFDR) through the two-stage
#' procedure specified by \code{mdfdr_control}. The mdFDR and the procedure are
#' defined in \code{\link{ANCOMBC-concepts}}.
#'
#' @param data the input data. The \code{data} parameter should be either a
#' \code{matrix}, \code{data.frame}, \code{phyloseq} or a \code{TreeSummarizedExperiment}
#' object. Both \code{phyloseq} and \code{TreeSummarizedExperiment} objects
#' consist of a feature table (microbial count table), a sample metadata table,
#' a taxonomy table (optional), and a phylogenetic tree (optional).
#' If a \code{matrix} or \code{data.frame} is provided, ensure that the row
#' names of the \code{metadata} match the sample names (column names if
#' \code{taxa_are_rows} is TRUE, and row names otherwise) in \code{data}.
#' if a \code{phyloseq} or a \code{TreeSummarizedExperiment} is used, this
#' standard has already been enforced. For detailed information, refer to
#' \code{?phyloseq::phyloseq} or
#' \code{?TreeSummarizedExperiment::TreeSummarizedExperiment}.
#' It is recommended to use low taxonomic levels, such as OTU or species level,
#' as the estimation of sampling fractions requires a large number of taxa.
#' @param taxa_are_rows logical. Whether taxa are positioned in the rows of the
#' feature table. Default is TRUE.
#' @param assay_name character. Name of the count table in the data object
#' (only applicable if data object is a \code{(Tree)SummarizedExperiment}).
#' Default is "counts".
#' See \code{?SummarizedExperiment::assay} for more details.
#' @param assay.type alias for \code{assay_name}.
#' @param tax_level character. The taxonomic or non taxonomic(rowData) level of interest. The input data
#' can be analyzed at any taxonomic or rowData level without prior agglomeration.
#' Note that \code{tax_level} must be a value from \code{taxonomyRanks} or \code{rowData}, which
#' includes "Kingdom", "Phylum" "Class", "Order", "Family" "Genus" "Species" etc.
#' See \code{?mia::taxonomyRanks} for more details.
#' Default is NULL, i.e., do not perform agglomeration, and the
#' ANCOM-BC2 analysis will be performed at the lowest taxonomic level of the
#' input \code{data}.
#' @param rank alias for \code{tax_level}.
#' @param aggregate_data The abundance data that has been aggregated to the desired
#' taxonomic level. This parameter is required only when the input data is in
#' \code{matrix} or \code{data.frame} format. For \code{phyloseq} or \code{TreeSummarizedExperiment}
#' data, aggregation is performed by specifying the \code{tax_level} parameter.
#' @param meta_data a \code{data.frame} containing sample metadata.
#' This parameter is mandatory when the input \code{data} is a generic
#' \code{matrix} or \code{data.frame}. Ensure that the row names of the \code{metadata} match the
#' sample names (column names if \code{taxa_are_rows} is TRUE, and row names
#' otherwise) in \code{data}.
#' @param fix_formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the fixed effects in metadata. When
#' specifying the \code{fix_formula}, make sure to include the \code{group}
#' variable in the formula if it is not NULL. This argument is required and has
#' no default.
#' @param rand_formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the random effects in metadata. ANCOM-BC2
#' follows the \code{lmerTest} package in formulating the random effects. See
#' \code{?lmerTest::lmer} for more details. Default is \code{NULL}.
#' @param p_adj_method character. method to adjust p-values. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{?stats::p.adjust} for more details.
#' @param pseudo numeric. A nonnegative value added to all counts before the log
#' transformation. Default is 0, in which case no pseudo-count is added and zero
#' counts are treated as missing on the log scale. A nonzero value replaces the
#' log of a zero count by the log of \code{pseudo}, so that every count
#' contributes to the estimation, at the cost of a result that depends on the
#' value chosen.
#' @param pseudo_sens logical. Whether to perform the sensitivity analysis to
#' the pseudo-count addition. Default is \code{TRUE}. While ANCOM-BC2 utilizes
#' complete data (nonzero counts) by default for its analysis, a comprehensive
#' evaluation of result robustness is performed by assessing how pseudo-count
#' addition to zeros may affect the outcomes. For a detailed discussion on this
#' sensitivity analysis, refer to the \code{Details} section.
#' @param conservative logical. Whether to add the pseudo-counts of the
#' sensitivity analysis before the estimation of sampling fractions. Default is
#' \code{TRUE}, which adds the pseudo-counts 0.1, 0.5, and 1 to the feature
#' table, re-runs the entire ANCOM-BC2 algorithm for each value, and evaluates
#' the agreement on the adjusted p-values. Setting \code{conservative = FALSE}
#' estimates the sampling fractions once on the complete data, adds the
#' pseudo-counts 0.01 to 0.5 in steps of 0.01 to the zero counts of the
#' bias-corrected data, re-estimates the fixed effects by ordinary least squares
#' for each value, and evaluates the agreement on the unadjusted p-values. Use
#' the default unless power is the primary concern, since
#' \code{conservative = FALSE} has higher power and a higher false positive
#' rate. This argument applies only when \code{pseudo_sens = TRUE}. The
#' \code{Details} section states the sensitivity score and the passing rule of
#' each setting.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' (the proportion of samples in which the taxon is present)
#' less than \code{prv_cut} will be excluded in the analysis. For example,
#' if there are 100 samples, and a taxon has nonzero counts present in less than
#' 100*prv_cut samples, it will not be considered in the analysis.
#' Default is 0.10.
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
#' @param group character. the name of the group variable in metadata.
#' The \code{group} parameter should be a character string representing the name
#' of the group variable in the metadata. The \code{group} variable should be
#' discrete, meaning it consists of categorical values. Specifying the
#' \code{group} variable is required if you are interested in detecting
#' structural zeros and performing multi-group comparisons (global
#' test, pairwise directional test, Dunnett's type of test, and trend test).
#' However, if these analyses are not of interest to you, you can leave the
#' \code{group} parameter as NULL. If the \code{group} variable of interest
#' contains only two categories, you can also leave the \code{group} parameter
#' as NULL. Default is NULL.
#' @param struc_zero logical. Whether to detect structural zeros based on
#' \code{group}. Default is FALSE. See \code{Details} for
#' a more comprehensive discussion on structural zeros.
#' @param neg_lb logical. Whether to classify a taxon as a structural zero using
#' its asymptotic lower bound. Default is FALSE.
#' @param alpha numeric. Level of significance. Default is 0.05.
#' @param n_cl numeric. The number of nodes to be forked. For details, see
#' \code{?parallel::makeCluster}. Default is 1 (no parallel computing).
#' @param verbose logical. Whether to generate verbose output during the
#' ANCOM-BC2 fitting process. Default is TRUE.
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
#' false discovery rate (mdFDR), including 1) \code{fwer_ctrl_method}: family
#' wise error (FWER) controlling procedure, such as "holm", "hochberg",
#' "bonferroni", etc (default is "holm") and 2) \code{B}: the number of
#' bootstrap samples (default is 100). Increase \code{B} will lead to a more
#' accurate p-values. See \code{Details} for a more comprehensive discussion on
#' mdFDR.
#' @param trend_control a named list of control parameters for the trend test,
#' including 1) \code{contrast}: the list of contrast matrices for
#' constructing inequalities, 2) \code{node}: the list of positions for the
#' nodal parameter, and 3) \code{B}: the number of bootstrap samples
#' (default is 100). Increase \code{B} will lead to a more accurate p-values.
#' See \code{vignette} for the corresponding trend test examples.
#'
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{feature_table}, a \code{data.frame} of pre-processed
#'         (based on \code{prv_cut} and \code{lib_cut}) microbial count table.}
#'         \item{ \code{bias_correct_log_table}, a \code{data.frame} of
#'         bias-corrected log abundance table.}
#'         \item{ \code{ss_tab}, a \code{data.frame} of sensitivity scores for
#'         pseudo-count addition to 0s.}
#'         \item{ \code{zero_ind}, a logical \code{data.frame} with TRUE
#'         indicating the taxon is detected to contain structural zeros in
#'         some specific groups.}
#'         \item{ \code{samp_frac}, a numeric vector of estimated sampling
#'         fractions in log scale (natural log).}
#'         \item{ \code{delta_em}, estimated sample-specific biases
#'         through E-M algorithm.}
#'         \item{ \code{delta_wls}, estimated sample-specific biases through
#'         weighted least squares (WLS) algorithm.}
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
#'         \item{ columns started with \code{passed_ss}: TRUE if the
#'         taxon passed the sensitivity analysis, i.e., adding different
#'         pseudo-counts to 0s would not change the results.}
#'         \item{ columns started with \code{diff_robust}: TRUE if the
#'         taxon is significant (has \code{q} less than \code{alpha}) and robust
#'         in the sensitivity analysis (\code{passed_ss} is \code{TRUE}).}
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
#'         \item{ \code{passed_ss}, A logical vector. TRUE if the taxon has
#'         passed the sensitivity analysis.}
#'         \item{ \code{diff_robust_abn}, A logical vector. TRUE if the taxon
#'         has \code{diff_abn} equal to TRUE and \code{passed_ss} equal to
#'         TRUE.}
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
#'         \item{ columns started with \code{passed_ss}: TRUE if the
#'         taxon has passed the sensitivity analysis.}
#'         \item{ columns started with \code{diff_robust}: TRUE if the
#'         taxon is significant (has \code{q} less than \code{alpha}) and has
#'         passed the sensitivity analysis.}
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
#'         \item{ columns started with \code{passed_ss}: TRUE if the
#'         taxon has passed the sensitivity analysis.}
#'         \item{ columns started with \code{diff_robust}: TRUE if the
#'         taxon is significant (has \code{q} less than \code{alpha}) and has
#'         passed the sensitivity analysis.}
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
#'         \item{ \code{passed_ss}, A logical vector. TRUE if the taxon has
#'         passed the sensitivity analysis. The sensitivity scores of the trend
#'         test are those of the global test. Under
#'         \code{conservative = FALSE}, these scores are compared with the
#'         p-value of the trend test.}
#'         \item{ \code{diff_robust_abn}, A logical vector. TRUE if the taxon
#'         has \code{diff_abn} equal to TRUE and \code{passed_ss} equal to
#'         TRUE.}
#'         }
#'         }
#'
#' @seealso \code{\link{ancom}} \code{\link{ancombc}}
#' \code{\link{data_sanity_check}} \code{\link{ANCOMBC-concepts}}
#'
#' @examples
#' library(ANCOMBC)
#' if (requireNamespace("microbiome", quietly = TRUE)) {
#'     data(dietswap, package = "microbiome")
#'
#'     # run ancombc2 function
#'     set.seed(123)
#'     # Setting max_iter = 1 and B = 1 reduces the running time of the example
#'     # The results printed below are for illustration only
#'     # Use the default or larger values for max_iter and B in an analysis
#'     out = ancombc2(data = dietswap, tax_level = "Family",
#'                    fix_formula = "nationality + timepoint + bmi_group",
#'                    rand_formula = NULL,
#'                    p_adj_method = "holm", pseudo_sens = TRUE,
#'                    prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#'                    group = "bmi_group", struc_zero = TRUE, neg_lb = TRUE,
#'                    alpha = 0.05, n_cl = 1, verbose = TRUE,
#'                    global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
#'                    iter_control = list(tol = 1e-2, max_iter = 1, verbose = TRUE),
#'                    em_control = list(tol = 1e-5, max_iter = 1),
#'                    lme_control = lme4::lmerControl(),
#'                    mdfdr_control = list(fwer_ctrl_method = "holm", B = 1),
#'                    trend_control = list(contrast =
#'                                             list(matrix(c(1, 0, -1, 1),
#'                                                         nrow = 2,
#'                                                         byrow = TRUE)),
#'                                         node = list(2),
#'                                         B = 1))
#'     res_prim = out$res
#'     res_global = out$res_global
#'     res_pair = out$res_pair
#'     res_dunn = out$res_dunn
#'     res_trend = out$res_trend
#' } else {
#'     message("The 'microbiome' package is not installed. Please install it to use this example.")
#' }
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
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
#' @importFrom multcomp glht mcp adjusted
#' @importFrom quadprog solve.QP
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar% %:% registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG %dorng%
#' @importFrom MASS ginv
#' @importFrom nloptr neldermead
#' @importFrom Rdpack reprompt
#'
#' @export
ancombc2 = function(data, taxa_are_rows = TRUE,
                    assay.type = assay_name, assay_name = "counts",
                    rank = tax_level, tax_level = NULL,
                    aggregate_data = NULL, meta_data = NULL,
                    fix_formula, rand_formula = NULL,
                    p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                    conservative = TRUE,
                    prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                    group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                    alpha = 0.05, n_cl = 1, verbose = TRUE,
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
                                         B = 100)){
    # Run with cluster
    if (n_cl > 1) {
      cl = parallel::makeCluster(n_cl)
      doParallel::registerDoParallel(cl)
    } else {
      foreach::registerDoSEQ()
    }

    # 1. Data pre-processing
    # Data sanity check
    check_results = data_sanity_check(data = data,
                                      taxa_are_rows = taxa_are_rows,
                                      assay.type = assay_name,
                                      assay_name = assay_name,
                                      rank = tax_level,
                                      tax_level = tax_level,
                                      aggregate_data = aggregate_data,
                                      meta_data = meta_data,
                                      fix_formula = fix_formula,
                                      group = group,
                                      struc_zero = struc_zero,
                                      global = global,
                                      pairwise = pairwise,
                                      dunnet = dunnet,
                                      mdfdr_control = mdfdr_control,
                                      trend = trend,
                                      trend_control = trend_control,
                                      verbose = verbose)
    feature_table = check_results$feature_table
    feature_table_aggregate = check_results$feature_table_aggregate
    meta_data = check_results$meta_data
    global = check_results$global
    pairwise = check_results$pairwise
    dunnet = check_results$dunnet
    trend = check_results$trend
    trend_control = check_results$trend_control

    # Identify taxa with structural zeros
    if (struc_zero) {
        zero_ind = .get_struc_zero(data = feature_table_aggregate,
                                   meta_data = meta_data,
                                   group = group, neg_lb = neg_lb)
        # Taxa with structural zeros will be removed from ANCOM-BC2 analyses
        tax_idx = apply(zero_ind[, -1], 1, function(x) all(x == FALSE))
        tax_keep = which(tax_idx)
    }else{
        zero_ind = NULL
        # Taxa with structural zeros will be removed from ANCOM-BC2 analyses
        tax_keep = seq(nrow(feature_table_aggregate))}

    # Filter data by prevalence and library size
    core1 = .data_core(data = feature_table,
                       meta_data = meta_data,
                       prv_cut = prv_cut, lib_cut = lib_cut,
                       tax_keep = NULL, samp_keep = NULL)
    O1 = core1$feature_table
    samp_keep = colnames(O1)

    core2 = .data_core(data = feature_table_aggregate,
                       meta_data = meta_data,
                       prv_cut = prv_cut, lib_cut = lib_cut,
                       tax_keep = tax_keep, samp_keep = samp_keep)
    O2 = core2$feature_table
    meta_data = core2$meta_data

    # 2. ANCOM-BC2 main analysis
    # The sensitivity analysis of the trend test uses the global test results
    global_req = global
    if (trend) global = TRUE

    res_main = .ancombc2_core(data = O1, aggregate_data = O2,
                              meta_data = meta_data, fix_formula = fix_formula,
                              rand_formula = rand_formula,
                              p_adj_method = p_adj_method,
                              pseudo = pseudo, s0_perc = s0_perc,
                              group = group, alpha = alpha, verbose = verbose,
                              global = global, pairwise = pairwise,
                              dunnet = dunnet, trend = trend,
                              iter_control = iter_control,
                              em_control = em_control,
                              lme_control = lme_control,
                              mdfdr_control = mdfdr_control,
                              trend_control = trend_control)

    # 3. Sensitivity analysis for pseudo-count addition to 0s
    if (!pseudo_sens) {
        warn_txt = paste("Sensitivity analysis is currently turned off.",
                         "Since sensitivity analysis is essential for reducing false positives,",
                         "it is highly recommended to enable it unless your primary focus is power.",
                         "Are you sure you want to proceed without it?", sep = "\n")
        message(warn_txt)
        ss_tab = NULL
        res = res_main$res
        res_global = res_main$res_global
        res_pair = res_main$res_pair
        res_dunn = res_main$res_dunn
        res_trend = res_main$res_trend
    }

    if (pseudo_sens) {
        message_txt = paste("Conducting sensitivity analysis for pseudo-count addition to 0s ...",
                            "For taxa that are significant but do not pass the sensitivity analysis,",
                            "they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.",
                            "For detailed instructions on performing sensitivity analysis, please refer to the package vignette.",
                            sep = "\n")
        message(message_txt)
    }

    # Pseudo-count addition before the estimation of sampling fractions
    if (pseudo_sens && conservative) {
        pseudo_list = c(0.1, 0.5, 1)
        iter_control$verbose = FALSE

        ss_list = lapply(pseudo_list, function(pseudo_count) {
            res_pseudo = .ancombc2_core(data = O1, aggregate_data = O2,
                                        meta_data = meta_data,
                                        fix_formula = fix_formula,
                                        rand_formula = rand_formula,
                                        p_adj_method = p_adj_method,
                                        pseudo = pseudo_count,
                                        s0_perc = s0_perc,
                                        group = group, alpha = alpha,
                                        verbose = FALSE,
                                        global = global, pairwise = pairwise,
                                        dunnet = dunnet, trend = FALSE,
                                        iter_control = iter_control,
                                        em_control = em_control,
                                        lme_control = lme_control,
                                        mdfdr_control = mdfdr_control)
            return(res_pseudo)
        })

        # Combine main results with sensitivity analysis results
        ss_list = c(list(res_main), ss_list)
        pseudo_list = c(0, pseudo_list)

        ## Primary results
        ss_list_prim = lapply(ss_list, function(res_pseudo)
            res_pseudo$res[, grepl("q_", colnames(res_pseudo$res))])
        ss_3d_prim = array(unlist(ss_list_prim),
                           c(dim(ss_list_prim[[1]]), length(ss_list_prim)))
        ss_tab_prim = apply(ss_3d_prim, c(1, 2), function(x) {
            sum(x > alpha)/length(pseudo_list)
        })
        colnames(ss_tab_prim) = gsub("^q_", "ss_prim_",
                                     colnames(ss_list_prim[[1]]))

        ss_tab_log = (ss_tab_prim == 0 | ss_tab_prim == 1)
        colnames(ss_tab_log) = gsub("ss_prim_", "passed_ss_",
                                    colnames(ss_tab_log))
        res = cbind(res_main$res, ss_tab_log)
        diff_cols = grep("^diff_", names(res), value = TRUE, perl = TRUE)
        suffixes = sub("^diff_", "", diff_cols)
        for (suffix in suffixes) {
            diff_col = paste0("diff_", suffix)
            passed_col = paste0("passed_ss_", suffix)
            new_col = paste0("diff_robust_", suffix)
            res[[new_col]] = res[[diff_col]] & res[[passed_col]]
        }

        ## Global and trend test results
        if (global) {
            ss_list_global = lapply(ss_list, function(res_pseudo)
                res_pseudo$res_global[, "q_val", drop = FALSE])
            ss_3d_global = array(unlist(ss_list_global),
                                 c(dim(ss_list_global[[1]]), length(ss_list_global)))
            ss_tab_global = apply(ss_3d_global, c(1, 2), function(x) {
                sum(x > alpha)/length(pseudo_list)
            })
            colnames(ss_tab_global) = "ss_global"

            ss_tab_log = (ss_tab_global == 0 | ss_tab_global == 1)
            colnames(ss_tab_log) = "passed_ss"
            res_global = cbind(res_main$res_global, ss_tab_log)
            res_global[["diff_robust_abn"]] = res_global[["diff_abn"]] &
                res_global[["passed_ss"]]
        } else { res_global = NULL }

        if (trend) {
            ss_tab_trend = ss_tab_global
            colnames(ss_tab_trend) = "ss_trend"

            ss_tab_log = (ss_tab_trend == 0 | ss_tab_trend == 1)
            colnames(ss_tab_log) = "passed_ss"
            res_trend = cbind(res_main$res_trend, ss_tab_log)
            res_trend[["diff_robust_abn"]] = res_trend[["diff_abn"]] &
                res_trend[["passed_ss"]]
        } else { res_trend = NULL }

        ## Pairwise test results
        if (pairwise) {
            ss_list_pair = lapply(ss_list, function(res_pseudo)
                res_pseudo$res_pair[, grepl("q_", colnames(res_pseudo$res_pair))])
            ss_3d_pair = array(unlist(ss_list_pair),
                               c(dim(ss_list_pair[[1]]), length(ss_list_pair)))
            ss_tab_pair = apply(ss_3d_pair, c(1, 2), function(x) {
                sum(x > alpha)/length(pseudo_list)
            })
            colnames(ss_tab_pair) = gsub("^q_", "ss_pair_",
                                         colnames(ss_list_pair[[1]]))

            ss_tab_log = (ss_tab_pair == 0 | ss_tab_pair == 1)
            colnames(ss_tab_log) = gsub("ss_pair_", "passed_ss_",
                                        colnames(ss_tab_log))
            res_pair = cbind(res_main$res_pair, ss_tab_log)
            diff_cols = grep("^diff_", names(res_pair), value = TRUE, perl = TRUE)
            suffixes = sub("^diff_", "", diff_cols)
            for (suffix in suffixes) {
                diff_col = paste0("diff_", suffix)
                passed_col = paste0("passed_ss_", suffix)
                new_col = paste0("diff_robust_", suffix)
                res_pair[[new_col]] = res_pair[[diff_col]] & res_pair[[passed_col]]
            }
        } else { res_pair = NULL }

        ## Dunnet's type of test results
        if (dunnet) {
            ss_list_dunn = lapply(ss_list, function(res_pseudo)
                res_pseudo$res_dunn[, grepl("q_", colnames(res_pseudo$res_dunn))])
            ss_3d_dunn = array(unlist(ss_list_dunn),
                               c(dim(ss_list_dunn[[1]]), length(ss_list_dunn)))
            ss_tab_dunn = apply(ss_3d_dunn, c(1, 2), function(x) {
                sum(x > alpha)/length(pseudo_list)
            })
            colnames(ss_tab_dunn) = gsub("^q_", "ss_dunn_",
                                         colnames(ss_list_dunn[[1]]))

            ss_tab_log = (ss_tab_dunn == 0 | ss_tab_dunn == 1)
            colnames(ss_tab_log) = gsub("ss_dunn_", "passed_ss_",
                                        colnames(ss_tab_log))
            res_dunn = cbind(res_main$res_dunn, ss_tab_log)
            diff_cols = grep("^diff_", names(res_dunn), value = TRUE, perl = TRUE)
            suffixes = sub("^diff_", "", diff_cols)
            for (suffix in suffixes) {
                diff_col = paste0("diff_", suffix)
                passed_col = paste0("passed_ss_", suffix)
                new_col = paste0("diff_robust_", suffix)
                res_dunn[[new_col]] = res_dunn[[diff_col]] & res_dunn[[passed_col]]
            }
        } else { res_dunn = NULL }

        ## Table of all sensitivity analysis results
        ss_tab_cols = list(taxon = rownames(O2), ss_tab_prim)
        if (global) ss_tab_cols$ss_tab_global = ss_tab_global
        if (pairwise) ss_tab_cols$ss_tab_pair = ss_tab_pair
        if (dunnet) ss_tab_cols$ss_tab_dunn = ss_tab_dunn
        if (trend) ss_tab_cols$ss_tab_trend = ss_tab_trend
        ss_tab = do.call(data.frame, c(ss_tab_cols, list(check.names = FALSE)))
        }

    # Pseudo-count addition to 0s of the bias-corrected data
    if (pseudo_sens && !conservative) {
        pseudo_list = seq(0.01, 0.5, 0.01)
        n_tax = nrow(O2)
        fix_eff = gsub("^q_", "",
                       grep("^q_", colnames(res_main$res), value = TRUE))
        if (is.null(group)) {
            group_eff = NULL
        } else {
            group_eff = fix_eff[grepl(group, fix_eff) & !grepl(":", fix_eff)]
        }

        fun_list = list(.ancombc2_sens_fit, .ancombc2_sens_p)

        pseudo_count = NULL
        ss_list = foreach(pseudo_count = pseudo_list) %dorng% {
            output = fun_list[[1]](data = O2, samp_frac = res_main$samp_frac,
                                   meta_data = meta_data,
                                   fix_formula = fix_formula,
                                   pseudo = pseudo_count,
                                   fix_eff = fix_eff, group = group,
                                   group_eff = group_eff,
                                   pairwise = pairwise, global = global,
                                   p_fun = fun_list[[2]])
        }

        # Proportion of the pseudo-count runs with a p-value above alpha
        ss_tab_fun = function(sens_col) {
            ss = vapply(seq_along(sens_col), function(j) {
                p_pseudo = vapply(ss_list, function(p) p[, sens_col[j]],
                                  FUN.VALUE = double(n_tax))
                rowMeans(p_pseudo > alpha)
            }, FUN.VALUE = double(n_tax))
            ss = matrix(ss, nrow = n_tax, ncol = length(sens_col))
            return(ss)
        }

        # Agreement between the main run and the pseudo-count runs
        passed_fun = function(ss_tab, p_main) {
            p_main = as.matrix(p_main)
            p_main[is.na(p_main)] = 1
            (ss_tab == 0 & p_main <= alpha) | (ss_tab == 1 & p_main > alpha)
        }

        # Flag the taxa that are robust to the pseudo-count addition
        flag_fun = function(res_tab, ss_tab, p_main) {
            ss_tab_log = passed_fun(ss_tab, p_main)
            colnames(ss_tab_log) = sub("^ss_(prim|pair|dunn)_", "passed_ss_",
                                       colnames(ss_tab))
            res_tab = cbind(res_tab, ss_tab_log)
            diff_cols = grep("^diff_", names(res_tab), value = TRUE, perl = TRUE)
            suffixes = sub("^diff_", "", diff_cols)
            for (suffix in suffixes) {
                diff_col = paste0("diff_", suffix)
                passed_col = paste0("passed_ss_", suffix)
                new_col = paste0("diff_robust_", suffix)
                res_tab[[new_col]] = res_tab[[diff_col]] & res_tab[[passed_col]]
            }
            return(res_tab)
        }

        ## Primary results
        ss_tab_prim = ss_tab_fun(fix_eff)
        colnames(ss_tab_prim) = paste0("ss_prim_", fix_eff)
        res = flag_fun(res_main$res, ss_tab_prim,
                       res_main$res[, paste0("p_", fix_eff)])

        ## Global and trend test results
        if (global) {
            ss_tab_global = ss_tab_fun("global")
            colnames(ss_tab_global) = "ss_global"

            ss_tab_log = passed_fun(ss_tab_global, res_main$res_global[, "p_val"])
            colnames(ss_tab_log) = "passed_ss"
            res_global = cbind(res_main$res_global, ss_tab_log)
            res_global[["diff_robust_abn"]] = res_global[["diff_abn"]] &
                res_global[["passed_ss"]]
        } else { res_global = NULL }

        if (trend) {
            ss_tab_trend = ss_tab_fun("global")
            colnames(ss_tab_trend) = "ss_trend"

            ss_tab_log = passed_fun(ss_tab_trend, res_main$res_trend[, "p_val"])
            colnames(ss_tab_log) = "passed_ss"
            res_trend = cbind(res_main$res_trend, ss_tab_log)
            res_trend[["diff_robust_abn"]] = res_trend[["diff_abn"]] &
                res_trend[["passed_ss"]]
        } else { res_trend = NULL }

        ## Pairwise test results
        if (pairwise) {
            pair_col = gsub("^q_", "",
                            grep("^q_", colnames(res_main$res_pair), value = TRUE))
            ss_tab_pair = ss_tab_fun(pair_col)
            colnames(ss_tab_pair) = paste0("ss_pair_", pair_col)
            res_pair = flag_fun(res_main$res_pair, ss_tab_pair,
                                res_main$res_pair[, paste0("p_", pair_col)])
        } else { res_pair = NULL }

        ## Dunnet's type of test results
        if (dunnet) {
            dunn_col = gsub("^q_", "",
                            grep("^q_", colnames(res_main$res_dunn), value = TRUE))
            ss_tab_dunn = ss_tab_fun(dunn_col)
            colnames(ss_tab_dunn) = paste0("ss_dunn_", dunn_col)
            res_dunn = flag_fun(res_main$res_dunn, ss_tab_dunn,
                                res_main$res_dunn[, paste0("p_", dunn_col)])
        } else { res_dunn = NULL }

        ## Table of all sensitivity analysis results
        ss_tab_cols = list(taxon = rownames(O2), ss_tab_prim)
        if (global) ss_tab_cols$ss_tab_global = ss_tab_global
        if (pairwise) ss_tab_cols$ss_tab_pair = ss_tab_pair
        if (dunnet) ss_tab_cols$ss_tab_dunn = ss_tab_dunn
        if (trend) ss_tab_cols$ss_tab_trend = ss_tab_trend
        ss_tab = do.call(data.frame, c(ss_tab_cols, list(check.names = FALSE)))
        }

    # 4. Outputs
    # The global test result is reported only when the global test is requested
    if (!global_req) res_global = NULL

    out = list(feature_table = O2,
               bias_correct_log_table = res_main$bias_correct_log_table,
               ss_tab = ss_tab,
               zero_ind = zero_ind,
               samp_frac = res_main$samp_frac,
               delta_em = res_main$delta_em,
               delta_wls = res_main$delta_wls,
               res = res,
               res_global = res_global,
               res_pair = res_pair,
               res_dunn = res_dunn,
               res_trend = res_trend)

    # Ensure cluster cleanup
    if (n_cl > 1) {
        parallel::stopCluster(cl)
        }

    return(out)
}



# Inference step of the sensitivity analysis for one pseudo-count
.ancombc2_sens_fit = function(data, samp_frac, meta_data, fix_formula, fix_eff,
                              pseudo, group, group_eff,
                              pairwise, global, p_fun) {
    O = as.matrix(data)
    O[O == 0] = pseudo
    o = log(O)
    y = o - rowMeans(o)
    y_bias_crt = t(t(y) - samp_frac)
    Y = data.frame(t(y_bias_crt), check.names = FALSE)

    p_list = lapply(Y, p_fun, meta_data = meta_data,
                    fix_formula = fix_formula, fix_eff = fix_eff,
                    group = group, group_eff = group_eff,
                    pairwise = pairwise, global = global)
    p_hat = do.call(rbind, p_list)
    p_hat[is.na(p_hat)] = 1
    return(p_hat)
}

# P-values of one taxon from the linear model fitted to the bias-corrected data
# Fixed effects that are not estimable are assigned a p-value of NA
.ancombc2_sens_p = function(y, meta_data, fix_formula, fix_eff, group,
                            group_eff, pairwise, global) {
    df = data.frame(y = y, meta_data)
    lm_fit = stats::lm(stats::formula(paste0("y ~ ", fix_formula)), data = df)
    coef_tab = summary(lm_fit)$coefficients
    p_val = rep(NA_real_, length(fix_eff))
    names(p_val) = fix_eff
    p_val[rownames(coef_tab)] = coef_tab[, "Pr(>|t|)"]

    if (pairwise) {
        beta_hat = stats::coef(lm_fit)[group_eff]
        vcov_hat = matrix(NA_real_, nrow = length(fix_eff),
                          ncol = length(fix_eff),
                          dimnames = list(fix_eff, fix_eff))
        vcov_fit = stats::vcov(lm_fit)
        vcov_hat[rownames(vcov_fit), colnames(vcov_fit)] = vcov_fit
        dof = lm_fit$df.residual
        combn_mat = utils::combn(group_eff, 2)
        pair_p_val = vapply(seq_len(ncol(combn_mat)), function(i) {
            id1 = combn_mat[2, i]
            id2 = combn_mat[1, i]
            beta_diff = beta_hat[id1] - beta_hat[id2]
            var_diff = vcov_hat[id1, id1] + vcov_hat[id2, id2] -
                2 * vcov_hat[id1, id2]
            2 * stats::pt(abs(beta_diff/sqrt(var_diff)), df = dof,
                          lower.tail = FALSE)
        }, FUN.VALUE = double(1))
        names(pair_p_val) = paste(combn_mat[2, ], combn_mat[1, ], sep = "_")
        p_val = c(p_val, pair_p_val)
    }

    if (global) {
        anova_tab = stats::anova(lm_fit)
        p_val = c(p_val,
                  global = anova_tab$`Pr(>F)`[rownames(anova_tab) == group])
    }

    return(p_val)
}
