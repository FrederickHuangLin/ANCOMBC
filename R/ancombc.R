#' @title Differential abundance (DA) analysis for
#' microbial absolute abundance data.
#'
#' @aliases ancom
#'
#' @description Determine taxa whose absolute abundances, per unit volume, of
#' the ecosystem (e.g. gut) are significantly different with changes in the
#' covariate of interest (e.g. the group effect). The current version of
#' \code{ancombc} function implements Analysis of Compositions of Microbiomes
#' with Bias Correction (ANCOM-BC) in cross-sectional data while allowing
#' the adjustment of covariates.
#'
#' @details The definition of structural zero can be found at
#' \href{https://doi.org/10.3389/fmicb.2017.02114}{ANCOM-II}.
#' Setting \code{neg_lb = TRUE} indicates that you are using both criteria
#' stated in section 3.2 of
#' \href{https://doi.org/10.3389/fmicb.2017.02114}{ANCOM-II}
#' to detect structural zeros; otherwise, the algorithm will only use the
#' equation 1 in section 3.2 for declaring structural zeros. Generally, it is
#' recommended to set \code{neg_lb = TRUE} when the sample size per group is
#' relatively large (e.g. > 30).
#'
#' @param phyloseq a phyloseq-class object, which consists of a feature table
#' (microbial observed abundance table), a sample metadata, a taxonomy table
#' (optional), and a phylogenetic tree (optional). The row names of the
#' metadata must match the sample names of the feature table, and the row names
#' of the taxonomy table must match the taxon (feature) names of the feature
#' table. See \code{\link[phyloseq]{phyloseq}} for more details.
#' @param formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the variables in metadata.
#' @param p_adj_method method to adjust p-values by. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{\link[stats]{p.adjust}} for more details.
#' @param zero_cut a numerical fraction between 0 and 1. Taxa with proportion of
#' zeroes greater than \code{zero_cut} will be excluded in the analysis. Default
#' is 0.90.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis.
#' @param group the name of the group variable in metadata. Specifying
#' \code{group} is required for detecting structural zeros and
#' performing global test.
#' @param struc_zero whether to detect structural zeros. Default is FALSE.
#' @param neg_lb whether to classify a taxon as a structural zero in the
#' corresponding study group using its asymptotic lower bound.
#' Default is FALSE.
#' @param tol the iteration convergence tolerance for the E-M algorithm.
#' Default is 1e-05.
#' @param max_iter the maximum number of iterations for the E-M algorithm.
#' Default is 100.
#' @param conserve whether to use a conservative variance estimate of
#' the test statistic. It is recommended if the sample size is small and/or
#' the number of differentially abundant taxa is believed to be large.
#' Default is FALSE.
#' @param alpha level of significance. Default is 0.05.
#' @param global whether to perform global test. Default is FALSE.
#'
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{feature_table}, a \code{data.frame} of pre-processed
#'         (based on \code{zero_cut} and \code{lib_cut}) microbial observed
#'         abundance table. }
#'         \item{ \code{zero_ind}, a logical \code{matrix} with TRUE indicating
#'         the taxon is identified as a structural zero for the specified
#'         \code{group} variable.}
#'         \item{ \code{samp_frac}, a numeric vector of estimated sampling
#'         fractions in log scale (natural log). Note that for each sample,
#'         if it contains missing values for any variable specified in the
#'         \code{formula}}, the corresponding sampling fraction estimate
#'         for this sample will return \code{NA} since the sampling fraction
#'         is not estimable with the presence of missing values.
#'         \item{ \code{resid}, a \code{matrix} of residuals from the ANCOM-BC
#'         log-linear (natural log) model.
#'         Rows are taxa and columns are samples.}
#'         \item{ \code{delta_em}, estimated bias terms through E-M algorithm. }
#'         \item{ \code{delta_wls}, estimated bias terms through weighted
#'         least squares (WLS) algorithm.}
#'         \item{ \code{res},  a \code{list} containing ANCOM-BC primary result,
#'         which consists of:}
#'         \itemize{
#'         \item{ \code{beta}, a \code{data.frame} of coefficients obtained
#'         from the ANCOM-BC log-linear (natural log) model. }
#'         \item{ \code{se}, a \code{data.frame} of standard errors (SEs) of
#'         \code{beta}. }
#'         \item{ \code{W}, a \code{data.frame} of test statistics.
#'         \code{W = beta/se}. }
#'         \item{ \code{p_val}, a \code{data.frame} of p-values. P-values are
#'         obtained from two-sided Z-test using the test statistic \code{W}. }
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
#' @examples
#' #================Build a Phyloseq-Class Object from Scratch==================
#' library(phyloseq)
#'
#' otu_mat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
#' rownames(otu_mat) = paste0("taxon", 1:nrow(otu_mat))
#' colnames(otu_mat) = paste0("sample", 1:ncol(otu_mat))
#'
#'
#' meta = data.frame(group = sample(LETTERS[1:4], size = 10, replace = TRUE),
#'                   row.names = paste0("sample", 1:ncol(otu_mat)),
#'                   stringsAsFactors = FALSE)
#'
#' tax_mat = matrix(sample(letters, 70, replace = TRUE),
#'                  nrow = nrow(otu_mat), ncol = 7)
#' rownames(tax_mat) = rownames(otu_mat)
#' colnames(tax_mat) = c("Kingdom", "Phylum", "Class", "Order",
#'                       "Family", "Genus", "Species")
#'
#' OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
#' META = sample_data(meta)
#' TAX = tax_table(tax_mat)
#' physeq = phyloseq(OTU, META, TAX)
#'
#' #========================Run ANCOMBC Using a Real Data=======================
#'
#' library(phyloseq)
#' library(microbiome)
#' library(tidyverse)
#' data(GlobalPatterns)
#'
#' # Aggregate to phylum level
#' phylum_data = aggregate_taxa(GlobalPatterns, "Phylum")
#' # The taxonomy table
#' tax_mat = as(tax_table(phylum_data), "matrix")
#'
#' # Run ancombc function
#' out = ancombc(phyloseq = phylum_data, formula = "SampleType",
#'               p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
#'               group = "SampleType", struc_zero = TRUE, neg_lb = FALSE,
#'               tol = 1e-5, max_iter = 100, conserve = TRUE,
#'               alpha = 0.05, global = TRUE)
#'
#' res = out$res
#' res_global = out$res_global
#'
#' @author Huang Lin
#'
#' @references
#' \insertRef{kaul2017analysis}{ANCOMBC}
#'
#' \insertRef{lin2020analysis}{ANCOMBC}
#'
#' @import stats
#' @import phyloseq
#' @import microbiome
#' @importFrom MASS ginv
#' @importFrom nloptr neldermead
#' @importFrom Rdpack reprompt
#'
#' @export
ancombc = function(phyloseq, formula, p_adj_method = "holm", zero_cut = 0.90,
                   lib_cut, group = NULL, struc_zero = FALSE, neg_lb = FALSE,
                   tol = 1e-05, max_iter = 100, conserve = FALSE, alpha = 0.05,
                   global = FALSE){
  # 1. Data pre-processing
  fiuo_prep = data_prep(phyloseq, group, zero_cut, lib_cut, global = global)
  feature_table = fiuo_prep$feature_table
  meta_data = fiuo_prep$meta_data
  global = fiuo_prep$global
  # samp_id = colnames(feature_table)
  # taxa_id = rownames(feature_table)
  # n_samp = ncol(feature_table)
  n_taxa = nrow(feature_table)
  if (n_taxa < 10) {
    warning_message = paste0("The number of taxa is too small.", "\n",
                             "ANCOMBC results would be unreliable since its methodology requires a relatively large number of taxa.", "\n",
                             "Typically, the number of taxa is >= 10, the current number of taxa is ", n_taxa, ".")
    warning(warning_message)
  }
  # Add pseudocount (1) and take logarithm.
  y = log(feature_table + 1)
  x = get_x(formula, meta_data)
  covariates = colnames(x)
  n_covariates = length(covariates)

  # 2. Identify taxa with structural zeros
  if (struc_zero) {
    if (is.null(group)) {
      stop("Please specify the group variable for detecting structural zeros.")
    }
    zero_ind = get_struc_zero(feature_table, meta_data, group, neg_lb)
  }else{ zero_ind = NULL }

  # 3. Estimation of parameters
  fiuo_para = para_est(y, meta_data, formula, tol, max_iter)
  beta = fiuo_para$beta
  d = fiuo_para$d
  e = fiuo_para$e
  var_cov_hat = fiuo_para$var_cov_hat
  var_hat = fiuo_para$var_hat

  # 4. Estimation of the between-sample bias
  fiuo_bias = bias_est(beta, var_hat, tol, max_iter, n_taxa)
  delta_em = fiuo_bias$delta_em
  delta_wls = fiuo_bias$delta_wls
  var_delta = fiuo_bias$var_delta

  # 5. Coefficients, standard error, and sampling fractions
  fiuo_fit = fit_summary(y, x, beta, var_hat, delta_em, var_delta, conserve)
  beta_hat = fiuo_fit$beta_hat
  se_hat = fiuo_fit$se_hat
  d_hat = fiuo_fit$d_hat

  # 6. Primary results
  W = beta_hat/se_hat
  p = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)
  q = apply(p, 2, function(x) p.adjust(x, method = p_adj_method))
  diff_abn = q < alpha & !is.na(q)
  res = list(beta = data.frame(beta_hat, check.names = FALSE),
             se = data.frame(se_hat, check.names = FALSE),
             W = data.frame(W, check.names = FALSE),
             p_val = data.frame(p, check.names = FALSE),
             q_val = data.frame(q, check.names = FALSE),
             diff_abn = data.frame(diff_abn, check.names = FALSE))

  # 7. Global test results
  if (global) {
    res_global = global_test(y, x, group, beta_hat, var_cov_hat,
                             p_adj_method, alpha)
  } else { res_global = NULL }

  # 8. Combine the information of structural zeros
  fiuo_out = res_combine_zero(x, group, struc_zero, zero_ind, alpha,
                              global, res, res_global)
  res = fiuo_out$res
  res_global = fiuo_out$res_global

  # 9. Outputs
  out = list(feature_table = feature_table, zero_ind = zero_ind,
             samp_frac = d_hat, resid = e,
             delta_em = delta_em, delta_wls = delta_wls,
             res = res, res_global = res_global)
  return(out)
}



