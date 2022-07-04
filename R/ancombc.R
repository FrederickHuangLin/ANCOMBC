#' @title Analysis of Compositions of Microbiomes with Bias Correction
#' (ANCOM-BC)
#'
#' @description Determine taxa whose absolute abundances, per unit volume, of
#' the ecosystem (e.g. gut) are significantly different with changes in the
#' covariate of interest (e.g. group). The current version of
#' \code{ancombc} function implements Analysis of Compositions of Microbiomes
#' with Bias Correction (ANCOM-BC) in cross-sectional data while allowing
#' for covariate adjustment.
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
#' table. See \code{?phyloseq::phyloseq} for more details.
#' @param formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the variables in metadata.
#' @param p_adj_method character. method to adjust p-values. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{?stats::p.adjust} for more details.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' less than \code{prv_cut} will be excluded in the analysis. Default
#' is 0.10.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis. Default is 0, i.e. do not discard any sample.
#' @param group character. the name of the group variable in metadata.
#' \code{group} should be discrete. Specifying \code{group} is required for
#' detecting structural zeros and performing global test.
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
#' @param global logical. whether to perform global test. Default is FALSE.
#' @param assay_name character. Name of the abundance table in the data object
#' (only applicable if data object is a (Tree)SummarizedExperiment).
#'
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{feature_table}, a \code{data.frame} of pre-processed
#'         (based on \code{prv_cut} and \code{lib_cut}) microbial observed
#'         abundance table.}
#'         \item{ \code{zero_ind}, a logical \code{matrix} with TRUE indicating
#'         the taxon is identified as a structural zero for the specified
#'         \code{group} variable.}
#'         \item{ \code{samp_frac}, a numeric vector of estimated sampling
#'         fractions in log scale (natural log). Note that for each sample,
#'         if it contains missing values for any variable specified in the
#'         \code{formula}, the corresponding sampling fraction estimate
#'         for this sample will be \code{NA} since the sampling fraction
#'         is not estimable with the presence of missing values.}
#'         \item{ \code{resid}, a \code{matrix} of residuals from the ANCOM-BC
#'         log-linear (natural log) model.
#'         Rows are taxa and columns are samples.}
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
#' @seealso \code{\link{ancom}}
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
#' library(microbiome)
#' library(tidyverse)
#' data(dietswap)
#'
#' # Subset to baseline
#' pseq = subset_samples(dietswap, timepoint == 1)
#' # Aggregate to family level
#' family_data = aggregate_taxa(pseq, "Family")
#'
#' # Run ancombc function
#' out = ancombc(phyloseq = family_data, formula = "bmi_group + nationality",
#'               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#'               group = "bmi_group", struc_zero = TRUE, neg_lb = FALSE,
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
#' @rawNamespace import(stats, except = filter)
#' @import phyloseq
#' @import microbiome
#' @importFrom MASS ginv
#' @importFrom nloptr neldermead
#' @importFrom Rdpack reprompt
#'
#' @export
ancombc = function(data, formula, p_adj_method = "holm", prv_cut = 0.10,
                   lib_cut = 0, group = NULL, struc_zero = FALSE,
                   neg_lb = FALSE, tol = 1e-05, max_iter = 100,
                   conserve = FALSE, alpha = 0.05, global = FALSE, 
                   assay_name = "counts"){
  # 1. Data pre-processing
  fiuo_core = data_core(data, prv_cut, lib_cut,
                        tax_keep = NULL, samp_keep = NULL, assay_name)
  feature_table = fiuo_core$feature_table
  meta_data = fiuo_core$meta_data
  n_taxa = nrow(feature_table)
  if (n_taxa < 10) {
    warn_txt = sprintf(paste("ANCOM-BC results would be unreliable when the number of taxa is too small (e.g. < 10)",
                             "The number of taxa in the current dataset is: ",
                             n_taxa, sep = "\n"))
    warning(warn_txt, call. = FALSE)
  }

  fiuo_qc = data_qc(meta_data, group, global)
  meta_data = fiuo_qc$meta_data
  global = fiuo_qc$global

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
      stop_txt = sprintf(paste("Please specify the group variable for detecting structural zeros.",
                               "Otherwise, set struc_zero = FALSE to proceed"))
      stop(stop_txt, call. = FALSE)
    }
    zero_ind = get_struc_zero(feature_table, meta_data, group, neg_lb)
  }else{ zero_ind = NULL }

  # 3. Estimation of parameters
  fiuo_para = para_est(x, y, meta_data, formula, tol, max_iter)
  beta = fiuo_para$beta
  d = fiuo_para$d
  e = fiuo_para$e
  var_cov_hat = fiuo_para$var_cov_hat
  var_hat = fiuo_para$var_hat

  # 4. Estimation of the sample-specific bias
  fiuo_bias = bias_est(beta, var_hat, tol, max_iter)
  delta_em = fiuo_bias$delta_em
  delta_wls = fiuo_bias$delta_wls
  var_delta = fiuo_bias$var_delta

  # 5. Obtain coefficients, standard errors, and sampling fractions
  beta_hat = beta
  beta_hat[, -1] = t(t(beta_hat[, -1]) - delta_em)

  if (conserve) {
    # Account for the variance of delta_hat
    se_hat = sqrt(sweep(var_hat, 2, c(0, var_delta), "+") +
                    2 * sqrt(sweep(var_hat, 2, c(0, var_delta), "*")))
  }else{ se_hat = sqrt(var_hat) }

  d_hat = matrix(NA, nrow = nrow(y), ncol = ncol(y))
  for (i in seq_len(n_taxa)) {
    d_hat[i, ] = y[i, ] - x %*% beta_hat[i, ]
  }
  d_hat = colMeans(d_hat, na.rm = TRUE)

  # Remove uninformative intercept column
  beta_hat = beta_hat[, setdiff(colnames(beta_hat), "(Intercept)"),
                      drop = FALSE]
  se_hat = se_hat[, setdiff(colnames(se_hat), "(Intercept)"),
                  drop = FALSE]

  # 6. Primary results
  W = beta_hat/se_hat
  p = 2 * pnorm(abs(W), mean = 0, sd = 1, lower.tail = FALSE)
  q = apply(p, 2, function(x) p.adjust(x, method = p_adj_method))
  diff_abn = q < alpha & !is.na(q)
  res = list(lfc = data.frame(beta_hat, check.names = FALSE),
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
  # Set p/q-values and SEs of structural zeros to be 0s.
  if (struc_zero) {
    group_ind = grepl(group, setdiff(covariates, "(Intercept)"))
    zero_mask = 1 - abs((zero_ind - zero_ind[, 1]))
    zero_mask = zero_mask[, -1, drop = FALSE]
    res$se[, group_ind] = res$se[, group_ind] * zero_mask
    res$p_val[, group_ind] = res$p_val[, group_ind] * zero_mask
    res$q_val[, group_ind] = res$q_val[, group_ind] * zero_mask
    res$diff_abn = data.frame(res$q_val < alpha & !is.na(res$q_val),
                              check.names = FALSE)

    # Global test
    if (global) {
      zero_mask = 1 - apply(zero_ind, 1, function(x)
        sum(x) > 0 & sum(x) < ncol(zero_ind))
      res_global[, "p_val"] = res_global[, "p_val"] * zero_mask
      res_global[, "q_val"] = res_global[, "q_val"] * zero_mask
      res_global[, "diff_abn"] = res_global[, "q_val"] < alpha &
        !is.na(res_global[, "q_val"])
    }
  }

  # 9. Outputs
  out = list(feature_table = feature_table, zero_ind = zero_ind,
             samp_frac = d_hat, resid = e,
             delta_em = delta_em, delta_wls = delta_wls,
             res = res, res_global = res_global)
  return(out)
}



