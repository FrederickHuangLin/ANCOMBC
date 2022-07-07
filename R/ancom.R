#' @title Analysis of Composition of Microbiomes (ANCOM)
#'
#' @description Determine taxa whose absolute abundances, per unit volume, of
#' the ecosystem (e.g. gut) are significantly different with changes in the
#' covariate of interest (e.g. group). The current version of
#' \code{ancom} function implements ANCOM in cross-sectional and longitudinal
#' datasets while allowing for covariate adjustment.
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
#' @param x a phyloseq or (Tree)SummarizedExperiment object, which consists of
#' a feature table (microbial observed abundance table), a sample metadata, a
#' taxonomy table (optional), and a phylogenetic tree (optional). The row names
#' of the metadata must match the sample names of the feature table, and the
#' row names of the taxonomy table must match the taxon (feature) names of the
#' feature table. See \code{?phyloseq::phyloseq} or
#' \code{?SummarizedExperiment::SummarizedExperiment} for more details.
#' @param p_adj_method character. method to adjust p-values. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{?stats::p.adjust} for more details.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' less than \code{prv_cut} will be excluded in the analysis. Default
#' is 0.10.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis. Default is 0, i.e. do not discard any sample.
#' @param main_var character. The name of the main variable of interest.
#' @param adj_formula  character string representing the formula for
#' covariate adjustment. Default is \code{NULL}.
#' @param rand_formula character string representing the formula for random
#' effects. For details, see \code{?nlme::lme}. Default is \code{NULL}.
#' @param lme_control a list specifying control values for \code{lme} fit.
#' For details, see ?\code{nlme::lmeControl}. Default is \code{NULL}.
#' @param struc_zero logical. whether to detect structural zeros based on
#' \code{main_var}. \code{main_var} should be discrete. Default is FALSE.
#' @param neg_lb logical. whether to classify a taxon as a structural zero using
#' its asymptotic lower bound. Default is FALSE.
#' @param alpha numeric. level of significance. Default is 0.05.
#' @param n_cl numeric. The number of nodes to be forked. For details, see
#' \code{?parallel::makeCluster}. Default is 1 (no parallel computing).
#' @param assay_name character. Name of the abundance table in the data object
#' (only applicable if data object is a (Tree)SummarizedExperiment).
#' @param phyloseq a phyloseq-class object. Will be deprecated.
#'
#' @return a \code{list} with components:
#'         \itemize{
#'         \item{ \code{res},  a \code{data.frame} containing ANCOM
#'         result for the variable specified in \code{main_var},
#'         each column is:}
#'         \itemize{
#'         \item{ \code{W}, test statistics.}
#'         \item{ \code{detected_0.9, detected_0.8, detected_0.7, detected_0.6},
#'         logical vectors representing whether a taxon is differentially
#'         abundant under a series of cutoffs. For example, TRUE in
#'         \code{detected_0.7} means the number of ALR transformed models where
#'         the taxon is differentially abundant with regard to the main variable
#'         outnumbers \code{0.7 * (n_taxa - 1)}. \code{detected_0.7} is commonly
#'         used. Choose \code{detected_0.8} or \code{detected_0.9} for more
#'         conservative results, or choose \code{detected_0.6} for more liberal
#'         results.}
#'         }
#'         \item{ \code{zero_ind}, a logical \code{matrix} with TRUE indicating
#'         the taxon is identified as a structural zero for the specified
#'         main variable.}
#'         \item{ \code{beta_data}, a numeric \code{matrix} containing pairwise
#'         coefficients for the main variable of interest in ALR transformed
#'         regression models.}
#'         \item{ \code{p_data}, a numeric \code{matrix} containing pairwise
#'         p-values for the main variable of interest in ALR transformed
#'         regression models.}
#'         \item{ \code{q_data}, a numeric \code{matrix} containing adjusted
#'         p-values by applying the \code{p_adj_method} to the \code{p_data}
#'         matrix.}
#'         }
#'
#' @seealso \code{\link{ancombc}}
#'
#' @examples
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
#' set.seed(123)
#' out = ancom(phyloseq = family_data,  p_adj_method = "holm",
#'             prv_cut = 0.10, lib_cut = 0, main_var = "nationality",
#'             adj_formula = "bmi_group",
#'             rand_formula = NULL, lme_control = NULL,
#'             struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 2)
#'
#' res = out$res
#'
#' @author Huang Lin
#'
#' @references
#' \insertRef{mandal2015analysis}{ANCOMBC}
#'
#' @rawNamespace import(stats, except = filter)
#' @import phyloseq
#' @import microbiome
#' @importFrom nlme lme
#' @importFrom dplyr mutate
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom Rdpack reprompt
#'
#' @export
ancom = function(x = phyloseq,  p_adj_method = "holm", prv_cut = 0.10, 
                 lib_cut = 0, main_var, adj_formula = NULL, rand_formula = NULL,
                 lme_control = NULL, struc_zero = FALSE, neg_lb = FALSE,
                 alpha = 0.05, n_cl = 1, assay_name = "counts", phyloseq){
  # 1. Data pre-processing
  fiuo_core = data_core(x, prv_cut, lib_cut,
                        tax_keep = NULL, samp_keep = NULL, assay_name)
  feature_table = fiuo_core$feature_table
  meta_data = fiuo_core$meta_data
  meta_data[] = lapply(meta_data, function(x)
    if(is.factor(x)) factor(x) else x)
  # Check the type of main variable
  main_val = meta_data[, main_var]
  main_class = class(main_val)
  if (main_class == "character") {
    if (length(unique(main_val)) == length(main_val)) {
      warn_txt = sprintf(paste("The class of main varible is:",
                               main_class,
                               "but it contains too many categories. Perhaps it should be numeric?",
                               sep = "\n"))
      warning(warn_txt, call. = FALSE)
    } else if (length(unique(main_val)) > 2) {
      main_cat = 1
    } else if (length(unique(main_val)) == 2) {
      main_cat = 0
    } else {
      stop_txt = sprintf(paste("The class of main varible is:",
                               main_class,
                               "but it contains < 2 categories",
                               sep = "\n"))
      stop(stop_txt, call. = FALSE)
    }
  } else if (main_class == "factor") {
    if (nlevels(main_val) == length(main_val)) {
      warn_txt = sprintf(paste("The class of main varible is:",
                               main_class,
                               "but it contains too many categories. Perhaps it should be numeric?",
                               sep = "\n"))
      warning(warn_txt, call. = FALSE)
    } else if (nlevels(main_val) > 2) {
      main_cat = 1
    } else if (nlevels(main_val) == 2) {
      main_cat = 0
    } else {
      stop_txt = sprintf(paste("The class of main varible is:",
                               main_class,
                               "but it contains < 2 categories",
                               sep = "\n"))
      stop(stop_txt, call. = FALSE)
    }
  } else {
    main_cat = 0
  }
  # Check the number of taxa
  n_taxa = nrow(feature_table)
  if (n_taxa < 10) {
    warn_txt = sprintf(paste("ANCOM results would be unreliable when the number of taxa is too small (e.g. < 10)",
                             "The number of taxa after filtering is: ",
                             n_taxa, sep = "\n"))
    warning(warn_txt, call. = FALSE)
  }

  # 2. Identify taxa with structural zeros
  if (struc_zero) {
    if (! main_class %in% c("character", "factor")) {
      stop_txt = sprintf(paste("The main variable should be discrete for detecting structural zeros.",
                               "Otherwise, set struc_zero = FALSE to proceed"))
      stop(stop_txt, call. = FALSE)
    }
    zero_ind = get_struc_zero(feature_table, meta_data, main_var, neg_lb)
    num_struc_zero = apply(zero_ind, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    zero_ind = NULL
    comp_table = feature_table
  }
  # Add pseudo-count (1) and take logarithm
  comp_table = log(as.matrix(comp_table) + 1)
  n_taxa = dim(comp_table)[1]
  taxon_id = rownames(comp_table)
  n_samp = dim(comp_table)[2]

  # 3. Determine the type of statistical test and its formula.
  if (is.null(rand_formula)) {
    # Model: linear model
    tfun = stats::lm
    # Formula
    if (is.null(adj_formula)) {
      tformula = formula(paste("x ~", main_var), sep = " ")
    }else {
      tformula = formula(paste("x ~", main_var, "+", adj_formula), sep = " ")
    }
  }else if (!is.null(rand_formula)) {
    # Model: linear mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var), sep = " ")
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }

  # 4. Computing pairwise p-values and effect sizes
  comb = function(...) {
    mapply("rbind", ..., SIMPLIFY = FALSE)
  }

  idx1 = NULL
  if (main_cat == 0) {
    cl = makeCluster(n_cl)
    registerDoParallel(cl)

    result = foreach(idx1 = seq_len(n_taxa), .combine = comb, .multicombine = TRUE) %dopar% {
      alr_data = apply(comp_table, 1, function(x) x - comp_table[idx1, ])
      alr_data = cbind(alr_data, meta_data)

      p_vec = rep(NA, n_taxa)
      beta_vec = rep(NA, n_taxa)

      idx2 = NULL
      if (is.null(rand_formula)) {
        for (idx2 in seq_len(n_taxa)) {
          test_data = data.frame(x = alr_data[, idx2],
                                 meta_data,
                                 check.names = FALSE)
          lm_fit = suppressWarnings(tfun(tformula, data = test_data))
          p_vec[idx2] = summary(lm_fit)$coef[2, "Pr(>|t|)"]
          beta_vec[idx2] = summary(lm_fit)$coef[2, "t value"]
        }
      } else {
        for (idx2 in seq_len(n_taxa)) {
          test_data = data.frame(x = alr_data[, idx2],
                                 meta_data,
                                 check.names = FALSE)
          lme_fit = try(tfun(fixed = tformula,
                             data = test_data,
                             random = formula(rand_formula),
                             na.action = na.omit,
                             control = lme_control),
                        silent = TRUE)
          if (inherits(lme_fit, "try-error")) {
            p_vec[idx2] = NA
            beta_vec[idx2] = NA
          } else {
            p_vec[idx2] = anova(lme_fit)[main_var, "p-value"]
            beta_vec[idx2] = lme_fit$coefficients$fixed[2]
          }
        }
      }

      list(p_vec, beta_vec)
    }
    stopCluster(cl)
  } else {
    cl = makeCluster(n_cl)
    registerDoParallel(cl)

    result = foreach(idx1 = seq_len(n_taxa), .combine = comb, .multicombine = TRUE) %dopar% {
      alr_data = apply(comp_table, 1, function(x) x - comp_table[idx1, ])
      alr_data = cbind(alr_data, meta_data)

      p_vec = rep(NA, n_taxa)
      beta_vec = rep(NA, n_taxa)

      idx2 = NULL
      if (is.null(rand_formula)) {
        for (idx2 in seq_len(n_taxa)) {
          test_data = data.frame(x = alr_data[, idx2],
                                 meta_data,
                                 check.names = FALSE)
          lm_fit = suppressWarnings(tfun(tformula, data = test_data))
          p_vec[idx2] = anova(lm_fit)[main_var, "Pr(>F)"]
          beta_vec[idx2] = anova(lm_fit)[main_var, "F value"]
        }
      } else {
        for (idx2 in seq_len(n_taxa)) {
          test_data = data.frame(x = alr_data[, idx2],
                                 meta_data,
                                 check.names = FALSE)
          lme_fit = try(tfun(fixed = tformula,
                             data = test_data,
                             random = formula(rand_formula),
                             na.action = na.omit,
                             control = lme_control),
                        silent = TRUE)
          if (inherits(lme_fit, "try-error")) {
            p_vec[idx2] = NA
            beta_vec[idx2] = NA
          } else {
            p_vec[idx2] = anova(lme_fit)[main_var, "p-value"]
            beta_vec[idx2] = anova(lme_fit)[main_var, "F-value"]
          }
        }
      }

      list(p_vec, beta_vec)
    }
    stopCluster(cl)
  }

  p_data = result[[1]]
  beta_data = result[[2]]

  colnames(p_data) = taxon_id
  rownames(p_data) = taxon_id
  colnames(beta_data) = taxon_id
  rownames(beta_data) = taxon_id

  diag(p_data) = 1
  p_data[is.na(p_data)] = 1
  diag(beta_data) = 0
  beta_data[is.na(beta_data)] = 0

  # 5. Primary results
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))
  W = apply(q_data, 2, function(x) sum(x < alpha))
  res_comp = data.frame(taxon_id, W, row.names = NULL, check.names = FALSE)
  res_comp = res_comp %>%
    mutate(detected_0.9 = ifelse(.data$W > 0.9 * (n_taxa - 1), TRUE, FALSE),
           detected_0.8 = ifelse(.data$W > 0.8 * (n_taxa - 1), TRUE, FALSE),
           detected_0.7 = ifelse(.data$W > 0.7 * (n_taxa - 1), TRUE, FALSE),
           detected_0.6 = ifelse(.data$W > 0.6 * (n_taxa - 1), TRUE, FALSE))

  # Combine the information of structural zeros
  if (struc_zero){
    res = data.frame(taxon_id = rownames(zero_ind), W = Inf,
                     detected_0.9 = TRUE, detected_0.8 = TRUE,
                     detected_0.7 = TRUE, detected_0.6 = TRUE,
                     row.names = NULL, check.names = FALSE)
    res[match(taxon_id, res$taxon_id), ] = res_comp
  }else{
    res = res_comp
  }

  out = list(res = res, zero_ind = zero_ind,
             beta_data = beta_data, p_data = p_data, q_data = q_data)
  return(out)
}






