#' @title Analysis of Composition of Microbiomes (ANCOM)
#'
#' @description Determine taxa whose absolute abundances, per unit volume, of
#' the ecosystem (e.g. gut) are significantly different with changes in the
#' covariate of interest (e.g. group). The current version of
#' \code{ancom} function implements ANCOM in cross-sectional and repeated
#' measurements data while allowing for covariate adjustment.
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
#' @param assay_name character. Name of the count table in the data object
#' (only applicable if data object is a \code{(Tree)SummarizedExperiment}).
#' Default is "counts".
#' See \code{?SummarizedExperiment::assay} for more details.
#' @param tax_level character. The taxonomic level of interest. The input data
#' can be agglomerated at different taxonomic levels based on your research
#' interest. Default is NULL, i.e., do not perform agglomeration, and the
#' ANCOM anlysis will be performed at the lowest taxonomic level of the
#' input \code{data}.
#' @param phyloseq a \code{phyloseq} object. Will be deprecated.
#' @param p_adj_method character. method to adjust p-values. Default is "holm".
#' Options include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". See \code{?stats::p.adjust} for more details.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' less than \code{prv_cut} will be excluded in the analysis. For instance,
#' suppose there are 100 samples, if a taxon has nonzero counts presented in
#' less than 10 samples, it will not be further analyzed. Default is 0.10.
#' @param lib_cut a numerical threshold for filtering samples based on library
#' sizes. Samples with library sizes less than \code{lib_cut} will be
#' excluded in the analysis. Default is 0, i.e. do not discard any sample.
#' @param main_var character. The name of the main variable of interest.
#' @param adj_formula  character string representing the formula for
#' covariate adjustment. Default is \code{NULL}.
#' @param rand_formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the random effects in metadata. ANCOM
#' follows the \code{lmerTest} package in formulating the random effects. See
#' \code{?lmerTest::lmer} for more details. Default is \code{NULL}.
#' @param lme_control a list of control parameters for mixed model fitting.
#' See \code{?lme4::lmerControl} for details.
#' @param struc_zero logical. whether to detect structural zeros based on
#' \code{main_var}. \code{main_var} should be discrete. Default is FALSE.
#' @param neg_lb logical. whether to classify a taxon as a structural zero using
#' its asymptotic lower bound. Default is FALSE.
#' @param alpha numeric. level of significance. Default is 0.05.
#' @param n_cl numeric. The number of nodes to be forked. For details, see
#' \code{?parallel::makeCluster}. Default is 1 (no parallel computing).
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
#'         outnumbers \code{0.7 * (n_tax - 1)}. \code{detected_0.7} is commonly
#'         used. Choose \code{detected_0.8} or \code{detected_0.9} for more
#'         conservative results, or choose \code{detected_0.6} for more liberal
#'         results.}
#'         }
#'         \item{ \code{zero_ind}, a logical \code{data.frame} with TRUE
#'         indicating the taxon is detected to contain structural zeros in
#'         some specific groups.}
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
#' @seealso \code{\link{ancombc}} \code{\link{ancombc2}}
#'
#' @examples
#' library(ANCOMBC)
#' library(mia)
#' data(atlas1006)
#'
#' # subset to baseline
#' tse = atlas1006[, atlas1006$time == 0]
#'
#' # run ancom function
#' set.seed(123)
#' out = ancom(data = tse, assay_name = "counts",
#'             tax_level = "Family", phyloseq = NULL,
#'             p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#'             main_var = "bmi_group", adj_formula = "age + nationality",
#'             rand_formula = NULL, lme_control = NULL,
#'             struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 1)
#'
#' res = out$res
#'
#' # to run ancom using the phyloseq object
#' tse_alt = agglomerateByRank(tse, "Family")
#' pseq = makePhyloseqFromTreeSummarizedExperiment(tse_alt)
#' set.seed(123)
#' out = ancom(data = NULL, assay_name = NULL,
#'             tax_level = "Family", phyloseq = pseq,
#'             p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#'             main_var = "bmi_group", adj_formula = "age + nationality",
#'             rand_formula = NULL, lme_control = NULL,
#'             struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 1)
#'
#' @author Huang Lin
#'
#' @references
#' \insertRef{mandal2015analysis}{ANCOMBC}
#'
#' \insertRef{kaul2017analysis}{ANCOMBC}
#'
#' @rawNamespace import(stats, except = filter)
#' @importFrom mia makeTreeSummarizedExperimentFromPhyloseq taxonomyRanks agglomerateByRank
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
#' @importFrom dplyr mutate
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom Rdpack reprompt
#'
#' @export
ancom = function(data = NULL, assay_name = "counts", tax_level = NULL,
                 phyloseq = NULL, p_adj_method = "holm", prv_cut = 0.10,
                 lib_cut = 0, main_var, adj_formula = NULL, rand_formula = NULL,
                 lme_control = lme4::lmerControl(), struc_zero = FALSE,
                 neg_lb = FALSE, alpha = 0.05, n_cl = 1){
    # 1. Data pre-processing
    # TSE data construction
    tse_obj = tse_construct(data = data, assay_name = assay_name,
                            tax_level = tax_level, phyloseq = phyloseq)
    tse = tse_obj$tse
    assay_name = tse_obj$assay_name

    # Filter data by prevalence and library size
    core = data_core(tse = tse, assay_name = assay_name, alt = TRUE,
                     prv_cut = prv_cut, lib_cut = lib_cut,
                     tax_keep = NULL, samp_keep = NULL)
    feature_table = core$feature_table
    tax_keep = core$tax_keep
    tax_name = rownames(feature_table)
    meta_data = core$meta_data
    meta_data[] = lapply(meta_data, function(x)
        if(is.factor(x)) factor(x) else x)
    # Check the type of main variable
    main_val = meta_data[, main_var]
    main_class = class(main_val)
    if (main_class == "character") {
        if (length(unique(main_val)) == length(main_val)) {
            warn_txt = sprintf(paste("The class of main varible is:",
                                     main_class,
                                     "but it contains too many categories.",
                                     "Perhaps it should be numeric?",
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
                                     "but it contains too many categories.",
                                     "Perhaps it should be numeric?",
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
    n_tax = nrow(feature_table)
    if (n_tax < 10) {
        warn_txt = sprintf(paste("ANCOM results would be unreliable when the number of taxa is too small (e.g. < 10)",
                                 "The number of taxa after filtering is: ",
                                 n_tax, sep = "\n"))
        warning(warn_txt, call. = FALSE)
    }

    # 2. Identify taxa with structural zeros
    if (struc_zero) {
        if (! main_class %in% c("character", "factor")) {
            stop_txt = sprintf(paste("The main variable should be discrete for detecting structural zeros.",
                                     "Otherwise, set struc_zero = FALSE to proceed"))
            stop(stop_txt, call. = FALSE)
        }
        zero_ind = get_struc_zero(tse = tse, assay_name = assay_name,
                                  alt = TRUE, group = main_var,
                                  neg_lb = neg_lb)
        zero_ind = zero_ind[tax_keep, ]
        rownames(zero_ind) = NULL
        num_struc_zero = apply(zero_ind[, -1], 1, sum)
        comp_idx = which(num_struc_zero == 0)
        comp_table = feature_table[comp_idx, ]
    }else{
        zero_ind = NULL
        comp_table = feature_table
    }
    # Add pseudo-count (1) and take logarithm
    comp_table = log(as.matrix(comp_table) + 1)
    n_tax = dim(comp_table)[1]
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
        tfun = lmerTest::lmer
        # Formula
        if (is.null(adj_formula)) {
            # Random intercept model
            tformula = formula(paste0("x ~ ", main_var, "+ ", rand_formula))
        }else {
            # Random coefficients/slope model
            tformula = formula(paste0("x ~ ", main_var, "+ ",
                                      adj_formula, "+ ", rand_formula))
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

        result = foreach(idx1 = seq_len(n_tax), .combine = comb, .multicombine = TRUE) %dopar% {
            alr_data = apply(comp_table, 1, function(x) x - comp_table[idx1, ])
            alr_data = cbind(alr_data, meta_data)

            p_vec = rep(NA, n_tax)
            beta_vec = rep(NA, n_tax)

            idx2 = NULL
            if (is.null(rand_formula)) {
                for (idx2 in seq_len(n_tax)) {
                    test_data = data.frame(x = alr_data[, idx2],
                                           meta_data,
                                           check.names = FALSE)
                    lm_fit = suppressWarnings(tfun(tformula, data = test_data))
                    p_vec[idx2] = summary(lm_fit)$coef[2, "Pr(>|t|)"]
                    beta_vec[idx2] = summary(lm_fit)$coef[2, "t value"]
                }
            } else {
                for (idx2 in seq_len(n_tax)) {
                    test_data = data.frame(x = alr_data[, idx2],
                                           meta_data,
                                           check.names = FALSE)
                    lme_fit = try(tfun(formula = tformula,
                                       data = test_data,
                                       na.action = na.omit,
                                       control = lme_control),
                                  silent = TRUE)
                    if (inherits(lme_fit, "try-error")) {
                        p_vec[idx2] = NA
                        beta_vec[idx2] = NA
                    } else {
                        summary_fit = summary(lme_fit)
                        p_vec[idx2] = summary_fit$coefficients[main_var, "Pr(>|t|)"]
                        beta_vec[idx2] = summary_fit$coefficients[main_var, "Estimate"]
                    }
                }
            }

            list(p_vec, beta_vec)
        }
        stopCluster(cl)
    } else {
        cl = makeCluster(n_cl)
        registerDoParallel(cl)

        result = foreach(idx1 = seq_len(n_tax), .combine = comb, .multicombine = TRUE) %dopar% {
            alr_data = apply(comp_table, 1, function(x) x - comp_table[idx1, ])
            alr_data = cbind(alr_data, meta_data)

            p_vec = rep(NA, n_tax)
            beta_vec = rep(NA, n_tax)

            idx2 = NULL
            if (is.null(rand_formula)) {
                for (idx2 in seq_len(n_tax)) {
                    test_data = data.frame(x = alr_data[, idx2],
                                           meta_data,
                                           check.names = FALSE)
                    lm_fit = suppressWarnings(tfun(tformula, data = test_data))
                    p_vec[idx2] = anova(lm_fit)[main_var, "Pr(>F)"]
                    beta_vec[idx2] = anova(lm_fit)[main_var, "F value"]
                }
            } else {
                for (idx2 in seq_len(n_tax)) {
                    test_data = data.frame(x = alr_data[, idx2],
                                           meta_data,
                                           check.names = FALSE)
                    lme_fit = try(tfun(formula = tformula,
                                       data = test_data,
                                       na.action = na.omit,
                                       control = lme_control),
                                  silent = TRUE)
                    if (inherits(lme_fit, "try-error")) {
                        p_vec[idx2] = NA
                        beta_vec[idx2] = NA
                    } else {
                        p_vec[idx2] = anova(lme_fit)[main_var, "Pr(>F)"]
                        beta_vec[idx2] = anova(lme_fit)[main_var, "F value"]
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
    res_comp = data.frame(taxon = taxon_id, W, row.names = NULL, check.names = FALSE)
    res_comp = res_comp %>%
        mutate(detected_0.9 = ifelse(.data$W > 0.9 * (n_tax - 1), TRUE, FALSE),
               detected_0.8 = ifelse(.data$W > 0.8 * (n_tax - 1), TRUE, FALSE),
               detected_0.7 = ifelse(.data$W > 0.7 * (n_tax - 1), TRUE, FALSE),
               detected_0.6 = ifelse(.data$W > 0.6 * (n_tax - 1), TRUE, FALSE))

    # Combine the information of structural zeros
    if (struc_zero){
        res = data.frame(taxon = tax_name, W = Inf,
                         detected_0.9 = TRUE, detected_0.8 = TRUE,
                         detected_0.7 = TRUE, detected_0.6 = TRUE,
                         row.names = NULL, check.names = FALSE)
        res[comp_idx, ] = res_comp
    }else{
        res = res_comp
    }

    out = list(res = res, zero_ind = zero_ind,
               beta_data = beta_data, p_data = p_data, q_data = q_data)
    return(out)
}






