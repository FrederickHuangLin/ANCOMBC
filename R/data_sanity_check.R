#' @title data sanity check
#'
#' @description Determine if the input data is in a correct format
#' The current version of \code{data_sanity_check} function implements data input check for Analysis of Compositions of Microbiomes
#' with Bias Correction (ANCOM-BC2) in cross-sectional and repeated measurements
#' data. 
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
#' @param data the input data. The \code{data} parameter should be either a
#' \code{phyloseq} or a \code{TreeSummarizedExperiment} object, which
#' consists of a feature table (microbial count table), a sample metadata table,
#' a taxonomy table (optional), and a phylogenetic tree (optional).
#' Ensure that the row names of the metadata table match the sample names in the
#' feature table, and the row names of the taxonomy table match the taxon
#' (feature) names in the feature table. For detailed information, refer to
#' \code{?phyloseq::phyloseq} or
#' \code{?TreeSummarizedExperiment::TreeSummarizedExperiment}.
#' It is recommended to use low taxonomic levels, such as OTU or species level,
#' as the estimation of sampling fractions requires a large number of taxa.
#' 
#' @param formula the character string expresses how the microbial absolute
#' abundances for each taxon depend on the fixed effects in metadata. When
#' specifying the \code{fix_formula}, make sure to include the \code{group}
#' variable in the formula if it is not NULL.\code{rand_formula} the character string expresses how the microbial absolute
#' abundances for each taxon depend on the random effects in metadata. ANCOM-BC2
#' follows the \code{lmerTest} package in formulating the random effects. See
#' \code{?lmerTest::lmer} for more details. Default is \code{NULL}.
#' @param group character. the name of the group variable in metadata.
#' The \code{group} parameter should be a character string representing the name
#' of the group variable in the metadata. The \code{group} variable should be
#' discrete, meaning it consists of categorical values. Specifying the
#' \code{group} variable is required if you are interested in detecting
#' structural zeros and performing performing multi-group comparisons (global
#' test, pairwise directional test, Dunnett's type of test, and trend test).
#' However, if these analyses are not of interest to you, you can leave the
#' \code{group} parameter as NULL. If the \code{group} variable of interest
#' contains only two categories, you can also leave the \code{group} parameter
#' as NULL. Default is NULL.
#' @param struc_zero logical. Whether to detect structural zeros based on
#' \code{group}. Default is FALSE. See \code{Details} for
#' a more comprehensive discussion on structural zeros.
#' @param global logical. Whether to perform the global test. Default is FALSE.
#' @param pairwise logical. Whether to perform the pairwise directional test.
#' Default is FALSE.
#' @param dunnet logical. Whether to perform the Dunnett's type of test.
#' Default is FALSE.
#' @param mdfdr_control a named list of control parameters for mixed directional
#' false discover rate (mdFDR), including 1) \code{fwer_ctrl_method}: family
#' wise error (FWER) controlling procedure, such as "holm", "hochberg",
#' "bonferroni", etc (default is "holm") and 2) \code{B}: the number of
#' bootstrap samples (default is 100). Increase \code{B} will lead to a more
#' accurate p-values. See \code{Details} for a more comprehensive discussion on
#' mdFDR.
#' @param trend logical. Whether to perform trend test. Default is FALSE.
#' @param trend_control a named list of control parameters for the trend test,
#' including 1) \code{contrast}: the list of contrast matrices for
#' constructing inequalities, 2) \code{node}: the list of positions for the
#' nodal parameter, 3) \code{solver}: a string indicating the solver to use
#' (default is "ECOS"), and 4) \code{B}: the number of bootstrap samples
#' (default is 100). Increase \code{B} will lead to a more accurate p-values.
#' See \code{vignette} for the corresponding trend test examples.
#'
#' @return a \code{Data} in correct format

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
#' # Can also contain non-taxonomic information, for instance
#' # colnames(tax_tab) = c("G1", "G2", "G3", "G4", "G5", "G6", "G7")
#' tax_tab = DataFrame(tax_tab)
#'
#' # create TSE
#' tse = TreeSummarizedExperiment(assays = assays,
#'                                colData = smd,
#'                                rowData = tax_tab)
#'
#'dat<-data_sanity_check(tse)
#'class(dat)
#'
#'
.data_sanity_check<-function(data){
   if(inherits(data,"phyloseq")){
     if(!requireNamespace("phyloseq", quietly = TRUE)){
       stop("The 'phyloseq' package is required but is not installed. Please install it to use this feature.")
     }
      # Process the phyloseq object
      data<-phyloseq::otu_table(data)
   }else if(inherits(data,"TreeSummarizedExperiment")){
  if(!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)){
    stop("The 'mia' package is required but is not installed. Please install it to use this feature.")
  }
   # Process the tse object
  data 
 } else if(is.data.frame(data)) {
    # Process the data.frame
      data <- data
  }else {stop("Unsupported data type. Please provide a data.frame, phyloseq, or tse object.")}
return(data)
  }

# Filter data by prevalence and library size
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
    
    

