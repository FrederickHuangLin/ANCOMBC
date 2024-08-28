#' @title Data Sanity and Integrity Check
#'
#' @description Determine if the input data is in a correct format
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
#' variable in the formula if it is not NULL.
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
#' @param verbose logical. Whether to display detailed progress messages.
#'
#' @return a \code{list} containing the outputs formatted appropriately for
#' downstream analysis.
#'
#' @examples
#' data(atlas1006, package = "microbiome")
#' check_results = data_sanity_check(data = atlas1006,
#'                                   tax_level = "Family",
#'                                   fix_formula = "age + sex + bmi_group",
#'                                   group = "bmi_group",
#'                                   struc_zero = TRUE,
#'                                   global = TRUE,
#'                                   verbose = TRUE)
#'
#' @author Huang Lin
#'
#' @importFrom microbiome abundances meta aggregate_taxa
#' @importFrom SummarizedExperiment assay colData
#'
#' @export
data_sanity_check = function(data, taxa_are_rows = TRUE,
                             assay.type = assay_name, assay_name = "counts",
                             rank = tax_level, tax_level = NULL,
                             aggregate_data = NULL, meta_data = NULL,
                             fix_formula, group = NULL, struc_zero = FALSE,
                             global = FALSE, pairwise = FALSE, dunnet = FALSE,
                             mdfdr_control = list(fwer_ctrl_method = "holm",
                                                  B = 100),
                             trend = FALSE,
                             trend_control = list(contrast = NULL,
                                                  node = NULL,
                                                  solver = "ECOS",
                                                  B = 100),
                             verbose = TRUE) {
    #=========== Check for aliases ===========
    if (!is.null(assay.type)) {
        assay_name = assay.type
    }

    if (!is.null(rank)) {
        tax_level = rank
    }

    #=========== Check input data type ===========
    if (verbose) {
        message("Checking the input data type ...")
        message(paste("The input data is of type:", class(data)[1]))
    }

    if (inherits(data, "phyloseq")) {
        if (!requireNamespace("microbiome", quietly = TRUE)){
            stop(paste("The 'microbiome' package is needed to process the imported data but is not installed.",
                       "Please install the package to continue."))
            }
        # Process the phyloseq object
        feature_table = microbiome::abundances(data)
        meta_data = microbiome::meta(data)
        if (!is.null(tax_level)) {
            aggregate_data = microbiome::aggregate_taxa(data, tax_level)
            feature_table_aggregate = microbiome::abundances(aggregate_data)
        } else {
            feature_table_aggregate = feature_table
        }
        }
    else if (inherits(data, "TreeSummarizedExperiment")) {
        if (!requireNamespace("mia", quietly = TRUE)) {
            stop(paste("The 'mia' package is needed to process the imported data but is not installed.",
                       "Please install the package to continue."))
            }
        # Process the tse object
        feature_table = SummarizedExperiment::assay(data, assay_name)
        meta_data = as.data.frame(SummarizedExperiment::colData(data))
        if (!is.null(tax_level)) {
            aggregate_data = .merge_features(data, tax_level)
            feature_table_aggregate = SummarizedExperiment::assay(aggregate_data, assay_name)
        } else {
            feature_table_aggregate = feature_table
            }
        }
    else if (inherits(data, c("data.frame", "matrix"))) {
        message("The imported data is in a generic 'matrix'/'data.frame' format.")

        if (is.null(aggregate_data)) {
            aggregate_data = data
        }

        if (is.null(meta_data)) {
            stop("Missing sample metadata. Please provide the sample metadata in 'data.frame' format.")
        }
        if (!taxa_are_rows) {
            feature_table = as.matrix(t(data))
            feature_table_aggregate = as.matrix(t(aggregate_data))
        } else {
            feature_table = as.matrix(data)
            feature_table_aggregate = as.matrix(aggregate_data)
        }

        non_numeric_columns = vapply(feature_table, function(col) !is.numeric(col), logical(1))

        if (any(non_numeric_columns)) {
            stop("All columns must be numeric. Non-numeric columns found: ",
                 paste(colnames(feature_table)[non_numeric_columns], collapse = ", "))
        }

        if (!all(colnames(feature_table) %in% rownames(meta_data))) {
            stop(paste("Sample names do not match between the feature table and sample metadata.",
                       "Please ensure the column names of the feature abundance matrix and the row names of the metadata data frame are consistent.",
                       sep = "\n"))
        }
        meta_data = meta_data[colnames(feature_table), ]
    }
    else {
        stop("Unsupported data type. Please provide a 'matrix', 'data.frame', 'phyloseq', or 'tse' object.")
    }

    if (verbose) {
        message("PASS")
    }

    #=========== Check sample metadata ===========
    if (verbose) {
        message("Checking the sample metadata ...")
    }

    # Drop unused levels
    meta_data[] = lapply(meta_data, function(x)
        if(is.factor(x)) factor(x) else x)

    # Check if all covariates specified in the fix_formula are columns in metadata
    fix_formula_non_interact = gsub("\\*", "+", fix_formula)
    vars = unlist(strsplit(fix_formula_non_interact, split = "\\s*\\+\\s*")) # \\s* matches zero or more spaces

    if (verbose) {
        message("The specified variables in the formula: ",
                paste(vars, collapse = ", "))
        message("The available variables in the sample metadata: ",
                paste(colnames(meta_data), collapse = ", "))
    }

    missing_vars = vars[!vars %in% colnames(meta_data)]
    if (length(missing_vars) > 0) {
        stop("The following variables specified are not in the meta data: ",
             paste(missing_vars, collapse = ", "))
    }

    if (verbose) {
        message("PASS")
    }

    #=========== Check other arguments ===========
    if (verbose) {
        message("Checking other arguments ...")
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
    } else if (! is.numeric(meta_data[, group])) {
        meta_data[, group] = as.factor(meta_data[, group])
        # Check the number of groups
        n_level = nlevels(meta_data[, group])

        if (verbose) {
            message(paste("The number of groups of interest is:", n_level))
        }

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

        if (verbose) {
            message_content = paste(names(size_per_group), size_per_group, sep = " = ", collapse = ", ")
            message("The sample size per group is: ", message_content)
        }

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
    } else {
        warn_txt = sprintf(paste("The group variable contains numerical values",
                                 "Are you sure this is the correct configuration? ",
                                 sep = "\n"))
        warning(warn_txt, call. = FALSE)
        cat("\n")
    }

    if (verbose) {
        message("PASS")
    }

    #=========== output ===========
    output = list(feature_table = feature_table,
                  feature_table_aggregate = feature_table_aggregate,
                  meta_data = meta_data,
                  global = global,
                  pairwise = pairwise,
                  dunnet = dunnet,
                  trend = trend,
                  trend_control = trend_control)
    return(output)
  }



