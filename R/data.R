#' @title Quantitative Microbiome Project data
#' @description The data containing quantitative microbiome count data of
#' dimension 106 samples/subjects (in rows) and 91 OTUs (in columns).
#' The raw dataset is pruned the taxa present less than 30% of samples and
#' final dataset contains only healthy subjects from two cohorts:
#' Study cohort and Disease cohort.
#' For details, see \url{https://doi.org/10.1038/nature24460}.
#' @name QMP
#' @return Loads the dataset in R.
#' @details The dataset is also available via the SPRING R package
#' \url{https://github.com/GraceYoon/SPRING} in \code{matrix} format.
#' @docType data
#' @author Huang Lin \email{huanglinfrederick@@gmail.com}
#' @references
#' Vanderputte et al. Nature. 551: 507-511, 2017.
#' \url{https://doi.org/10.1038/nature24460}
#' @usage data(QMP)
#' @format The dataset in \code{matrix} format.
#' @keywords data
NULL

