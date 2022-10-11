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


#' @title Diet Swap Data
#' @description The diet swap dataset represents a study with African
#' and African American groups undergoing a two-week diet swap.
#' For details, see
#' \url{https://doi.org/10.1038/ncomms7342}.
#' @name dietswap
#' @return Loads the dataset in R.
#' @details The dataset is also available via the microbiome R package
#' \url{http://microbiome.github.com/microbiome} in \code{phyloseq} format.
#' @docType data
#' @author Huang Lin \email{huanglinfrederick@@gmail.com}
#' @references
#' O'Keefe et al. Nature Communications 6:6342, 2015.
#' \url{https://doi.org/10.1038/ncomms7342}
#' @usage data(dietswap)
#' @format The dataset in \code{TreeSummarizedExperiment} format.
#' @keywords data
NULL



#' @title HITChip Atlas with 1006 Western Adults
#' @description This dataset contains genus-level microbiota profiling with
#' HITChip for 1006 western adults with no reported health complications,
#' reported in Lahti et al. (2014) \url{https://doi.org/10.1038/ncomms5344}.
#' @name atlas1006
#' @details The dataset is also available via the microbiome R package
#' \url{http://microbiome.github.com/microbiome} in \code{phyloseq} format.
#' @docType data
#' @author Huang Lin \email{huanglinfrederick@@gmail.com}
#' @return Loads the dataset in R.
#' @references
#' Lahti et al. Nature Communications 5:4344, 2014.
#' \url{https://doi.org/10.1038/ncomms5344}
#' @usage data(atlas1006)
#' @format The dataset in \code{TreeSummarizedExperiment} format.
#' @keywords data
NULL
