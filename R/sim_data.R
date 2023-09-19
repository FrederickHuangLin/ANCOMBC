#' @title Simulate Microbial Absolute Abundance Data by Poisson lognormal
#' (PLN) model Based on a Real Dataset
#'
#' @description Generate microbial absolute abundances using the Poisson
#' lognormal (PLN) model based on the mechanism described in the
#' \href{https://doi.org/10.1093/bioinformatics/btaa260}{LDM}
#' paper (supplementary text S2).
#'
#' @details The PLN model relates the abundance vector with a Gaussian latent
#' vector. Because of the presence of a latent layer, the PLN model displays a
#' larger variance than the Poisson model (over-dispersion). Also, the
#' covariance (correlation) between abundances has the same sign as the
#' covariance (correlation) between the corresponding latent variables.
#' This property gives enormous flexibility in modeling the variance-covariance
#' structure of microbial abundances since it is easy to specify different
#' variance-covariance matrices in the multivariate Gaussian distribution.
#'
#' However, instead of manually specifying the variance-covariance matrix, we
#' choose to estimate the variance-covariance matrix from a real dataset,
#' which will make the simulated data more resemble real data.
#'
#' @param abn_table the input microbial count table. It is used to obtain
#' the estimated variance-covariance matrix, can be in either \code{matrix}
#' or \code{data.frame} format.
#' @param taxa_are_rows logical. TRUE if the input dataset has rows
#' represent taxa. Default is TRUE.
#' @param prv_cut a numerical fraction between 0 and 1. Taxa with prevalences
#' less than \code{prv_cut} will be excluded in the analysis. For instance,
#' suppose there are 100 samples, if a taxon has nonzero counts presented in
#' less than 10 samples, it will not be further analyzed. Default is 0.10.
#' @param n numeric. The desired sample size for the simulated data.
#' @param lib_mean numeric. Mean of the library size. Library sizes are
#' generated from the negative binomial distribution with parameters
#' \code{lib_mean} and \code{disp}. For details, see \code{?rnbinom}.
#' @param disp numeric. The dispersion parameter for the library size.
#' For details, see \code{?rnbinom}.
#' @return a \code{matrix} of microbial absolute abundances, where taxa are in
#' rows and samples are in columns.
#'
#' @examples
#' library(ANCOMBC)
#' data(QMP)
#' abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05,
#'                     n = 100, lib_mean = 1e8, disp = 0.5)
#' rownames(abn_data) = paste0("Taxon", seq_len(nrow(abn_data)))
#' colnames(abn_data) = paste0("Sample", seq_len(ncol(abn_data)))
#'
#' @author Huang Lin
#'
#' @references
#' \insertRef{hu2020testing}{ANCOMBC}
#'
#' @rawNamespace import(stats, except = filter)
#' @importFrom Matrix nearPD
#'
#' @export
sim_plnm = function(abn_table, taxa_are_rows = TRUE,
                    prv_cut = 0.1, n, lib_mean, disp) {
    abn_table = as.matrix(abn_table)
    if (!taxa_are_rows) {
        abn_table = t(abn_table)
    }
    prevalence = apply(abn_table, 1, function(x)
        sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
    tax_keep = which(prevalence >= prv_cut)
    txt = paste0("The number of taxa after filtering is: ", length(tax_keep))
    message(txt)

    if (length(tax_keep) > 0) {
        abn_table = abn_table[tax_keep, , drop = FALSE]
    } else {
        stop("No taxa remain under the current cutoff", call. = FALSE)
    }
    rel_table = t(t(abn_table)/colSums(abn_table))

    n_taxa = nrow(rel_table)
    mean_rel = rowMeans(rel_table)
    cov_rel = cov(t(rel_table))
    overdisp = cov_rel - diag(mean_rel) + outer(mean_rel, mean_rel)
    exp_cov = matrix(1, nrow = n_taxa, ncol = n_taxa) +
        diag(1/mean_rel) %*% overdisp %*% diag(1/mean_rel)

    ev = eigen(exp_cov)
    Lambda = ev$values
    if (all(Lambda <= 0)) {
        top_txt = paste0("All eigenvalues are nonpositive \n",
                         "Please use a different dataset")
        stop(top_txt, call. = FALSE)
    }
    Q = ev$vectors
    pos_idx = which(Lambda > 0)
    if (length(pos_idx) == 1) {
        log_Lambda = log(Lambda[pos_idx])
    } else {
        log_Lambda = diag(log(Lambda[pos_idx]))
    }
    cov_est = Q[, pos_idx, drop = FALSE] %*% log_Lambda %*%
        t(Q[, pos_idx, drop = FALSE])

    if (!all(eigen(cov_est)$values > 0)) {
      cov_est = as.matrix(Matrix::nearPD(cov_est)$mat)
    }

    mean_est = log(mean_rel) - 0.5 * diag(cov_est)

    N = rnbinom(n = n, mu = lib_mean, size = disp)
    sim_data = .rplnm(mu = mean_est, sigma = sqrt(cov_est),
                      n = n, N = N)
    return(t(sim_data))
}



































