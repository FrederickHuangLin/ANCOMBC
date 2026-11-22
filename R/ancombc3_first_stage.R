# =============================================================================
# ANCOM-BC3 first stage: nuisance estimation with the ZIPLN-VAE
# =============================================================================
#
# ancombc3_fitted_vae_offsets() carries out the five steps of the first stage:
#   1. fit the working ZIPLN-VAE (zipln_vae.R);
#   2. decode the candidate offsets with the tested covariates at their
#      reference values, the sample means;
#   3. abundance offset: the double-centered log decoded at-risk mean,
#      projected onto the orthogonal complement of the residualized tested
#      covariates and double-centered again;
#   4. structural offset: the decoded output of the structural network,
#      projected in the same way, without centering;
#   5. decode at the observed covariates: at-risk probabilities,
#      sampling-zero probabilities, and posterior probabilities of structural
#      absence.
# The outputs are inputs of the second-stage DA fit (fit_ancombc3_da()), the
# DA eligibility (ancombc3_da_eligibility()), the DP eligibility
# (dp_eligibility()), and the DP test (fit_ancombc3_dp_zero_pattern()).

# Double centering of a matrix, dc(M)_ij = M_ij - Mbar_i. - Mbar_.j + Mbar_..:
# row means and column means are removed and the grand mean is added back.
.ac3_double_center <- function(mat) {
  mat <- as.matrix(mat)
  row_center <- rowMeans(mat, na.rm = TRUE)
  col_center <- colMeans(mat, na.rm = TRUE)
  grand <- mean(mat, na.rm = TRUE)
  out <- sweep(mat, 1L, row_center, "-")
  out <- sweep(out, 2L, col_center, "-")
  out <- out + grand
  dimnames(out) <- dimnames(mat)
  out
}

# Projection of each column of an n x p offset matrix onto the orthogonal
# complement of the residualized tested covariates, (I - Q Q') offset, with Q
# from .ac3_da_residualized_tested_basis(), followed by double centering when
# double_center is TRUE.
ancombc3_project_offset <- function(offset,
                                    x,
                                    test_terms = colnames(x),
                                    double_center = TRUE) {
  offset <- as.matrix(offset)
  x <- as.matrix(x)
  if (nrow(offset) != nrow(x)) {
    stop("offset and x must have the same number of samples.")
  }
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  test_terms <- intersect(test_terms, colnames(x))
  if (length(test_terms) == 0L) return(offset)
  Q <- .ac3_da_residualized_tested_basis(x, test_terms)
  if (ncol(Q) == 0L) return(offset)
  out <- offset - Q %*% crossprod(Q, offset)
  dimnames(out) <- dimnames(offset)
  if (isTRUE(double_center)) out <- .ac3_double_center(out)
  out
}

## First-stage nuisance quantities from one ZIPLN-VAE fit to the n x p count
## matrix y with the n x q_x abundance design x and the n x q_r prevalence
## design r (ancombc3_make_design(), without intercept columns). The library
## size of sample i is the observed total N_i = sum_j y_ij, used in training
## and in both decodings. The VAE is trained at the observed covariates, and
## its structural network excludes the tested columns of r.
## Returns
##   abundance_offset: n x p abundance offset O^A (steps 2 and 3);
##   structural_offset: n x p structural offset O^S (steps 2 and 4);
##   sigma_vae: length-p log-scale standard deviations sigma_j^VAE (step 1);
## and, from the decoding at the observed covariates (step 5),
##   completed_zero: n x p completed zeros v_ij mu_ij;
##   sampling_zero_posterior: n x p posterior probabilities v_ij that a zero
##     is a sampling zero, for every cell;
##   at_risk_prob: n x p at-risk probabilities pi_ij^full;
##   structural_absence_prob: 1 - v_ij, the posterior probability of
##     structural absence w_ij, for every cell;
##   structural_absence_weight: s_ij^el = 1{y_ij = 0} w_ij, the input of
##     dp_eligibility();
##   sampling_zero_prob: n x p sampling-zero probabilities f0_ij^full.
## Every n x p matrix carries the dimnames of y (samples, taxa), and sigma_vae
## carries its column names, so that the downstream functions match the taxa
## by name; when y has no dimnames, the outputs have none and are matched by
## position (taxa T1, ..., Tp). The decoder computes in single precision.
## structural_absence_prob and structural_absence_weight are computed from
## sampling_zero_posterior in double precision. With repeated measures,
## cluster is the subject of each sample, and the samples withheld in the
## training of the VAE are those of whole subjects
## (.ac3_cluster_validation_split()). progress is an optional per-epoch
## function passed to train_zipln_vae(). The result also holds
## training_history, the per-epoch losses of train_zipln_vae(), and
## n_validation, the number of withheld samples.
ancombc3_fitted_vae_offsets <- function(y,
                                        x,
                                        r,
                                        tested_terms = colnames(x),
                                        epochs = 100L,
                                        device = "cpu",
                                        seed = 1L,
                                        cluster = NULL,
                                        progress = NULL) {
  # Lower bound of the decoded at-risk mean before the log transform.
  abundance_floor <- 1e-8
  tested_terms <- intersect(tested_terms, colnames(x))
  if (length(tested_terms) == 0L) stop("No tested_terms are present in x.")
  y <- as.matrix(y)
  x <- as.matrix(x)
  r <- as.matrix(r)
  if (is.null(colnames(x))) colnames(x) <- paste0("x", seq_len(ncol(x)))
  if (is.null(colnames(r))) colnames(r) <- paste0("r", seq_len(ncol(r)))
  x_reference <- reference_test_columns(x, tested_terms)
  r_reference <- reference_test_columns(r, tested_terms)
  library_size <- as.numeric(rowSums(y))
  fit <- train_zipln_vae(
    y = y, x = x, r = r, N = library_size,
    k_abundance = min(6, max(2, ncol(x) + 3)),
    k_structural = min(4, max(2, ncol(r) + 2)),
    hidden = 96, offset_hidden = 64,
    epochs = as.integer(epochs), batch_size = min(64, nrow(y)),
    lr = 1e-3, weight_decay = 5e-3,
    kl_warmup_epochs = 40L,
    val_frac = 0.1,
    positive_ll_weight = 0.5,
    structural_init_fraction = 0.5,
    structural_offset_exclude_r_columns = tested_terms,
    seed = seed, device = device,
    cluster = cluster, progress = progress
  )
  reference_completion <- zipln_vae_completion(
    fit$model, y, x, r, library_size,
    x_decode = x_reference,
    r_decode = r_reference,
    device = device
  )
  observed_completion <- zipln_vae_completion(
    fit$model, y, x, r, library_size,
    device = device
  )
  abundance_offset <- .ac3_double_center(log(pmax(reference_completion$mu, abundance_floor)))
  abundance_offset <- ancombc3_project_offset(
    abundance_offset,
    x = x,
    test_terms = tested_terms,
    double_center = TRUE
  )
  ## Structural offset: the output of the structural network decoded with the
  ## tested covariates at their reference values, projected onto the
  ## orthogonal complement of the residualized tested covariates.
  structural_offset <- ancombc3_project_offset(
    as.matrix(reference_completion$h_structural),
    x = x,
    test_terms = tested_terms,
    double_center = FALSE
  )
  dimnames(abundance_offset) <- dimnames(y)
  dimnames(structural_offset) <- dimnames(y)
  ## Poisson-log-normal sampling-zero probability f(0; mu, sigma) of the
  ## decoding at the observed covariates.
  sampling_zero_prob <- observed_completion$f0
  dimnames(sampling_zero_prob) <- dimnames(y)
  ## Posterior probability of structural absence, w_ij = 1 - v_ij, and the
  ## structural-absence weight s_ij^el, zero for a positive count.
  sampling_zero_posterior <- observed_completion$v
  dimnames(sampling_zero_posterior) <- dimnames(y)
  structural_absence_prob <- 1 - sampling_zero_posterior
  structural_absence_weight <- ifelse(y == 0, structural_absence_prob, 0)
  dimnames(structural_absence_weight) <- dimnames(y)
  completed_zero <- observed_completion$completed_zero
  dimnames(completed_zero) <- dimnames(y)
  at_risk_prob <- observed_completion$pi
  dimnames(at_risk_prob) <- dimnames(y)
  sigma_vae <- as.numeric(observed_completion$sigma)
  names(sigma_vae) <- colnames(y)
  list(
    abundance_offset = abundance_offset,
    structural_offset = structural_offset,
    sigma_vae = sigma_vae,
    completed_zero = completed_zero,
    sampling_zero_posterior = sampling_zero_posterior,
    at_risk_prob = at_risk_prob,
    structural_absence_prob = structural_absence_prob,
    structural_absence_weight = structural_absence_weight,
    sampling_zero_prob = sampling_zero_prob,
    training_history = fit$history,
    n_validation = fit$n_validation
  )
}
