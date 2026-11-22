# =============================================================================
# ANCOM-BC3 ZIPLN-VAE (first stage, step 1)
# =============================================================================
# Working zero-inflated Poisson-log-normal variational autoencoder of the first
# stage (ancombc3_fitted_vae_offsets()). For sample i and taxon j:
#   latent scores z_i = (z_i^A, z_i^S), with the variational distribution
#     q_phi(z_i | y_i); the encoder input is the library-size-normalized log
#     counts, and the covariates are inputs of the decoder only;
#   abundance part: eta_ij = alpha_j^VAE + h_j^A(z_i^A, x_i),
#     rho_ij = exp(eta_ij) / sum_k exp(eta_ik), and the at-risk mean
#     mu_ij = N_i rho_ij;
#   structural part: m_ij^S = a_j^VAE + r_i' gamma_j^VAE + h_j^S(z_i^S, z_i^A, r_i,N)
#     and the structural-presence probability pi_ij = expit(m_ij^S), where the
#     input r_i,N of h^S holds the prevalence covariates with each tested
#     column replaced by its sample mean;
#   working distribution of y_ij: (1 - pi_ij) 1{y = 0} + pi_ij f(y; mu_ij, sigma_j),
#     with f the Poisson-log-normal mass function whose log intensity is normal
#     with mean c_ij = log mu_ij - sigma_j^2 / 2 and standard deviation sigma_j;
#   sampling-zero posterior v_ij = Pr(S_ij = 1 | Y_ij = 0)
#     = pi_ij f0_ij / {(1 - pi_ij) + pi_ij f0_ij}, with f0_ij = f(0; mu_ij, sigma_j);
#   completed zero v_ij mu_ij.
# Code names: eta (eta^dec), rho, mu, presence_logit (m^S,dec), h_structural
# (the output of h^S), linear_presence (a^VAE + r' gamma^VAE), pi (pi^dec),
# sigma, center (c), v, and completed_zero; the parameters abundance_intercept
# (alpha^VAE), presence_intercept (a^VAE), gamma (gamma^VAE), and log_sigma
# (log sigma^VAE); the networks abundance_network with abundance_head (h^A)
# and structural_network (h^S).
#
# Dimension convention:
#   y: n x p observed count matrix, samples by taxa.
#   x: n x q_x abundance covariates; r: n x q_r prevalence covariates.
#   N: length-n library sizes, internally reshaped to n x 1.
#   z_abundance: n x k_abundance; z_structural: n x k_structural.
#   eta, rho, mu, presence_logit, h_structural, pi, v, completed_zero, f0: n x p
#   decoder outputs.
#   gamma: q_r x p linear prevalence coefficients (positive = higher prevalence).
#   sigma: length-p taxon log-scale standard deviations.
# =============================================================================

as_numeric_matrix <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.")
  }
  x
}

# Gauss-Hermite rule for weight exp(-t^2) (Golub-Welsch): nodes t_k and the log
# of w_k exp(t_k^2), so that sum_k exp(lw_k) g(t_k) approximates the integral of g.
.zipln_gh_rule <- function(n_nodes) {
  i <- seq_len(n_nodes - 1L)
  J <- matrix(0, n_nodes, n_nodes)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  nodes <- e$values
  w <- sqrt(pi) * e$vectors[1L, ]^2
  ord <- order(nodes)
  list(t = nodes[ord], lw = log(w[ord]) + nodes[ord]^2)
}

.zipln_gh <- .zipln_gh_rule(20L)

pln_log_prob <- function(y, center, sigma, newton_steps = 30L) {
  # Elementwise Poisson-log-normal log mass function log f(y; mu, sigma), where
  # f(y; mu, sigma) = int Poisson(y; e^u) N(u; c, sigma^2) du and
  # c = log mu - sigma^2 / 2 (center) is the normal mean of the log intensity.
  # y, center, and sigma are tensors of the same shape, usually n x p.
  # Adaptive Gauss-Hermite quadrature: the nodes are centered at the mode u* of
  # the log integrand L(u) = y u - e^u - (u - c)^2 / (2 sigma^2) and scaled by
  # its curvature, s = (e^{u*} + sigma^{-2})^{-1/2}. The center and scale of
  # the nodes are computed without gradients; the gradient with respect to
  # (c, sigma) is that of the quadrature sum through the normal density at the
  # nodes. Newton's method for L'(u) = 0 converges monotonically from
  # u0 = max(c, log y), an upper bound of u*, because L' is decreasing and
  # concave; it stops after newton_steps steps, or earlier, after at least
  # four steps, when the largest step is below 1e-14 max(1, |u|). The Poisson
  # term y u - e^u - log y! of the log integrand does not
  # depend on (c, sigma); its terms are of size y log y, so it is evaluated at
  # the nodes in double precision and then converted to the precision of y.
  torch::with_no_grad({
    yd <- y$to(dtype = torch::torch_double())
    cd <- center$detach()$to(dtype = torch::torch_double())
    s2 <- sigma$detach()$to(dtype = torch::torch_double())$pow(2)
    u <- torch::torch_maximum(cd, torch::torch_log(yd$clamp_min(1e-300)))
    for (k in seq_len(newton_steps)) {
      eu <- torch::torch_exp(u)
      step <- (yd - eu - (u - cd) / s2) / (-eu - 1 / s2)
      u <- u - step
      if (k >= 4L && torch::torch_max(torch::torch_abs(step) / torch::torch_abs(u)$clamp_min(1))$item() < 1e-14) break
    }
    scale <- torch::torch_rsqrt(torch::torch_exp(u) + 1 / s2)
    gh_nodes <- torch::torch_tensor(.zipln_gh$t, dtype = torch::torch_double(), device = y$device)
    lw <- torch::torch_tensor(.zipln_gh$lw, dtype = torch::torch_double(), device = y$device)
    nodes_d <- u$unsqueeze(-1L) + sqrt(2) * scale$unsqueeze(-1L) * gh_nodes
    ydk <- yd$unsqueeze(-1L)
    log_poisson <- (ydk * nodes_d - torch::torch_exp(nodes_d) - torch::torch_lgamma(ydk + 1))$to(dtype = y$dtype)
    nodes <- nodes_d$to(dtype = y$dtype)
    log_q <- (lw + torch::torch_log(sqrt(2) * scale)$unsqueeze(-1L))$to(dtype = y$dtype)
  })
  ck <- center$unsqueeze(-1L)
  sk <- sigma$unsqueeze(-1L)
  log_integrand <- log_poisson -
    0.5 * log(2 * pi) - torch::torch_log(sk) - (nodes - ck)$pow(2) / (2 * sk$pow(2))
  torch::torch_logsumexp(log_integrand + log_q, dim = -1L)
}

zipln_log_prob <- function(y, log_f_y, log_f0, presence_logit) {
  # Elementwise ZIPLN log-likelihood from the Poisson-log-normal log mass at y
  # (log_f_y) and at zero (log_f0). presence_logit is the logit of the
  # structural-presence probability pi, so log_presence = log pi and
  # log_absence = log(1 - pi). A zero is a structural zero, with probability
  # 1 - pi, or a sampling zero, with probability pi f0; a positive count arises
  # from the at-risk component only.
  log_presence <- -torch::nnf_softplus(-presence_logit)
  log_absence <- -torch::nnf_softplus(presence_logit)
  log_zero <- torch::torch_logsumexp(
    torch::torch_stack(list(log_absence, log_presence + log_f0), dim = 1L),
    dim = 1L
  )
  log_positive <- log_presence + log_f_y
  torch::torch_where(y$eq(0), log_zero, log_positive)
}

.zipln_sampling_zero_posterior <- function(log_f0, presence_logit) {
  # Sampling-zero posterior v_ij = Pr(S_ij = 1 | Y_ij = 0, z_i), the posterior
  # probability that an observed zero is a sampling zero of an at-risk taxon.
  # It is computed for every n x p cell; only the zero cells are used.
  log_presence <- -torch::nnf_softplus(-presence_logit)
  log_absence <- -torch::nnf_softplus(presence_logit)
  denom <- torch::torch_logsumexp(
    torch::torch_stack(list(log_absence, log_presence + log_f0), dim = 1L),
    dim = 1L
  )
  torch::torch_exp(log_presence + log_f0 - denom)
}

# Bounds of the taxon log-scale standard deviation in the decoder.
.zipln_sigma_bounds <- c(0.02, 4)

# self is the module instance inside the methods of the module generator.
utils::globalVariables("self")

# The module generator is created at the first call and kept, because torch is
# not available when the namespace of the package is loaded.
.zipln_vae_state <- new.env(parent = emptyenv())
.zipln_vae_module <- function() {
  if (is.null(.zipln_vae_state$module)) .zipln_vae_state$module <- .zipln_vae_generator()
  .zipln_vae_state$module
}

.zipln_vae_generator <- function() torch::nn_module(
  "ZIPLNVAE",

  initialize = function(p, q_x, q_r, k_abundance = 8, k_structural = 4,
                        hidden = 128, offset_hidden = 64,
                        abundance_intercept_init = NULL,
                        presence_intercept_init = NULL,
                        offset_r_mask = NULL,
                        offset_r_fill = NULL,
                        encoder_count_scale = 1e4) {
    self$p <- as.integer(p)
    self$q_x <- as.integer(q_x)
    self$q_r <- as.integer(q_r)
    self$k_abundance <- as.integer(k_abundance)
    self$k_structural <- as.integer(k_structural)
    self$encoder_count_scale <- encoder_count_scale
    # Mask and fill of the prevalence covariates used as inputs of the
    # structural network: column k is r_k m_k + f_k (1 - m_k), so a column with
    # mask 0 is replaced by its fill value.
    self$offset_r_mask_values <- if (is.null(offset_r_mask)) rep(1, self$q_r) else as.numeric(offset_r_mask)
    self$offset_r_fill_values <- if (is.null(offset_r_fill)) rep(0, self$q_r) else as.numeric(offset_r_fill)

    # The encoder reads the transformed counts y_i only; the covariates are
    # inputs of the decoder only.
    encoder_input <- self$p
    self$encoder <- torch::nn_sequential(
      torch::nn_linear(encoder_input, hidden),
      torch::nn_relu(),
      torch::nn_linear(hidden, hidden),
      torch::nn_relu()
    )
    self$enc_mean_abundance <- torch::nn_linear(hidden, self$k_abundance)
    self$enc_logvar_abundance <- torch::nn_linear(hidden, self$k_abundance)
    self$enc_mean_structural <- torch::nn_linear(hidden, self$k_structural)
    self$enc_logvar_structural <- torch::nn_linear(hidden, self$k_structural)

    # Abundance network h^A: hidden layers on (z^A, x) and a linear output
    # layer without bias.
    self$abundance_network <- torch::nn_sequential(
      torch::nn_linear(self$k_abundance + self$q_x, hidden),
      torch::nn_relu(),
      torch::nn_linear(hidden, hidden),
      torch::nn_relu()
    )
    self$abundance_head <- torch::nn_linear(hidden, self$p, bias = FALSE)
    if (is.null(abundance_intercept_init)) {
      abundance_intercept_init <- torch::torch_zeros(self$p)
    }
    # Taxon-specific intercepts alpha_j^VAE of the abundance logits eta; rho is
    # normalized by the softmax over taxa.
    self$abundance_intercept <- torch::nn_parameter(abundance_intercept_init$clone()$detach())

    # Structural network h^S: an n x p term of the presence logit. The
    # prevalence covariates also contribute linearly through r gamma.
    self$structural_network <- torch::nn_sequential(
      torch::nn_linear(self$k_structural + self$k_abundance + self$q_r, offset_hidden),
      torch::nn_relu(),
      torch::nn_linear(offset_hidden, offset_hidden),
      torch::nn_relu(),
      torch::nn_linear(offset_hidden, self$p)
    )

    if (is.null(presence_intercept_init)) {
      presence_intercept_init <- torch::torch_zeros(self$p)
    }
    self$presence_intercept <- torch::nn_parameter(presence_intercept_init$clone()$detach())
    self$gamma <- torch::nn_parameter(torch::torch_zeros(self$q_r, self$p))
    # log sigma_j, initialized at 0 (sigma_j = 1); the weight decay of the
    # optimizer shrinks it toward 0.
    self$log_sigma <- torch::nn_parameter(torch::torch_zeros(self$p))
  },

  sigma = function() {
    torch::torch_exp(self$log_sigma)$clamp(min = .zipln_sigma_bounds[1], max = .zipln_sigma_bounds[2])
  },

  encode = function(y_log) {
    # y_log: n x p transformed counts. Returns the means and log variances of
    # z^A (n x k_abundance) and z^S (n x k_structural).
    h <- self$encoder(y_log)
    list(
      mean_abundance = self$enc_mean_abundance(h),
      logvar_abundance = self$enc_logvar_abundance(h),
      mean_structural = self$enc_mean_structural(h),
      logvar_structural = self$enc_logvar_structural(h)
    )
  },

  reparameterize = function(mean, logvar) {
    mean + torch::torch_exp(0.5 * logvar) * torch::torch_randn_like(mean)
  },

  encoder_transform = function(y, N) {
    # Library-size-normalized log counts log(1 + scale y / max(N, 1)), the
    # encoder input; the likelihood uses the counts y.
    torch::torch_log1p(y / N$clamp_min(1) * self$encoder_count_scale)
  },

  offset_transform_r = function(r) {
    if (self$q_r == 0L) return(r)
    mask <- torch::torch_tensor(self$offset_r_mask_values, dtype = torch::torch_float(), device = r$device)$unsqueeze(1L)
    fill <- torch::torch_tensor(self$offset_r_fill_values, dtype = torch::torch_float(), device = r$device)$unsqueeze(1L)
    r * mask + fill * (1 - mask)
  },

  decode = function(z_abundance, z_structural, x, r, N) {
    # Decoded sample-by-taxon quantities, each n x p. N is n x 1, so
    # mu = rho N is on the observed count scale. The covariates x and r are
    # those supplied by the caller: the observed covariates, or covariates
    # whose tested columns are set to their reference values (the sample
    # means) for the offsets.
    r_structural <- self$offset_transform_r(r)
    h_abundance <- self$abundance_network(torch::torch_cat(list(z_abundance, x), dim = 2L))
    eta <- self$abundance_head(h_abundance) +
      self$abundance_intercept$unsqueeze(1L)$expand(c(x$size(1), self$p))
    rho <- torch::nnf_softmax(eta, dim = 2L)
    mu <- rho * N

    h_structural <- self$structural_network(torch::torch_cat(list(z_structural, z_abundance, r_structural), dim = 2L))
    linear_presence <- torch::torch_matmul(r, self$gamma) +
      self$presence_intercept$unsqueeze(1L)$expand(c(r$size(1), self$p))
    presence_logit <- linear_presence + h_structural
    pi_dec <- torch::torch_sigmoid(presence_logit)

    # Normal mean c of the log intensity, chosen so that the at-risk mean is mu.
    sigma <- self$sigma()
    sigma_b <- sigma$unsqueeze(1L)$expand(c(r$size(1), self$p))
    center <- torch::torch_log(mu$clamp_min(1e-8)) - sigma_b$pow(2) / 2
    log_f0 <- pln_log_prob(torch::torch_zeros_like(mu), center, sigma_b)
    v <- .zipln_sampling_zero_posterior(log_f0, presence_logit)
    completed_zero <- v * mu

    list(
      eta = eta,
      rho = rho,
      mu = mu,
      h_structural = h_structural,
      linear_presence = linear_presence,
      presence_logit = presence_logit,
      pi = pi_dec,
      sigma = sigma,
      sigma_b = sigma_b,
      center = center,
      log_f0 = log_f0,
      v = v,
      completed_zero = completed_zero
    )
  },

  forward = function(y, x, r, N, kl_weight = 1, positive_ll_weight = 0) {
    y_log <- self$encoder_transform(y, N)
    enc <- self$encode(y_log)
    z_abundance <- self$reparameterize(enc$mean_abundance, enc$logvar_abundance)
    z_structural <- self$reparameterize(enc$mean_structural, enc$logvar_structural)
    dec <- self$decode(z_abundance, z_structural, x, r, N)

    log_f_y <- pln_log_prob(y, dec$center, dec$sigma_b)
    logp <- zipln_log_prob(y, log_f_y, dec$log_f0, dec$presence_logit)
    # ll is the sum over taxa averaged over samples. positive_ll is the
    # Poisson-log-normal log-likelihood averaged over the positive cells, the
    # third term of the training loss.
    ll <- logp$sum(dim = 2L)$mean()
    is_positive <- y$gt(0)
    positive_count <- is_positive$sum()$clamp_min(1)
    positive_ll <- (
      torch::torch_where(is_positive, log_f_y, torch::torch_zeros_like(y))$sum() /
        positive_count
    )

    kl_abundance <- (-0.5 * (
      1 + enc$logvar_abundance -
        enc$mean_abundance$pow(2) -
        torch::torch_exp(enc$logvar_abundance)
    ))$sum(dim = 2L)$mean()

    kl_structural <- (-0.5 * (
      1 + enc$logvar_structural -
        enc$mean_structural$pow(2) -
        torch::torch_exp(enc$logvar_structural)
    ))$sum(dim = 2L)$mean()

    kl <- kl_abundance + kl_structural
    loss <- -ll + kl_weight * kl - positive_ll_weight * positive_ll

    c(
      list(
        loss = loss,
        ll = ll,
        positive_ll = positive_ll,
        kl = kl,
        kl_abundance = kl_abundance,
        kl_structural = kl_structural,
        z_abundance = z_abundance,
        z_structural = z_structural
      ),
      enc,
      dec
    )
  }
)

resolve_matrix_columns <- function(x, columns, name) {
  if (is.null(columns)) return(integer(0L))
  if (is.character(columns)) {
    idx <- match(columns, colnames(x))
  } else {
    idx <- as.integer(columns)
  }
  if (any(is.na(idx)) || any(idx < 1L) || any(idx > ncol(x))) {
    stop(name, " must identify columns of the supplied matrix.")
  }
  unique(idx)
}

prepare_zipln_vae_data <- function(y, x, r, N, device = "cpu") {
  y <- as_numeric_matrix(y, "y")
  x <- as_numeric_matrix(x, "x")
  r <- as_numeric_matrix(r, "r")
  N <- as.numeric(N)

  if (nrow(y) != nrow(x) || nrow(y) != nrow(r) || nrow(y) != length(N)) {
    stop("y, x, r, and N must have the same sample dimension.")
  }
  if (any(y < 0) || any(N <= 0)) {
    stop("y must be nonnegative and N must be positive.")
  }

  list(
    y = y,
    x = x,
    r = r,
    N = N,
    y_t = torch::torch_tensor(y, dtype = torch::torch_float(), device = device),
    x_t = torch::torch_tensor(x, dtype = torch::torch_float(), device = device),
    r_t = torch::torch_tensor(r, dtype = torch::torch_float(), device = device),
    # N_t is n x 1, enabling row-wise broadcasting across p taxa.
    N_t = torch::torch_tensor(N, dtype = torch::torch_float(), device = device)$unsqueeze(2L)
  )
}

# Training of the ZIPLN-VAE. The default settings are those of the
# Supplementary Table of settings of the ZIPLN-VAE: latent dimensions
# K_A = min(6, max(2, q_x + 3)) and K_S = min(4, max(2, q_r + 2)), hidden
# layers of 96 units in the encoder and the abundance network h^A and of 64
# units in the structural network h^S, 100 epochs, mini-batches of min(64, n)
# samples, AdamW with learning rate 1e-3 and weight decay 5e-3, the
# Kullback-Leibler weight w_t = min(1, t / 40), withheld fraction 0.1, weight
# 0.5 of the positive-count term, and structural share 0.5 of the zeros in the
# initial presence intercepts. The result holds the model and history, a data
# frame with one row per epoch: epoch, train_loss (the mean loss of the
# mini-batches of the epoch), val_loss (the loss of the withheld samples after
# the epoch; the training loss when no sample is withheld), and kl_weight.
# progress is an optional function
# called after every epoch as progress(epoch, epochs, train_loss, val_loss,
# kl_weight); it receives these values only.
train_zipln_vae <- function(y, x, r, N,
                            k_abundance = min(6, max(2, NCOL(x) + 3)),
                            k_structural = min(4, max(2, NCOL(r) + 2)),
                            hidden = 96,
                            offset_hidden = 64,
                            epochs = 100,
                            batch_size = min(64, NROW(y)),
                            lr = 1e-3,
                            weight_decay = 5e-3,
                            kl_warmup_epochs = 40,
                            val_frac = 0.1,
                            positive_ll_weight = 0.5,
                            structural_offset_exclude_r_columns = NULL,
                            structural_init_fraction = 0.5,
                            seed = 123,
                            device = "cpu",
                            cluster = NULL,
                            progress = NULL) {
  # The R random-number stream supplies the validation split and the order of
  # the mini-batches; the state of the caller is restored on exit.
  withr::local_seed(seed, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                     .rng_sample_kind = "Rejection")
  torch::torch_manual_seed(seed)

  dat <- prepare_zipln_vae_data(y, x, r, N, device = device)
  n <- nrow(dat$y)
  p <- ncol(dat$y)
  # Initial values:
  #   abundance_intercept_init: length p, the log composition of the pooled
  #   counts.
  #   presence_intercept_init: length p, the presence logit when the fraction
  #   structural_init_fraction of the observed zeros is structural.
  taxon_mass <- colSums(dat$y) + 0.5
  abundance_intercept_init <- log(taxon_mass / sum(taxon_mass))
  zero_rate <- colMeans(dat$y == 0)
  presence_intercept_init <- qlogis(pmin(pmax(1 - structural_init_fraction * zero_rate, 0.02), 0.999))

  # The decoder is trained at the observed covariates. The offsets are decoded
  # with the tested covariates at their reference values (zipln_vae_completion()
  # with x_decode and r_decode), so one fit serves every tested covariate. The
  # structural network excludes the tested prevalence covariates in training
  # and in evaluation (structural_offset_exclude_r_columns), so that their
  # association with presence is represented only by the linear term.
  offset_r_mask <- rep(1, ncol(dat$r))
  offset_r_fill <- colMeans(dat$r)
  r_offset_exclude_idx <- resolve_matrix_columns(
    dat$r,
    structural_offset_exclude_r_columns,
    "structural_offset_exclude_r_columns"
  )
  if (length(r_offset_exclude_idx) > 0L) offset_r_mask[r_offset_exclude_idx] <- 0

  model <- .zipln_vae_module()(
    p = p,
    q_x = ncol(dat$x),
    q_r = ncol(dat$r),
    k_abundance = k_abundance,
    k_structural = k_structural,
    hidden = hidden,
    offset_hidden = offset_hidden,
    abundance_intercept_init = torch::torch_tensor(abundance_intercept_init, dtype = torch::torch_float(), device = device),
    presence_intercept_init = torch::torch_tensor(presence_intercept_init, dtype = torch::torch_float(), device = device),
    offset_r_mask = offset_r_mask,
    offset_r_fill = offset_r_fill
  )
  model$to(device = device)
  optimizer <- torch::optim_adamw(model$parameters, lr = lr, weight_decay = weight_decay)

  n_val <- if (n > 10L) floor(n * val_frac) else 0L
  all_idx <- seq_len(n)
  if (n_val > 0L) {
    # With repeated measures (cluster, the subject of each sample), the
    # withheld samples are those of whole subjects.
    val_idx <- if (is.null(cluster)) {
      sort(sample.int(n, n_val))
    } else {
      .ac3_cluster_validation_split(cluster, n_val)
    }
    train_idx <- setdiff(all_idx, val_idx)
  } else {
    val_idx <- integer(0L)
    train_idx <- all_idx
  }

  history_train <- history_val <- history_kl <- rep(NA_real_, epochs)
  for (epoch in seq_len(epochs)) {
    model$train()
    kl_weight <- min(1, epoch / max(1, kl_warmup_epochs))
    perm <- sample(train_idx, length(train_idx))
    starts <- seq.int(1L, length(perm), by = batch_size)

    batch_loss <- numeric(length(starts))
    for (bb in seq_along(starts)) {
      # Each mini-batch holds a subset of the samples and all p taxa.
      idx <- perm[starts[bb]:min(length(perm), starts[bb] + batch_size - 1L)]

      out <- model(
        dat$y_t[idx, , drop = FALSE],
        dat$x_t[idx, , drop = FALSE],
        dat$r_t[idx, , drop = FALSE],
        dat$N_t[idx, , drop = FALSE],
        kl_weight = kl_weight,
        positive_ll_weight = positive_ll_weight
      )
      total_loss <- out$loss

      optimizer$zero_grad()
      total_loss$backward()
      # The Euclidean norm of the gradient is bounded by 5.
      torch::nn_utils_clip_grad_norm_(model$parameters, max_norm = 5)
      optimizer$step()

      batch_loss[bb] <- total_loss$item()
    }

    # Loss on the withheld samples after each epoch. The forward pass draws
    # the latent scores from the torch random-number stream, so the training
    # batches of later epochs depend on it.
    if (length(val_idx) > 0L) {
      model$eval()
      torch::with_no_grad({
        val_out <- model(dat$y_t[val_idx, , drop = FALSE], dat$x_t[val_idx, , drop = FALSE],
                         dat$r_t[val_idx, , drop = FALSE], dat$N_t[val_idx, , drop = FALSE],
                         kl_weight = kl_weight,
                         positive_ll_weight = positive_ll_weight)
      })
      val_loss <- val_out$loss$item()
    } else {
      val_loss <- mean(batch_loss)
    }

    history_train[epoch] <- mean(batch_loss)
    history_val[epoch] <- val_loss
    history_kl[epoch] <- kl_weight
    if (!is.null(progress)) progress(epoch, epochs, history_train[epoch], val_loss, kl_weight)
  }

  list(
    model = model,
    history = data.frame(epoch = seq_len(epochs), train_loss = history_train, val_loss = history_val,
                         kl_weight = history_kl),
    n_validation = length(val_idx)
  )
}

zipln_vae_completion <- function(model, y, x, r, N,
                                 x_decode = NULL,
                                 r_decode = NULL,
                                 device = "cpu") {
  dat <- prepare_zipln_vae_data(y, x, r, N, device = device)
  model$to(device = device)
  model$eval()

  # The encoder reads the observed counts. The decoder is evaluated at the
  # covariates x_decode and r_decode (in ancombc3_fitted_vae_offsets(), the
  # covariates with the tested columns at their reference values); the
  # defaults are the observed covariates.
  x_decode_mat <- if (is.null(x_decode)) dat$x else as_numeric_matrix(x_decode, "x_decode")
  r_decode_mat <- if (is.null(r_decode)) dat$r else as_numeric_matrix(r_decode, "r_decode")
  x_dec_t <- torch::torch_tensor(x_decode_mat, dtype = torch::torch_float(), device = device)
  r_dec_t <- torch::torch_tensor(r_decode_mat, dtype = torch::torch_float(), device = device)

  torch::with_no_grad({
    enc <- model$encode(model$encoder_transform(dat$y_t, dat$N_t))

    # The decoder is evaluated once at the variational posterior means. All
    # decoder outputs below are n x p.
    dec <- model$decode(enc$mean_abundance, enc$mean_structural, x_dec_t, r_dec_t, dat$N_t)
    completed_zero <- dec$completed_zero
    mu <- dec$mu
    pi_dec <- dec$pi
    v <- dec$v
    h_structural <- dec$h_structural
  })

  list(
    # Plain R matrices on the CPU, n x p, except sigma (length p).
    # completed_zero is v mu; f0 is the Poisson-log-normal sampling-zero
    # probability f(0; mu, sigma) = E exp(-Lambda), which is at most 1. When
    # mu is close to 0, the exact f0 lies within mu of 1, and the
    # single-precision quadrature of log f0 can exceed 0 by rounding error
    # (by at most about 1e-5 over the range of sigma). f0 is therefore
    # bounded above by 1, its limit as mu tends to 0; the bound changes only
    # the cells whose computed f0 exceeds 1.
    completed_zero = torch::as_array(completed_zero$to(device = "cpu")),
    mu = torch::as_array(mu$to(device = "cpu")),
    pi = torch::as_array(pi_dec$to(device = "cpu")),
    v = torch::as_array(v$to(device = "cpu")),
    f0 = torch::as_array(torch::torch_exp(dec$log_f0)$clamp_max(1)$to(device = "cpu")),
    h_structural = torch::as_array(h_structural$to(device = "cpu")),
    sigma = torch::as_array(dec$sigma$detach()$to(device = "cpu"))
  )
}
