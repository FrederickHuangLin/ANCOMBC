# ---------------------------------------------------------------------------
# workloads.R
#
# Reproducibility / timing / memory workloads for ANCOMBC.
#
# Each element of `workloads` is a list with
#   scale : "vignette" (an evaluated vignette chunk, reproduced verbatim) or
#           "large"    (a larger configuration reflecting real use)
#   fn    : the ANCOMBC entry point exercised
#   prep  : function() returning the data objects used by the call. Run once,
#           outside the timed region.
#   call  : function(d) returning the result object. This is the timed region.
#   note  : free text, used for the large workloads to record sizing choices.
#
# This file does NOT attach ANCOMBC. run.R attaches it from a specific library
# path so that the same workload definitions can be run against several
# installations.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(dplyr)
    library(phyloseq)
})

# --- shared data preparation ------------------------------------------------

# atlas1006, prepared exactly as in ANCOM.Rmd / ANCOMBC.Rmd / ANCOMBC2.Rmd /
# SECOM.Rmd (the same block is repeated in all four vignettes).
.prep_atlas = local({
    cached = NULL
    function() {
        if (!is.null(cached)) return(cached)
        data(atlas1006, package = "microbiome", envir = environment())

        # Subset to baseline
        pseq = phyloseq::subset_samples(atlas1006, time == 0)

        # Re-code the bmi group
        meta_data = microbiome::meta(pseq)
        meta_data$bmi = recode(meta_data$bmi_group,
                               obese = "obese",
                               severeobese = "obese",
                               morbidobese = "obese")
        meta_data$bmi = factor(meta_data$bmi,
                               levels = c("obese", "overweight", "lean"))

        # Create the region variable
        meta_data$region = recode(as.character(meta_data$nationality),
                                  Scandinavia = "NE", UKIE = "NE",
                                  SouthEurope = "SE", CentralEurope = "CE",
                                  EasternEurope = "EE", .missing = "unknown")

        phyloseq::sample_data(pseq) = meta_data

        pseq = phyloseq::subset_samples(pseq,
                                        bmi %in% c("lean", "overweight", "obese"))
        pseq = phyloseq::subset_samples(pseq, ! region %in% c("EE", "unknown"))
        cached <<- pseq
        pseq
    }
})

# dietswap, as in ANCOM.Rmd / ANCOMBC2.Rmd.
.prep_dietswap = local({
    cached = NULL
    function() {
        if (!is.null(cached)) return(cached)
        data(dietswap, package = "microbiome", envir = environment())
        cached <<- dietswap
        dietswap
    }
})

# Permuted-bmi copy of atlas1006, as in ANCOMBC2.Rmd section 4.
.prep_atlas_perm = function() {
    pseq = .prep_atlas()
    set.seed(123)
    pseq_perm = pseq
    meta_data_perm = microbiome::meta(pseq_perm)
    meta_data_perm$bmi = sample(meta_data_perm$bmi)
    phyloseq::sample_data(pseq_perm) = meta_data_perm
    pseq_perm
}

# Poisson log-normal simulation seeded by QMP, as in ANCOMBC2.Rmd section 3.
# `n` is the sample size; the vignette uses n = 150.
.prep_qmp_sim = function(n = 150) {
    data(QMP, package = "ANCOMBC", envir = environment())
    set.seed(123)
    d = ncol(QMP)
    diff_prop = 0.1
    lfc_cont = 1
    lfc_cat2_vs_1 = -2
    lfc_cat3_vs_1 = 1

    # Generate the true abundances
    abn_data = sim_plnm(abn_table = QMP, taxa_are_rows = FALSE, prv_cut = 0.05,
                        n = n, lib_mean = 1e8, disp = 0.5)
    log_abn_data = log(abn_data + 1e-5)
    rownames(log_abn_data) = paste0("T", seq_len(d))
    colnames(log_abn_data) = paste0("S", seq_len(n))

    # Generate the sample and feature meta data
    smd = data.frame(samp_frac = log(c(runif(n/3, min = 1e-4, max = 1e-3),
                                       runif(n/3, min = 1e-3, max = 1e-2),
                                       runif(n/3, min = 1e-2, max = 1e-1))),
                     cont_cov = rnorm(n),
                     cat_cov = as.factor(rep(seq_len(3), each = n/3)))
    rownames(smd) = paste0("S", seq_len(n))

    fmd = data.frame(taxon = paste0("T", seq_len(d)),
                     seq_eff = log(runif(d, min = 0.1, max = 1)),
                     lfc_cont = sample(c(0, lfc_cont),
                                       size = d, replace = TRUE,
                                       prob = c(1 - diff_prop, diff_prop)),
                     lfc_cat2_vs_1 = sample(c(0, lfc_cat2_vs_1),
                                            size = d, replace = TRUE,
                                            prob = c(1 - diff_prop, diff_prop)),
                     lfc_cat3_vs_1 = sample(c(0, lfc_cat3_vs_1),
                                            size = d, replace = TRUE,
                                            prob = c(1 - diff_prop, diff_prop))) %>%
        mutate(lfc_cat3_vs_2 = lfc_cat3_vs_1 - lfc_cat2_vs_1)

    # Add effect sizes of covariates to the true abundances
    smd_dmy = model.matrix(~ 0 + cont_cov + cat_cov, data = smd)
    log_abn_data = log_abn_data + outer(fmd$lfc_cont, smd_dmy[, "cont_cov"])
    log_abn_data = log_abn_data + outer(fmd$lfc_cat2_vs_1, smd_dmy[, "cat_cov2"])
    log_abn_data = log_abn_data + outer(fmd$lfc_cat3_vs_1, smd_dmy[, "cat_cov3"])

    # Add sample- and taxon-specific biases
    log_otu_data = t(t(log_abn_data) + smd$samp_frac)
    log_otu_data = log_otu_data + fmd$seq_eff
    otu_data = round(exp(log_otu_data))

    list(otu_data = otu_data, smd = smd, fmd = fmd)
}

# Two-group atlas1006 split, as in SECOM.Rmd section 4.
.prep_atlas_2groups = function() {
    pseq = .prep_atlas()
    pseq1 = phyloseq::subset_samples(pseq, region == "CE")
    pseq2 = phyloseq::subset_samples(pseq, region == "NE")
    phyloseq::sample_names(pseq1) =
        paste0("Sample-", seq_len(phyloseq::nsamples(pseq1)))
    phyloseq::sample_names(pseq2) =
        paste0("Sample-", seq_len(phyloseq::nsamples(pseq2)))
    list(pseq1 = pseq1, pseq2 = pseq2)
}

# --- workloads --------------------------------------------------------------

workloads = list()

## ---- vignette scale: ANCOM.Rmd --------------------------------------------

workloads$ancom_atlas = list(
    scale = "vignette", fn = "ancom",
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        ancom(data = d$pseq, tax_level = "Family", meta_data = NULL,
              p_adj_method = "holm", prv_cut = 0.10,
              lib_cut = 1000, main_var = "bmi", adj_formula = "age + region",
              rand_formula = NULL, lme_control = NULL, struc_zero = TRUE,
              neg_lb = TRUE, alpha = 0.05, n_cl = 2, verbose = TRUE)
    })

workloads$ancom_dietswap = list(
    scale = "vignette", fn = "ancom",
    prep = function() list(dietswap = .prep_dietswap()),
    call = function(d) {
        ancom(data = d$dietswap, tax_level = "Family",
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
              main_var = "group",
              adj_formula = "nationality + timepoint",
              rand_formula = "(timepoint | subject)",
              lme_control = lme4::lmerControl(),
              struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 2)
    })

## ---- vignette scale: ANCOMBC.Rmd ------------------------------------------

workloads$ancombc_atlas = list(
    scale = "vignette", fn = "ancombc",
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        ancombc(data = d$pseq, tax_level = "Family",
                formula = "age + region + bmi",
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
                group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                n_cl = 1, verbose = TRUE)
    })

## ---- vignette scale: ANCOMBC2.Rmd -----------------------------------------

workloads$ancombc2_qmp = list(
    scale = "vignette", fn = "ancombc2",
    prep = function() .prep_qmp_sim(n = 150),
    call = function(d) {
        ancombc2(data = d$otu_data, meta_data = d$smd,
                 fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                 group = "cat_cov", struc_zero = FALSE, neg_lb = FALSE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = FALSE, pairwise = TRUE,
                 dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20,
                                     verbose = FALSE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = NULL,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                 trend_control = NULL)
    })

workloads$ancombc2_atlas_perm = list(
    scale = "vignette", fn = "ancombc2",
    prep = function() list(pseq_perm = .prep_atlas_perm()),
    call = function(d) {
        ancombc2(data = d$pseq_perm, tax_level = "Genus",
                 fix_formula = "bmi", rand_formula = NULL,
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
                 group = "bmi", struc_zero = TRUE, neg_lb = TRUE)
    })

workloads$ancombc2_atlas_family = list(
    scale = "vignette", fn = "ancombc2",
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        ancombc2(data = d$pseq, tax_level = "Family",
                 fix_formula = "age + region + bmi", rand_formula = NULL,
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                 group = "bmi", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                 iter_control = list(tol = 1e-2, max_iter = 20,
                                     verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                            nrow = 2,
                                                            byrow = TRUE),
                                                     matrix(c(-1, 0, 1, -1),
                                                            nrow = 2,
                                                            byrow = TRUE),
                                                     matrix(c(1, 0, 1, -1),
                                                            nrow = 2,
                                                            byrow = TRUE)),
                                      node = list(2, 2, 1),
                                      B = 10))
    })

workloads$ancombc2_atlas_interact = list(
    scale = "vignette", fn = "ancombc2",
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        ancombc2(data = d$pseq, tax_level = "Family",
                 fix_formula = "age * bmi + region", rand_formula = NULL,
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                 group = "bmi", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE)
    })

workloads$ancombc2_dietswap = list(
    scale = "vignette", fn = "ancombc2",
    prep = function() list(dietswap = .prep_dietswap()),
    call = function(d) {
        ancombc2(data = d$dietswap, tax_level = "Family",
                 fix_formula = "nationality + timepoint + group",
                 rand_formula = "(timepoint | subject)",
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                 group = "group", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                 iter_control = list(tol = 1e-2, max_iter = 20,
                                     verbose = TRUE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = lme4::lmerControl(),
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                            nrow = 2,
                                                            byrow = TRUE)),
                                      node = list(2),
                                      B = 10))
    })

## ---- vignette scale: SECOM.Rmd --------------------------------------------

workloads$secom_linear_atlas = list(
    scale = "vignette", fn = "secom_linear",
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        secom_linear(data = list(d$pseq), taxa_are_rows = TRUE,
                     tax_level = "Phylum",
                     aggregate_data = NULL, meta_data = NULL, pseudo = 0,
                     prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
                     wins_quant = c(0.05, 0.95), method = "pearson",
                     soft = FALSE, thresh_len = 20, n_cv = 10,
                     thresh_hard = 0.3, max_p = 0.005, n_cl = 2)
    })

workloads$secom_dist_atlas = list(
    scale = "vignette", fn = "secom_dist",
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        secom_dist(data = list(d$pseq), taxa_are_rows = TRUE,
                   tax_level = "Phylum",
                   aggregate_data = NULL, meta_data = NULL, pseudo = 0,
                   prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
                   wins_quant = c(0.05, 0.95), R = 1000,
                   thresh_hard = 0.3, max_p = 0.005, n_cl = 2)
    })

workloads$secom_linear_2groups = list(
    scale = "vignette", fn = "secom_linear",
    prep = function() .prep_atlas_2groups(),
    call = function(d) {
        secom_linear(data = list(CE = d$pseq1, NE = d$pseq2),
                     taxa_are_rows = TRUE,
                     tax_level = c("Phylum", "Phylum"),
                     aggregate_data = NULL, meta_data = NULL, pseudo = 0,
                     prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
                     wins_quant = c(0.05, 0.95), method = "pearson",
                     soft = FALSE, thresh_len = 20, n_cv = 10,
                     thresh_hard = 0.3, max_p = 0.005, n_cl = 2)
    })

workloads$secom_dist_2groups = list(
    scale = "vignette", fn = "secom_dist",
    prep = function() .prep_atlas_2groups(),
    call = function(d) {
        secom_dist(data = list(CE = d$pseq1, NE = d$pseq2),
                   taxa_are_rows = TRUE,
                   tax_level = c("Phylum", "Phylum"),
                   aggregate_data = NULL, meta_data = NULL, pseudo = 0,
                   prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
                   wins_quant = c(0.05, 0.95), R = 1000,
                   thresh_hard = 0.3, max_p = 0.005, n_cl = 2)
    })

## ---- large scale ----------------------------------------------------------
# Sizing rule: each large workload is kept under roughly five minutes on
# lib_base. Where the direct scale-up exceeded that, the configuration was
# reduced and the reduction is recorded in `note`.

workloads$ancombc2_atlas_genus_large = list(
    scale = "large", fn = "ancombc2",
    note = paste("atlas1006 at Genus level (114 taxa after prv_cut, 873",
                 "samples) with the same fix_formula as the ANCOMBC2",
                 "vignette, and the global, pairwise, Dunnett and trend",
                 "tests. The trend bootstrap uses B = 10, as in the vignette;",
                 "the recommended B = 100 at this taxon count exceeds the",
                 "five-minute budget."),
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        ancombc2(data = d$pseq, tax_level = "Genus",
                 fix_formula = "age + region + bmi", rand_formula = NULL,
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                 group = "bmi", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                 iter_control = list(tol = 1e-2, max_iter = 20,
                                     verbose = FALSE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = NULL,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                 trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                            nrow = 2,
                                                            byrow = TRUE),
                                                     matrix(c(-1, 0, 1, -1),
                                                            nrow = 2,
                                                            byrow = TRUE),
                                                     matrix(c(1, 0, 1, -1),
                                                            nrow = 2,
                                                            byrow = TRUE)),
                                      node = list(2, 2, 1),
                                      B = 10))
    })

workloads$ancombc2_qmp_otu_large = list(
    scale = "large", fn = "ancombc2",
    note = paste("QMP-seeded Poisson log-normal simulation at the OTU level",
                 "(91 taxa) with n = 900 samples, six times the vignette's",
                 "n = 150, plus pseudo_sens = TRUE, structural zeros and the",
                 "global, pairwise and Dunnett screens."),
    prep = function() .prep_qmp_sim(n = 900),
    call = function(d) {
        ancombc2(data = d$otu_data, meta_data = d$smd,
                 fix_formula = "cont_cov + cat_cov", rand_formula = NULL,
                 p_adj_method = "holm", pseudo_sens = TRUE,
                 prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                 group = "cat_cov", struc_zero = TRUE, neg_lb = TRUE,
                 alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = TRUE, pairwise = TRUE,
                 dunnet = TRUE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20,
                                     verbose = FALSE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = NULL,
                 mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                 trend_control = NULL)
    })

workloads$ancombc_atlas_genus_large = list(
    scale = "large", fn = "ancombc",
    note = "atlas1006 at Genus level, same formula as the ANCOMBC vignette.",
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        ancombc(data = d$pseq, tax_level = "Genus",
                formula = "age + region + bmi",
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
                group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                n_cl = 1, verbose = TRUE)
    })

workloads$ancom_atlas_genus_large = list(
    scale = "large", fn = "ancom",
    note = paste("atlas1006 at Genus level with the vignette's prv_cut =",
                 "0.10, giving 114 taxa and therefore 114^2 ALR models,",
                 "against 22 taxa at the vignette's Family level."),
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        ancom(data = d$pseq, tax_level = "Genus", meta_data = NULL,
              p_adj_method = "holm", prv_cut = 0.10,
              lib_cut = 1000, main_var = "bmi", adj_formula = "age + region",
              rand_formula = NULL, lme_control = NULL, struc_zero = TRUE,
              neg_lb = TRUE, alpha = 0.05, n_cl = 2, verbose = FALSE)
    })

workloads$secom_linear_atlas_genus_large = list(
    scale = "large", fn = "secom_linear",
    note = paste("atlas1006 at Genus level (114 taxa) with prv_cut = 0.10 and",
                 "a threshold grid of 100 points, against 6 taxa and 20 grid",
                 "points at the vignette's Phylum level. secom_dist is not",
                 "included at Genus level: its permutation test costs O(d^2 R)",
                 "distance correlations, far beyond the five-minute budget."),
    prep = function() list(pseq = .prep_atlas()),
    call = function(d) {
        secom_linear(data = list(d$pseq), taxa_are_rows = TRUE,
                     tax_level = "Genus",
                     aggregate_data = NULL, meta_data = NULL, pseudo = 0,
                     prv_cut = 0.10, lib_cut = 1000, corr_cut = 0.5,
                     wins_quant = c(0.05, 0.95), method = "pearson",
                     soft = FALSE, thresh_len = 100, n_cv = 10,
                     thresh_hard = 0.3, max_p = 0.005, n_cl = 2)
    })
