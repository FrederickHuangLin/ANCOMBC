context("Testing bias estimation function")
library(ANCOMBC)
library(testthat)
library(tidyverse)
library(microbiome)

data(atlas1006)
# Subset to baseline
pseq = subset_samples(atlas1006, time == 0)
# Re-code the bmi group
sample_data(pseq)$bmi_group = recode(sample_data(pseq)$bmi_group,
                                     underweight = "lean",
                                     lean = "lean",
                                     overweight = "overweight",
                                     obese = "obese",
                                     severeobese = "obese",
                                     morbidobese = "obese")
# Re-code the nationality group
sample_data(pseq)$nation = recode(sample_data(pseq)$nationality,
                                  Scandinavia = "NE",
                                  UKIE = "NE",
                                  SouthEurope = "SE",
                                  CentralEurope = "CE",
                                  EasternEurope = "EE")

# Aggregate to phylum level
phylum_data = aggregate_taxa(pseq, "Phylum")

# Test
test_that("`bias_test` function provides expected values for bias terms", {
   phyloseq = phylum_data; group = "nation"; zero_cut = 0.90; lib_cut = 1000
   global = TRUE; formula = "age + nation + bmi_group"
   tol = 1e-5; max_iter = 100

   # Data pre-processing
   fiuo_prep = data_prep(phyloseq, group, zero_cut, lib_cut, global)
   feature_table = fiuo_prep$feature_table
   meta_data = fiuo_prep$meta_data
   y = log(feature_table + 1)

   # Parameters estimation
   fiuo_para = para_est(y, meta_data, formula, tol, max_iter)
   beta = fiuo_para$beta; d = fiuo_para$d; e = fiuo_para$e
   var_cov_hat = fiuo_para$var_cov_hat; var_hat = fiuo_para$var_hat

   # Bias estimation
   fiuo_bias = bias_est(beta, var_hat, tol, max_iter, n_taxa = nrow(y))
   delta_em = fiuo_bias$delta_em; delta_wls = fiuo_bias$delta_wls
   var_delta = fiuo_bias$var_delta

   # Testing
   test_output = signif(c(delta_em[1], delta_wls[2], var_delta[3]), 2)

   expect_equal(test_output, c(0.00330, 0.43000, 0.00029))
})


