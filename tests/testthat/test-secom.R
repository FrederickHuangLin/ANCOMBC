context("Testing secom functions")
library(ANCOMBC)
library(testthat)
library(tidyverse)
library(microbiome)

data(atlas1006)
# Subset to baseline
pseq = subset_samples(atlas1006, time == 0)

# Re-code the bmi group
sample_data(pseq)$bmi_group = recode(sample_data(pseq)$bmi_group,
                                     lean = "lean",
                                     overweight = "overweight",
                                     obese = "obese",
                                     severeobese = "obese",
                                     morbidobese = "obese")
# Subset to lean and obese subjects
pseq = subset_samples(pseq, bmi_group %in% c("lean", "obese"))

# Create the region variable
sample_data(pseq)$region = recode(sample_data(pseq)$nationality,
                                  Scandinavia = "NE",
                                  UKIE = "NE",
                                  SouthEurope = "SE",
                                  CentralEurope = "CE",
                                  EasternEurope = "EE")
# Discard "EE" as it contains only 1 subject
pseq = subset_samples(pseq, region != "EE")

# Genus level data
genus_data = pseq
# Phylum level data
phylum_data = aggregate_taxa(pseq, "Phylum")

# Test 1
test_that("`secom_linear` function provides expected results", {
   pseqs = list(c(genus_data, phylum_data)); pseudo = 0
   prv_cut = 0.5; lib_cut = 1000; corr_cut = 0.5
   wins_quant = c(0.05, 0.95); method = "pearson"
   soft = FALSE; thresh_len = 20; n_cv = 10
   thresh_hard = 0.3; max_p = 0.005; n_cl = 1

   set.seed(123)
   res_linear = secom_linear(pseqs, pseudo, prv_cut, lib_cut, corr_cut,
                             wins_quant, method, soft, thresh_len, n_cv,
                             thresh_hard, max_p, n_cl)
   # Testing
   test_output = round(c(res_linear$corr_th[2, 1],
                         res_linear$corr_p[2, 1]), 2)

   expect_equal(test_output, c(-0.47, 0.00))
})

# Test 2
test_that("`secom_dist` function provides expected results", {
   pseqs = list(c(genus_data, phylum_data)); pseudo = 0
   prv_cut = 0.5; lib_cut = 1000; corr_cut = 0.5
   wins_quant = c(0.05, 0.95); R = 100
   thresh_hard = 0.3; max_p = 0.05; n_cl = 1

   set.seed(123)
   res_dist = secom_dist(pseqs, pseudo, prv_cut, lib_cut, corr_cut,
                         wins_quant, R, thresh_hard, max_p, n_cl)
   # Testing
   test_output = round(c(res_dist$dcorr_fl[2, 1],
                         res_dist$dcorr_p[2, 1]), 2)

   expect_equal(test_output, c(0.49, 0.01))
})


