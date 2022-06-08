context("Testing ancom function")
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

# Aggregate to family level
family_data = aggregate_taxa(pseq, "Family")

# Test
test_that("`ancom` function provides expected results", {
   phyloseq = family_data; p_adj_method = "holm"
   prv_cut = 0.10; lib_cut = 1000; main_var = "bmi_group"
   adj_formula = "age + region"; rand_formula = NULL
   lme_control = NULL; struc_zero = TRUE; neg_lb = TRUE
   alpha = 0.05; n_cl = 1

   set.seed(123)
   out = ancom(phyloseq, p_adj_method, prv_cut, lib_cut, main_var,
               adj_formula, rand_formula, lme_control, struc_zero, neg_lb,
               alpha, n_cl)
   res = out$res

   # Testing
   test_output = res$W[1]

   expect_equal(test_output, 3)
})


