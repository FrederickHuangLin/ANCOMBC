context("Testing ancombc function")
library(ANCOMBC); library(testthat); library(tidyverse); library(microbiome)

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
test_that("`ancombc` function provides expected results", {
      feature_table = abundances(phylum_data); meta_data = meta(phylum_data)
      meta_data = meta_data %>% rownames_to_column("sample_id")
      sample_id = "sample_id"; formula = "age + nation + bmi_group"
      p_adj_method = "holm"; zero_cut = 0.90; lib_cut = 1000; group = "nation"
      struc_zero = TRUE; neg_lb = TRUE; tol = 1e-5; max_iter = 100
      conserve = TRUE; alpha = 0.05; global = TRUE

      out = ancombc(feature_table, meta_data, sample_id, formula, p_adj_method,
                    zero_cut, lib_cut, group, struc_zero, neg_lb,
                    tol, max_iter, conserve, alpha, global)

      res = out$res
      res_global = out$res_global


      # Testing
      test_output = signif(c(res$W[1, 1], res_global$W[2]), 2)

      expect_equal(test_output, c(48, 75))
})


