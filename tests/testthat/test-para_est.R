context("Testing parameter estimation function")
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
test_that("`para_test` function provides expected values prior to the
          E-M algorithm", {
      feature_table = abundances(phylum_data); meta_data = meta(phylum_data)
      meta_data = meta_data %>% rownames_to_column("sample_id")
      sample_id = "sample_id"; group = "nation"; zero_cut = 0.90; lib_cut = 1000
      global = TRUE; formula = "age + nation + bmi_group"
      tol = 1e-5; max_iter = 100

      # Data pre-processing
      fiuo_prep = data_prep(feature_table, meta_data, sample_id,
                            group, zero_cut, lib_cut, global)
      feature_table = fiuo_prep$feature_table; meta_data = fiuo_prep$meta_data
      y = log(feature_table + 1)

      # Parameters estimation
      fiuo_para = para_est(y, meta_data, formula, tol, max_iter)
      beta = fiuo_para$beta; d = fiuo_para$d; e = fiuo_para$e
      var_cov_hat = fiuo_para$var_cov_hat; var_hat = fiuo_para$var_hat

      # Testing
      test_output = signif(c(beta[1, 1], var_hat[2, 1]), 2)

      expect_equal(test_output, c(5.80, 0.02))
})


