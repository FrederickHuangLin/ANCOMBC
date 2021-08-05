context("Testing data preprocessing function")
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

# Test 1
test_that("`data_prep` function returns error when the group variable is
          missing for the global test", {
    phyloseq = phylum_data; group = NULL; zero_cut = 0.90
    lib_cut = 1000; global = TRUE

    expect_error(data_prep(phyloseq, group, zero_cut, lib_cut, global),
                 "Please specify the group variable for the global test")
})

# Test 2
test_that("`data_prep` function returns error when the number of categories for
          the group variable is less than 2", {
    # Subset data
    phylum_data2 = subset_samples(phylum_data, nation == "NE")

    # Testing
    phyloseq = phylum_data2; group = "nation"; zero_cut = 0.90
    lib_cut = 1000; global = TRUE

    expect_error(data_prep(phyloseq, group, zero_cut, lib_cut, global),
                 "The group variable should have >= 2 categories")
})

# Test 3
test_that("`data_prep` function should only activate the global test when the
          number of categories for the group variable is greater than or equal
          to 3", {
      # Subset data
      phylum_data3 = subset_samples(phylum_data, nation %in% c("NE", "SE"))

      # Testing
      phyloseq = phylum_data3; group = "nation"; zero_cut = 0.90
      lib_cut = 1000; global = TRUE

      fiuo_prep = data_prep(phyloseq, group, zero_cut, lib_cut, global)
      expect_equal(fiuo_prep$global, FALSE)
})



