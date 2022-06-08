context("Testing data QC function")
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

# Test 1
test_that("`data_qc` function returns error when the group variable is
          missing for the global test", {
    meta_data = meta(family_data); group = NULL; global = TRUE

    expect_error(data_qc(meta_data, group, global),
                 "Please specify the group variable for the global test")
})

# Test 2
test_that("`data_qc` function returns error when the number of categories for
          the group variable is less than 2", {
    # Subset data
    family_data2 = subset_samples(family_data, region == "NE")

    # Testing
    meta_data = meta(family_data2); group = "region"; global = TRUE

    expect_error(data_qc(meta_data, group, global),
                 "The group variable should have >= 2 categories")
})


