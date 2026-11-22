# Agreement of ancombc3() with the reference fit of the method code.
#
# fixtures/ancombc3_reference_toy.rds holds the result tables of the method
# code at its tag (element method_tag) for the toy data and the call recorded
# in its element call. The package code performs the same operations, so the
# two agree up to the differences between platforms in the single-precision
# arithmetic of the first stage.

test_that("ancombc3() agrees with the reference fit of the method code", {
  skip_if_no_torch()
  torch::torch_set_num_threads(1L)
  ref <- readRDS(test_path("fixtures", "ancombc3_reference_toy.rds"))
  toy <- ac3_toy_data()
  fit <- ancombc3(data = toy$counts, meta_data = toy$meta, fix_formula = ~ group + bmi_z,
                  tested = c("group", "bmi_z"), alpha = 0.10, seed = 1, engine = "compiled",
                  verbose = FALSE)
  res <- fit$res
  expect_identical(res[c("taxon", "covariate", "covariate_type")],
                   ref$res[c("taxon", "covariate", "covariate_type")])
  expect_identical(res$da_eligible, ref$res$da_eligible)
  expect_identical(res$dp_eligible, ref$res$dp_eligible)
  tol <- ac3_reference_tolerance
  for (nm in c("da_estimate", "da_se", "dp_estimate")) {
    expect_equal(res[[nm]], ref$res[[nm]], tolerance = tol, label = nm)
  }
  expect_equal(unname(fit$sample_bias), unname(ref$sample_bias), tolerance = tol)
  expect_equal(unname(fit$sigma), unname(ref$sigma), tolerance = tol)
  expect_identical(res$da_discovery, ref$res$da_discovery)
  expect_identical(res$dp_discovery, ref$res$dp_discovery)
  expect_identical(res$omnibus_discovery, ref$res$omnibus_discovery)
})
