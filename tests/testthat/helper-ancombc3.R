# Shared helpers of the ancombc3() tests.

# Skips the test when the package torch or its LibTorch library is absent.
skip_if_no_torch <- function() {
  testthat::skip_if_not_installed("torch")
  testthat::skip_if_not(torch::torch_is_installed(), "LibTorch is not installed")
}

# Skips the test unless the environment variable ANCOMBC_EXTENDED_TESTS is
# "true".
skip_if_not_extended <- function() {
  testthat::skip_if_not(identical(Sys.getenv("ANCOMBC_EXTENDED_TESTS"), "true"),
                        "extended tests run only with ANCOMBC_EXTENDED_TESTS=true")
}

# Toy data installed with the package: a count table with taxa in rows and
# the metadata of its samples (group, a factor with levels control and
# treated, and bmi_z).
ac3_toy_data <- function() {
  counts <- as.matrix(utils::read.csv(system.file("extdata", "ancombc3_toy_counts.csv", package = "ANCOMBC"),
                                      row.names = 1, check.names = FALSE))
  meta <- utils::read.csv(system.file("extdata", "ancombc3_toy_metadata.csv", package = "ANCOMBC"),
                          row.names = 1, check.names = FALSE)
  meta$group <- factor(meta$group, levels = c("control", "treated"))
  list(counts = counts, meta = meta)
}

# TRUE when the evaluation of expr stops with an error whose message matches
# pattern.
ac3_errs <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))
  grepl(pattern, msg)
}

# Function that evaluates compute() at its first call and returns the stored
# value at later calls.
ac3_cached <- function(compute) {
  value <- NULL
  done <- FALSE
  function() {
    if (!done) {
      value <<- compute()
      done <<- TRUE
    }
    value
  }
}

# Relative tolerance of the comparison with the reference fit of the method
# code (test-ancombc3-reference.R).
ac3_reference_tolerance <- 1e-6
