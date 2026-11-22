test_that("the compiled engine does not contract floating-point expressions", {
    # (1 + 2^-27) * (1 - 2^-27) = 1 - 2^-54, which rounds to 1. The sum with
    # -1 is 0 when the product is rounded first and -2^-54 under a fused
    # multiply-add.
    expect_identical(ANCOMBC:::.ac3_mul_add(1 + 2^-27, 1 - 2^-27, -1), 0)
    expect_identical(ANCOMBC:::.ac3_mul_add(2, 3, 4), 10)
})

test_that("an absent torch gives the installation commands", {
    skip_if(ANCOMBC:::.ac3_torch_available(), "torch is available")
    expect_error(ANCOMBC:::.ac3_require_torch(), "install")
})

test_that("an available torch passes the check", {
    skip_if_not_installed("torch")
    skip_if_not(torch::torch_is_installed(), "LibTorch is not installed")
    expect_true(ANCOMBC:::.ac3_torch_available())
    expect_true(ANCOMBC:::.ac3_require_torch())
})
