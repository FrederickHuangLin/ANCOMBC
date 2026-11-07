test_that("variance search escapes a zero boundary", {
    # Equal known variances have analytic optimum max(mean(residual^2)-v, 0).
    expect_equal(.ancombc_variance(1, 0.99, 1, 0, 0), 0.01,
                 tolerance = 1e-7)
    expect_equal(.ancombc_variance(c(-2, -1, 1, 2), rep(0.5, 4),
                                  rep(1, 4), 0, 0), 2, tolerance = 1e-7)
    expect_equal(.ancombc_variance(c(-1, 1), c(2, 2), c(1, 1), 0, 3), 0)
    expect_equal(.ancombc_variance(c(2, 2), c(1, 2), c(1, 1), 2, 1), 0)
    expect_equal(.ancombc_variance(c(1, 2), c(1, 2), c(0, 0), 0, 3), 3)
})

test_that("variance updates preserve the current conditional likelihood", {
    beta = c(-4, -1, 0, 2, 5)
    known = c(0.05, 0.1, 1, 0.2, 2)
    weights = c(0.01, 0.1, 0.8, 0.05, 0.04)
    score = function(t) -sum(weights * dnorm(beta, 0, sqrt(known + t), log = TRUE))
    for (current in c(0, 0.01, 0.1, 1, 5, 30)) {
        updated = .ancombc_variance(beta, known, weights, 0, current)
        expect_true(is.finite(updated) && updated >= 0)
        expect_lte(score(updated), score(current) + 1e-12)
    }
    # log=TRUE keeps an extreme tail finite instead of replacing log(0) by 0.
    expect_true(is.finite(.ancombc_variance(c(0, 1000), c(1, 1),
                                          c(0.999, 0.001), 0, 0)))
})
