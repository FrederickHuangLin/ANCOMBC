# Infrastructure of ancombc3(): the optional dependency torch and the compiled
# engine.

#' @useDynLib ANCOMBC, .registration = TRUE
NULL

# TRUE when the package torch can be loaded and its LibTorch library is
# installed.
.ac3_torch_available <- function() {
    requireNamespace("torch", quietly = TRUE) &&
        isTRUE(tryCatch(torch::torch_is_installed(),
                        error = function(e) FALSE))
}

# Stops with the installation commands when torch or its LibTorch library is
# absent.
.ac3_require_torch <- function() {
    if (!requireNamespace("torch", quietly = TRUE)) {
        stop("ancombc3() requires the package 'torch', which is not ",
             "installed. Install it with install.packages(\"torch\") and ",
             "then run torch::install_torch() once to download the LibTorch ",
             "library.", call. = FALSE)
    }
    installed <- tryCatch(torch::torch_is_installed(),
                          error = function(e) FALSE)
    if (!isTRUE(installed)) {
        stop("ancombc3() requires the LibTorch library of the package ",
             "'torch', which is not installed. Run torch::install_torch() ",
             "once.", call. = FALSE)
    }
    invisible(TRUE)
}

# a * b + c evaluated by the compiled engine, with the product rounded before
# the addition.
.ac3_mul_add <- function(a, b, c) {
    .Call(ac3_mul_add_c, as.double(a), as.double(b), as.double(c))
}
