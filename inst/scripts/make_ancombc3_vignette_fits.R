# Stored results of the vignette "ANCOM-BC3 Tutorial"
# (inst/extdata/ancombc3_vignette_fits.rds).
#
# The vignette fits the models when the package torch and its LibTorch library
# are available, and loads these stored results otherwise. The calls below are
# those of the vignette. Run from the package directory, with ANCOMBC, torch,
# phyloseq, and microbiome installed:
#   Rscript inst/scripts/make_ancombc3_vignette_fits.R

library(ANCOMBC)

read_example = function(name) {
    counts = read.csv(system.file("extdata",
                                  paste0("ancombc3_", name, "_counts.csv"),
                                  package = "ANCOMBC"),
                      row.names = 1, check.names = FALSE)
    meta_data = read.csv(system.file("extdata",
                                     paste0("ancombc3_", name, "_metadata.csv"),
                                     package = "ANCOMBC"),
                         row.names = 1)
    list(counts = counts, meta_data = meta_data)
}

# The sample-by-taxon matrices are kept for the first analysis only, which is
# the one whose matrices the vignette displays.
slim = function(fit) {
    fit$bias_correct_log_table = NULL
    fit$struc_absence_prob = NULL
    fit
}

fits = list()

# Primary analysis
toy = read_example("toy")
toy$meta_data$group = factor(toy$meta_data$group,
                             levels = c("control", "treated"))
fits$toy = ancombc3(data = toy$counts, meta_data = toy$meta_data,
                    fix_formula = ~ group + bmi_z,
                    tested = c("group", "bmi_z"),
                    alpha = 0.10, seed = 1, verbose = FALSE)

# Multigroup analysis
mg = read_example("multigroup")
mg$meta_data$dose = factor(mg$meta_data$dose,
                           levels = c("low", "medium", "high"))
fits$multigroup = slim(ancombc3(
    data = mg$counts, meta_data = mg$meta_data,
    fix_formula = ~ dose + age_z, tested = "age_z", group = "dose",
    global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
    alpha = 0.10, seed = 1, verbose = FALSE))

# Trend test with user-specified patterns
patterns = list(
    increasing = matrix(c(1, 0,
                          -1, 1), nrow = 2, byrow = TRUE),
    decreasing = matrix(c(-1, 0,
                          1, -1), nrow = 2, byrow = TRUE),
    umbrella = matrix(c(1, 0,
                        1, -1), nrow = 2, byrow = TRUE))
fits$trend = slim(ancombc3(
    data = mg$counts, meta_data = mg$meta_data,
    fix_formula = ~ dose + age_z, tested = "age_z", group = "dose",
    trend = TRUE,
    trend_control = list(contrast = patterns, node = list(2, 2, 1), B = 1000),
    alpha = 0.10, seed = 1, verbose = FALSE))

# Repeated measures
long = read_example("longitudinal")
long$meta_data$treatment = factor(long$meta_data$treatment,
                                  levels = c("control", "treated"))
fits$longitudinal = slim(ancombc3(
    data = long$counts, meta_data = long$meta_data,
    fix_formula = ~ treatment + visit, tested = c("treatment", "visit"),
    cluster = "subject", alpha = 0.10, seed = 1, verbose = FALSE))

# Real data
data(dietswap, package = "microbiome")
fits$dietswap = slim(ancombc3(
    data = dietswap,
    fix_formula = ~ nationality + timepoint + bmi_group,
    tested = c("nationality", "timepoint"), group = "bmi_group",
    cluster = "subject", global = TRUE, pairwise = TRUE,
    alpha = 0.10, seed = 1, verbose = FALSE))

# The call is not needed by the vignette.
fits = lapply(fits, function(fit) {
    fit$settings$call = NULL
    fit
})

saveRDS(fits, file.path("inst", "extdata", "ancombc3_vignette_fits.rds"),
        compress = "xz")
