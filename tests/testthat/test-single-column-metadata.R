test_that("generic input preserves a single factor metadata column", {
    counts = matrix(seq_len(16), nrow = 4,
                    dimnames = list(paste0("taxon", 1:4), paste0("sample", 1:4)))
    metadata = data.frame(condition = factor(c("reference", "reference", "contrast", "contrast"),
                                             levels = c("reference", "contrast", "unused")),
                          row.names = colnames(counts))
    # Reordering must match samples without dropping the data.frame dimension.
    metadata = metadata[4:1, , drop = FALSE]
    for (taxa_in_rows in c(TRUE, FALSE)) {
        input = if (taxa_in_rows) counts else t(counts)
        checked = suppressWarnings(data_sanity_check(
            data = input, taxa_are_rows = taxa_in_rows, meta_data = metadata,
            fix_formula = "condition", group = "condition", verbose = FALSE))
        expect_s3_class(checked$meta_data, "data.frame")
        expect_identical(dim(checked$meta_data), c(4L, 1L))
        expect_identical(rownames(checked$meta_data), colnames(counts))
        expect_identical(as.character(checked$meta_data$condition),
                         c("reference", "reference", "contrast", "contrast"))
        expect_false("unused" %in% levels(checked$meta_data$condition))
        expect_equal(checked$feature_table, counts)
    }
})

test_that("generic input preserves a single numeric metadata column", {
    counts = matrix(seq_len(16), nrow = 4,
                    dimnames = list(paste0("taxon", 1:4), paste0("sample", 1:4)))
    metadata = data.frame(age = c(20, 30, 40, 50), row.names = colnames(counts))
    checked = data_sanity_check(data = as.data.frame(counts),
                                meta_data = metadata[4:1, , drop = FALSE],
                                fix_formula = "age", verbose = FALSE)
    expect_s3_class(checked$meta_data, "data.frame")
    expect_identical(checked$meta_data, metadata)
})
