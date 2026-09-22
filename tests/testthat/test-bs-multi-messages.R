## BS_Multi() reporting when a component fails.

make_minimal_input <- function(pheno) {
    sample_ids <- pheno$PatientID
    expression <- matrix(
        seq_along(sample_ids),
        nrow = 1L,
        dimnames = list("GENE1", sample_ids)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = expression),
        colData = S4Vectors::DataFrame(pheno, row.names = sample_ids)
    )
    list(se_NC = se, se_SSP = se)
}

mock_ssp_calls <- function(se_obj, call) {
    sample_ids <- colnames(SummarizedExperiment::assay(se_obj))
    list(cl = matrix(
        rep(call, length(sample_ids)),
        ncol = 1L,
        dimnames = list(sample_ids, "Subtype")
    ))
}

test_that("a PCAPAM50 failure warning carries the original error message", {
    old_options <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old_options), add = TRUE)

    pheno <- data.frame(
        PatientID = c("P1", "P2", "P3"),
        ER = c("ER+", "ER-", "ER+"),
        HER2 = c("HER2-", "HER2-", "HER2+"),
        stringsAsFactors = FALSE
    )
    input <- make_minimal_input(pheno)

    testthat::local_mocked_bindings(
        BS_PCAPAM50 = function(...) stop("synthetic PCAPAM50 failure"),
        BS_AIMS = function(se_obj, ...) mock_ssp_calls(se_obj, "LumA"),
        BS_sspbc = function(se_obj, ...) mock_ssp_calls(se_obj, "LumA"),
        .package = "BreastSubtypeR"
    )

    expect_warning(
        result <- suppressMessages(BreastSubtypeR::BS_Multi(
            input,
            methods = c("PCAPAM50", "AIMS", "sspbc")
        )),
        "PCAPAM50 failed in this iteration: synthetic PCAPAM50 failure",
        fixed = TRUE
    )
    expect_true(all(is.na(result$res_subtypes$PCAPAM50)))
    expect_identical(result$res_subtypes$AIMS, rep("LumA", 3L))
})
