## Regression tests for the AUTO-mode ssBC / ssBC.v2 sample subsets returned by
## get_methods(). Samples whose ER (or HER2) status is missing must not appear
## in the subsets: a missing value used to produce an NA sample name, which
## either broke the SummarizedExperiment subsetting in BS_Multi() or, when the
## padded vector was as long as the cohort, silently disabled the subsetting.

make_pheno <- function(er, her2) {
    pheno <- data.frame(
        PatientID = sprintf("S%03d", seq_along(er)),
        ER = er,
        HER2 = her2,
        stringsAsFactors = FALSE
    )
    rownames(pheno) <- pheno$PatientID
    pheno
}

expect_valid_subset <- function(subset, pheno, expected_ids) {
    expect_false(anyNA(subset))
    expect_true(all(subset %in% rownames(pheno)))
    expect_setequal(subset, expected_ids)
}

test_that("missing ER values do not leak NA sample names into the ssBC subsets", {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)

    pheno <- make_pheno(
        er = c(rep("ER+", 50), rep("ER-", 45), rep(NA_character_, 5)),
        her2 = rep(c("HER2-", "HER2+"), 50)
    )
    out <- BreastSubtypeR:::get_methods(pheno)

    expect_true(all(c("ssBC", "ssBC.v2") %in% out$methods))
    known_er <- rownames(pheno)[!is.na(pheno$ER)]
    expect_valid_subset(out$samples_ER.icd, pheno, known_er)
    expect_valid_subset(out$samples_ERHER2.icd, pheno, known_er)
    expect_length(out$samples_ER.icd, 95L)
})

test_that("missing HER2 values do not leak NA sample names into the ssBC.v2 subset", {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)

    pheno <- make_pheno(
        er = c(rep("ER+", 10), rep("ER-", 80), rep(NA_character_, 10)),
        her2 = c(rep("HER2-", 45), rep("HER2+", 45), rep(NA_character_, 10))
    )
    out <- BreastSubtypeR:::get_methods(pheno)

    expect_identical(out$cohort.select, "ERneg")
    expect_identical(out$methods, c("ssBC", "ssBC.v2", "AIMS", "sspbc"))
    ## ssBC: only the ER- group reaches its threshold (ER+ has 10 < 15).
    er_neg <- rownames(pheno)[which(pheno$ER == "ER-")]
    expect_valid_subset(out$samples_ER.icd, pheno, er_neg)
    ## ssBC.v2: every ER/HER2 subgroup with known status is large enough
    ## (ER+/HER2- 10 >= 8, ER-/HER2+ 45 >= 9, ER-/HER2- 35 >= 9).
    known_both <- rownames(pheno)[!is.na(pheno$ER) & !is.na(pheno$HER2)]
    expect_valid_subset(out$samples_ERHER2.icd, pheno, known_both)
    expect_length(out$samples_ERHER2.icd, 90L)
    ## The subset is shorter than the cohort, so BS_Multi() will subset the
    ## SummarizedExperiment with it; every name must be a real sample.
    expect_lt(length(out$samples_ER.icd), nrow(pheno))
})
