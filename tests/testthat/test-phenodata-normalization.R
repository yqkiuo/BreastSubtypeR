## Regression tests for .normalize_er_her2_tn() with factor columns. The
## ifelse() fallback used to return the integer codes of unmatched factor
## levels (for example "Unknown" -> "3"), so a phenotype table with factor
## ER/HER2/TN columns was normalized differently from the same table with
## character columns.

test_that("factor and character phenodata normalize identically", {
    chr <- data.frame(
        PatientID = paste0("P", 1:5),
        ER = c("ER+", "positive", "ER-", "Unknown", NA),
        HER2 = c("HER2+", "2+", "neg", "HER2-", "equivocal"),
        TN = c("TN", "no", "maybe", "nonTN", "1"),
        stringsAsFactors = FALSE
    )
    fct <- chr
    fct$ER <- factor(fct$ER)
    fct$HER2 <- factor(fct$HER2)
    fct$TN <- factor(fct$TN)

    from_chr <- suppressWarnings(BreastSubtypeR:::.normalize_er_her2_tn(chr))
    from_fct <- suppressWarnings(BreastSubtypeR:::.normalize_er_her2_tn(fct))

    expect_identical(from_fct, from_chr)
    expect_identical(from_fct$ER, c("ER+", "ER+", "ER-", "Unknown", NA))
    expect_identical(from_fct$HER2, c("HER2+", "2+", "HER2-", "HER2-", "equivocal"))
    expect_identical(from_fct$TN, c("TN", "nonTN", "maybe", "nonTN", "TN"))
})

test_that("already canonical factor columns produce no coercion warning", {
    fct <- data.frame(
        PatientID = paste0("P", 1:4),
        ER = factor(c("ER+", "ER-", "ER+", "ER-")),
        HER2 = factor(c("HER2+", "HER2-", "HER2-", "HER2+"))
    )
    expect_no_warning(res <- BreastSubtypeR:::.normalize_er_her2_tn(fct))
    expect_identical(res$ER, c("ER+", "ER-", "ER+", "ER-"))
    expect_identical(res$HER2, c("HER2+", "HER2-", "HER2-", "HER2+"))
})

test_that("get_methods() gives the same AUTO decision for factor phenodata", {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)

    chr <- data.frame(
        PatientID = sprintf("S%03d", 1:100),
        ER = c(rep("ER+", 55), rep("ER-", 45)),
        HER2 = rep(c("HER2-", "HER2+"), 50),
        stringsAsFactors = FALSE
    )
    rownames(chr) <- chr$PatientID
    fct <- chr
    fct$ER <- factor(fct$ER)
    fct$HER2 <- factor(fct$HER2)

    out_chr <- BreastSubtypeR:::get_methods(chr)
    out_fct <- BreastSubtypeR:::get_methods(fct)
    expect_identical(out_fct, out_chr)
})
