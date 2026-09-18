## AUTO routing rules of get_methods(): size gating of the ER/HER2-defined
## cohorts and the handling of cohorts without HER2 information.


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

auto <- function(pheno) {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)
    BreastSubtypeR:::get_methods(pheno)
}

test_that("ER/HER2-defined cohorts are size-gated in AUTO", {
    ## ER+/HER2-: needs ER+ >= 15 and HER2- >= 8
    small <- auto(make_pheno(rep("ER+", 12), rep("HER2-", 12)))
    expect_identical(small$methods, c("AIMS", "sspbc"))
    expect_null(small$samples_ERHER2.icd)

    large <- auto(make_pheno(rep("ER+", 20), rep("HER2-", 20)))
    expect_identical(large$methods, c("ssBC.v2", "AIMS", "sspbc"))
    expect_identical(large$cohort.select, "ERpos")
    expect_length(large$samples_ERHER2.icd, 20L)

    ## ER-/HER2-: needs ER- >= 18 and HER2- >= 9
    small_neg <- auto(make_pheno(rep("ER-", 12), rep("HER2-", 12)))
    expect_identical(small_neg$methods, c("AIMS", "sspbc"))

    large_neg <- auto(make_pheno(rep("ER-", 20), rep("HER2-", 20)))
    expect_identical(large_neg$methods, c("ssBC.v2", "AIMS", "sspbc"))
    expect_identical(large_neg$cohort.select, "ERneg")

    ## HER2+ cohorts: subgroup minimums 8 (ER+) and 9 (ER-)
    small_her2 <- auto(make_pheno(c(rep("ER+", 5), rep("ER-", 5)), rep("HER2+", 10)))
    expect_identical(small_her2$methods, c("AIMS", "sspbc"))
    expect_identical(small_her2$cohort.select, "HER2pos")

    large_her2 <- auto(make_pheno(c(rep("ER+", 10), rep("ER-", 10)), rep("HER2+", 20)))
    expect_identical(large_her2$methods, c("ssBC.v2", "AIMS", "sspbc"))
    expect_identical(large_her2$cohort.select, "HER2pos")
})

test_that("cohorts without HER2 information use the ER-based AUTO rules", {
    ## Mixed ER cohort, HER2 entirely missing: previously routed to the HER2+
    ## branch (AIMS and sspbc only); now the balanced mixed panel.
    no_her2 <- auto(make_pheno(c(rep("ER+", 60), rep("ER-", 40)), rep(NA_character_, 100)))
    expect_identical(no_her2$cohort.select, "mixed")
    expect_identical(
        no_her2$methods,
        c("parker.original", "genefu.scale", "genefu.robust", "ssBC", "ssBC.v2",
          "cIHC", "cIHC.itr", "PCAPAM50", "AIMS", "sspbc")
    )
    expect_length(no_her2$samples_ER.icd, 100L)
    expect_null(no_her2$samples_ERHER2.icd)

    ## Unrecognized HER2 codes count as missing.
    odd_her2 <- auto(make_pheno(c(rep("ER+", 60), rep("ER-", 40)), rep("equivocal", 100)))
    expect_identical(odd_her2$cohort.select, "mixed")

    ## ER+ only without HER2: the ER+ sub-rules need HER2 subgroups, so the
    ## announced fallback applies.
    erp_no_her2 <- auto(make_pheno(rep("ER+", 20), rep(NA_character_, 20)))
    expect_identical(erp_no_her2$cohort.select, "ERpos")
    expect_identical(erp_no_her2$methods, c("AIMS", "sspbc"))
})

test_that("HER2+ cohorts are still detected with partially missing HER2", {
    her2_pos <- auto(make_pheno(
        c(rep("ER+", 10), rep("ER-", 10), "ER+", "ER-"),
        c(rep("HER2+", 20), NA_character_, NA_character_)
    ))
    expect_identical(her2_pos$cohort.select, "HER2pos")
    expect_identical(her2_pos$methods, c("ssBC.v2", "AIMS", "sspbc"))

    ## An equivocal HER2 code is not evaluable either, and is treated the
    ## same way as a missing value.
    her2_pos_equivocal <- auto(make_pheno(
        c(rep("ER+", 10), rep("ER-", 10), "ER+"),
        c(rep("HER2+", 20), "2+")
    ))
    expect_identical(her2_pos_equivocal$cohort.select, "HER2pos")
})

test_that("one HER2- sample ends the HER2+ classification whatever its ER", {
    ## The HER2+ test reads the HER2 column alone, so a HER2- sample counts
    ## whether or not its ER value is known. Before this rule the joint
    ## ER/HER2 counts silently ignored a HER2- sample with missing ER and the
    ## cohort was still routed as HER2+.
    her2_neg_known_er <- auto(make_pheno(
        c(rep("ER+", 10), rep("ER-", 10), "ER+"),
        c(rep("HER2+", 20), "HER2-")
    ))
    expect_false(identical(her2_neg_known_er$cohort.select, "HER2pos"))

    her2_neg_na_er <- auto(make_pheno(
        c(rep("ER+", 10), rep("ER-", 10), NA_character_),
        c(rep("HER2+", 20), "HER2-")
    ))
    expect_false(identical(her2_neg_na_er$cohort.select, "HER2pos"))

    ## Both routes agree, which is the point of the rule.
    expect_identical(
        her2_neg_na_er$cohort.select,
        her2_neg_known_er$cohort.select
    )
})

test_that("AUTO reports the samples that took no part in the HER2+ decision", {
    pheno <- make_pheno(
        c(rep("ER+", 10), rep("ER-", 10), "ER+", "ER-"),
        c(rep("HER2+", 20), NA_character_, NA_character_)
    )
    expect_message(
        BreastSubtypeR:::get_methods(pheno),
        "2 of 22 samples have no evaluable HER2 value"
    )

    ## Nothing is reported when every HER2 value is evaluable.
    complete <- make_pheno(
        c(rep("ER+", 10), rep("ER-", 10)),
        rep("HER2+", 20)
    )
    msgs <- testthat::capture_messages(BreastSubtypeR:::get_methods(complete))
    expect_false(any(grepl("took no part in this decision", msgs, fixed = TRUE)))
})
