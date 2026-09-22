## Regression tests for the ROR risk-group columns produced by RORgroup().
## The synthetic 'out' object isolates the combined (ROR-C) and combined +
## proliferation (ROR-PC) models: only the LumA correlation and the
## proliferation genes vary, so the expected risk groups follow directly from
## the published coefficients and the thresholds hard-coded in RORgroup()
## (ROR-C: low < -0.1, high > 0.2; ROR-PC: low < -0.2, high > 0.2).

make_ror_input <- function() {
    ids <- c("S1", "S2", "S3")
    genes <- c("ESR1", "ERBB2", "CCNB1", "UBE2C", "MKI67")
    prolif <- c(-1, 1, 0)
    testData <- rbind(
        ESR1 = c(0.5, 0.4, 0.3),
        ERBB2 = c(-0.5, -0.4, -0.3),
        CCNB1 = prolif,
        UBE2C = prolif,
        MKI67 = prolif
    )
    colnames(testData) <- ids
    centroids <- matrix(
        c(
            -1.0, 0.2, 1.2, 0.8, 0.1,
            0.1, 1.1, -0.4, 0.2, -0.2,
            0.9, 0.4, -0.8, 0.6, -0.1,
            0.8, 0.3, -0.6, 0.9, 0.0,
            0.7, 0.5, -0.7, 0.7, 0.2
        ),
        nrow = length(genes),
        byrow = TRUE,
        dimnames = list(genes, c("Basal", "Her2", "LumA", "LumB", "Normal"))
    )
    ## Basal, Her2 and LumB correlations are zero, so the combined score is
    ## -0.2608388 * LumA and the combined + proliferation score is
    ## -0.090436516 * LumA + 0.327259375 * prolif.
    distances <- matrix(
        0,
        nrow = length(ids),
        ncol = 5,
        dimnames = list(ids, colnames(centroids))
    )
    distances[, "LumA"] <- c(-0.9, 0.9, 0)
    distances[, "Normal"] <- c(0.1, 0.2, 0.3)
    predictions <- setNames(c("Basal", "LumA", "Normal"), ids)
    out <- list(
        predictions = predictions,
        testData = testData,
        distances = distances,
        dist.RORSubtype = distances[, c("Basal", "Her2", "LumA", "LumB")],
        centroids = centroids
    )
    df.cln <- data.frame(
        PatientID = ids,
        TSIZE = c(0, 0, 0),
        NODE = c(0, 0, 0),
        stringsAsFactors = FALSE
    )
    list(out = out, df.cln = df.cln)
}

test_that("ROR-C group follows the ROR-C score, not the ROR-PC score", {
    input <- make_ror_input()
    res <- suppressMessages(BreastSubtypeR:::RORgroup(
        input$out,
        input$df.cln,
        Subtype = FALSE,
        hasClinical = TRUE
    ))

    ## S1: combined = 0.2348 (high), combined + prolif = -0.2459 (low)
    ## S2: combined = -0.2348 (low), combined + prolif = 0.2459 (high)
    ## S3: both scores are 0 (med)
    expect_identical(
        unname(as.character(res[["ROR-C Group (Subtype + Clinic)"]])),
        c("high", "low", "med")
    )
    expect_identical(
        unname(as.character(res[["ROR-PC Group (Subtype + Clinic + Prolif)"]])),
        c("low", "high", "med")
    )

    ## The group must be reproducible from the reported ROR-C score. The score
    ## is reported on the 100 * (x + 0.35) / 0.85 scale used for all ROR columns.
    combined <- res[["ROR-C (Subtype + Clinic)"]] * 0.85 / 100 - 0.35
    expected <- ifelse(combined > 0.2, "high", ifelse(combined < -0.1, "low", "med"))
    expect_identical(
        unname(as.character(res[["ROR-C Group (Subtype + Clinic)"]])),
        unname(expected)
    )
})

test_that("ROR-C and ROR-PC groups are unchanged without clinical data", {
    input <- make_ror_input()
    res <- suppressMessages(BreastSubtypeR:::RORgroup(
        input$out,
        input$df.cln,
        Subtype = FALSE,
        hasClinical = FALSE
    ))
    expect_false("ROR-C Group (Subtype + Clinic)" %in% colnames(res))
    expect_identical(
        unname(as.character(res[["ROR-S Group (Subtype Only)"]])),
        c("high", "low", "med")
    )
})
