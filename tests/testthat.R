library(testthat)
test_that("BS_Multi performs subtyping correctly", {
    # Load required dataset
    data("OSLO2EMIT0obj")
    
    # Define methods to use for subtyping
    methods <- c("parker.original", "PCAPAM50", "ssBC","AIMS")
    
    # Perform subtyping
    res.test <- BS_Multi(
        data_input = OSLO2EMIT0obj$data_input,
        pheno = OSLO2EMIT0obj$clinic.oslo,
        methods = methods,
        Subtype = FALSE,
        hasClinical = FALSE
    )
    
    # Validate results
    expect_true(res.test$res_subtypes[,methods], OSLO2EMIT0obj$res$res_subtypes[,methods] )  # Validate results
})


