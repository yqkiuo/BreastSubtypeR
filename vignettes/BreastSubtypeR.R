## ----eval = FALSE-------------------------------------------------------------
# if (!require("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# 
# BiocManager::install("BreastSubtypeR")

## ----eval = FALSE-------------------------------------------------------------
# # Install devtools package if you haven't already
# install.packages("devtools")
# 
# # Install BreastSubtypeR from GitHub
# devtools::install_github("yqkiuo/BreastSubtypeR")

## -----------------------------------------------------------------------------
library(BreastSubtypeR)

# Load example data
data("BreastSubtypeRobj")
data("OSLO2EMIT0obj")

# Perform gene mapping before subtyping
data_input <- Mapping(OSLO2EMIT0obj$se_obj, method = "max", impute = TRUE, verbose = FALSE)

# Perform multi-method subtyping
methods <- c("parker.original", "PCAPAM50", "sspbc")
result <- BS_Multi(
    data_input = data_input,
    methods = methods,
    Subtype = FALSE,
    hasClinical = FALSE
)

# View the results
head(result$res_subtypes[, 1:min(5, ncol(result$res_subtypes))], 5)

# Visualize results
plot <- Vis_Multi(result$res_subtypes)
plot(plot)

