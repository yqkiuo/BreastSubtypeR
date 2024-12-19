# BreastSubtypeR <a href='https://github.com/yqkiuo/BreastSubtypeR.git'><img src='inst/ShinyBreastSubtypeR/logo.svg' align="right" height="110" /></a>

<!-- badges: start -->
<!-- badges: end -->

## Overview

**BreastSubtypeR** is an R package designed to unify and streamline intrinsic molecular subtyping methods for breast cancer (BC). It integrates both nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches, along with an innovative "AUTO" mode feature (described below). The package utilizes standardized input and output formats, providing a cohesive framework that is fully compatible with other R packages in the gene expression profiling field. Additionally, its core functions are accessible through an **interactive Shiny app**, making it user-friendly for researchers and clinicians with limited R programming experience. 

## Features
- **Comprehensive Intrinsic Subtyping for Breast Cancer**: Integrates multiple published intrinsic subtyping methods, including NC-based approaches like the original PAM50 (Parker et al., J Clin Oncol, 2009) and SSP-based methods.
- **AUTO Mode Feature**: The "AUTO" mode streamlines subtyping by analyzing the distribution of ER and HER2 status and selecting methods that align with the test cohort's characteristics. It ensures that only methods with assumptions compatible with the data are applied, enhancing the accuracy and reliability of the results.
- **Optimized Gene Mapping**: Provides gene mapping through Entrez IDs to maximize the number of genes included in each method.
- **Flexible Input and Output**: Features standardized input/output formats, facilitating seamless integration with other gene expression profiling tools.
- **Shiny App Interface**: Provides a user-friendly, web-based interface designed for users who prefer not to use R directly.

## Methods included

| Methods | Group | Citation |
|-----------------|-----------------|-----------------|
| parker.orginal   | NC-based   | [Parker et al., 2009](https://doi.org/10.1200/JCO.2008.18.1370)   | 
| genefu (scale and robust)  | NC-based   |  [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693)  |
| cIHC-based iterative ER subset (cIHC.itr)    | NC-based    | [Curtis  et al., 2012](https://doi.org/10.1038/nature10983)  |
| conventional IHC (cIHC)    | NC-based    | [Ciriello et al., 2015](https://doi.org/10.1016/j.cell.2015.09.033)   | 
| PCA-PAM50   | NC-based    | [Raj-Kumar et al., 2019](https://doi.org/10.1038/s41598-019-44339-4)    |
| ssBC    | NC-based    | [Zhao et al., 2015](https://doi.org/10.1186/s13058-015-0520-4) |
| ssBC.v2   | NC-based    | [Fernandez-Martinez  et al., 2020](https://doi.org/10.1200/JCO.20.01276)    |
| AIMS   | SSP-based    | [Paquet et al., 2015](https://doi.org/10.1093/jnci/dju357)    |
| sspbc        | SSP-based    | [Staaf et al., 2022](https://doi.org/10.1038/s41523-022-00465-3)   |



## Installation

You can install **BreastSubtypeR** from GitHub:

```R
# Install devtools package if you haven't already
install.packages("devtools")

# Install BreastSubtypeR from GitHub
devtools::install_github("yqkiuo/BreastSubtypeR")
```

## Getting Started

Here is an example of how to use **BreastSubtypeR** for breast cancer subtyping using multiple methods:
```R
library(BreastSubtypeR)

# Example data input: gene expression and clinical data
data("OSLO2MEIT0obj")

## do mapping before subtyping
data = OSLO2EMIT0.103.genematrix_noNeg.subset
data_input = Mapping(gene_expression_matrix = data, featuredata = anno_feature.subset, impute = TRUE, verbose = TRUE )
# Run the subtyping
methods = c("parker.median", "PCAPAM50", "sspbc")
result = BS_Multi(data_input = data_input, phenodata = clinic.oslo, methods = methods, Subtype = TRUE)

# View the results
head(result$res_subtypes)

## visualization
plot = Vis_Multi(res$res_subtypes)
plot(plot)

```

Here is an example of how to use **BreastSubtypeR** with **AUTO** mode feature for breast cancer subtyping. AUTO mode automatically selects methods based on the ER/HER2 distribution of the test cohort:
```R
library(BreastSubtypeR)

# Example data input: gene expression and clinical data
data("OSLO2MEIT0obj")

## do mapping before subtyping
data = OSLO2EMIT0.103.genematrix_noNeg.subset
data_input = Mapping(gene_expression_matrix = data, featuredata = anno_feature.subset, impute = TRUE, verbose = TRUE )
# Run the subtyping with AUTO mode
result = BS_Multi(data_input = data_input, phenodata = clinic.oslo, methods = "AUTO")

# View the results
head(result$res_subtypes)

## visualization
plot = Vis_Multi(res$res_subtypes)
plot(plot)

```

## Method usage

| Methods | Usage |
|-----------------|-----------------|
| parker.orginal   |  BS_parker(calibration = "Internal", internal = "medianCtr" ,...) |
| genefu (scale and robust)  | BS_parker(calibration = "Internal", internal = "meanCtr" ,...); BS_parker(calibration = "Internal", internal = "qCtr",...) |
| cIHC-based iterative ER subset (cIHC.itr)    |  BS_cIHC.itr(...) |
| conventional IHC (cIHC) | BS_cIHC(...) |
| PCA-PAM50   |  BS_PCAPAM50(...) |
| ssBC    |  BS_ssBC(s = "ER", ...) |
| ssBC.v2   | BS_ssBC(s = "ER.v2", ...) |
| AIMS   | BS_AIMS(...) |
| sspbc | BS_Multi(...) |


## Shiny App
For users unfamiliar with R, we provide an easy-to-use Shiny app that allows you to perform molecular subtyping interactively.

### Launch the Shiny App:
```R
library(BreastSubtypeR)
iBreastSubtypeR()
```

The Shiny app allows you to:
- Upload gene expression, clinical, and annotation data.
- Perform subtyping using a preferred method.
- Visualize the results in real-time.
- Download the results locally


## Contributing
We welcome contributions to the package. If you find any bugs or have feature requests, feel free to open an issue or submit a pull request.

## Citation


