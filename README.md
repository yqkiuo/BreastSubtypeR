# BreastSubtypeR <a href='https://github.com/yqkiuo/BreastSubtypeR.git'><img src='inst/ShinyBreastSubtypeR/logo.svg' align="right" height="110" /></a>

<!-- badges: start -->
<!-- badges: end -->

## Overview

**BreastSubtypeR** is an R package designed to unify and streamline intrinsic molecular subtyping methods for breast cancer (BC). It integrates both nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches, along with an innovative "AUTO" mode feature (described below). The package utilizes standardized input and output formats, providing a cohesive framework that is fully compatible with other R packages in the gene expression profiling field. Additionally, its core functions are accessible through an **interactive Shiny app**, making it user-friendly for researchers and clinicians with limited R programming experience. 

## Features
- **Comprehensive Intrinsic Subtyping for Breast Cancer**: Integrates multiple published intrinsic subtyping methods, including NC-based approaches like the original PAM50 (Parker et al., J Clin Oncol, 2009) and SSP-based methods like AIMS (Paquet et al., J Natl Cancer Inst, 2015).
- **Multi-Method Subtyping Functionality**: Simultaneously predicts breast cancer intrinsic subtypes using a variety of validated methods for comparative analysis.
- **AUTO Mode Feature**: Evaluates the distribution of ER and HER2 status in the test cohort to automatically select subtyping methods that align with the cohort's characteristics, ensuring compatibility with method-specific assumptions for greater accuracy and reliability.
- **Optimized Gene Mapping**: Optimizes gene mapping using Entrez IDs to maximize the inclusion of genes across subtyping methods.
- **Streamlined Input and Output**: Provides standardized input/output formats to ensure smooth integration with other gene expression analysis tools.
- **User-Friendly Shiny App Interface**: Features a web-based GUI that runs entirely locally, ensuring data privacy with no online sharing, ideal for users who prefer a visual interface over R scripting.  


## Single-Method Subtyping Approaches

| **Approach**        | **Description**         | **Group**  | **Citation**   |
|---------------------|-----------------------------------------------------------------------------|----------------|----------------------------------------|
| `parker.original`   | Original PAM50 by Parker et al., 2009                                                          | NC-based   | [Parker et al., 2009](https://doi.org/10.1200/JCO.2008.18.1370)                                  |
| `genefu.scale`      | PAM50 implementation as in the genefu R package (scaled version)                                | NC-based   | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693)                             |
| `genefu.robust`     | PAM50 implementation as in the genefu R package (robust version)                                | NC-based   | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693)                             |
| `cIHC`              | Conventional ER-balancing using immunohistochemistry (IHC)                                    | NC-based   | [Ciriello et al., 2015](https://doi.org/10.1016/j.cell.2015.09.033)                              |
| `cIHC.itr`          | Iterative cIHC                                                                                 | NC-based   | [Curtis et al., 2012](https://doi.org/10.1038/nature10983)                                       |
| `PCAPAM50`          | PCA-based iterative PAM50 (ER-balancing using ESR1 gene expression)                            | NC-based   | [Raj-Kumar et al., 2019](https://doi.org/10.1038/s41598-019-44339-4)                             |
| `ssBC`              | Subgroup-specific gene-centering PAM50                                                        | NC-based   | [Zhao et al., 2015](https://doi.org/10.1186/s13058-015-0520-4)                                   |
| `ssBC.v2`           | Updated subgroup-specific gene-centering PAM50 (with refined quantiles)                        | NC-based   | [Fernandez-Martinez et al., 2020](https://doi.org/10.1200/JCO.20.01276)                          |
| `AIMS`              | Absolute Intrinsic Molecular Subtyping (AIMS) method                                           | SSP-based  | [Paquet et al., 2015](https://doi.org/10.1093/jnci/dju357)                                       |
| `sspbc`             | Single-Sample Predictors for Breast Cancer (an adaptation of the AIMS method)                  | SSP-based  | [Staaf et al., 2022](https://doi.org/10.1038/s41523-022-00465-3)                                 |

## Multi-Method Subtyping Functionality

| **Approach**                       | **Description**                                                                                              |
|------------------------------------|------------------------------------------------------------------------------------------------------------|
| **User-defined Multi-Method**      | Intrinsic subtyping that integrates multiple user-selected approaches                                        |
| **AUTO Mode Multi-Method**         | Intrinsic subtyping with methods automatically selected based on the ER/HER2 distribution of the test cohort |


## Installation

You can install **BreastSubtypeR** from GitHub:

```R
# Install devtools package if you haven't already
install.packages("devtools")

# Install BreastSubtypeR from GitHub
devtools::install_github("yqkiuo/BreastSubtypeR")
```

## Getting Started

**Example: User-defined Multi-Method Subtyping**
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

**Example: AUTO Mode Multi-Method Subtyping**
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

### Usage

#### Single-Method Subtyping

| **Approach**       | **Usage**                                                               |
|---------------------|-------------------------------------------------------------------------|
| `parker.original`   | `BS_parker(calibration = "Internal", internal = "medianCtr", ...)`     |
| `genefu.scale`      | `BS_parker(calibration = "Internal", internal = "meanCtr", ...)`       |
| `genefu.robust`     | `BS_parker(calibration = "Internal", internal = "qCtr", ...)`          |
| `cIHC`              | `BS_cIHC(...)`                                                        |
| `cIHC.itr`          | `BS_cIHC.itr(...)`                                                    |
| `PCAPAM50`          | `BS_PCAPAM50(...)`                                                    |
| `ssBC`              | `BS_ssBC(s = "ER", ...)`                                              |
| `ssBC.v2`           | `BS_ssBC(s = "ER.v2", ...)`                                           |
| `AIMS`              | `BS_AIMS(...)`                                                        |
| `sspbc`             | `BS_sspbc(...)`                                                       |

#### Multi-Method Subtyping

| **Mode**                         | **Usage**                                                               |
|---------------------|-------------------------------------------------------------------------|
| User-defined                     | `BS_Multi(methods = c("parker.original", "ssBC.v2", "sspbc", ...), ...)`|
| AUTO Mode                        | `BS_Multi(methods = "AUTO", ...)`                                       |


## Shiny App
For users who prefer a graphical interface, **BreastSubtypeR** provides a Shiny app for interactive subtyping.

### Launch the Shiny App
```R
library(BreastSubtypeR)
iBreastSubtypeR()
```

The Shiny app allows you to:

- Upload gene expression, clinical, and annotation data.    
- Perform subtyping using a preferred method.   
- Visualize results in real-time.    
- Download results directly to your local machine.   


## Contributing
We welcome contributions! If you encounter bugs or have feature suggestions, please open an issue [here](https://github.com/yqkiuo/BreastSubtypeR/issues).

## Citation
If you use **BreastSubtypeR** in your work, please cite:

- Yang, Q. [aut] & Sifakis, E. G. [cre], *BreastSubtypeR: Streamlining Molecular Intrinsic Subtyping for Breast Cancer*. Available at: [https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR](https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR).
- Additional relevant citations based on the methods you use (refer to the specific methods section for details).

