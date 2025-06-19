# BreastSubtypeR <a href='https://github.com/yqkiuo/BreastSubtypeR.git'><img src="inst/ShinyBreastSubtypeR/logo.svg" align="right" height="110"/></a>

<!-- badges: start -->

<!-- badges: end -->

## Overview

**BreastSubtypeR** is an R package designed to unify and streamline intrinsic molecular subtyping methods for breast cancer (BC). It integrates both nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches, along with an innovative **AUTO mode** feature (described below). The package utilizes standardized input and output formats, providing a cohesive framework that is fully compatible with other R packages in the gene expression profiling field. Additionally, its core functions are accessible through an **interactive Shiny app**, making it user-friendly for researchers and clinicians with limited R programming experience.

## Features

-   **Comprehensive Intrinsic Subtyping for Breast Cancer**: Integrates multiple published intrinsic subtyping methods, including NC-based approaches like the original PAM50 (Parker et al., J Clin Oncol, 2009) and SSP-based methods like AIMS (Paquet et al., J Natl Cancer Inst, 2015).
-   **Multi-Method Subtyping Functionality**: Simultaneously predicts breast cancer intrinsic subtypes using a variety of validated methods for comparative analysis.
-   **AUTO Mode**: Automatically selects subtyping methods based on the ER/HER2 distribution of the test cohort, ensuring compatibility with the method-specific assumptions and improving accuracy.
-   **Optimized Gene Mapping**: Uses Entrez IDs for gene mapping to ensure the maximum inclusion of genes across subtyping methods.
-   **Streamlined Input/Output**: Standardized input/output formats to ensure smooth integration with other gene expression analysis tools.
-   **Shiny App Interface**: An intuitive web-based graphical user interface (GUI) for local, single-method subtyping analysis, ensuring privacy and data security.

### Single-Method Subtyping Approaches

| **Approach** | **Description** | **Group** | **Citation** |
|-----------------|----------------------|-----------------|-----------------|
| `parker.original` | Original PAM50 by Parker et al., 2009 | NC-based | [Parker et al., 2009](https://doi.org/10.1200/JCO.2008.18.1370) |
| `genefu.scale` | PAM50 implementation as in the genefu R package (scaled version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `genefu.robust` | PAM50 implementation as in the genefu R package (robust version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `cIHC` | Conventional estrogen receptor (ER)-balancing via immunohistochemistry (cIHC) | NC-based | [Ciriello et al., 2015](https://doi.org/10.1016/j.cell.2015.09.033) |
| `cIHC.itr` | Iterative version of cIHC | NC-based | [Curtis et al., 2012](https://doi.org/10.1038/nature10983) |
| `PCAPAM50` | PCA-based iterative PAM50 (ER-balancing using ESR1 gene expression) | NC-based | [Raj-Kumar et al., 2019](https://doi.org/10.1038/s41598-019-44339-4) |
| `ssBC` | Subgroup-specific gene-centering PAM50 | NC-based | [Zhao et al., 2015](https://doi.org/10.1186/s13058-015-0520-4) |
| `ssBC.v2` | Updated subgroup-specific gene-centering PAM50 with refined quantiles | NC-based | [Fernandez-Martinez et al., 2020](https://doi.org/10.1200/JCO.20.01276) |
| `AIMS` | Absolute Intrinsic Molecular Subtyping (AIMS) method | SSP-based | [Paquet & Hallett, 2015](https://doi.org/10.1093/jnci/dju357) |
| `sspbc` | Single-Sample Predictors for Breast Cancer (AIMS adaptation) | SSP-based | [Staaf et al., 2022](https://doi.org/10.1038/s41523-022-00465-3) |

### Multi-Method Subtyping Functionality

| **Approach** | **Description** |
|----------------------|--------------------------------------------------|
| **User-defined Multi-Method** | Allows users to select multiple subtyping methods for comparative analysis. |
| **AUTO Mode Multi-Method** | Automatically selects subtyping methods based on the ER/HER2 distribution of the test cohort. |

## Installation

To install **BreastSubtypeR** from Biocondunctor, run:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BreastSubtypeR")
```

To install **BreastSubtypeR** from GitHub, run:

``` r
# Install devtools package if you haven't already
install.packages("devtools")

# Install BreastSubtypeR from GitHub
devtools::install_github("yqkiuo/BreastSubtypeR")
```

## Getting Started

**Example: User-defined Multi-Method Subtyping**

Here's an example of how to use **BreastSubtypeR** for multi-method breast cancer subtyping. The user manually selects the methods to be used:

``` r
library(BreastSubtypeR)

# Load example data
data("BreastSubtypeRobj")
data("OSLO2EMIT0obj")

# Perform gene mapping before subtyping
data_input <- Mapping( OSLO2EMIT0obj$se_obj, method = "max", impute = TRUE, verbose = FALSE )

# Perform multi-method subtyping
methods <- c("parker.original", "PCAPAM50", "sspbc")
result <- BS_Multi(
    data_input = data_input,
    methods = methods,
    Subtype = FALSE,
    hasClinical = FALSE)

# View the results
head(result$res_subtypes[, 1:min(5, ncol(result$res_subtypes))], 5)

# Visualize results
plot <- Vis_Multi(result$res_subtypes)
plot(plot)
```

**Example: AUTO Mode Multi-Method Subtyping**

Here’s how to use **BreastSubtypeR** for multi-method subtyping with **AUTO** mode. AUTO mode automatically selects methods based on the ER/HER2 distribution of the test cohort:

``` r
library(BreastSubtypeR)

# Load example data
data("BreastSubtypeRobj")
data("OSLO2EMIT0obj")

# Perform gene mapping before subtyping
data_input <- Mapping( OSLO2EMIT0obj$se_obj, method = "max", impute = TRUE, verbose = FALSE )

# Run subtyping with AUTO mode
result <- BS_Multi(
  data_input = data_input,
  methods = "AUTO",
  Subtype = FALSE,
  hasClinical = FALSE
)

# View the results
head(result$res_subtypes[, 1:min(5, ncol(result$res_subtypes))], 5)

# Visualize results
plot <- Vis_Multi(result$res_subtypes)
plot(plot)
```

### Usage

#### Single-Method Subtyping

| **Approach** | **Usage** |
|--------------------|----------------------------------------------------|
| `parker.original` | `BS_parker(calibration = "Internal", internal = "medianCtr", ...)` |
| `genefu.scale` | `BS_parker(calibration = "Internal", internal = "meanCtr", ...)` |
| `genefu.robust` | `BS_parker(calibration = "Internal", internal = "qCtr", ...)` |
| `cIHC` | `BS_cIHC(...)` |
| `cIHC.itr` | `BS_cIHC.itr(...)` |
| `PCAPAM50` | `BS_PCAPAM50(...)` |
| `ssBC` | `BS_ssBC(s = "ER", ...)` |
| `ssBC.v2` | `BS_ssBC(s = "ER.v2", ...)` |
| `AIMS` | `BS_AIMS(...)` |
| `sspbc` | `BS_sspbc(...)` |

#### Multi-Method Subtyping

| **Mode** | **Usage** |
|--------------------|----------------------------------------------------|
| User-defined | `BS_Multi(methods = c("parker.original", "ssBC.v2", "sspbc", ...), ...)` |
| AUTO Mode | `BS_Multi(methods = "AUTO", ...)` |

## Shiny App

For users new to R, we offer an intuitive Shiny app for interactive molecular subtyping.

### Launch the Shiny App

To run iBreastSubtypeR locally with your data, first install and load the package as described above. Afterward, you can interactively access the Shiny app to visualize and analyze your dataset. Here’s an example of how to launch it:

``` r
# Launch iBreastSubtypeR for interactive analysis
library(BreastSubtypeR)
library(tidyverse)
library(shiny)
library(bslib)
iBreastSubtypeR()
```

The Shiny app allows you to:

-   Upload gene expression, clinical, and annotation data.\
-   Perform subtyping using a preferred method.\
-   Visualize results in real-time.\
-   Download results directly to your local machine.

## Contributing

We welcome contributions to the package. If you find any bugs or have feature requests, feel free to open an issue [here](https://github.com/yqkiuo/BreastSubtypeR/issues).

## Citation

If you use **BreastSubtypeR** in your work, please cite:

-   Yang, Q. [aut] & Sifakis, E. G. [cre], *BreastSubtypeR: A Unified R Package for Comprehensive Intrinsic Molecular Subtyping in Breast Cancer Research*. Available at: <https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR>.
-   Additional relevant citations based on the methods you use (refer to the specific methods section for details).
