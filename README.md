# BreastSubtypeR <a href='https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR'><img src="inst/ShinyBreastSubtypeR/logo.svg" align="right" height="110"/></a>

<!-- badges: start -->
[![Bioconductor Release](https://bioconductor.org/shields/years-in-bioc/BreastSubtypeR.svg)](https://bioconductor.org/packages/BreastSubtypeR)
[![Bioconductor Devel](https://bioconductor.org/shields/build/devel/bioc/BreastSubtypeR.svg)](https://bioconductor.org/packages/devel/bioc/html/BreastSubtypeR.html)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
<!-- badges: end -->

**BreastSubtypeR** — assumption-aware, multi-method intrinsic molecular subtyping for breast cancer (R / Bioconductor)

**Authors:** Qiao Yang, Emmanouil G. Sifakis\
**Affiliation:** Department of Oncology-Pathology, Karolinska Institutet (Stockholm, Sweden)\
**Paper:** Yang Q., Hartman J., Sifakis E.G. *BreastSubtypeR: A Unified R/Bioconductor Package for Intrinsic Molecular Subtyping in Breast Cancer Research*. **NAR Genomics and Bioinformatics** (2025). **Editor’s Choice**.

------------------------------------------------------------------------

## Overview

**BreastSubtypeR** consolidates established gene-expression–based intrinsic subtyping methods into a single, reproducible R/Bioconductor package and provides a local Shiny app (`iBreastSubtypeR`) for users without programming experience.

Key goals: - reduce method misapplication across heterogeneous research cohorts, - enable direct cross-method benchmarking and exploration of discordant calls, - provide method-specific preprocessing and robust probe-to-gene mapping, - support privacy-preserving local analyses via a Shiny GUI.

------------------------------------------------------------------------

## Main features

-   **Comprehensive Intrinsic Subtyping:** Integrates multiple published intrinsic subtyping algorithms (NC- and SSP-based), including PAM50 variants, AIMS, ssBC, sspbc, and others.
-   **Unified Multi-Method Interface (`BS_Multi`)**: Run many classifiers from one consistent API and compare results side-by-side.
-   **AUTO Mode (cohort-aware selection):** Evaluates cohort diagnostics (e.g., receptor-status distribution, subtype purity, subgroup sizes) and programmatically disables classifiers whose assumptions are likely violated—reducing misclassification in skewed or small cohorts.
-   **Standardised Input & Method-Specific Normalisation:** Supports raw RNA-seq counts, precomputed FPKM, and log₂-normalised microarray/nCounter matrices with automated, method-appropriate transformations.
-   **Optimised Probe/Gene Mapping:** Entrez ID–based mapping and conflict resolution to maximise marker coverage across platforms.
-   **Interactive Shiny App (`iBreastSubtypeR`):** Local GUI that replicates core workflows for non-programmers and preserves data privacy.
-   **Bioconductor distribution & reproducibility:** Unit tests, vignettes and SummarizedExperiment compatibility to support reproducible deployment.

------------------------------------------------------------------------

### Methods included (single-method implementations)

| **Method id** | **Short description** | **Group** | **Reference** |
|------------------|-------------------|------------------|------------------|
| `parker.original` | Original PAM50 by Parker et al., 2009 | NC-based | [Parker et al., 2009](https://doi.org/10.1200/JCO.2008.18.1370) |
| `genefu.scale` | PAM50 implementation as in the genefu R package (scaled version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `genefu.robust` | PAM50 implementation as in the genefu R package (robust version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `cIHC` | Conventional ER-balancing using immunohistochemistry (IHC) | NC-based | [Ciriello et al., 2015](https://doi.org/10.1016/j.cell.2015.09.033) |
| `cIHC.itr` | Iterative version of cIHC | NC-based | [Curtis et al., 2012](https://doi.org/10.1038/nature10983) |
| `PCAPAM50` | Selects IHC-defined ER subsets, then uses Principal component analysis (PCA) to create ESR1 expression-based ER-balancing | NC-based | [Raj-Kumar et al., 2019](https://doi.org/10.1038/s41598-019-44339-4) |
| `ssBC` | Subgroup-specific gene-centering PAM50 | NC-based | [Zhao et al., 2015](https://doi.org/10.1186/s13058-015-0520-4) |
| `ssBC.v2` | Updated subgroup-specific gene-centering PAM50 with refined quantiles | NC-based | [Fernandez-Martinez et al., 2020](https://doi.org/10.1200/JCO.20.01276) |
| `AIMS` | Absolute Intrinsic Molecular Subtyping (AIMS) method | SSP-based | [Paquet & Hallett, 2015](https://doi.org/10.1093/jnci/dju357) |
| `sspbc` | Single-Sample Predictors for Breast Cancer (AIMS adaptation) | SSP-based | [Staaf et al., 2022](https://doi.org/10.1038/s41523-022-00465-3) |

(See the package vignette for implementation details.)

------------------------------------------------------------------------

## Installation

Install the released version from Bioconductor:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BreastSubtypeR")
```

Or install the development version from GitHub:

``` r
# Install devtools package if you haven't already
install.packages("devtools")

# Install BreastSubtypeR from GitHub
devtools::install_github("yqkiuo/BreastSubtypeR")
```

## Quick start examples

Note: this README uses example datasets included in the package. Replace BreastSubtypeRobj / OSLO2EMIT0obj with your own `SummarizedExperiment` object and clinical metadata..

**1) Preprocessing & mapping**

``` r
library(BreastSubtypeR)

# Example data shipped with the package
data("BreastSubtypeRobj") # package-specific object
data("OSLO2EMIT0obj")

# Map probes/ids to Entrez
data_input <- Mapping( OSLO2EMIT0obj$se_obj, RawCounts = FALSE, method = "max", impute = TRUE, verbose = FALSE )
```

**2) Multi-method subtyping (user-defined methods)**

``` r
methods <- c("parker.original", "PCAPAM50", "sspbc")
result <- BS_Multi(
    data_input = data_input,
    methods = methods,
    Subtype = FALSE,
    hasClinical = FALSE
    )

# View per-sample subtype calls (methods x samples)
head(result$res_subtypes[, 1:min(5, ncol(result$res_subtypes))], 5)
```

**3) AUTO mode (cohort-aware selection)**

``` r
result_auto <- BS_Multi(
  data_input = data_input,
  methods = "AUTO",
  Subtype = FALSE,
  hasClinical = FALSE
)

# Visualize subtype calls and inter-method concordance
Vis_Multi(result_auto$res_subtypes)
```

**4) Launch the local Shiny app**

``` r
library(BreastSubtypeR)
library(tidyverse)
library(shiny)
library(bslib)

iBreastSubtypeR() # interactive GUI (local)
```

## Vignette & Documentation

A comprehensive usage guide is included as a vignette with the package—browse it through your R help system or find the rendered version in your documentation files.

For specific functions (like `BS_Multi`, `Mapping`, or `iBreastSubtypeR`), see their help pages (e.g., `?BS_Multi`). The package manual also lists accepted input formats and parameter descriptions for ease of reference.


## Contributing & issues

Contributions and issue reports are welcome. Please open issues or pull requests on the GitHub repository: <https://github.com/yqkiuo/BreastSubtypeR/issues>.

## Citation

When using **BreastSubtypeR** in publications, please cite the package and the paper:

Yang Q., Hartman J., Sifakis E. G. (2025) *BreastSubtypeR: A Unified R/Bioconductor Package for Intrinsic Molecular Subtyping in Breast Cancer Research*. **NAR Genomics and Bioinformatics**. (Editor’s Choice).

You can also use citation("BreastSubtypeR") after installing the package to retrieve the canonical citation(s).

## License

This project is released under the GPL-3 license.
