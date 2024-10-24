# BreastSubtypeR

**BreastSubtypeR** is an R package designed to integrate multiple intrinsic molecular subtyping methods for breast cancer (BC), including nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches, along with a new consensus approach. It employs standardized input and output formats, providing a unified framework that is highly compatible with other R packages in the gene expression profiling field. The core functions of the package can also be used through an **interactive Shiny app**, making it user-friendly for those unfamiliar with R. 

## Features
- **Comprehensive Intrinsic Subtyping for Breast Cancer**: Integrates multiple published intrinsic subtyping methods, including NC-based approaches like the original PAM50 (Parker et al., J Clin Oncol, 2009) and SSP-based methods.
- **Consensus Intrinsic Subtyping Approach**: Introduces a new consensus approach, running multiple intrinsic subtyping methods and using a majority voting scheme to determine the final subtype.
- **Optimized Gene Mapping**: Provides gene mapping through Entrez IDs to maximize the number of genes included in each method.
- **Flexible Input and Output**: Features standardized input/output formats, facilitating seamless integration with other gene expression profiling tools.
- **Shiny App Interface**: Provides a user-friendly, web-based interface designed for users who prefer not to use R directly.

## Methods included

| Methods | Group | Citation |
|-----------------|-----------------|-----------------|
| parker.orginal   | NC-based   | [Parker et al., 2009](https://doi.org/10.1200/JCO.2008.18.1370)   |
| genefu           | NC-based   |  [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693)  |
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

Here is an example of how to use **BreastSubtypeR** for breast cancer subtyping:
```R
library(BreastSubtypeR)

# Example data input: gene expression and clinical data
data("BreastSubtypeR")
data("OSLO2MEITOobj")

## do mapping before subtyping
data = OSLO2EMIT0.103.genematrix_noNeg[,clinic.oslo$PatientID]
data_input = Mapping(gene_expression_matrix = data, featuredata = anno_feature, impute = TRUE, verbose = TRUE )
# Run the subtyping
methods = c("parker.median", "PCAPAM50", "sspbc")
result = BS_Multi(data_input = data_input, phenodata = clinic_data, methods = methods, Subtype = TRUE)

# View the results
head(result$res_subtypes)

## visualization
plot = Vis_consensus(res$res_subtypes)
plot(plot)

```

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



