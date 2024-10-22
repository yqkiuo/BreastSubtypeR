# BreastSubtypeR

**BreastSubtypeR** is an R package designed to integrate intrinsic molecular subtyping methods for breast cancer, including nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches. It employs standardized input and output formats, offering a unified framework that is highly compatible with other R packages in the gene expression profiling field. 

The package aims to facilitate intrinsic molecular subtyping of breast cancer tumors by providing a comprehensive and user-friendly interface for researchers and clinicians, even for those unfamiliar with R code. The package is accompanied by an **interactive Shiny app** for ease of use.

## Features
- **Nearest-Centroid (NC-based) Methods**: Implements standard NC-based approaches such as PAM50.
- **Single-Sample Predictor (SSP-based) Methods**: Provides SSP-based methods, making it adaptable for various breast cancer subtyping tasks.
- **Flexible Input and Output**: Standardized input/output formats allow for seamless integration with other gene expression profiling packages.
- **Shiny App Interface**: A web-based interface designed for users who prefer not to use R directly.

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
result = BS_Check(data_input = data_input, phenodata = clinic_data, methods = methods, Subtype = TRUE)

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
- Select subtyping methods.
- Visualize the results in real-time.


## Contributing
We welcome contributions to the package. If you find any bugs or have feature requests, feel free to open an issue or submit a pull request.



