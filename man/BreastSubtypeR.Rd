% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BreastSubtypeR_methods.R
\docType{package}
\name{BreastSubtypeR}
\alias{BreastSubtypeR}
\alias{BreastSubtypeR-package}
\title{BreastSubtypeR: A Unified R Package for Intrinsic Molecular Subtyping in Breast Cancer Research}
\description{
\strong{BreastSubtypeR} is an R package designed to unify and
streamline intrinsic molecular subtyping methods for breast cancer (BC).

It integrates both nearest-centroid (NC-based) and single-sample predictor
(SSP-based) approaches, along with an innovative \strong{AUTO mode} feature
(described below).The package utilizes standardized input and output
formats, providing a cohesive framework that is fully compatible with other
R packages in the gene expression profiling field. Additionally, its core
functions are accessible through an \strong{interactive Shiny app}, making it
user-friendly for researchers and clinicians with limited R programming
experience.
\subsection{\strong{Workflow}}{
\enumerate{
\item \strong{Data Input}: Load example data or supply your own gene expression dataset as a SummarizedExperiment object.
\item \strong{Gene Mapping}: Prepare your dataset for subtyping using the \code{\link{Mapping}} function.
\item \strong{Subtyping}: Run multiple subtyping methods (or leverage AUTO mode) with the \code{\link{BS_Multi}} function.
\item \strong{Visualization}: Explore and interpret the subtyping results using the \code{\link{Vis_Multi}} function.
}
}

\subsection{\strong{Key Functions}}{
\itemize{
\item \code{\link{Mapping}}: Prepares gene expression data for subtyping.
\item \code{\link{BS_Multi}}: Executes multiple subtyping methods simultaneously, including an \strong{AUTO} mode for method selection based on cohort characteristics.
\item \code{\link{Vis_Multi}}: Generates visualizations to facilitate interpretation of the subtyping outcomes.
}
}
}
\seealso{
\code{\link{Mapping}}, \code{\link{BS_Multi}}, \code{\link{Vis_Multi}}
}
\author{
\strong{Maintainer}: Qiao Yang \email{yq.kiuo@gmail.com} (\href{https://orcid.org/0000-0002-4098-3246}{ORCID}) [funder]

Authors:
\itemize{
  \item Emmanouil Sifakis \email{emmanouil.sifakis@ki.se}
}

}
