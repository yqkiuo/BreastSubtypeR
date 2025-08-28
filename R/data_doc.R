#' @title OSLO2EMIT0obj: Example Dataset for the OSLO2EMIT0 Cohort
#'
#' @description
#' This example dataset is based on the OSLO2EMIT0 cohort
#' described in Staaf et al., 2022. It contains subsetted data for gene
#' expression, clinical information, feature annotations, and example outputs
#' from the `Mapping` and `BS_Multi` functions.
#'
#' @docType data
#' @usage data("OSLO2EMIT0obj")
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{\code{se_obj}}{A SummarizedExperiment object containing a subset of
#'   the log2-transformed, normalized gene expression matrix, clinical information, and gene feature
#'   annotations for the OSLO2EMIT0 cohort.}
#'   \item{\code{data_input}}{Example output from the \code{Mapping} function.}
#'   \item{\code{res}}{Example output from the \code{BS_Multi} function with *AUTO* mode.}
#' }
#'
#'
#' @references
#' - Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I,
#' et al. *RNA sequencing-based single sample predictors of molecular subtype and
#' risk of recurrence for clinical assessment of early-stage breast cancer*. NPJ
#' Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#'
#'
#' @examples
#' library(BreastSubtypeR)
#' data("OSLO2EMIT0obj")
#'
"OSLO2EMIT0obj"


#' @title BreastSubtypeRobj: Data for NC-based Methods
#'
#' @description
#' A list object containing the data required for nearest-centroid (NC)-based
#' molecular subtyping methods. This includes platform medians, centroids, gene
#' signatures, and subgroup quantiles, as well as information from the UNC232
#' training cohort.
#'
#' @docType data
#' @usage data("BreastSubtypeRobj")
#'
#' @format A list with the following elements:
#' \describe{
#'   \item{\code{medians}}{A matrix of medians prepared for nine sequencing platforms.}
#'   \item{\code{centroid}}{The centroids provided by the \code{parker.original} method.}
#'   \item{\code{genes.sig50}}{A data frame of 50 genes used in NC-based methods, including their proliferation information.}
#'   \item{\code{ssBC.subgroupQuantile}}{Subgroup medians prepared by the \code{ssBC} method.}
#'   \item{\code{genes.signature}}{A collection of genes used in both NC-based and SSP-based methods.}
#'   \item{\code{UNC232}}{Data from the UNC232 training cohort.}
#'   \item{\code{platform.UNC232}}{The sequencing platform used for the UNC232 training cohort.}
#' }
#'
#'
#' @references
#' -Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al.
#' *Supervised risk predictor of breast cancer based on intrinsic subtypes*.
#' Journal of Clinical Oncology. 2009;27(8).
#' https://doi.org/10.1200/JCO.2008.18.1370
#' - Zhao X, Rodland EA, Tibshirani R, Plevritis S. *Molecular
#' subtyping for clinically defined breast cancer subgroups.* Breast Cancer
#' Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4
#' - Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L,
#' et al. *Survival, pathologic response, and genomics in CALGB 40601
#' (Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or
#' without lapatinib in HER2-positive breast cancer.* Journal of Clinical
#' Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
#' - Picornell AC, Echavarria I, Alvarez E, López-Tarruella S, Jerez Y, Hoadley K,
#' et al. Breast cancer PAM50 signature: Correlation and concordance between
#' RNA-Seq and digital multiplexed gene expression technologies in a triple
#' negative breast cancer series. BMC Genomics. 2019;20(1).
#' https://doi.org/10.1186/s12864-019-5849-0
#'
#'
#' @examples
#' library(BreastSubtypeR)
#' data("BreastSubtypeRobj")
#'
"BreastSubtypeRobj"


#' @title The AIMS model
#'
#' @description
#' This is the model definition for AIMS. It contains the
#' naive bayes classifier composed of the 100 rules described in Paquet et al.
#' "Absolute assignment of breast cancer intrinsic molecular subtype" (under
#' review at JNCI).
#'
#'
#' @details This is the AIMS model define using 100 simple rules of the form
#' gene A < gene B and combine within a naive bayes classifier within e1071.
#' (Paquet et al. under review JNCI).
#'
#' Briefly, using a suitably large training dataset(~5000 gene breast cancer
#' gene expression profiles), the approach identifies a small set of simple
#' binary rules (~20) that examine the raw expression measurements for pairs of
#' genes from a single breast cancer patient, and only that patient. The binary
#' rules are of the form "if the expression of gene x is greater than gene y,
#' then tend to assign subtype z for that patient". Subtypes could be : Basal,
#' Her2, LumA, LumB, or Normal. The collection of binary rules is combined for a
#' single estimation of a patient subtype via a single probabilistic model using
#' naiveBayes in e1071. In this way, since only expression levels of genes with
#' a single patient is considered, the method represents a promising approach to
#' ablate the instability caused by relativistic approaches (Paquet et al. in
#' review at JNCI).
#'
#' @docType data
#' @usage data("AIMSmodel")
#'
#' @return
#' \item{all.pairs}{The 100 rules in AIMS in the form EntrezID gene A < EntrezID gene B}
#' \item{k}{The selected number of optimal rules. For AIMS we have shown it is 20}
#' \item{one.vs.all.tsp}{The Naive bayes classifier used in combination with the 100 rules}
#' \item{selected.pairs.list}{The list of rules sorted from the best discriminating rule to the least discriminating rules subdivided by subtype}
#'
#'
#' @references
#' - Paquet ER, Hallett MT. *Absolute assignment of breast cancer
#' intrinsic molecular subtype.* J Natl Cancer Inst. 2015;107(1).
#' https://doi.org/10.1093/jnci/dju357
#'
#' @examples
#' library(BreastSubtypeR)
#' data("AIMSmodel")
#'
"AIMSmodel"


#' @title The sspbc models (short name) for the 11 developed predictors
#'
#' @description
#' List with the 11 SSP models from Staaf J. et al. medRxiv 2021.12.03.21267116.
#'
#' @docType data
#' @usage data("sspbc.models")
#'
#' @details
#' List elements are named with short name for respective ssp model.
#' The ssp models in list sspbc.models are the same as in list sspbc.models.fullname.
#'
#' @return
#' \item{sspbc.models}{The collection of 11 ssp models used by sspbc.}
#'
#'
#'
#' #' @references
#' - Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I,
#' et al. *RNA sequencing-based single sample predictors of molecular subtype and
#' risk of recurrence for clinical assessment of early-stage breast cancer*. NPJ
#' Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' ## Load the sspbc.models
#' data("sspbc.models")
#'
"sspbc.models"



#' @title The sspbc models (full name) for the 11 developed predictors
#'
#' @description
#' List with the 11 SSP models from Staaf J. et al. medRxiv 2021.12.03.21267116.
#'
#' @docType data
#' @usage data("sspbc.models.fullname")
#'
#' @details
#' List elements are named with full name for respective ssp model.
#' Note that ssp models in list sspbc.models.fullname are the same as in list sspbc.models.
#'
#' @return
#' \item{sspbc.models.fullname}{The collection of 11 ssp models used by sspbc.}
#'
#'
#' @references
#' - Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I,
#' et al. *RNA sequencing-based single sample predictors of molecular subtype and
#' risk of recurrence for clinical assessment of early-stage breast cancer*. NPJ
#' Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' ## Load the sspbc.models
#' data(sspbc.models.fullname)
#'
"sspbc.models.fullname"


#' @title Gene annotation table.
#'
#' @description Annotation table for the GENCODE Human Release 27 genes
#' (Gene.ID) included in the StringTie target when summarizing gene expression.
#' The annotation are from from GENCODE Human Release 27 metadata files and
#' include HGNC, EntrezGene and RefSeq.
#'
#'
#' @docType data
#' @usage data("Gene.ID.ann")
#'
#' @details Annotation table used by the applySSP to translate gene identifiers
#' as needed before classification with provided ssp models.
#'
#' @return \item{Gene.ID.ann}{Annotation table for GENCODE Human Release 27
#' genes.}
#'
#'
#'
#' @references
#' - Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I,
#' et al. *RNA sequencing-based single sample predictors of molecular subtype and
#' risk of recurrence for clinical assessment of early-stage breast cancer*. NPJ
#' Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' ## Load the Gene.ID.ann
#' data(Gene.ID.ann)
#'
"Gene.ID.ann"
