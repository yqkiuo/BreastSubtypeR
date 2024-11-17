#' Collection of breast cancer intrinsic subtyping methods
#' 
#' @title Collection of breast cancer intrinsic subtyping methods
#' @description BreastSubtypeR is an R package that grant the access to intrinsic subtyping methods. 
#' 
#' @name BreastSubtypeR
#' 
#' @docType package
#' 
#' @import AIMS
#' @import stringr
#' @importFrom data.table data.table
#' @importFrom data.table set
#' @importFrom rlang call2
#' @importFrom rlang dots_list
#' @noRd 
#' 
NULL

#' Gene ID Mapping
#'
#' @name Mapping
#' @description
#' A function to map gene IDs and preprocess gene expression data for downstream analyses.
#'
#' @param gene_expression_matrix Gene expression matrix with probes in rows and samples in columns.
#' @param featuredata Feature data provided by the user. The table must contain at least three columns: `probe` (ProbeID or TranscriptID), `SYMBOL`, and `ENTREZID`.
#' @param method Method for deduplicating probes in microarray or RNA-seq data. Choose from:
#'   - `"IQR"` for Affymetrix arrays,
#'   - `"mean"` for Agilent arrays,
#'   - `"max"` for RNA-seq data.
#' @param impute Logical. Specify whether to perform K-Nearest Neighbors (KNN) imputation on missing data (`NA` values).
#' @param verbose Logical. If `TRUE`, progress messages will be displayed during execution.
#' 
#' @return A list containing three gene expression matrices:
#'   - `"x_NC.log"`: Preprocessed matrix for nearest-centroid (NC)-based methods.
#'   - `"x_SSP"`: Preprocessed matrix for single-sample predictor (SSP)-based methods.
#'
#' @details
#' If the rows of the gene expression matrix/table represent gene symbols, an additional `SYMBOL` column must be added to the feature table and renamed as `probe`.
#'
#' @examples
#' data("OSLO2EMITOobj")
#' data <- OSLO2EMITO.103.genematrix_noNeg[, clinic.oslo$PatientID]
#' data_input <- Mapping(
#'   gene_expression_matrix = data,
#'   featuredata = anno_feature,
#'   method = "max",
#'   impute = TRUE,
#'   verbose = TRUE
#' )
#'
#' @export

Mapping = function(gene_expression_matrix, featuredata, method = "max", impute = TRUE, verbose = TRUE, ...){
  
  arguments = rlang::dots_list(
    x = gene_expression_matrix,
    y = featuredata,
    method = method,
    impute = impute,
    verbose = verbose, ..., .homonyms = "last"
  )
  
  call = rlang::call2(domapping, !!!arguments)
  res = eval(call)
  
  return(res)
  
}


#' Original Parker Intrinsic Subtyping (BS_parker)
#'
#' @name BS_parker
#' @description
#' This function performs intrinsic subtyping of breast cancer samples using Parker-based methods.
#'
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must be named `"PatientID"`.
#' @param calibration The calibration method to use. Options include:
#'   - `"None"`: No calibration is applied.
#'   - `"Internal"`: Calibration is performed using internal strategies. See the `internal` parameter for details.
#'   - `"External"`: Calibration is performed using external medians. See the `external` parameter for details.
#' @param internal The internal calibration strategy to apply when `calibration = "Internal"`. Options are:
#'   - `"medianCtr"` (default): Median-centered calibration.
#'   - `"meanCtr"`: Mean-centered calibration (aligned with `genefu.scale`).
#'   - `"qCtr"`: Quantile-centered calibration (aligned with `genefu.robust`).
#' @param external Specify the platform name (i.e., column name) for external medians calculated from the training cohort. To use user-provided medians, set this parameter to `"Given.mdns"` and provide the values via the `medians` parameter.
#' @param medians A matrix or table of user-provided medians. Required if `external = "Given.mdns"`. The first column should list 50 genes, and the second column should provide the corresponding median values.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from the `phenodata` table. Required columns include:
#'   - `"T"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#'
#' @return A list containing the intrinsic subtypes assigned using the Parker-based method.
#'
#' @references
#' - Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al. *Supervised risk predictor of breast cancer based on intrinsic subtypes*. Journal of Clinical Oncology. 2009;27(8). https://doi.org/10.1200/JCO.2008.18.1370
#' - Gendoo DMA, Ratanasirigulchai N, Schr"der MS, Par' L, Parker JS, Prat A, et al. *Genefu: An R/Bioconductor package for computation of gene expression-based signatures in breast cancer*. Bioinformatics. 2016;32(7). https://doi.org/10.1093/bioinformatics/btv693
#'
#' @examples
#' data("OSLO2EMITOobj")
#' res <- BS_parker(
#'   gene_expression_matrix = data_input$x_NC.log[],
#'   phenodata = NA,
#'   calibration = "Internal",
#'   internal = "medianCtr",
#'   Subtype = FALSE,
#'   hasClinical = FALSE
#' )
#'
#' @export

BS_parker = function(gene_expression_matrix, phenodata = NA, calibration = "None", internal = NA, external=NA, medians = NA, Subtype = FALSE, hasClinical = FALSE, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    calibration = calibration,
    internal = internal,
    external = external,
    medians = medians,
    Subtype = Subtype,
    hasClinical = hasClinical, ...,.homonyms = "last"
  )
  
  call = rlang::call2(makeCalls.parker, !!!arguments)
  res_parker = eval(call)

  return(res_parker)
  
}

#' Conventional IHC Intrinsic Subtyping (BS_cIHC)
#'
#' @name BS_cIHC
#' @description
#' This function predicts breast cancer intrinsic subtypes using the conventional IHC (cIHC) method, balancing estrogen receptor (ER) subsets for gene centering.
#'
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must contain sample or patient names and be named `"PatientID"`. Additionally, the table must include an `"ER"` column, where estrogen receptor (ER) status is recorded as `"ER+"` or `"ER-"`.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from the `phenodata` table. Required columns include:
#'   - `"T"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#'
#' @return A data frame containing the intrinsic subtypes estimated using the conventional IHC (cIHC) method.
#'
#' @references
#' Ciriello G, Gatza ML, Beck AH, Wilkerson MD, Rhie SK, Pastore A, et al. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer*. Cell. 2015;163(2). https://doi.org/10.1016/j.cell.2015.09.033
#'
#' @examples
#' data("OSLO2EMITOobj")
#' res <- BS_cIHC(
#'   gene_expression_matrix = data_input$x_NC.log,
#'   phenodata = clinic.oslo,
#'   Subtype = FALSE,
#'   hasClinical = FALSE
#' )
#'
#' @export

BS_cIHC = function(gene_expression_matrix, phenodata, Subtype = FALSE ,  hasClinical = FALSE, ...){

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    Subtype = Subtype, 
    hasClinical = hasClinical,..., .homonyms = "last"
  )
  
  
  call = rlang::call2(makeCalls.ihc, !!!arguments)
  res_IHC= eval(call)
  
  
  return(res_IHC)
}


#' Iterative conventional IHC Intrinsic Subtyping (BS_cIHC.itr)
#'
#' @name BS_cIHC.itr
#' @description
#' Predicts breast cancer intrinsic subtypes using an iterative ER-balanced procedure for gene centering. Supports customization of ER subset ratios.
#'
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must contain sample or patient names and be named `"PatientID"`. Additionally, the table must include an `"ER"` column, where estrogen receptor (ER) status is recorded as `"ER+"` or `"ER-"`.
#' @param iterative Integer. Number of iterations for the ER-balanced procedure with a specified ratio.
#' @param ratio Numeric. Specifies the ER+ to ER- ratio for balancing. Options are `1:1` or `54:64` (default). The latter corresponds to the ER ratio used in the UNC230 training cohort.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from the `phenodata` table. Required columns include:
#'   - `"T"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#'
#' @return A list containing:
#' - Intrinsic subtype predictions.
#' - Confidence levels for each subtype.
#' - Percentages of ER+ and ER- subsets across iterations.
#'
#' @references
#' Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, et al. *The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups*. Nature. 2012;486(7403). https://doi.org/10.1038/nature10983
#'
#' @examples
#' data("OSLO2EMITOobj")
#' res <- BS_cIHC.itr(
#'   gene_expression_matrix = data_input$x_NC.log,
#'   phenodata = clinic.oslo,
#'   Subtype = FALSE,
#'   hasClinical = FALSE
#' )
#'
#' @export

BS_cIHC.itr = function(gene_expression_matrix, phenodata, iteration = 100, ratio = 54/64, Subtype = FALSE, hasClinical = FALSE,seed=118, ...){
  
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    iteration = iteration,
    ratio = ratio, 
    Subtype = Subtype, 
    hasClinical = hasClinical,
    seed = seed, .homonyms = "last"
    )
  
  
  call = rlang::call2(makeCalls.ihc.iterative, !!!arguments)
  res_IHC.itr = eval(call)

  return(res_IHC.itr)
}


#' PCA-PAM50 Intrinsic Subtyping (BS_PCAPAM50)
#'
#' @name BS_PCAPAM50
#' @description
#' Performs intrinsic subtyping of breast cancer using the PCA-PAM50 method. This approach integrates Principal Component Analysis (PCA) to enhance the consistency of PAM50-based subtyping.
#'
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must contain sample or patient names and be named `"PatientID"`. Additionally, the table must include an `"ER"` column, where estrogen receptor (ER) status is recorded as `"ER+"` or `"ER-"`.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from the `phenodata` table. Required columns include:
#'   - `"T"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#'
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated by the PCA-PAM50 method.
#'
#' @references
#' Raj-Kumar PK, Liu J, Hooke JA, Kovatich AJ, Kvecher L, Shriver CD, et al. *PCA-PAM50 improves consistency between breast cancer intrinsic and clinical subtyping, reclassifying a subset of luminal A tumors as luminal B.* Sci Rep. 2019;9(1). https://doi.org/10.1038/s41598-019-44339-4
#'
#' @examples
#' data("OSLO2EMITOobj")
#' res <- BS_PCAPAM50(
#'   gene_expression_matrix = data_input$x_NC.log,
#'   phenodata = clinic.oslo,
#'   Subtype = FALSE,
#'   hasClinical = FALSE,
#'   seed = 118
#' )
#'
#' @export

BS_PCAPAM50 = function(gene_expression_matrix, phenodata, Subtype = FALSE, hasClinical =FALSE,seed=118){

  samples = phenodata$PatientID

  if ( "ER"  %in% colnames(phenodata) ) {
    
    ## create IHC column for PCAPAM50
    phenodata$IHC = case_when(
      phenodata$ER == "ER+" ~ "LA",
      phenodata$ER == "ER-" ~ "TN",
      .default = NA
    )
    
    phenodata$ER_status = NA
    
    phenodata$ER_status[which(phenodata$ER == "ER+")] = "pos"
    phenodata$ER_status[which(phenodata$ER == "ER-")] = "neg"
    phenodata = phenodata[order(phenodata$ER_status,decreasing=T),]
  } else {
    stop("Please prepare ER status in clinical table")
  }
  
  gene_expression_matrix = as.matrix( gene_expression_matrix[,phenodata$PatientID])

  ## first step
  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    Subtype = Subtype,
    hasClinical = hasClinical, .homonyms = "last"
    )
  
  call = rlang::call2(makeCalls.PC1ihc, !!!arguments)
  res_PC1IHC = eval(call)
  
  
  ## second step
  if (hasClinical) {
    df.pc1pam = data.frame(PatientID=res_PC1IHC$BS.all$PatientID, PAM50=res_PC1IHC$BS.all$BS,
                           T = phenodata[res_PC1IHC$BS.all$PatientID,]$T, NODE= phenodata[res_PC1IHC$BS.all$PatientID,]$NODE,
                           stringsAsFactors=F)
  } else {
    df.pc1pam = data.frame(PatientID=res_PC1IHC$BS.all$PatientID, PAM50=res_PC1IHC$BS.all$BS,
                           stringsAsFactors=F)
  }
  

  arguments2 = rlang::dots_list(
    mat = gene_expression_matrix,
    df.pam = df.pc1pam,
    Subtype = Subtype,
    hasClinical = hasClinical,
    seed=118, .homonyms = "last"
  )
  
  call = rlang::call2(makeCalls.v1PAM, !!!arguments2)
  res_PCAPAM50 = eval(call)
  
  ## reorder
  res_PCAPAM50$BS.all = res_PCAPAM50$BS.all[na.omit(match( samples, res_PCAPAM50$BS.all$PatientID)),]

  return(res_PCAPAM50)
  
}


 
#' ssBC Intrinsic Subtyping (BS_ssBC)
#'
#' @name BS_ssBC
#' @description
#' Performs intrinsic subtyping of breast cancer using the ssBC (single-sample Breast Cancer) method. This method allows for quantile selection tailored to specific breast cancer cohorts and published guidelines.
#'
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data should be log-transformed.
#' @param phenodata A clinical information table. The first column must contain sample or patient names, named `"PatientID"`. 
#'   - When `"s"` is set as `"ER"`, the table must include an `"ER"` column, where estrogen receptor (ER) status is recorded as `"ER+"` or `"ER-"`.
#'   - When `"s"` is set as `"ER.v2"`, the table must include an `"ER"` column, where estrogen receptor (ER) status is recorded as `"ER+"` or `"ER-"`, and a `"HER2"` column, where human epidermal growth factor receptor 2 (HER2) status is recorded as `"HER2+"` or `"HER2-"`.
#'   - When `"s"` is set as `"TN"` or `"TNBC"`, the table must include a `"TN"` column, recording triple-negative samples as `"TN"`.
#' @param s Character. Options are `"ER"`, `"TN"`, `"ER.v2"`, or `"TNBC"`. Specifies the quantiles to use:
#'   - `"ER"` and `"TN"`: Original quantiles published in *Breast Cancer Research* (2015).
#'   - `"ER.v2"` and `"TNBC"`: Updated quantiles published in *Journal of Clinical Oncology* (2024).
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from the `phenodata` table. Required columns include:
#'   - `"T"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#'
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated by the ssBC method.
#'
#' @references
#' Zhao X, Rodland EA, Tibshirani R, Plevritis S. *Molecular subtyping for clinically defined breast cancer subgroups.* Breast Cancer Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4  
#' Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L, et al. *Survival, pathologic response, and genomics in CALGB 40601 (Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or without lapatinib in HER2-positive breast cancer.* Journal of Clinical Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
#'
#' @examples
#' data("OSLO2EMITOobj")
#' res <- BS_ssBC(
#'   gene_expression_matrix = data_input$x_NC.log,
#'   phenodata = clinic.oslo,
#'   s = "ER.v2",
#'   Subtype = FALSE,
#'   hasClinical = FALSE
#' )
#'
#' @export

BS_ssBC = function(gene_expression_matrix, phenodata, s , Subtype = FALSE, hasClinical =FALSE, ...) {

  arguments = rlang::dots_list(
    mat = gene_expression_matrix,
    df.cln = phenodata,
    s = s,
    Subtype = Subtype,
    hasClinical = hasClinical, ..., .homonyms = "last"
  )
  
  call = rlang::call2(makeCalls.ssBC, !!!arguments)
  res_ssBC = eval(call)
  
  return(res_ssBC)
  
}


#' AIMS Intrinsic Subtyping (BS_AIMS)
#'
#' @name BS_AIMS
#' @description
#' Performs intrinsic subtyping of breast cancer using the AIMS (Absolute assignment of Intrinsic Molecular Subtype) method. This method provides a robust assignment of breast cancer molecular subtypes based on Entrez gene IDs.
#'
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data must be log-transformed.
#' @param EntrezID A vector of Entrez gene IDs corresponding to the genes in the gene expression matrix.
#'
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated by the AIMS method.
#'
#' @references
#' Paquet ER, Hallett MT. *Absolute assignment of breast cancer intrinsic molecular subtype.* J Natl Cancer Inst. 2015;107(1). https://doi.org/10.1093/jnci/dju357
#'
#' @examples
#' # Load required datasets
#' data("BreastSubtypeR")
#' data("OSLO2EMITOobj")
#'
#' # Extract AIMS-specific genes
#' genes <- as.character(
#'   BreastSubtypeR$genes.signature$EntrezGene.ID[
#'     which(BreastSubtypeR$genes.signature$AIMS == "Yes")
#'   ]
#' )
#'
#' # Perform subtyping
#' res <- BS_AIMS(
#'   gene_expression_matrix = data_input$x_SSP[genes, ],
#'   EntrezID = rownames(data_input$x_SSP[genes, ])
#' )
#'
#' @export

BS_AIMS = function(gene_expression_matrix, EntrezID, ...){
  
  ## check dependencies
  suppressMessages(suppressWarnings({
    require(AIMS,quietly = TRUE)
  }))
  
  arguments = rlang::dots_list(
    eset =as.matrix(gene_expression_matrix),
    EntrezID = EntrezID 
  )
  
  call = rlang::call2(AIMS::applyAIMS, !!!arguments)
  
  res_AIMS = eval(call)

}


#' Intrinsic Subtyping using SSPBC (BS_sspbc)
#'
#' @name BS_sspbc
#' @description
#' Performs intrinsic subtyping of breast cancer using the SSPBC (Single Sample Predictor for Breast Cancer) method. This method supports predicting subtypes based on RNA sequencing data and offers flexibility in selecting the prediction model.
#'
#' @param gene_expression_matrix A gene expression matrix with genes in rows and samples in columns. The data must be log-transformed.
#' @param ssp.name Specifies the model to use. Options are "ssp.pam50" (for PAM50-based predictions) or "ssp.subtype" (for predicting four subtypes by excluding the Normal-like subtype).
#'
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated by the SSPBC method.
#'
#' @references
#' Staaf J, H"kkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al. *RNA sequencing-based single sample predictors of molecular subtype and risk of recurrence for clinical assessment of early-stage breast cancer.* NPJ Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' # Load required dataset
#' data("OSLO2EMITOobj")
#'
#' # Perform subtyping with the SSPBC method
#' res <- BS_sspbc(
#'   gene_expression_matrix = as.matrix(data_input$x_SSP),
#'   ssp.name = "ssp.pam50"
#' )
#'
#' @export

BS_sspbc = function(gene_expression_matrix, ssp.name= "ssp.pam50" ,...){
  
  ## check dependencies
  if (!requireNamespace("sspbc", quietly = TRUE)) {
    sspbc_1.0.tar.gz = system.file("extdata", "sspbc","sspbc_1.0.tar.gz", package = "BreastSubtypeR", mustWork = TRUE)
    install.packages(sspbc_1.0.tar.gz, repos = NULL, type="source")
    if (!requireNamespace("sspbc", quietly = TRUE)) {
      stop("Package 'sspbc' should be installed manually.")
    }
  }
  
  #data(sspbc.models)
  suppressMessages(suppressWarnings({
    require(sspbc,quietly = TRUE)
  }))
  
  arguments = rlang::dots_list(
    gex = gene_expression_matrix,
    id = rownames(gene_expression_matrix),
    id.type = "EntrezGene",
    ssp.name= ssp.name, ..., .homonyms = "last"
  )
  
  call = rlang::call2(sspbc::applySSP, !!!arguments)
  res_sspbc = eval(call)

}



#' Consensus Intrinsic Subtyping with Multiple Methods (BS_Multi)
#'
#' @name BS_Multi
#' @description
#' Performs consensus intrinsic subtyping of breast cancer using multiple methods, allowing users to specify or automatically select subtyping approaches based on the test cohort's biomarker distribution.
#'
#' @param data_input The output from the `Mapping()` function, containing processed gene expression data prepared for subtyping analysis.
#' @param phenodata A clinical information table. The first column must be named "PatientID".
#' @param methods A character vector specifying the subtyping methods to be used. Options include:
#'   - "parker.original"
#'   - "genefu.scale"
#'   - "genefu.robust"
#'   - "ssBC"
#'   - "ssBC.v2"
#'   - "cIHC"
#'   - "cIHC.itr"
#'   - "PCAPAM50"
#'   - "AIMS"
#'   - "sspbc"
#'   - "AUTO" (automatically selects methods based on the biomarker distribution of the test cohort).
#' @param Subtype Logical. If `TRUE`, predicts four subtypes by excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, includes clinical information in the analysis (e.g., tumor size in the "T" column and lymph node status in the "NODE" column).
#'
#' @return A list of intrinsic subtypes estimated by the selected methods.
#'
#' @examples
#' # Load required dataset
#' data("OSLO2EMITOobj")
#'
#' # Define methods to use for consensus subtyping
#' methods <- c("parker.original", "genefu.scale", "genefu.robust")
#'
#' # Perform subtyping
#' res.test <- BS_Multi(
#'   data_input = data_input,
#'   phenodata = clinic.oslo,
#'   methods = methods,
#'   Subtype = FALSE,
#'   hasClinical = FALSE
#' )
#'
#' @export

BS_Multi = function(data_input, phenodata, methods = NA, Subtype = FALSE, hasClinical = FALSE,... ){

  ## minor control
  if(length(methods) < 2 ){
    stop("Please select two methods at least")
  } 
  
  if (length(methods[str_detect(methods, pattern =  "parker.original|genefu.scale|genefu.robust|ssBC|ssBC.v2|cIHC|cIHC.itr|PCAPAM50|AIMS|sspbc")] ) < length(methods)){
    stop("Please provide right method names.")
  }

  ## check ER and if methods are feasible
  if(  !("ER" %in% colnames(phenodata)) & ( length(methods[str_detect(methods, pattern =  "ssBC|ssBC.v2|cIHC|cIHC.itr")] ) > 0 ) ) {
    stop("Please do not select any of ssBC, ssBC.v2, cIHC and cIHC.itr")
  } 
  
  if(  !("HER2" %in% colnames(phenodata)) & ( "ssBC.v2" %in% methods  ) ) {
    stop("ssBC.v2 is not supported")
  }

  ## run each method
  results = sapply(methods, function(method) {
    
    ## try NC-based
    if( method == "parker.original") {
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "medianCtr", Subtype = Subtype ,  hasClinical = hasClinical))
    }
    if( method == "genefu.scale"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "meanCtr", Subtype = Subtype , hasClinical = hasClinical))
    }
    if( method == "genefu.robust"){
      print(paste0(method," is running!"))
      return(BS_parker(data_input$x_NC.log, phenodata, calibration = "Internal", internal = "qCtr", Subtype = Subtype , hasClinical = hasClinical))
    }
    if(method == "cIHC"){
      ## try conventional IHC
      print(paste0(method," is running!"))
      return(BS_cIHC(data_input$x_NC.log, phenodata,Subtype = Subtype ,  hasClinical = hasClinical))
    }
    
    if(method == "cIHC.itr"){
      ## try iterative IHC
      print(paste0(method," is running!"))
      return( BS_cIHC.itr(data_input$x_NC.log, phenodata,Subtype = Subtype , hasClinical = hasClinical))
    }
    
    if(method == "PCAPAM50"){
      print(paste0(method," is running!"))
      return(BS_PCAPAM50(data_input$x_NC.log, phenodata, Subtype = Subtype , hasClinical = hasClinical))
    }
    
    if(method == "ssBC"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_NC.log, phenodata, s= "ER",Subtype = Subtype, hasClinical = hasClinical ))
    }
    
    if(method == "ssBC.v2"){
      print(paste0(method," is running!"))
      return(BS_ssBC( data_input$x_NC.log, phenodata, s= "ER.v2", Subtype = Subtype, hasClinical = hasClinical ))
    }
    
    if(method == "AIMS"){
      print(paste0(method," is running!"))
      data("BreastSubtypeR")
      genes.signature = BreastSubtypeR$genes.signature
      
      ## loading library first or model
      genes = as.character( genes.signature$EntrezGene.ID[which( genes.signature$AIMS == "Yes" )])
      res_AIMS = BS_AIMS(data_input$x_SSP[genes,], rownames(data_input$x_SSP[genes,]) )
      res_AIMS$BS.all = data.frame( PatientID = rownames(res_AIMS$cl) ,
                                    BS = res_AIMS$cl[,1])
      
      if( Subtype){
        res_AIMS$BS.all$BS.Subtype = res_AIMS$BS.all$BS
      }
      return(res_AIMS)
    }
    
    
    if(method == "sspbc" ){
      print(paste0(method," is running!"))
      
      res_sspbc = BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_SSP), ssp.name= "ssp.pam50"  )
      
      BS.all = data.frame( PatinetID = rownames(res_sspbc),
                           BS = res_sspbc[,1],
                           row.names = rownames(res_sspbc) )

      if(Subtype) {
        res_sspbc.Subtype = BS_sspbc( gene_expression_matrix = as.matrix(data_input$x_SSP), ssp.name= "ssp.subtype"  )
        BS.all$BS.Subtype = res_sspbc.Subtype[,1]
      } 
      
      return(list(BS.all = BS.all ) )
    }
  
    
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  names(results) = paste0(names(results))
  

  res_subtypes = data.table(row_id = colnames(data_input$x_NC))
  if (Subtype) {
    res_subtypes.Subtype = data.table(row_id = colnames(data_input$x_NC))
  }

  for (method in methods) {
    if (!is.null(results[[method]])) {
      set(res_subtypes, j = method, value = results[[method]]$BS.all$BS)
      if (Subtype) {
        set(res_subtypes.Subtype, j = method, value = results[[method]]$BS.all$BS.Subtype)
      }
    }
  }
  
  ## entropy index
  res_subtypes = as.data.frame(res_subtypes)
  rownames(res_subtypes) = colnames(data_input$x_NC); res_subtypes[,1]= NULL
  
  entropy = apply(res_subtypes, 1, get_entropy  )
  res_subtypes$entropy = entropy

  if(Subtype){
    res_subtypes.Subtype = as.data.frame(res_subtypes.Subtype)
    rownames(res_subtypes.Subtype) = colnames(data_input$x_NC); res_subtypes.Subtype[,1]= NULL
    
    entropy = apply(res_subtypes.Subtype, 1, get_entropy  )
    res_subtypes.Subtype$entropy = entropy
  }
  
  if (Subtype) {
    res = list(res_subtypes = res_subtypes, res_subtypes.Subtype = res_subtypes.Subtype, results = results)
  } else {
    res = list(res_subtypes = res_subtypes, results = results)
  }

}



