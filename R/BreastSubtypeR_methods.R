#' Collection of Breast Cancer Intrinsic Subtyping methods
#'
#' @title Collection of Breast Cancer Intrinsic Subtyping methods
#' @description is an R package designed to unify and streamline intrinsic
#'   molecular subtyping approaches for breast cancer.
#'
#' @name BreastSubtypeR
#'
#' @docType _PACKAGE
#'
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
#' @description A function to map gene IDs and preprocess gene expression data
#' for subsequent analyses.
#'
#' @param gene_expr Gene expression matrix with probes in rows and
#'   samples in columns.
#' @param featuredata Feature data provided by the user. The table must contain
#'   at least three columns: `probe` (ProbeID or TranscriptID), `SYMBOL`, and
#'   `ENTREZID`.
#' @param method Method for deduplicating probes in microarray or RNA-seq data.
#'   Choose one of the following options:
#'   - `"IQR"` for Affymetrix arrays,
#'   - `"mean"` for Agilent/Illumina arrays,
#'   - `"max"` for RNA-seq data.
#' @param impute Logical. Specify whether to perform K-Nearest Neighbors (KNN)
#'   imputation on missing data (`NA` values).
#' @param verbose Logical. If `TRUE`, progress messages will be displayed during
#'   execution.
#'
#' @return A list containing three gene expression matrices:
#'   - `"x_NC.log"`:
#'   Preprocessed matrix for input into nearest-centroid (NC)-based methods.
#'   - `"x_SSP"`: Preprocessed matrix for input into single-sample predictor
#'   (SSP)-based methods.
#'
#' @details If the rows of the gene expression matrix/table represent gene
#' symbols, an additional `SYMBOL` column must be added to the feature table and
#' renamed as `probe`.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' data <- OSLO2EMITO.103.genematrix_noNeg.subset
#' data_input <- Mapping(
#'     gene_expr = data,
#'     featuredata = anno_feature.subset,
#'     method = "max",
#'     impute = TRUE,
#'     verbose = TRUE
#' )
#'
#' @export

Mapping <- function(gene_expr,
    featuredata,
    method = "max",
    impute = TRUE,
    verbose = TRUE) {
    arguments <- rlang::dots_list(
        x = gene_expr,
        y = featuredata,
        method = method,
        impute = impute,
        verbose = verbose,
        .homonyms = "last"
    )

    call <- rlang::call2(domapping, !!!arguments)
    res <- eval(call)

    return(res)
}


#' Original Parker Intrinsic Subtyping (BS_parker)
#'
#' @name BS_parker
#' @description This function predicts breast cancer intrinsic subtypes using
#'   the original Parker et al., 2019 method, as well as some variations of the
#'   original approach.
#'
#' @param gene_expr A gene expression matrix with genes in rows and
#'   samples in columns. The data should be log-transformed.
#' @param pheno A clinical information table. The first column must be named
#'   `"PatientID"`. The default is NA.
#' @param calibration The calibration method to use. Options include:
#'   - `"None"`: No calibration is applied.
#'   - `"Internal"`: Calibration is performed using internal strategies.
#'   See the `internal` parameter for details.
#'   - `"External"`: Calibration is performed using external medians.
#'   See the `external` parameter for details.
#' @param internal The internal calibration strategy to apply when `calibration
#'   = "Internal"`. Options are:
#'   - `"medianCtr"` (default): Median-centered calibration.
#'   - `"meanCtr"`:
#'   Mean-centered calibration (aligned with `genefu.scale` function).
#'   - `"qCtr"`:
#'   Quantile-based calibration (aligned with `genefu.robust` function).
#' @param external Specify the platform name (i.e., column name) for external
#'   medians calculated from the training cohort. To use user-provided medians,
#'   set this parameter to `"Given.mdns"` and provide the values via the
#'   `medians` parameter.
#' @param medians A matrix or table of user-provided medians. Required if
#'   `external = "Given.mdns"`. The first column should list 50 genes, and the
#'   second column should provide the corresponding median values.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by
#'   excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from
#'   the `pheno` table. Required columns include:
#'   - `"TSIZE"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#'
#' @return A list containing the intrinsic subtypes assigned using the
#'   Parker-based method.
#'
#' @references
#' - Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al.
#' *Supervised risk predictor of breast cancer based on intrinsic subtypes*.
#' Journal of Clinical Oncology. 2009;27(8).
#' https://doi.org/10.1200/JCO.2008.18.1370
#'
#' - Gendoo DMA, Ratanasirigulchai N, Schröder MS, Paré L, Parker JS, Prat A, et
#' al. *Genefu: An R/Bioconductor package for computation of gene
#' expression-based signatures in breast cancer*. Bioinformatics. 2016;32(7).
#' https://doi.org/10.1093/bioinformatics/btv693
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_parker(
#'     gene_expr = data_input$x_NC.log,
#'     pheno = NA,
#'     calibration = "Internal",
#'     internal = "medianCtr",
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_parker <- function(gene_expr,
    pheno = NA,
    calibration = "None",
    internal = NA,
    external = NA,
    medians = NA,
    Subtype = FALSE,
    hasClinical = FALSE) {
    if (!is.null(pheno)) {
        rownames(pheno) <- pheno$PatientID
    } else {
        stop("Please provide pheno table as required.")
    }

    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        calibration = calibration,
        internal = internal,
        external = external,
        medians = medians,
        Subtype = Subtype,
        hasClinical = hasClinical,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.parker, !!!arguments)
    res_parker <- eval(call)

    return(res_parker)
}

#' Conventional IHC Intrinsic Subtyping (BS_cIHC)
#'
#' @name BS_cIHC
#' @description This function predicts breast cancer intrinsic subtypes using
#' the conventional IHC (cIHC) method, which employs ER-balanced subsets for
#' gene centering.
#'
#' @param gene_expr A gene expression matrix with genes in rows and
#'   samples in columns. The data should be log-transformed.
#' @param pheno A clinical information table. The first column must contain
#'   sample or patient names and be named `"PatientID"`. Additionally, the table
#'   must include an `"ER"` column, where estrogen receptor (ER) status is
#'   recorded as `"ER+"` or `"ER-"`.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by
#'   excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from
#'   the `pheno` table. Required columns include:
#'   - `"TSIZE"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#' @param seed An integer value is used to set the random seed.
#' @return A data frame containing the intrinsic subtypes estimated using the
#'   conventional IHC (cIHC) method.
#'
#' @references
#' - Ciriello G, Gatza ML, Beck AH, Wilkerson MD, Rhie SK, Pastore A,
#' et al. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer*.
#' Cell. 2015;163(2). https://doi.org/10.1016/j.cell.2015.09.033
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_cIHC(
#'     gene_expr = data_input$x_NC.log,
#'     pheno = clinic.oslo,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_cIHC <- function(gene_expr,
    pheno,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    if (!is.null(pheno)) {
        rownames(pheno) <- pheno$PatientID
    } else {
        stop("Please provide pheno table as required.")
    }

    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )


    call <- rlang::call2(makeCalls_ihc, !!!arguments)
    res_IHC <- eval(call)


    return(res_IHC)
}


#' Iterative conventional IHC Intrinsic Subtyping (BS_cIHC.itr)
#'
#' @name BS_cIHC.itr
#' @description This function predicts breast cancer intrinsic subtypes using an
#' iterative ER-balanced procedure for gene centering. It supports customization
#' of the ER+/ER- ratio.
#'
#' @param gene_expr A gene expression matrix with genes in rows and
#'   samples in columns. The data should be log-transformed.
#' @param pheno A clinical information table. The first column must contain
#'   sample or patient names and be named `"PatientID"`. Additionally, the table
#'   must include an `"ER"` column, where estrogen receptor (ER) status is
#'   recorded as `"ER+"` or `"ER-"`.
#' @param iterative Integer. Number of iterations for the ER-balanced procedure
#'   with the specified ratio.
#' @param ratio Numeric. Specifies the ER+/ER??? ratio for balancing. Options
#'   are `1:1` or `54:64` (default). The latter corresponds to the ER+/ER- ratio
#'   used in the UNC232 training cohort.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by
#'   excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from
#'   the `pheno` table. Required columns include:
#'   - `"TSIZE"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#' @param seed An integer value is used to set the random seed.
#' @return A list containing:
#' - Intrinsic subtype predictions.
#' - Confidence levels for each subtype.
#' - Percentages of ER+ and ER??? subsets across iterations.
#'
#' @references
#' - Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ,
#' et al. *The genomic and transcriptomic architecture of 2,000 breast tumours
#' reveals novel subgroups*. Nature. 2012;486(7403).
#' https://doi.org/10.1038/nature10983
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_cIHC.itr(
#'     gene_expr = data_input$x_NC.log,
#'     pheno = clinic.oslo,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_cIHC.itr <- function(gene_expr,
    pheno,
    iteration = 100,
    ratio = 54 / 64,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    if (!is.null(pheno)) {
        rownames(pheno) <- pheno$PatientID
    } else {
        stop("Please provide pheno table as required.")
    }


    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        iteration = iteration,
        ratio = ratio,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )


    call <- rlang::call2(makeCalls_ihc.iterative, !!!arguments)
    res_IHC.itr <- eval(call)

    return(res_IHC.itr)
}


#' PCA-PAM50 Intrinsic Subtyping (BS_PCAPAM50)
#'
#' @name BS_PCAPAM50
#' @description This function predicts breast cancer intrinsic subtypes using
#' the PCA-PAM50 method. This approach integrates Principal Component Analysis
#' (PCA) to enhance the consistency of PAM50-based subtyping.
#'
#' @param gene_expr A gene expression matrix with genes in rows and
#'   samples in columns. The data should be log-transformed.
#' @param pheno A clinical information table. The first column must contain
#'   sample or patient names and be named `"PatientID"`. Additionally, the table
#'   must include an `"ER"` column, where estrogen receptor (ER) status is
#'   recorded as `"ER+"` or `"ER-"`.
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by
#'   excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from
#'   the `pheno` table. Required columns include:
#'   - `"TSIZE"`: Tumor size.
#'   - `"NODE"`: Lymph node status.
#' @param seed An integer value is used to set the random seed.
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated
#'   by the PCA-PAM50 method.
#'
#' @references
#' - Raj-Kumar PK, Liu J, Hooke JA, Kovatich AJ, Kvecher L, Shriver
#' CD, et al. *PCA-PAM50 improves consistency between breast cancer intrinsic
#' and clinical subtyping, reclassifying a subset of luminal A tumors as luminal
#' B.* Sci Rep. 2019;9(1). https://doi.org/10.1038/s41598-019-44339-4
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_PCAPAM50(
#'     gene_expr = data_input$x_NC.log,
#'     pheno = clinic.oslo,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_PCAPAM50 <- function(gene_expr,
    pheno,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    if (!is.null(gene_expr)) {
        gene_expr <- as.matrix(gene_expr)
    } else {
        stop("Please provide pheno table as required.")
    }

    if (!is.null(pheno)) {
        rownames(pheno) <- pheno$PatientID
    } else {
        stop("Please provide pheno table as required.")
    }

    samples <- pheno$PatientID

    if ("ER" %in% colnames(pheno)) {
        ## create IHC column for PCAPAM50
        pheno$IHC <- dplyr::case_when(pheno$ER == "ER+" ~ "Luminal",
            pheno$ER == "ER-" ~ "non-Luminal",
            .default = NA
        )

        pheno$ER_status <- NA

        pheno$ER_status[which(pheno$ER == "ER+")] <- "pos"
        pheno$ER_status[which(pheno$ER == "ER-")] <- "neg"
        pheno <- pheno[order(pheno$ER_status, decreasing = TRUE), ]
    } else {
        stop("Please prepare ER status in clinical table")
    }


    gene_expr <- as.matrix(gene_expr[, pheno$PatientID])

    ## first step
    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.PC1ihc, !!!arguments)
    res_PC1IHC <- eval(call)


    ## second step
    if (hasClinical) {
        df.pc1pam <- data.frame(
            PatientID = res_PC1IHC$BS.all$PatientID,
            PAM50 = res_PC1IHC$BS.all$BS,
            TSIZE = pheno[res_PC1IHC$BS.all$PatientID, ]$TSIZE,
            NODE = pheno[res_PC1IHC$BS.all$PatientID, ]$NODE,
            stringsAsFactors = FALSE
        )
    } else {
        df.pc1pam <- data.frame(
            PatientID = res_PC1IHC$BS.all$PatientID,
            PAM50 = res_PC1IHC$BS.all$BS,
            stringsAsFactors = FALSE
        )
    }


    arguments2 <- rlang::dots_list(
        mat = gene_expr,
        df.pam = df.pc1pam,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.v1PAM, !!!arguments2)
    res_PCAPAM50 <- eval(call)

    ## reorder
    idx <- na.omit(match(samples, res_PCAPAM50$BS.all$PatientID))
    res_PCAPAM50$BS.all <- res_PCAPAM50$BS.all[idx, ]

    return(res_PCAPAM50)
}



#' ssBC Intrinsic Subtyping (BS_ssBC)
#'
#' @name BS_ssBC
#' @description This function predicts breast cancer intrinsic subtypes using
#'   the subgroup-specific (ssBC) method. This method applies a
#'   subgroup-specific gene-centering approach to a test cohort with a skewed
#'   distribution of clinicopathological characteristics compared to the
#'   original training cohort (e.g., an ER+ selected cohort).
#'
#' @param gene_expr A gene expression matrix with genes in rows and
#'   samples in columns. The data should be log-transformed.
#' @param pheno A clinical information table. The first column must contain
#'   sample or patient names, named `"PatientID"`.
#'   - When `"s"` is set as `"ER"`,
#'   the table must include an `"ER"` column, where estrogen receptor (ER)
#'   status is recorded as `"ER+"` or `"ER-"`.
#'   - When `"s"` is set as `"ER.v2"`,
#'   the table must include an `"ER"` column, where estrogen receptor (ER)
#'   status is recorded as `"ER+"` or `"ER-"`, and a `"HER2"` column, where
#'   human epidermal growth factor receptor 2 (HER2) status is recorded as
#'   `"HER2+"` or `"HER2-"`.
#'   - When `"s"` is set as `"TN"` or `"TN.v2"`, the table
#'   must include a `"TN"` column, recording triple-negative samples as `"TN"`.
#' @param s Character. Options are `"ER"`, `"TN"`, `"ER.v2"`, or `"TN.v2"`.
#'   These specifies which subgroup-specific quantiles to use:
#'   - `"ER"` and `"TN"`: Original subgroup-specific quantiles published in
#'   *Breast Cancer Research* (2015).
#'   - `"ER.v2"` and `"TN.v2"`: Updated
#'   subgroup-specific quantiles published in *Journal of Clinical Oncology*
#'   (2020).
#' @param Subtype Logical. If `TRUE`, the function predicts four subtypes by
#'   excluding the Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from
#'   the `pheno` table. Required columns include:
#'   - `"T"`: Tumor size (0 for
#'   size <= 2cm or 1 for size > 2cm).
#'   - `"NODE"`: Lymph node status (0 for Lymph
#'   node negative or 1 for Lymph node positive).
#'
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated
#'   by the ssBC method.
#'
#' @references
#' - Zhao X, Rodland EA, Tibshirani R, Plevritis S. *Molecular
#'   subtyping for clinically defined breast cancer subgroups.* Breast Cancer
#'   Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4
#' - Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L,
#'   et al. *Survival, pathologic response, and genomics in CALGB 40601
#'   (Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or
#'   without lapatinib in HER2-positive breast cancer.* Journal of Clinical
#'   Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
#'
#' @examples
#' ## ssBC.v2
#' data("OSLO2EMIT0obj")
#' res <- BS_ssBC(
#'     gene_expr = data_input$x_NC.log,
#'     pheno = clinic.oslo,
#'     s = "ER.v2",
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_ssBC <- function(gene_expr,
    pheno,
    s,
    Subtype = FALSE,
    hasClinical = FALSE) {
    if (!is.null(pheno)) {
        rownames(pheno) <- pheno$PatientID
    } else {
        stop("Please provide pheno table as required.")
    }

    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        s = s,
        Subtype = Subtype,
        hasClinical = hasClinical,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.ssBC, !!!arguments)
    res_ssBC <- eval(call)

    return(res_ssBC)
}


#' AIMS Intrinsic Subtyping (BS_AIMS)
#'
#' @name BS_AIMS
#' @description This function predicts breast cancer intrinsic subtypes using
#' the AIMS (Absolute assignment of Intrinsic Molecular Subtype) method.
#'
#' @param gene_expr A gene expression matrix with genes in rows and
#'   samples in columns. Note: The data must be **unlogged** before using this
#'   function.
#' @param EntrezID A vector of Entrez gene IDs corresponding to the genes in the
#'   gene expression matrix.
#'
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated
#'   by the AIMS method.
#'
#' @references
#' - Paquet ER, Hallett MT. *Absolute assignment of breast cancer
#' intrinsic molecular subtype.* J Natl Cancer Inst. 2015;107(1).
#' https://doi.org/10.1093/jnci/dju357
#'
#' @examples
#' # Load required datasets
#' data("OSLO2EMIT0obj")
#' data("BreastSubtypeRobj")
#'
#' # Extract AIMS-specific genes
#' genes <- as.character(
#'     BreastSubtypeRobj$genes.signature$EntrezGene.ID[
#'         which(BreastSubtypeRobj$genes.signature$AIMS == "Yes")
#'     ]
#' )
#'
#' # Perform subtyping
#' res <- BS_AIMS(
#'     gene_expr = data_input$x_SSP[genes, ],
#'     EntrezID = rownames(data_input$x_SSP[genes, ])
#' )
#'
#' @export

BS_AIMS <- function(gene_expr, EntrezID) {
    ## loading datasets
    data("AIMSmodel")
    data("BreastSubtypeRobj")

    arguments <- rlang::dots_list(
        eset = as.matrix(gene_expr),
        EntrezID = EntrezID
    )

    call <- rlang::call2(applyAIMS_AIMS, !!!arguments)

    res_AIMS <- eval(call)
}


#' Intrinsic Subtyping using SSPBC (BS_sspbc)
#'
#' @name BS_sspbc
#' @description This function predicts breast cancer intrinsic subtypes using
#' the SSPBC (Single Sample Predictor for Breast Cancer) method. This method is
#' based on a refined version of AIMS' original methodology, using a large,
#' uniformly accrued population-based cohort (SCAN-B) for training. Performs
#' intrinsic subtyping of breast cancer using the SSPBC (Single Sample Predictor
#' for Breast Cancer) method. This method supports predicting subtypes based on
#' RNA sequencing data and offers flexibility in selecting the prediction model.
#'
#' @param gene_expr A gene expression matrix with genes in rows and
#'   samples in columns. Note: The data must be **unlogged** before using this
#'   function.
#' @param ssp.name Specifies the model to use. Options are "ssp.pam50" (for
#'   PAM50-based predictions) or "ssp.subtype" (for predicting four subtypes by
#'   excluding the Normal-like subtype).
#'
#' @return A vector of intrinsic subtypes assigned to the samples, as estimated
#'   by the SSPBC method.
#'
#' @references
#' - Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I,
#' et al. *RNA sequencing-based single sample predictors of molecular subtype and
#' risk of recurrence for clinical assessment of early-stage breast cancer*. NPJ
#' Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#'
#' # Load required dataset
#' data("OSLO2EMIT0obj")
#'
#' # Perform subtyping with the SSPBC method
#' res <- BS_sspbc(
#'     gene_expr = as.matrix(data_input$x_SSP),
#'     ssp.name = "ssp.pam50"
#' )
#'
#' @export

BS_sspbc <- function(gene_expr, ssp.name = "ssp.pam50") {
    data("sspbc.models")
    data("sspbc.models.fullname")

    arguments <- rlang::dots_list(
        gex = gene_expr,
        id = rownames(gene_expr),
        id.type = "EntrezGene",
        ssp.name = ssp.name,
        .homonyms = "last"
    )

    call <- rlang::call2(applySSP, !!!arguments)
    res_sspbc <- eval(call)
}

#' Intrinsic Subtyping with Multiple Approaches (BS_Multi)
#'
#' @name BS_Multi
#' @description This function predicts breast cancer intrinsic subtypes using
#'   multiple methods, allowing users to either specify the subtyping approaches
#'   directly or enable automatic selection based on the ER/HER2 distribution of
#'   the test cohort.
#'
#' @param data_input The output from the `Mapping()` function, containing
#'   processed gene expression data prepared for subtyping analysis.
#' @param pheno A clinical information table. The first column must be named
#'   "PatientID".
#' @param methods A character vector specifying the subtyping methods to be
#'   used. Options include:
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
#'   - "AUTO" (Automatically selects methods based on the ER/HER2 distribution of the test cohort).
#'
#'   If "AUTO" is selected, it must be the sole value in the vector. Otherwise,
#'   specify at least two other methods for subtyping. An error will occur if
#'   fewer than two methods (excluding "AUTO") are provided.
#' @param Subtype Logical. If `TRUE`, predicts four subtypes by excluding the
#'   Normal-like subtype.
#' @param hasClinical Logical. If `TRUE`, the function uses clinical data from
#'   the `pheno` table. Required columns include:
#'   - `"T"`: Tumor size (0 for
#'   size <= 2cm or 1 for size > 2cm).
#'   - `"NODE"`: Lymph node status (0 for Lymph
#'   node negative or 1 for Lymph node positive).
#'
#' @return A list of intrinsic subtypes estimated by the selected methods.
#'
#' @examples
#' # Load required dataset
#' data("OSLO2EMIT0obj")
#'
#' # Define methods to use for consensus subtyping
#' methods <- c("parker.original", "genefu.scale", "genefu.robust")
#'
#' # Perform subtyping
#' res.test <- BS_Multi(
#'     data_input = data_input,
#'     pheno = clinic.oslo,
#'     methods = methods,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_Multi <- function(data_input,
    pheno,
    methods = NA,
    Subtype = FALSE,
    hasClinical = FALSE) {
    if (length(methods) == 1 && methods[1] == "AUTO") {
        message("Running AUTO mode for subtyping.")
    } else if (length(methods) < 2) {
        stop("Select at least two methods or set method to 'AUTO'.")
    }

    valid_methods <- c(
        "parker.original", "genefu.scale", "genefu.robust",
        "ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50",
        "AIMS", "sspbc", "AUTO"
    )
    pattern <- paste(valid_methods, collapse = "|")
    methods_detected <- methods[str_detect(methods, pattern = pattern)]
    if (length(methods_detected) < length(methods)) {
        stop("Please specify the correct method names.")
    }

    ## check ER and if methods are feasible
    valid_methods <- c("ssBC", "ssBC.v2", "cIHC", "cIHC.itr")
    pattern <- paste(valid_methods, collapse = "|")
    methods_detected <- methods[str_detect(methods, pattern = pattern)]
    if (!("ER" %in% colnames(pheno)) &
        (length(methods_detected) > 0)) {
        stop(
            "Please avoid selecting any of the following methods:
            ssBC, ssBC.v2, cIHC, and cIHC.itr."
        )
    }

    if (!("HER2" %in% colnames(pheno)) &
        ("ssBC.v2" %in% methods)) {
        stop("The method ssBC.v2 is not supported")
    }

    #### AUTO mode
    samples_ER.icd <- NULL
    samples_ERHER2.icd <- NULL
    if (length(methods) == 1 && methods[1] == "AUTO") {
        ## first check ER and HER2 status
        if (!("ER" %in% colnames(pheno)) ||
            !("HER2" %in% colnames(pheno))) {
            stop("Please provide the ER and HER2 status for the \"AUTO\" mode.")
        }

        ## first check sample size
        n_ERpos <- length(pheno$ER[which(pheno$ER == "ER+")])
        n_ERneg <- length(pheno$ER[which(pheno$ER == "ER-")])
        n_ERnegHER2pos <- length(pheno$HER2[which(pheno$HER2 == "HER2+" &
            pheno$ER == "ER-")])
        n_ERnegHER2neg <- length(pheno$HER2[which(pheno$HER2 == "HER2-" &
            pheno$ER == "ER-")])
        n_ERposHER2pos <- length(pheno$HER2[which(pheno$HER2 == "HER2+" &
            pheno$ER == "ER+")])
        n_ERposHER2neg <- length(pheno$HER2[which(pheno$HER2 == "HER2-" &
            pheno$ER == "ER+")])

        ## setting cutoff
        n_ER <- 10
        n_ERHER2 <- 5
        per_ratio <- 0.2
        upper_ratio <- 54 / 64 + (54 / 64) * per_ratio
        lower_ratio <- 54 / 64 - (54 / 64) * per_ratio

        ## main panel
        if (n_ERposHER2neg == 0 &&
            n_ERnegHER2neg == 0) {
            message("A HER2+ cohort has been detected.")

            if (n_ERposHER2pos < n_ERHER2 && n_ERnegHER2pos < n_ERHER2) {
                message("A small HER2+ cohort has been detected.")
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("AIMS", "sspbc")
            } else if (n_ERposHER2pos >= n_ERHER2) {
                message("A ER+/HER2+ cohort has been detected.")
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERnegHER2pos >= n_ERHER2) {
                message("A ER-/HER2+ cohort has been detected.")
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERnegHER2pos == 0 &&
            n_ERposHER2pos == 0 && n_ERposHER2neg == 0) {
            message("A TNBC cohort has been detected.")

            if (!("TN" %in% colnames(pheno))) {
                stop("Provide \"TN\" in pheno for: ssBC(TN) & ssBC.v2 (TN)")
            }
            #### ??? later
            message("Running methods: ssBC (TN), ssBC.v2 (TN), AIMS & sspbc")
            methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
        } else if (n_ERpos < n_ER && n_ERneg < n_ER) {
            message("A small number of ER-/ER+ samples has been detected.")
            message("Running methods: AIMS & sspbc")
            methods <- c("AIMS", "sspbc")
        } else if (n_ERpos > n_ER && n_ERneg < n_ER) {
            if (n_ERposHER2pos > n_ERHER2 && n_ERposHER2neg > n_ERHER2) {
                message("Running methods for ER+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2pos > n_ERHER2) {
                message("Running methods for ER+/HER2+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2neg > n_ERHER2) {
                message("Running methods for ER+/HER2- samples:
                ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERpos < n_ER && n_ERneg > n_ER) {
            if (n_ERnegHER2pos > n_ERHER2 && n_ERnegHER2neg > n_ERHER2) {
                message("Running methods for ER- samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERnegHER2pos > n_ERHER2) {
                message("Running methods for ER-/HER2+ samples:
                ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERnegHER2neg > n_ERHER2) {
                message("Running methods for ER-/HER2- samples:
                ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERpos > n_ER && n_ERneg > n_ER) {
            ## for other NC-based methods
            ratio_ER <- n_ERpos / n_ERneg

            if (ratio_ER > lower_ratio && ratio_ER < upper_ratio) {
                message(
                    "Running methods:
                    parker.original, genefu.scale, genefu.robust,
                    ssBC, ssBC.v2, cIHC, cIHC.itr, PCAPAM50, AIMS & sspbc"
                )
                methods <- c(
                    "parker.original",
                    "genefu.scale",
                    "genefu.robust",
                    "ssBC",
                    "ssBC.v2",
                    "cIHC",
                    "cIHC.itr",
                    "PCAPAM50",
                    "AIMS",
                    "sspbc"
                )
            } else {
                message(
                    "The ER+/ER- ratio in the current dataset
                    differs from that observed in the UNC232 training cohort."
                )
                message("Running methods:
                        ssBC, ssBC.v2, cIHC, cIHC.itr, PCAPAM50, AIMS & sspbc")
                methods <- c(
                    "ssBC",
                    "ssBC.v2",
                    "cIHC",
                    "cIHC.itr",
                    "PCAPAM50",
                    "AIMS",
                    "sspbc"
                )
            }
        }

        ## subsetting samples for ssBC and ssBC.v2
        # Handle ssBC & ssBC.v2 method for imbalanced subtypes
        ERHER2_counts <- c(
            n_ERpos,
            n_ERneg,
            n_ERnegHER2pos,
            n_ERnegHER2neg,
            n_ERposHER2pos,
            n_ERposHER2neg
        )
        names(ERHER2_counts) <- c(
            "ERpos",
            "ERneg",
            "ERnegHER2pos",
            "ERnegHER2neg",
            "ERposHER2pos",
            "ERposHER2neg"
        )

        er_idx <- ERHER2_counts[seq(1, 2)] > n_ER
        samples_ER <- names(ERHER2_counts)[seq(1, 2)][er_idx]
        erher2_idx <- ERHER2_counts[seq(3, 6)] > n_ERHER2
        samples_ERHER2 <- names(ERHER2_counts)[seq(3, 6)][erher2_idx]

        if (length(samples_ER) > 0) {
            message("ssBC for samples: ", paste(samples_ER, collapse = ", "))
            samples_ER.icd <- unlist(lapply(samples_ER, function(subtype) {
                subtype <- str_replace_all(subtype, "pos", "+")
                subtype <- str_replace_all(subtype, "neg", "-")
                ER_status <- subtype
                rownames(pheno)[pheno$ER == ER_status]
            }))
        }


        if (length(samples_ERHER2) > 0) {
            message(
                "ssBC.v2 for samples: ",
                paste(samples_ERHER2, collapse = ", ")
            )

            samples_ERHER2.icd <- unlist(
                lapply(samples_ERHER2, function(subtype) {
                    subtype <- subtype |>
                        str_replace_all("pos", "+") |>
                        str_replace_all("neg", "-")
                    ER_sts <- substr(subtype, 1, 3)
                    HER2_sts <- substr(subtype, 4, 8)
                    rownames(pheno)[pheno$ER == ER_sts & pheno$HER2 == HER2_sts]
                })
            )
        }
    }

    ## run each method
    results <- lapply(methods, function(method) {
        ## try NC-based
        if (method == "parker.original") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$x_NC.log,
                    pheno,
                    calibration = "Internal",
                    internal = "medianCtr",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }
        if (method == "genefu.scale") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$x_NC.log,
                    pheno,
                    calibration = "Internal",
                    internal = "meanCtr",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }
        if (method == "genefu.robust") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$x_NC.log,
                    pheno,
                    calibration = "Internal",
                    internal = "qCtr",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }
        if (method == "cIHC") {
            ## try conventional IHC
            message(method, " is running!")
            return(
                BS_cIHC(
                    data_input$x_NC.log,
                    pheno,
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }

        if (method == "cIHC.itr") {
            ## try iterative IHC
            message(method, " is running!")
            return(
                BS_cIHC.itr(
                    data_input$x_NC.log,
                    pheno,
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }

        if (method == "PCAPAM50") {
            message(method, " is running!")
            return(
                BS_PCAPAM50(
                    data_input$x_NC.log,
                    pheno,
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }

        if (method == "ssBC") {
            message(method, " is running!")

            if ("TN" %in% colnames(pheno)) {
                res_ssBC <- BS_ssBC(
                    data_input$x_NC.log,
                    pheno,
                    s = "ER",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            } else {
                res_ssBC <- BS_ssBC(
                    data_input$x_NC.log,
                    pheno,
                    s = "ER",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            }

            ## removing results for AUTO mode
            if (!is.null(samples_ER.icd) &&
                length(samples_ER.icd) < nrow(pheno)) {
                idx <- !(res_ssBC$BS.all$PatientID %in% samples_ER.icd)
                res_ssBC$BS.all[idx, 2:ncol(res_ssBC$BS.all)] <- NA
            }


            return(res_ssBC)
        }

        if (method == "ssBC.v2") {
            message(method, " is running!")

            if ("TN" %in% colnames(pheno)) {
                res_ssBC.v2 <- BS_ssBC(
                    data_input$x_NC.log,
                    pheno,
                    s = "TN.v2",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            } else {
                res_ssBC.v2 <- BS_ssBC(
                    data_input$x_NC.log,
                    pheno,
                    s = "ER.v2",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            }

            ## removing results for AUTO mode
            if (!is.null(samples_ERHER2.icd) &&
                length(samples_ERHER2.icd) < nrow(pheno)) {
                idx <- !(res_ssBC.v2$BS.all$PatientID %in% samples_ERHER2.icd)
                res_ssBC.v2$BS.all[idx, 2:ncol(res_ssBC.v2$BS.all)] <- NA
            }

            return(res_ssBC.v2)
        }


        if (method == "AIMS") {
            message(method, " is running!")
            data("BreastSubtypeRobj")
            genes.s <- BreastSubtypeRobj$genes.signature

            ## loading library first or model
            genes <- genes.s$EntrezGene.ID[genes.s$AIMS == "Yes"] %>%
                as.character()
            res_AIMS <- BS_AIMS(
                data_input$x_SSP[genes, ],
                rownames(data_input$x_SSP[genes, ])
            )
            res_AIMS$BS.all <- data.frame(
                PatientID = rownames(res_AIMS$cl),
                BS = res_AIMS$cl[, 1]
            )

            if (Subtype) {
                res_AIMS$BS.all$BS.Subtype <- res_AIMS$BS.all$BS
            }
            return(res_AIMS)
        }


        if (method == "sspbc") {
            message(method, " is running!")

            res_sspbc <- BS_sspbc(
                gene_expr = as.matrix(data_input$x_SSP),
                ssp.name = "ssp.pam50"
            )

            BS.all <- data.frame(
                PatinetID = rownames(res_sspbc),
                BS = res_sspbc[, 1],
                row.names = rownames(res_sspbc)
            )

            if (Subtype) {
                res_sspbc.Subtype <- BS_sspbc(
                    gene_expr = as.matrix(data_input$x_SSP),
                    ssp.name = "ssp.subtype"
                )
                BS.all$BS.Subtype <- res_sspbc.Subtype[, 1]
            }

            return(list(BS.all = BS.all))
        }

        stop("Unknown method: ", method)
    })

    names(results) <- methods

    res_subtypes <- data.table(row_id = colnames(data_input$x_NC))
    if (Subtype) {
        res_subtypes.Subtype <- data.table(row_id = colnames(data_input$x_NC))
    }

    for (method in methods) {
        if (!is.null(results[[method]])) {
            set(res_subtypes, j = method, value = results[[method]]$BS.all$BS)
            if (Subtype) {
                set(res_subtypes.Subtype,
                    j = method,
                    value = results[[method]]$BS.all$BS.Subtype
                )
            }
        }
    }

    ## entropy index
    res_subtypes <- as.data.frame(res_subtypes)
    rownames(res_subtypes) <- colnames(data_input$x_NC)
    res_subtypes[, 1] <- NULL

    entropy <- apply(res_subtypes, 1, get_entropy)
    res_subtypes$entropy <- entropy

    if (Subtype) {
        res_subtypes.Subtype <- as.data.frame(res_subtypes.Subtype)
        rownames(res_subtypes.Subtype) <- colnames(data_input$x_NC)
        res_subtypes.Subtype[, 1] <- NULL

        entropy <- apply(res_subtypes.Subtype, 1, get_entropy)
        res_subtypes.Subtype$entropy <- entropy
    }

    if (Subtype) {
        res <- list(
            res_subtypes = res_subtypes,
            res_subtypes.Subtype = res_subtypes.Subtype,
            results = results
        )
    } else {
        res <- list(res_subtypes = res_subtypes, results = results)
    }
}
