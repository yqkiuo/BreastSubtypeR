#' @import stringr
#' @import e1071
#' @import Biobase
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom data.table data.table
#' @importFrom data.table set
#' @importFrom rlang call2
#' @importFrom rlang dots_list
#' @importFrom utils data installed.packages read.delim write.table
#' @importFrom graphics barplot mtext par
#' @importFrom grDevices dev.off pdf
#'
NULL

#' BreastSubtypeR: A Unified R Package for Intrinsic Molecular Subtyping in Breast Cancer Research
#'
#'
#' @name BreastSubtypeR
#' @aliases NULL
#' @docType package
#'
#' @description **BreastSubtypeR** is an R package designed to unify and
#'   streamline intrinsic molecular subtyping methods for breast cancer (BC).
#'
#'
#'   It integrates both nearest-centroid (NC-based) and single-sample predictor
#'   (SSP-based) approaches, along with an innovative **AUTO mode** feature
#'   (described below).The package utilizes standardized input and output
#'   formats, providing a cohesive framework that is fully compatible with other
#'   R packages in the gene expression profiling field. Additionally, its core
#'   functions are accessible through an **interactive Shiny app**, making it
#'   user-friendly for researchers and clinicians with limited R programming
#'   experience.
#'
#' ## **Workflow**
#' 1. **Data Input**: Load example data or supply your own gene expression dataset as a SummarizedExperiment object.
#' 2. **Gene Mapping**: Prepare your dataset for subtyping using the \code{\link{Mapping}} function.
#' 3. **Subtyping**: Run multiple subtyping methods (or leverage AUTO mode) with the \code{\link{BS_Multi}} function.
#' 4. **Visualization**: Explore and interpret the subtyping results using the \code{\link{Vis_Multi}} function.
#'
#'
#' ## **Key Functions**
#' - \code{\link{Mapping}}: Prepares gene expression data for subtyping.
#' - \code{\link{BS_Multi}}: Executes multiple subtyping methods simultaneously, including an **AUTO** mode for method selection based on cohort characteristics.
#' - \code{\link{Vis_Multi}}: Generates visualizations to facilitate interpretation of the subtyping outcomes.
#'
#'
#' @seealso \code{\link{Mapping}}, \code{\link{BS_Multi}}, \code{\link{Vis_Multi}}
#'
#'
"_PACKAGE"


#' Gene ID Mapping
#'
#' @name Mapping
#' @description Preprocesses and maps gene expression input to prepare for
#'   downstream subtyping workflows (both NC-based and SSP-based).
#'
#' @param se_obj A `SummarizedExperiment` object with:
#'   - **Assay data**: A gene expression matrix.
#'     - If `RawCounts = FALSE`: `assay()` must contain log₂-transformed, normalized expression (e.g., pre-normalized microarray, nCounter, or FPKM RNAseq data).
#'     - If `RawCounts = TRUE`: `assay()` contains raw RNA-seq read counts; see `RawCounts` parameter.
#'   - **Row metadata**: Must include:
#'     - `"probe"`: feature identifiers (e.g., Gene Symbols or Probe IDs)
#'     - `"ENTREZID"`: matching Entrez Gene IDs
#'     - When using gene symbols as row names, the feature data should have an additional `SYMBOL` column and renamed as`probe`
#'   - **Column metadata** (optional): Sample descriptions via `colData()`.
#'   
#' @param RawCounts Logical. If `TRUE`, indicates that `assay()` holds raw RNA-seq counts.
#'   You must supply gene lengths (e.g., in base pairs) in `rowData()` (column name e.g. `"Length"`),
#'   which will be used to compute:
#'   - NC-based methods: log₂ CPM using upper-quartile normalization
#'   - SSP-based methods: linear FPKM (not log-transformed)
#'
#' @param method A string specifying the method for resolving duplicate probes
#'   in microarray or RNA-seq data. Options include:
#'   - `"iqr"`: Selects the probe with the highest interquartile range (IQR), typically used for short-oligo arrays (e.g., Affymetrix).
#'   - `"mean"`: Chooses the probe with the highest average expression, commonly used for long-oligo arrays (e.g., Agilent, Illumina).
#'   - `"max"`: Retains the probe with the highest expression value, often used for RNA-seq data.
#'   - `"stdev"`: Selects the probe with the highest standard deviation.
#'   - `"median"`: Chooses the probe with the highest median expression value.
#'   
#' @param impute Logical. If `TRUE`, applies K-Nearest Neighbors (KNN) imputation for missing values.
#'   
#' @param verbose Logical. If `TRUE`, prints progress messages during execution.
#'
#' @return A named list with:
#'   \item{x_NC}{`SummarizedExperiment` with log₂-transformed data ready for NC-based methods, plus clinical metadata.}
#'   \item{x_SSP}{`SummarizedExperiment` with exponential-transformed data ready for SSP-based methods, plus clinical metadata.}
#'
#' @details
#' `Mapping()` can handle multiple input types seamlessly:
#' - **Raw RNA-seq counts** (need gene lengths) → normalized to CPM or FPKM as appropriate
#' - **Pre-computed log₂ FPKM** → used directly for NC, back-transformed for SSP
#' - **Microarray or nCounter data** → expected to be log₂-normalized already, ingested directly or back-transformed for SSP
#' This enables users to supply a single input type and let BreastSubtypeR manage the normalization pipeline automatically.
#' 
#'
#' @examples
#' \donttest{
#' library(BreastSubtypeR)
#' # Using raw counts (with gene lengths in rowData)
#' se_obj <- SummarizedExperiment(
#'   assays = list(counts = raw_counts_mat),
#'   rowData = DataFrame(
#'     probe = rownames(raw_counts_mat),
#'     ENTREZID = entrez_ids,
#'     Length = gene_lengths
#'   )
#' )
#' res <- Mapping(se_obj, RawCounts = TRUE)
#'
#' # Using pre-normalized log2 FPKM
#' se_obj_fpkm <- SummarizedExperiment(
#'   assays = list(expr = log2_fpkm_mat),
#'   rowData = DataFrame(probe = rownames(log2_fpkm_mat), ENTREZID = entrez_ids)
#' )
#' res <- Mapping(se_obj_fpkm, RawCounts = FALSE)
#' }
#' @export

Mapping <- function(
        se_obj,
        RawCounts = FALSE,
        method = c("max", "mean", "median", "iqr", "stdev"),
        impute = TRUE,
        verbose = TRUE) {
    method <- match.arg(method)

    arguments <- rlang::dots_list(
        se_obj = se_obj,
        RawCounts = RawCounts,
        method = method,
        impute = impute,
        verbose = verbose,
        .homonyms = "last"
    )

    call <- rlang::call2(domapping, !!!arguments)
    res <- eval(call)

    # Extract results for NC-based and SSP-based methods
    x_NC <- res$x_NC
    x_SSP <- res$x_SSP

    ## create se_obj for NC_based methods
    se_NC <- SummarizedExperiment(
        assays = list(counts = x_NC),
        colData = colData(se_obj)
    )
    ## create se_obj for SSP_based methods
    se_SSP <- SummarizedExperiment(
        assays = list(counts = x_SSP),
        colData = colData(se_obj)
    )
    # Return both objects as a list
    return(list(se_NC = se_NC, se_SSP = se_SSP))
}


#' Original Parker Intrinsic Subtyping (BS_parker)
#'
#' @name BS_parker
#' @description This function predicts breast cancer intrinsic subtypes using
#'   the original Parker et al. (2019) method, along with variations of the original approach.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log-transformed, normalized gene expression matrix with genes (Gene Symbols) as rows
#'   and samples as columns.
#'   - **Column metadata** (`colData`): Optional clinical information.
#' @param calibration Specifies the calibration method to apply. Options include:
#'   - `"None"`: No calibration is applied.
#'   - `"Internal"`: Uses internal calibration strategies (see `internal` argument).
#'   - `"External"`: Uses external medians (see `external` argument).
#' @param internal Specifies the internal calibration strategy when `calibration
#'   = "Internal"`. Options include:
#'   - `"-1"` (default): Median-centered calibration.
#'   - `"meanCtr"`: Mean-centered calibration (aligned with `genefu.scale`).
#'   - `"qCtr"`: Quantile-based calibration (aligned with `genefu.robust`).
#' @param external Specifies the platform name (i.e., column name) for external
#'   medians derived from the training cohort.
#'   - To use user-provided medians, set `"external = "Given.mdns"` and provide values via the
#'   `medians` argument.
#' @param medians A matrix or table of user-provided median values, required if
#'   `external = "Given.mdns"`.
#'   - The first column should contain 50 genes.
#'   - The second column should contain the corresponding median values.
#' @param Subtype Logical (`TRUE` or `FALSE`). If `TRUE`, the function predicts four subtypes,
#'   **excluding** the Normal-like subtype.
#' @param hasClinical Logical (`TRUE` or `FALSE`). If `TRUE`, the function incorporates clinical data from
#'   the phenotype (`pheno`) table. Required columns:
#'   - `"TSIZE"`: Tumor size (`0` for <= 2cm, `1` for > 2cm).
#'   - `"NODE"`: Lymph node status (`0` for negative, `1` or higher for positive nodes; this column must be numeric).
#'
#' @return Returns a list containing intrinsic subtypes assigned using the Parker-based method, or its variations.
#'
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
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     calibration = "Internal",
#'     internal = "-1",
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_parker <- function(
        se_obj,
        calibration = "None",
        internal = NA,
        external = NA,
        medians = NA,
        Subtype = FALSE,
        hasClinical = FALSE) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## maintain parameter
    if (internal == "-1") {
        internal <- "medianCtr"
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    pheno <- colData(se_obj) %>% data.frame()

    # Handle clinical metadata if required
    if (ncol(pheno) == 0) {
        pheno <- NULL
    } else {
        if (!"PatientID" %in% colnames(pheno)) {
            stop("The `colData` of `se_obj` must include a `PatientID` column when `hasClinical = TRUE`.")
        }
        rownames(pheno) <- pheno$PatientID
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
#' the conventional estrogen receptor (ER)-balancing via immunohistochemistry (cIHC).
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log-transformed, normalized gene expression matrix with genes (Gene Symbols)
#'   as rows and samples as columns.
#'   - **Column metadata** (`colData`): A clinical information table, which must include:
#'     - `"PatientID"`: Unique sample or patient identifiers.
#'     - `"ER"`: Estrogen receptor (ER) status, recorded as `"ER+"` or `"ER-"`.
#'
#' @param Subtype Logical (`TRUE` or `FALSE`). If `TRUE`, the function predicts four subtypes,
#'   **excluding** the Normal-like subtype.
#' @param hasClinical Logical (`TRUE` or `FALSE`). If `TRUE`, the function incorporates clinical data from
#'   the phenotype (`pheno`) table. Required columns:
#'   - `"TSIZE"`: Tumor size (`0` for <= 2cm, `1` for > 2cm).
#'   - `"NODE"`: Lymph node status (`0` for negative, `1` or higher for positive nodes; this column must be numeric).
#'
#' @param seed An integer used to set the random seed for reproducibility.
#' @return Returns a data frame containing intrinsic subtypes estimated using the
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
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_cIHC <- function(se_obj,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    # Extract clinical metadata if required
    pheno <- colData(se_obj) %>% data.frame()

    if (!all(c("PatientID", "ER") %in% colnames(pheno))) {
        stop("The 'colData' of 'se_obj' must include 'PatientID' and 'ER' columns when 'hasClinical = TRUE'.")
    }
    rownames(pheno) <- pheno$PatientID

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
#' **iterative** version of conventional estrogen receptor (ER)-balancing via immunohistochemistry (cIHC).
#' It allows customization of the ER+/ER- ratio to refine subtype classification..
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log-transformed, normalized gene expression matrix with genes (Gene Symbol)
#'   in rows and samples in columns.
#'   - **Column metadata** (`colData`): Clinical information table.
#    The column metadata must include:
#'     - `"PatientID"`: Unique sample or patient identifiers.
#'     - `"ER"`: Estrogen receptor (ER) status recorded as `"ER+"` or `"ER-"`.
#'
#' @param iteration Integer. The number of iterations for the ER-balanced procedure
#'   with the specified ratio. Default: 100.
#' @param ratio Numeric. Specifies the ER+/ER- ratio for balancing. Options:
#'   - `1:1`: Equal balancing.
#'   - `54:64`: Default, based on the ER+/ER- ratio in the UNC232 training cohort.
#'
#' @param Subtype Logical (`TRUE` or `FALSE`). If `TRUE`, the function predicts four subtypes,
#'   **excluding** the Normal-like subtype.
#' @param hasClinical Logical (`TRUE` or `FALSE`). If `TRUE`, the function incorporates clinical data from
#'   the phenotype (`pheno`) table. Required columns:
#'   - `"TSIZE"`: Tumor size (`0` for <= 2cm, `1` for > 2cm).
#'   - `"NODE"`: Lymph node status (`0` for negative, `1` or higher for positive nodes; this column must be numeric).
#'
#' @param seed An integer used to set the random seed for reproducibility.
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
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     iteration = 10, ## For final analysis, set iteration = 100
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_cIHC.itr <- function(se_obj,
    iteration = 100,
    ratio = 54 / 64,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    # Extract clinical metadata if required
    pheno <- colData(se_obj) %>% data.frame()

    if (!all(c("PatientID", "ER") %in% colnames(pheno))) {
        stop("The 'colData' of 'se_obj' must include 'PatientID' and 'ER' columns when 'hasClinical = TRUE'.")
    }
    rownames(pheno) <- pheno$PatientID

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
#' the PCA-PAM50 method. This approach integrates **Principal Component Analysis (PCA)**
#'  to perform estrogen receptor (ER) balancing based on ESR1 gene expression.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log-transformed, normalized gene expression matrix with genes (Gene Symbols)
#'   as rows and samples as columns.
#'   - **Column metadata** (`colData`): Clinical information table.
#    The column metadata must include:
#'     - `"PatientID"`: Unique sample or patient identifiers.
#'     - `"ER"`: Estrogen receptor (ER) status, recorded as `"ER+"` or `"ER-"`.
#'
#' @param Subtype Logical (`TRUE` or `FALSE`). If `TRUE`, the function predicts four subtypes,
#'   **excluding** the Normal-like subtype.
#' @param hasClinical Logical (`TRUE` or `FALSE`). If `TRUE`, the function incorporates clinical data from
#'   the phenotype (`pheno`) table. Required columns:
#'   - `"TSIZE"`: Tumor size (`0` for <= 2cm, `1` for > 2cm).
#'   - `"NODE"`: Lymph node status (`0` for negative, `1` or higher for positive nodes; this column must be numeric).
#'
#' @param seed An integer used to set the random seed for reproducibility.
#' @return Returns a vector of intrinsic subtypes assigned to the samples, as estimated
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
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_PCAPAM50 <- function(se_obj,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    # Extract clinical metadata if required
    pheno <- colData(se_obj) %>% data.frame()

    if (!all(c("PatientID", "ER") %in% colnames(pheno))) {
        stop("The 'colData' of 'se_obj' must include 'PatientID' and 'ER' columns when 'hasClinical = TRUE'.")
    }
    rownames(pheno) <- pheno$PatientID

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


#' Subgroup-specific gene-centering Intrinsic Subtyping (BS_ssBC)
#'
#' @name BS_ssBC
#' @description This function predicts breast cancer intrinsic subtypes using
#'   the **subgroup-specific (ssBC)** method. The ssBC method applies a
#'   subgroup-specific gene-centering approach to cohorts with a skewed
#'   distribution of clinicopathological characteristics compared to the
#'   original training cohort (e.g., an ER+ selected cohort).
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log-transformed, normalized gene expression matrix with genes (Gene Symbols)
#'   as rows and samples as columns.
#'   - **Column metadata** (`colData`): A clinical information table. If `hasClinical = TRUE`,
#'     this table must include:
#'     - `"PatientID"`: Unique identifiers for patients or samples.
#'     - Additional columns depending on the `s` parameter:
#'       - `"ER"`: Estrogen receptor status (`"ER+"` or `"ER-"`) if `s = "ER"`.
#'       - `"HER2"`: HER2 status (`"HER2+"` or `"HER2-"`) if `s = "ER.v2"`.
#'       - `"TN"`: Triple-negative status (`"TN"` or `"nonTN"`) if `s = "TN"` or `s = "TN.v2"`, indicating a triple-negative cohort.
#' @param s Character.  Specifies which subgroup-specific quantiles to use:
#'   - `"ER"` and `"TN"`: Original subgroup-specific quantiles published in *Breast Cancer Research* (2015).
#'   - `"ER.v2"` and `"TN.v2"`: Updated subgroup-specific quantiles published in *Journal of Clinical Oncology* (2020).
#' @param Subtype Logical (`TRUE` or `FALSE`). If `TRUE`, the function predicts four subtypes,
#'   **excluding** the Normal-like subtype.
#' @param hasClinical Logical (`TRUE` or `FALSE`). If `TRUE`, the function incorporates clinical data from
#'   the phenotype (`pheno`) table. Required columns:
#'   - `"TSIZE"`: Tumor size (`0` for <= 2cm, `1` for > 2cm).
#'   - `"NODE"`: Lymph node status (`0` for negative, `1` or higher for positive nodes; this column must be numeric).
#'
#' @return Returns a vector of intrinsic subtypes assigned to the samples, as estimated
#'   by the ssBC method.
#'
#' @references
#' - Zhao X, Rodland EA, Tibshirani R, Plevritis S. *Molecular
#' subtyping for clinically defined breast cancer subgroups.* Breast Cancer
#' Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4
#' - Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L,
#' et al. *Survival, pathologic response, and genomics in CALGB 40601
#' (Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or
#' without lapatinib in HER2-positive breast cancer.* Journal of Clinical
#' Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
#'
#' @examples
#' ## ssBC.v2
#' data("OSLO2EMIT0obj")
#' res <- BS_ssBC(
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     s = "ER.v2",
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_ssBC <- function(se_obj,
    s,
    Subtype = FALSE,
    hasClinical = FALSE) {
    # Check that input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    # Extract gene expression matrix
    gene_expr <- assay(se_obj)

    # Extract clinical metadata if hasClinical is TRUE
    pheno <- colData(se_obj) %>% data.frame()
    if (!"PatientID" %in% colnames(pheno)) {
        stop("The 'colData' of 'se_obj' must include a 'PatientID' column when 'hasClinical = TRUE'.")
    }
    rownames(pheno) <- pheno$PatientID

    # Additional checks based on `s`
    required_columns <- switch(s,
        "ER" = c("ER"),
        "ER.v2" = c("ER", "HER2"),
        "TN" = c("TN"),
        "TN.v2" = c("TN"),
        stop("Invalid value for 's'. Must be one of 'ER', 'ER.v2', 'TN', or 'TN.v2'.")
    )

    missing_columns <- setdiff(required_columns, colnames(pheno))
    if (length(missing_columns) > 0) {
        stop(paste("The following required columns are missing from 'colData':", paste(missing_columns, collapse = ", ")))
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
#' the **AIMS (Absolute assignment of Intrinsic Molecular Subtype)** method.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A gene expression matrix with genes (EntrezID) as rows and samples as columns.
#'     Important: The gene expression values should not be gene-centered.
#'      All expression values must be **positive**.
#'
#' @return Returns a vector of intrinsic subtypes assigned to the samples, as estimated
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
#'
#'
#' # Perform subtyping
#' res <- BS_AIMS(
#'     se_obj = OSLO2EMIT0obj$data_input$se_SSP
#' )
#'
#' @export

BS_AIMS <- function(se_obj) {
    ## loading datasets

    data_env <- new.env(parent = emptyenv())
    data("AIMSmodel", envir = data_env, package = "BreastSubtypeR")
    AIMSmodel <- data_env[["AIMSmodel"]]
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    # Check that input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    # Extract gene expression matrix
    gene_expr <- assay(se_obj)

    # Extract AIMS-specific genes
    genes <- as.character(
        BreastSubtypeRobj$genes.signature$EntrezGene.ID[
            which(BreastSubtypeRobj$genes.signature$AIMS == "Yes")
        ]
    )
    gene_expr <- gene_expr[rownames(gene_expr) %in% genes, ]

    arguments <- rlang::dots_list(
        eset = as.matrix(gene_expr),
        EntrezID = rownames(gene_expr)
    )

    call <- rlang::call2(applyAIMS_AIMS, !!!arguments)

    res_AIMS <- eval(call)
}


#' Intrinsic Subtyping using SSPBC (BS_sspbc)
#'
#' @name BS_sspbc
#' @description This function predicts breast cancer intrinsic subtypes using
#' SSPBC (Single Sample Predictor for Breast Cancer). SSPBC is
#' based on a refined version of the original AIMS methodology, utilizing a large,
#' uniformly accrued population-based cohort (SCAN-B) for training. This method supports
#' RNA sequencing data and provides flexibility in selecting the prediction model.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A gene expression matrix with genes (EntrezID) as rows and samples as columns.
#'     Important: The gene expression values should not be gene-centered.
#'      All expression values must be **positive**.
#'
#' @param ssp.name Specifies the model to use. Options are:
#' - "ssp.pam50": For PAM50-based predictions.
#' - "ssp.subtype": For predicting Prosigna-like subtypes (four subtypes, **excluding** the Normal-like subtype).
#'
#' @return Returns a vector of intrinsic subtypes assigned to the samples, as estimated
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
#'     se_obj = OSLO2EMIT0obj$data_input$se_SSP,
#'     ssp.name = "ssp.pam50"
#' )
#'
#' @export

BS_sspbc <- function(se_obj, ssp.name = "ssp.pam50") {
    data_env <- new.env(parent = emptyenv())
    data("sspbc.models", envir = data_env, package = "BreastSubtypeR")
    sspbc.models <- data_env[["sspbc.models"]]
    data("sspbc.models.fullname", envir = data_env, package = "BreastSubtypeR")
    sspbc.models.fullname <- data_env[["sspbc.models.fullname"]]

    # Check that input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    # Extract gene expression matrix
    gene_expr <- assay(se_obj)

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
#' @description This function predicts breast cancer intrinsic subtypes using multiple methods.
#' Users can either specify the subtyping approaches directly or enable automatic selection ("AUTO")
#' based on the ER/HER2 distribution of the test cohort.
#'
#' @param data_input The output from the `Mapping()` function, containing
#'   processed gene expression data prepared for subtyping analysis.
#' @param methods A character vector specifying the subtyping methods to be
#'   used. Available options:
#'   - "parker.original": Original PAM50 by Parker et al., 2009 (Parker et al., 2009)
#'   - "genefu.scale": PAM50 implementation as in the genefu R package (scaled version) (Gendoo et al., 2016)
#'   - "genefu.robust": PAM50 implementation as in the genefu R package (robust version) (Gendoo et al., 2016)
#'   - "cIHC": Conventional estrogen receptor (ER)-balancing using immunohistochemistry (cIHC) (Ciriello et al., 2015)
#'   - "cIHC.itr": Iterative version of cIHC (Curtis et al., 2012)
#'   - "PCAPAM50": PCA-based iterative PAM50 (ER-balancing using ESR1 gene expression) (Raj-Kumar et al., 2019)
#'   - "ssBC": Subgroup-specific gene-centering PAM50 (Zhao et al., 2015)
#'   - "ssBC.v2": Updated subgroup-specific gene-centering PAM50 with refined quantiles (Fernandez-Martinez et al., 2020)
#'   - "AIMS": Absolute Intrinsic Molecular Subtyping (AIMS) method (Paquet & Hallett, 2015)
#'   - "sspbc": Single-Sample Predictors for Breast Cancer (AIMS adaptation) (Staaf et al., 2022)
#'   - "AUTO": Automatically selects subtyping methods based on the ER/HER2 distribution of the test cohort.
#'
#'   Notes:
#'   - If "AUTO" is selected, it must be the sole value in the vector.
#'   - If "AUTO" is not selected, at least **two** methods must be specified; otherwise, an error will occur.
#' @param Subtype Logical (`TRUE` or `FALSE`). If `TRUE`, the function predicts four subtypes,
#'   **excluding** the Normal-like subtype.
#' @param hasClinical Logical (`TRUE` or `FALSE`). If `TRUE`, the function incorporates clinical data from
#'   the phenotype (`pheno`) table. Required columns:
#'   - `"TSIZE"`: Tumor size (`0` for <= 2cm, `1` for > 2cm).
#'   - `"NODE"`: Lymph node status (`0` for negative, `1` or higher for positive nodes; this column must be numeric).
#'
#' @return Returns a list of intrinsic subtypes estimated by the selected methods.
#'
#'
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
#' - Ciriello G, Gatza ML, Beck AH, Wilkerson MD, Rhie SK, Pastore A,
#' et al. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer*.
#' Cell. 2015;163(2). https://doi.org/10.1016/j.cell.2015.09.033
#'
#' - Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ,
#' et al. *The genomic and transcriptomic architecture of 2,000 breast tumours
#' reveals novel subgroups*. Nature. 2012;486(7403).
#' https://doi.org/10.1038/nature10983
#'
#' - Zhao X, Rodland EA, Tibshirani R, Plevritis S. *Molecular
#' subtyping for clinically defined breast cancer subgroups.* Breast Cancer
#' Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4
#'
#' - Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L,
#' et al. *Survival, pathologic response, and genomics in CALGB 40601
#' (Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or
#' without lapatinib in HER2-positive breast cancer.* Journal of Clinical
#' Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
#'
#' - Zhao X, Rodland EA, Tibshirani R, Plevritis S. *Molecular
#' subtyping for clinically defined breast cancer subgroups.* Breast Cancer
#' Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4
#'
#' - Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L,
#' et al. *Survival, pathologic response, and genomics in CALGB 40601
#' (Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or
#' without lapatinib in HER2-positive breast cancer.* Journal of Clinical
#' Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
#'
#' - Paquet ER, Hallett MT. *Absolute assignment of breast cancer
#' intrinsic molecular subtype.* J Natl Cancer Inst. 2015;107(1).
#' https://doi.org/10.1093/jnci/dju357
#'
#' - Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I,
#' et al. *RNA sequencing-based single sample predictors of molecular subtype
#' and risk of recurrence for clinical assessment of early-stage breast cancer*.
#' NPJ Breast Cancer. 2022;8(1). https://doi.org/10.1038/s41523-022-00465-3
#'
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
#'     data_input = OSLO2EMIT0obj$data_input,
#'     methods = methods,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_Multi <- function(
        data_input,
        methods = "AUTO",
        Subtype = FALSE,
        hasClinical = FALSE) {
    valid_methods <- c(
        "parker.original", "genefu.scale", "genefu.robust",
        "ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50",
        "AIMS", "sspbc", "AUTO"
    )

    invalid_methods <- methods[!methods %in% valid_methods]

    if (length(invalid_methods) > 0) {
        stop("Invalid method(s) specified: ", paste(invalid_methods, collapse = ", "))
    }

    if (length(methods) == 1 && methods[1] == "AUTO") {
        message("Running AUTO mode for subtyping.")
    } else if (length(methods) < 2) {
        stop("Select at least two methods or set method to 'AUTO'.")
    }

    ## extract pheno table
    pheno <- colData(data_input$se_NC) %>% data.frame()
    rownames(pheno) <- pheno$PatientID

    # Check ER and HER2 columns in pheno
    if (!("ER" %in% colnames(pheno)) && any(methods %in% c("ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50"))) {
        stop("The 'ER' column is required for selected methods.")
    }
    if (!("HER2" %in% colnames(pheno)) && "ssBC.v2" %in% methods) {
        stop("The 'HER2' column is required for the 'ssBC.v2' method.")
    }

    ## AUTO mode
    # methods = "AUTO"
    cohort.select <- "ERpos"
    samples_ER.icd <- NULL
    samples_ERHER2.icd <- NULL
    if (length(methods) == 1 && methods[1] == "AUTO") {
        AUTO.output <- get_methods(pheno)
        samples_ER.icd <- AUTO.output$samples_ER.icd
        samples_ERHER2.icd <- AUTO.output$samples_ERHER2.icd
        methods <- AUTO.output$methods
        cohort.select <- AUTO.output$cohort.select
    }

    ## run each method
    results <- lapply(methods, function(method) {
        ## try NC-based
        if (method == "parker.original") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$se_NC,
                    calibration = "Internal",
                    internal = "-1",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }
        if (method == "genefu.scale") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$se_NC,
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
                    data_input$se_NC,
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
                    data_input$se_NC,
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
                    data_input$se_NC,
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }

        if (method == "PCAPAM50") {
            message(method, " is running!")
            result <- tryCatch(
                {
                    # Attempt to run PCAPAM50
                    BS_PCAPAM50(
                        data_input$se_NC,
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                },
                error = function(e) {
                    # Error handling
                    warning("PCAPAM50 failed in this iteration. Error: ")
                    return(NULL) # Return NULL or a dummy tibble with NAs
                }
            )

            # Return result (or NULL if failed)
            return(result)
        }

        if (method == "ssBC") {
            message(method, " is running!")

            if ("TN" %in% colnames(pheno) && cohort.select == "TNBC") {
                res_ssBC <- BS_ssBC(
                    data_input$se_NC,
                    s = "TN",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            } else {
                if (!is.null(samples_ER.icd) & length(samples_ERHER2.icd) < nrow(pheno)) {
                    res_ssBC <- BS_ssBC(
                        data_input$se_NC[, samples_ER.icd],
                        s = "ER",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                } else {
                    res_ssBC <- BS_ssBC(
                        data_input$se_NC,
                        s = "ER",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                }
            }

            ## Keep patients' results for AUTO mode
            if (!is.null(samples_ER.icd) &&
                length(samples_ER.icd) < nrow(pheno)) {
                unprocessed_patients <- base::setdiff(pheno$PatientID, samples_ER.icd)

                # Create NA-filled dataframe for unprocessed patients with matching structure
                na_df <- data.frame(
                    PatientID = unprocessed_patients,
                    matrix(
                        nrow = length(unprocessed_patients),
                        ncol = ncol(res_ssBC$BS.all) - 1, # Subtract 1 for PatientID column
                        dimnames = list(NULL, colnames(res_ssBC$BS.all)[-1])
                    ),
                    row.names = unprocessed_patients
                )

                res_ssBC$BS.all <- rbind(res_ssBC$BS.all, na_df)
                res_ssBC$BS.all <- res_ssBC$BS.all[pheno$PatientID, ]
            }

            return(res_ssBC)
        }

        if (method == "ssBC.v2") {
            message(method, " is running!")

            if ("TN" %in% colnames(pheno) && cohort.select == "TNBC") {
                res_ssBC.v2 <- BS_ssBC(
                    data_input$se_NC,
                    s = "TN.v2",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            } else {
                if (!is.null(samples_ERHER2.icd) & length(samples_ERHER2.icd) < nrow(pheno)) {
                    res_ssBC.v2 <- BS_ssBC(
                        data_input$se_NC[, samples_ERHER2.icd],
                        s = "ER.v2",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                } else {
                    res_ssBC.v2 <- BS_ssBC(
                        data_input$se_NC,
                        s = "ER.v2",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                }
            }

            ## Keep patients' results for AUTO mode
            if (!is.null(samples_ERHER2.icd) &&
                length(samples_ERHER2.icd) < nrow(pheno)) {
                unprocessed_patients <- base::setdiff(pheno$PatientID, samples_ERHER2.icd)

                # Create NA-filled dataframe for unprocessed patients with matching structure
                na_df <- data.frame(
                    PatientID = unprocessed_patients,
                    matrix(
                        nrow = length(unprocessed_patients),
                        ncol = ncol(res_ssBC.v2$BS.all) - 1, # Subtract 1 for PatientID column
                        dimnames = list(NULL, colnames(res_ssBC.v2$BS.all)[-1])
                    ),
                    row.names = unprocessed_patients
                )

                res_ssBC.v2$BS.all <- rbind(res_ssBC.v2$BS.all, na_df)
                res_ssBC.v2$BS.all <- res_ssBC.v2$BS.all[pheno$PatientID, ]
            }


            return(res_ssBC.v2)
        }

        if (method == "AIMS") {
            message(method, " is running!")
            data_env <- new.env(parent = emptyenv())
            data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
            BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


            res_AIMS <- BS_AIMS(data_input$se_SSP)

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
                se_obj = data_input$se_SSP,
                ssp.name = "ssp.pam50"
            )

            BS.all <- data.frame(
                PatientID = rownames(res_sspbc),
                BS = res_sspbc[, 1],
                row.names = rownames(res_sspbc)
            )

            if (Subtype) {
                res_sspbc.Subtype <- BS_sspbc(
                    se_obj = data_input$se_SSP,
                    ssp.name = "ssp.subtype"
                )
                BS.all$BS.Subtype <- res_sspbc.Subtype[, 1]
            }

            return(list(BS.all = BS.all))
        }

        stop("Unknown method: ", method)
    })

    names(results) <- methods
    samples <- colnames(assay(data_input$se_NC))

    res_subtypes <- data.table(row_id = samples)
    if (Subtype) {
        res_subtypes.Subtype <- data.table(row_id = samples)
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
    rownames(res_subtypes) <- samples
    res_subtypes[, 1] <- NULL

    entropy <- apply(res_subtypes, 1, get_entropy)
    res_subtypes$entropy <- entropy

    if (Subtype) {
        res_subtypes.Subtype <- as.data.frame(res_subtypes.Subtype)
        rownames(res_subtypes.Subtype) <- samples
        res_subtypes.Subtype[, 1] <- NULL

        ## removing Normal-like in AIMS (to be comparable with other methods)
        if ("AIMS" %in% methods) {
            res_subtypes.Subtype$AIMS[which(res_subtypes.Subtype$AIMS == "Normal")] <- NA
        }
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
