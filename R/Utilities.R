#'
#' Functions for BreastSubtypeR package
#' @import ggplot2
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import circlize
#' @import ggrepel
#' @import magrittr
#' @import impute
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#'
#' @noRd
NULL


#' Map Gene IDs and Handle missing data
#'
#' @name Mapping
#' @description A function to map gene IDs and preprocess gene expression data
#' for subsequent analyses.
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log2-transformed gene expression matrix where rows represent probes (e.g., ProbeID, TranscriptID, or Gene Symbol)
#'     and columns represent samples. This should be stored in the `assay()` slot.
#'   - **Row metadata**: A data frame with annotations for probes, including at least the following columns:
#'     - `"probe"`: Unique identifiers for the probes (e.g., ProbeID, TranscriptID or gene symbol).
#'     - `"ENTREZID"`: Entrez gene IDs corresponding to the probes.
#' @param method Method for deduplicating probes in microarray or RNA-seq data.
#'   Choose one of the following options:
#'   - `"iqr"` for Affymetrix arrays,
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
#' @noRd

domapping <- function(
        se_obj,
        method = "max",
        impute = TRUE,
        verbose = TRUE) {
    data("BreastSubtypeRobj")

    ## Extract data from SummarizedExperiment
    x <- assay(se_obj)
    y <- rowData(se_obj)

    if (!is.data.frame(y)) {
        y <- as.data.frame(y)
    }

    samplenames <- colnames(x)

    ## check feature data
    if (length(y) == 0) {
        stop("Please provide feature annotation to do probeID mapping ")
    }

    ## loading genes.signature
    genes.s <- BreastSubtypeRobj$genes.signature

    ## filter by ENTREZID
    y <- y[y$ENTREZID %in% genes.s$EntrezGene.ID, ]
    x <- x[y$probe, ]

    ## first step
    ## for empty cells. imput or not ?
    if (sum(apply(x, 2, is.na)) > 0 & impute) {
        if (verbose) {
            probeid_NA <- rownames(x)[rowSums(is.na(x))]
            sample_NA <- colnames(x)[colSums(is.na(x))]
        }

        x <- impute::impute.knn(x)
        x <- x$data
    }


    #
    ## second step
    ## if mapping (microarray or transcript )

    probeid <- rownames(x)
    entrezid <- y$ENTREZID
    names(entrezid) <- y$probe

    ## process probeid in input data
    entrezid <- entrezid[probeid]
    ## remove NA
    entrezid <- entrezid[!(is.na(entrezid))]
    x <- x[names(entrezid), ]
    entrezid <- factor(entrezid, levels = unique(entrezid))
    ## names are unique probeid and content are redundant entrezid


    ## This is for probeID or transcriptID
    ## split expression matrix
    split_mat <- split(as.data.frame(x), entrezid, drop = FALSE)

    # function to calculate the desired statistic
    calculate_stat <- function(mat, method) {
        switch(method,
            "mean" = apply(mat, 2, mean, na.rm = TRUE),
            "median" = apply(mat, 2, median, na.rm = TRUE),
            "max" = mat[which.max(rowSums(mat)), ],
            "stdev" = {
                temp <- mat
                stdevs <- apply(mat, 1, sd, na.rm = TRUE)
                vals <- temp[match(max(stdevs), stdevs), ]
                return(vals)
            },
            "iqr" = {
                temp <- mat
                iqrs <- apply(mat, 1, function(x) {
                    quantile(x, 0.75, na.rm = TRUE) -
                        quantile(x, 0.25, na.rm = TRUE)
                })
                vals <- temp[match(max(iqrs), iqrs), ]
                return(vals)
            }
        )
    }

    ## keep processed x
    x <- mapply(
        calculate_stat,
        split_mat,
        MoreArgs = list(method = method),
        SIMPLIFY = TRUE,
        USE.NAMES = TRUE
    )
    x <- apply(x, 1, unlist)

    ## print necessary information
    ## Parker
    missing_ID_parker <- setdiff(
        BreastSubtypeRobj$genes.sig50$EntrezGene.ID,
        rownames(x)
    )
    if (length(missing_ID_parker) == 0 & verbose) {
        message("Genes used in NC-based methods are covered.")
    } else if (verbose) {
        message("These genes are missing for NC-based methods:")
        message(paste(missing_ID_parker, collapse = "\n"))
    }

    ## AIMS
    missing_ID_AIMS <- setdiff(
        genes.s[genes.s$SSP_based == "Yes", ]$EntrezGene.ID,
        rownames(x)
    )
    if (length(missing_ID_AIMS) == 0 & verbose) {
        message("Genes used in SSP-based methods are covered.")
    } else if (verbose) {
        message("These genes are missing for SSP-based methods:")
        message(paste(missing_ID_AIMS, collapse = "\n"))
    }

    ## get matrix for NC (symbol as rows, sample as col)
    genes_nc <- BreastSubtypeRobj$genes.sig50$EntrezGene.ID
    x_NC <- x[na.omit(match(genes_nc, rownames(x))), ]
    rownames(x_NC) <- BreastSubtypeRobj$genes.sig50$Symbol[match(
        rownames(x_NC),
        genes_nc
    )]
    x_NC <- data.frame(x_NC)
    colnames(x_NC) <- samplenames

    ## get matrix for AIMS (entrezID as colnames)
    genes_ssp <- as.character(genes.s$EntrezGene.ID[
        genes.s$SSP_based == "Yes"
    ])
    x_SSP <- x[na.omit(match(genes_ssp, rownames(x))), ]
    colnames(x_SSP) <- samplenames

    ## exponential transformation
    x_SSP <- 2^x_SSP

    result <- list(
        x_NC = x_NC,
        x_SSP = x_SSP
    )

    return(result)
}

#' Function for consensus subtype
#' @noRd
get_methods <- function(pheno) {
    #### AUTO mode
    samples_ER.icd <- NULL
    samples_ERHER2.icd <- NULL
    cohort.select <- "ERpos"
    if (ncol(pheno) == 0) {
        message("The pheno table has not been detected.")
        message("Running methods: AIMS, & sspbc")
        methods <- c("AIMS", "sspbc")
    } else if (!all(c("ER", "HER2") %in% colnames(pheno))) {
        stop("The 'AUTO' mode requires both 'ER' and 'HER2' columns in the 'pheno' dataframe.")
    } else {
        # Calculate sample sizes
        sample_counts <- with(pheno, table(ER, HER2))
        n_ERpos <- sum(pheno$ER == "ER+")
        n_ERneg <- sum(pheno$ER == "ER-")

        n_ERnegHER2pos <- sum(pheno$ER == "ER-" & pheno$HER2 == "HER2+")
        n_ERnegHER2neg <- sum(pheno$ER == "ER-" & pheno$HER2 == "HER2-")
        n_ERposHER2pos <- sum(pheno$ER == "ER+" & pheno$HER2 == "HER2+")
        n_ERposHER2neg <- sum(pheno$ER == "ER+" & pheno$HER2 == "HER2-")

        # Set thresholds
        n_ER_threshold <- 10
        n_ERHER2_threshold <- 5
        per_ratio <- 0.2
        upper_ratio <- 54 / 64 + (54 / 64) * per_ratio
        lower_ratio <- 54 / 64 - (54 / 64) * per_ratio

        ## main panel
        if (n_ERposHER2neg == 0 && n_ERnegHER2neg == 0) {
            message("A HER2+ cohort has been detected.")
            cohort.select <- "HER2pos"

            if (n_ERposHER2pos < n_ERHER2_threshold && n_ERnegHER2pos < n_ERHER2_threshold) {
                message("A small HER2+ cohort has been detected.")
                message("Running methods: AIMS, & sspbc")
                methods <- c("AIMS", "sspbc")
            } else if (n_ERposHER2pos >= n_ERHER2_threshold && n_ERnegHER2pos < n_ERHER2_threshold) {
                message("A ER+/HER2+ cohort has been detected.")
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2pos < n_ERHER2_threshold && n_ERnegHER2pos >= n_ERHER2_threshold) {
                message("A ER-/HER2+ cohort has been detected.")
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            } else {
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERnegHER2pos == 0 && n_ERposHER2pos == 0 && n_ERposHER2neg == 0) {
            message("A TNBC cohort has been detected.")
            cohort.select <- "TNBC"

            if (!("TN" %in% colnames(pheno))) {
                stop("Provide \"TN\" in pheno for: ssBC(TN) & ssBC.v2 (TN)")
            }

            message("Running methods: ssBC (TN), ssBC.v2 (TN), AIMS & sspbc")
            methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
        } else if (n_ERpos < n_ER_threshold && n_ERneg < n_ER_threshold) {
            message("A small number of ER-/ER+ samples has been detected.")
            message("Running methods: AIMS & sspbc")
            methods <- c("AIMS", "sspbc")
        } else if (n_ERpos >= n_ER_threshold && n_ERneg < n_ER_threshold) {
            if (n_ERposHER2pos >= n_ERHER2_threshold && n_ERposHER2neg >= n_ERHER2_threshold) {
                message("Running methods for ER+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2pos >= n_ERHER2_threshold && n_ERposHER2neg < n_ERHER2_threshold) {
                message("Running methods for ER+/HER2+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2pos < n_ERHER2_threshold && n_ERposHER2neg >= n_ERHER2_threshold) {
                message("Running methods for ER+/HER2- samples:
                        ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERpos < n_ER_threshold && n_ERneg >= n_ER_threshold) {
            if (n_ERnegHER2pos >= n_ERHER2_threshold && n_ERnegHER2neg >= n_ERHER2_threshold) {
                message("Running methods for ER- samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERnegHER2pos >= n_ERHER2_threshold && n_ERnegHER2neg < n_ERHER2_threshold) {
                message("Running methods for ER-/HER2+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERnegHER2pos < n_ERHER2_threshold && n_ERnegHER2neg >= n_ERHER2_threshold) {
                message("Running methods for ER-/HER2- samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERpos >= n_ER_threshold && n_ERneg >= n_ER_threshold) {
            ## for other NC-based methods
            ratio_ER <- n_ERpos / n_ERneg

            if (ratio_ER > lower_ratio && ratio_ER < upper_ratio) {
                message(
                    "Running methods:
                    parker.original, genefu.scale, genefu.robust, ssBC, ssBC.v2, cIHC, cIHC.itr, PCAPAM50, AIMS & sspbc"
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
                    "The ER+/ER- ratio in the current dataset differs from that observed in the UNC232 training cohort."
                )
                message("Running methods:
                        ssBC, ssBC.v2, cIHC, cIHC.itr, PCAPAM50, AIMS & sspbc")
                methods <- c(
                    "genefu.robust",
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

        if (cohort.select != "TNBC") {
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

            er_idx <- ERHER2_counts[seq(1, 2)] > n_ER_threshold
            samples_ER <- names(ERHER2_counts)[seq(1, 2)][er_idx]
            erher2_idx <- ERHER2_counts[seq(3, 6)] > n_ERHER2_threshold
            samples_ERHER2 <- names(ERHER2_counts)[seq(3, 6)][erher2_idx]

            if (cohort.select != "HER2pos" && cohort.select != "TNBC") {
                if (length(samples_ER) > 0) {
                    message("ssBC for samples: ", paste(samples_ER, collapse = ", "))
                    samples_ER.icd <- unlist(lapply(samples_ER, function(subtype) {
                        subtype <- str_replace_all(subtype, "pos", "+")
                        subtype <- str_replace_all(subtype, "neg", "-")
                        ER_status <- subtype
                        rownames(pheno)[pheno$ER == ER_status]
                    }))
                }
            }

            if (cohort.select != "TNBC") {
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
        }
    }

    return(list(samples_ER.icd = samples_ER.icd, samples_ERHER2.icd = samples_ERHER2.icd, methods = methods))
}

#' Function for consensus subtype
#' @noRd
get_consensus_subtype <- function(patient_row) {
    patient_row <- unlist(patient_row, use.names = FALSE)
    counts <- table(patient_row)
    max_subtype <- names(counts)[which.max(counts)]
    return(max_subtype)
}

#' Function for entropy calculation
#' @noRd
get_entropy <- function(patient_row) {
    freq <- table(patient_row)
    prob <- freq / sum(freq)
    entropy <- -sum(prob * log2(prob))
    return(entropy)
}


#' Function to get the average correlation and ROR
#' @noRd
get_average_subtype <- function(res_ihc_iterative, consensus_subtypes) {
    ## correlation and ROR to be averaged
    ## if hasclini ?? need to be added later
    sum_colnames <- c("Basal", "Her2", "LumA", "LumB", "Normal")

    all_patients <- names(res_ihc_iterative[[1]]$predictions)

    sum_cols_list <- mapply(
        function(res_ihc) {
            # res_ihc = as.data.frame( res_ihc_iterative[[1]])

            res_ihc$distances <- as.data.frame(res_ihc$distances)

            ## if FALSE, make the cell as NULL
            keep <- res_ihc$predictions == consensus_subtypes
            res_ihc$distances[!keep, ] <- as.list(rep(NA, 5))

            res <- dplyr::mutate_at(
                res_ihc$distances,
                vars(everything()),
                ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA)
            )

            return(res)
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )


    ## count_na
    count_weight_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
        x[!is.na(x)] <- 1
        x[is.na(x)] <- 0
        return(x)
    }))

    ## sum all for each cell
    sum_cols_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
        ## change NA cell to 0 cell
        x[is.na(x)] <- 0

        return(x)
    }))

    ## get the mean for each cell
    ## only when subtype is supported by consensus_subtypes
    ## for each iteration and each patient
    mean_cols_save <- sum_cols_save / count_weight_save

    ## get the mean testdata
    sum_cols_list.testdata <- mapply(
        function(res_ihc) {
            res_ihc$testData
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )

    ## sum all for each cell
    mean_cols_save.testdata <- Reduce(
        `+`,
        sum_cols_list.testdata
    ) / length(sum_cols_list.testdata)


    ## distances.Subtype for ROR
    sum_cols_list.Subtype <- mapply(
        function(res_ihc) {
            # res_ihc = res_ihc_iterative[[1]]

            res_ihc$distances.Subtype <- data.frame(res_ihc$dist.Subtype)

            ## if FALSE, make the cell as NULL
            keep <- res_ihc$predictions == consensus_subtypes
            res_ihc$distances.Subtype[!keep, ] <- as.list(rep(NA, 4))

            res <- dplyr::mutate_at(
                res_ihc$distances.Subtype,
                vars(everything()),
                ~ ifelse(!is.na(.), as.numeric(as.character(.)), NA)
            )

            return(res)
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )

    ## count_na
    count_weight_save.Subtype <- Reduce(
        `+`,
        lapply(sum_cols_list.Subtype, function(x) {
            x[!is.na(x)] <- 1
            x[is.na(x)] <- 0
            return(x)
        })
    )


    ## sum all for each cell
    sum_cols_save.Subtype <- Reduce(
        `+`,
        lapply(sum_cols_list.Subtype, function(x) {
            ## change NA cell to 0 cell
            x[is.na(x)] <- 0
            return(x)
        })
    )

    ## get the mean for each cell
    ## only when subtype is supported by consensus_subtypes
    ## for each iteration and each patient
    sum_cols_save.Subtype <- sum_cols_save.Subtype / count_weight_save.Subtype

    res <- list(
        mean_distance = mean_cols_save,
        mean_distance.Subtype = sum_cols_save.Subtype,
        testdata = mean_cols_save.testdata
    )


    return(res)
}

#### Visualization  ####

#' Boxplot of Correlation per Subtype
#'
#' @name Vis_boxplot
#' @description This function generates a boxplot to visualize the correlation
#' distribution between different subtypes of breast cancer, based on the
#' provided correlation table and subtype information.
#'
#' @param out A data frame containing the columns `"PatientID"` and `"Subtype"`.
#'   The `"PatientID"` column should have unique identifiers for each patient,
#'   and the `"Subtype"` column should specify the assigned subtype for each
#'   patient.
#' @param correlations A data frame or matrix containing the correlation values
#'   computed from NC-based methods.
#'
#' @return A `ggplot` object representing the boxplot visualization of the
#'   correlation distributions across the different subtypes.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Subtype information and correlation matrix
#' out <- data.frame(
#'     PatientID = res$results$parker.original$BS.all$PatientID,
#'     Subtype = res$results$parker.original$BS.all$BS
#' )
#' correlations <- res$results$parker.original$outList$distances
#'
#' # Generate the boxplot
#' p <- Vis_boxplot(out, correlations)
#' plot(p)
#'
#' @export

Vis_boxplot <- function(out, correlations) {
    df <- data.frame(
        predictions = out$Subtype,
        cor = apply(correlations, 1, max)
    )

    plot <- ggplot(df, aes(x = predictions, y = cor)) +
        geom_boxplot() +
        labs(x = "", y = "Correlation") +
        theme_classic() +
        theme(axis.text = element_text(size = 12))

    return(plot)
}


#' Heatmap Visualization of Gene Expression by Subtype
#'
#' @name Vis_heatmap
#' @description This function generates a heatmap to visualize gene expression
#' patterns across breast cancer subtypes, based on the provided gene expression
#' matrix and subtype information.
#'
#' @param x A gene expression matrix, where genes are rows and samples are
#'   columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#'
#' @return A `ggplot` or `heatmap` object (depending on implementation)
#'   representing the heatmap of gene expression across different subtypes.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Gene expression matrix and subtype information
#' x <- assay(OSLO2EMIT0obj$data_input$se_NC)
#' out <- data.frame(
#'     PatientID = res$results$parker.original$BS.all$PatientID,
#'     Subtype = res$results$parker.original$BS.all$BS
#' )
#'
#' # Generate the heatmap
#' p <- Vis_heatmap(x, out)
#' plot(p)
#'
#' @export

Vis_heatmap <- function(x, out) {
    scaled_mat <- t(scale(t(x[, out$PatientID])))

    ## color
    scaled_dat <- c(min(scaled_mat), 0, max(scaled_mat))
    col_fun <- circlize::colorRamp2(scaled_dat, c("green", "black", "red"))

    ## column annotation
    col_anno <- data.frame(
        row.names = out$PatientID,
        Subtype = out$Subtype
    )
    anno_col <- ComplexHeatmap::HeatmapAnnotation(
        df = col_anno,
        show_legend = FALSE,
        col = list(
            Subtype = c(
                "Basal" = "red",
                "Her2" = "hotpink",
                "LumA" = "darkblue",
                "LumB" = "skyblue",
                "Normal" = "green"
            )
        )
    )

    heatmap <- ComplexHeatmap::Heatmap(
        scaled_mat,
        name = "Expr",
        col = col_fun,
        ## annotation
        top_annotation = anno_col,

        ## clustering
        ## as original heatmap plot
        clustering_distance_rows = "spearman",
        clustering_method_rows = "average",
        show_row_dend = FALSE,
        cluster_column_slices = TRUE,
        column_split = col_anno$Subtype,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "average",

        ## general
        show_column_names = FALSE,
        show_heatmap_legend = TRUE,
        row_names_gp = grid::gpar(fontsize = 6)
    )

    return(heatmap)
}




#' PCA Plot Visualization of Gene Expression by Subtype
#'
#' @name Vis_PCA
#' @description This function generates a PCA plot to visualize the principal
#' components of gene expression data, colored by the assigned subtypes.
#' Optionally, it can display a scree plot of eigenvalues to evaluate the
#' explained variance.
#'
#' @param x A gene expression matrix, where genes are rows and samples are
#'   columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#' @param Eigen Logical. If `TRUE`, the function will display a scree plot
#'   showing the eigenvalues of the principal components.
#'
#' @return A `ggplot` object representing the PCA plot, colored by subtype. If
#'   `Eigen` is set to `TRUE`, a scree plot of the eigenvalues is also included.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Gene expression matrix and subtype information
#' x <- assay(OSLO2EMIT0obj$data_input$se_NC)
#' out <- data.frame(
#'     PatientID = res$results$parker.original$BS.all$PatientID,
#'     Subtype = res$results$parker.original$BS.all$BS
#' )
#'
#' # Generate the PCA plot
#' p <- Vis_PCA(x = x, out = out)
#' plot(p)
#'
#' # Generate PCA plot with scree plot of eigenvalues
#' p_with_eigen <- Vis_PCA(x = x, out = out, Eigen = TRUE)
#' plot(p_with_eigen)
#'
#' @export

Vis_PCA <- function(x, out, Eigen = FALSE) {
    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )
    pca <- prcomp(t(x), center = TRUE, scale. = TRUE)

    if (Eigen) {
        # Scree plot
        variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
        scree_data <- data.frame(PC = seq_along(variance), Variance = variance)

        screeplot <- ggplot(
            scree_data[seq(1, 10), ],
            aes(x = PC, y = Variance)
        ) +
            geom_bar(
                stat = "identity",
                width = 0.5,
                fill = "steelblue"
            ) +
            geom_line() +
            geom_point(size = 2) +
            scale_x_continuous(breaks = seq(1, 10, 1)) +
            geom_text(
                aes(
                    x = PC,
                    label = paste0(round(Variance, 2), "%")
                ),
                nudge_y = 1
            ) +
            labs(
                x = "Principal Component",
                y = "Percentage of variance Explained"
            ) +
            theme_classic() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)
            )

        return(screeplot)
    } else {
        ## PCA plot
        scores <- as.data.frame(pca$x)
        scores$PatientID <- rownames(scores)
        scores <- dplyr::left_join(scores, out, by = "PatientID")
        rownames(scores) <- scores$PatientID

        pc1_var <- round(100 * summary(pca)$importance[2, 1], 2)
        pc2_var <- round(100 * summary(pca)$importance[2, 2], 2)

        pcaplot <- ggplot(
            data = scores,
            aes(x = PC1, y = PC2, color = Subtype)
        ) +
            geom_point(size = 2) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_color_manual(name = "Subtype", values = Subtype.color) +
            labs(
                x = paste0("PC1 (", pc1_var, "% variance)"),
                y = paste0("PC2 (", pc2_var, "% variance)")
            ) +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12)
            )

        return(pcaplot)
    }
}

#' Pie Chart Visualization of Subtype Distribution
#'
#' @name Vis_pie
#' @description This function generates a pie chart to visualize the
#' distribution of breast cancer subtypes in a cohort, based on the provided
#' `Subtype` data.
#'
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#'
#' @return A `ggplot` object representing a pie chart showing the proportion of
#'   each subtype in the dataset.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Subtype information
#' out <- data.frame(
#'     PatientID = res$results$parker.original$BS.all$PatientID,
#'     Subtype = res$results$parker.original$BS.all$BS
#' )
#'
#' # Generate the pie chart
#' p <- Vis_pie(out = out)
#' plot(p)
#'
#' @export

Vis_pie <- function(out) {
    data <- data.frame(table(out$Subtype))
    data <- data %>%
        dplyr::mutate(
            perc = round(`Freq` / sum(`Freq`) * 100, 2),
            csum = rev(cumsum(rev(perc))),
            pos = perc / 2 + dplyr::lead(csum, 1),
            pos = dplyr::if_else(is.na(pos), perc / 2, pos)
        )


    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )

    plot <- ggplot(data, aes(x = "", y = perc, fill = Var1)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(name = "Subtype", values = Subtype.color) +
        geom_text_repel(
            data = data,
            aes(y = pos, label = paste0(Freq, " (", perc, "%", ")")),
            size = 4.5,
            nudge_x = 0.6,
            show.legend = FALSE,
            segment.color = NA
        ) +
        coord_polar("y", direction = -1) +
        theme_void() +
        theme(
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
        )


    return(plot)
}


#' Multi-Method Subtype Heatmap Visualization
#'
#' @name Vis_Multi
#' @description This function generates a heatmap visualization of breast cancer
#' subtypes based on multiple subtyping methods. It allows users to visualize
#' how different methods classify samples into subtypes.
#'
#' @param data A data frame or matrix that includes the subtypes predicted by
#'   different methods. The rows should represent individual samples (e.g.,
#'   `PatientID`), and columns should correspond to the subtypes predicted by
#'   different methods.
#'
#' @return A heatmap visualization of multi-method subtype results.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#'
#' # Assuming `OSLO2EMIT0obj$res$res_subtypes` contains multi-method subtype results
#' p <- Vis_Multi(OSLO2EMIT0obj$res$res_subtypes)
#' plot(p)
#'
#' @export

Vis_Multi <- function(data) {
    ## order by subtype
    data <- data[order(data[, ncol(data) - 1]), ]
    ## order by entropy
    data <- data[order(data[, ncol(data)], decreasing = FALSE), ]
    mat <- data[, -ncol(data)]

    Labels <- unique(as.vector(as.matrix(mat)))

    ## preset
    ctg <- data.frame(
        Cat = rep(c("NC-based", "SSP-based"), c(8, 2)),
        row.names = c(
            "parker.original",
            "genefu.scale",
            "genefu.robust",
            "cIHC",
            "cIHC.itr",
            "PCAPAM50",
            "ssBC",
            "ssBC.v2",
            "AIMS",
            "sspbc"
        )
    )
    ctg$Cat <- factor(ctg$Cat, levels = c("NC-based", "SSP-based"))

    ## color
    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )
    col_cat <- setNames(c("#fb9a99", "#a6cee3"), c("NC-based", "SSP-based"))

    ## make row annotation
    row_anno <- data.frame(
        Category = ctg[colnames(mat), ],
        row.names = colnames(mat)
    )
    row_anno <- HeatmapAnnotation(
        df = row_anno,
        which = c("row"),
        col = list(Category = col_cat),
        annotation_legend_param = list(
            title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
            gap = unit(2, "points"),
            labels_gp = grid::gpar(fontsize = 10),
            border = "white"
        )
    )
    # ## column annotation
    col_anno <- HeatmapAnnotation(
        which = c("column"),
        Entropy = anno_barplot(data[rownames(data), ncol(data)], bar_width = 1)
    )

    Subtype.color <- Subtype.color[names(Subtype.color) %in% Labels]
    p <- ComplexHeatmap::Heatmap(
        t(as.matrix(mat)),
        name = "Subtypes",
        col = Subtype.color,
        row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        right_annotation = row_anno,
        top_annotation = col_anno,
        show_column_names = FALSE,
        heatmap_legend_param = list(
            title = "Intrinsic Subtype",
            labels = Labels[match(names(Subtype.color), Labels)],
            title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
            gap = unit(2, "points"),
            labels_gp = grid::gpar(fontsize = 10),
            border = "white"
        )
    )

    return(p)
}
