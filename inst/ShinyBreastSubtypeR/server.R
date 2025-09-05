# Define server logic

server <- function(input, output) {
    output$logo <- renderImage(
        {
            list(src = "logo.svg", height = "100%", inline = FALSE)
        },
        deleteFile = FALSE
    )

#### upload data ####

reactive_files <- shiny::reactiveValues(GEX = NULL, clinic = NULL, anno = NULL,
                                        RawCounts = NULL, data_input = NULL, output_res = NULL)

shiny::observeEvent(input$map, {
  req(input$GEX, input$clinic, input$anno) # Ensure files are uploaded
  
  withProgress(message = "Running Mapping...", value = 0, {
    incProgress(0.2, detail = "Loading gene expression data...")
    
    # --- helper to read csv/tsv with safe defaults ---
    .read_table <- function(inFile) {
      ext <- tools::file_ext(inFile$name)
      sep <- if (tolower(ext) == "csv") "," else "\t"
      utils::read.table(inFile$datapath, header = TRUE, sep = sep,
                        row.names = NULL, check.names = FALSE,
                        quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
    }
    
    # GEX
    gex_df <- .read_table(input$GEX)
    
    # If the first column is gene IDs (no header for that column), set rownames from it
    # Case 1: you wrote a header with ONLY sample IDs → the first data column is gene IDs
    if (!any(grepl("^probe$", names(gex_df), ignore.case = TRUE))) {
      # treat column 1 as rownames if it looks like gene IDs (non-numeric)
      if (!is.numeric(gex_df[[1]])) {
        rn <- gex_df[[1]]; gex_df[[1]] <- NULL
        rownames(gex_df) <- rn
      }
    } else {
      # If a 'probe' column exists in GEX by mistake, use it as rownames
      rownames(gex_df) <- gex_df$probe
      gex_df$probe <- NULL
    }
    # Coerce to numeric matrix
    gex_mat <- as.matrix(sapply(gex_df, function(x) suppressWarnings(as.numeric(x))))
    rownames(gex_mat) <- rownames(gex_df)
    if (anyNA(gex_mat)) {
      # NA can appear if there were stray characters; warn but continue
      showNotification("Warning: NAs introduced while coercing GEX to numeric.", type = "warning")
    }
    reactive_files$GEX <- gex_mat
    
    # RawCounts flag from UI
    reactive_files$RawCounts <- isTRUE(input$is_raw_counts)
    if (reactive_files$RawCounts) {
      showNotification("Raw counts mode selected.", type = "message")
    } else {
      showNotification("Normalized data mode selected.", type = "message")
    }
    
    # Clinical
    incProgress(0.45, detail = "Loading clinical data...")
    clinic_df <- .read_table(input$clinic)
    if (!"PatientID" %in% names(clinic_df)) {
      stop("Clinical file must contain a 'PatientID' column.")
    }
    rownames(clinic_df) <- clinic_df$PatientID
    reactive_files$clinic <- clinic_df
    
    # Annotation
    incProgress(0.65, detail = "Loading feature annotation...")
    anno_df <- .read_table(input$anno)
    if (!all(c("probe","ENTREZID") %in% names(anno_df))) {
      stop("Annotation file must contain 'probe' and 'ENTREZID' columns.")
    }
    if (reactive_files$RawCounts && !"Length" %in% names(anno_df)) {
      showModal(modalDialog(
        title = "Warning: Gene lengths not found",
        "Raw counts mode is selected but the annotation data does not contain a 'Length' column. Gene lengths are required for proper normalization.",
        easyClose = TRUE,
        footer = tagList(modalButton("Continue anyway"),
                         actionButton("switch_to_normalized", "Switch to normalized data mode"))
      ))
    }
    rownames(anno_df) <- anno_df$probe
    reactive_files$anno <- anno_df
    
    observeEvent(input$switch_to_normalized, {
      reactive_files$RawCounts <- FALSE
      removeModal()
      showNotification("Switched to normalized data mode", type = "warning")
    })
    
    ## Align identities and validate overlaps
    incProgress(0.8, detail = "Aligning samples and features...")
    
    samples <- intersect(colnames(reactive_files$GEX), reactive_files$clinic$PatientID)
    if (length(samples) == 0) {
      stop("No overlapping sample IDs between GEX columns and clinical 'PatientID'. Check headers.")
    }
    probeID <- intersect(rownames(reactive_files$GEX), reactive_files$anno$probe)
    if (length(probeID) == 0) {
      stop("No overlapping probes between GEX rownames and annotation 'probe'.")
    }
    
    # Build SE — IMPORTANT: use drop=FALSE to avoid vector collapse (dim(X) error)
    se_obj <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = reactive_files$GEX[probeID, samples, drop = FALSE]),
      rowData = S4Vectors::DataFrame(reactive_files$anno[probeID, , drop = FALSE]),
      colData = S4Vectors::DataFrame(reactive_files$clinic[samples, , drop = FALSE])
    )
    
    # Run Mapping with the correct RawCounts flag
    tryCatch({
      output_text <- capture.output({
        data_input <- Mapping(se_obj = se_obj,
                              RawCounts = reactive_files$RawCounts,
                              impute = TRUE, verbose = TRUE)
        cat("Step 1 is complete.\nYou may now proceed to Step 2.\n")
      })
      reactive_files$data_input <- data_input
      showNotification(HTML(paste(output_text, collapse = "<br>")),
                       type = "message", duration = NULL)
    }, error = function(e) {
      showNotification(paste("Error in Mapping:", e$message), type = "error", duration = NULL)
    })
    
    incProgress(1, detail = "Completed.")
  })
})
    

    #### perform analysis ####

    shiny::observeEvent(input$run, {
        # Check if data has been mapped
        req(reactive_files$data_input)

        withProgress(message = "Performing analysis...", value = 0, {
            incProgress(0.2, detail = "Initializing the method...")

            # Execute the selected method
            if (input$BSmethod == "PAM50.parker") {

                args <- list(
                    se_obj = reactive_files$data_input$se_NC,
                    calibration = input$calibration,
                    hasClinical = input$hasClinical
                )

                # Modify the list of arguments based on conditions
                if (input$calibration == "Internal") {
                    args$internal <- input$internal
                } else if (input$calibration == "External" & input$external == "Given.mdns") {
                    req(input$medians)
                    inFile <- input$medians
                    fileExt <- tools::file_ext(inFile$name)

                    # Read the file based on its extension
                    medians <- if (fileExt == "csv") {
                        read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",")
                    } else if (fileExt == "txt") {
                        read.table(inFile$datapath, header = TRUE, row.names = NULL, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
                    } else {
                        stop("Unsupported file type")
                    }

                    args$external <- input$external
                    args$medians <- medians
                } else if (input$calibration == "External" & input$external != "Given.mdns") {
                    args$external <- input$external
                }


                incProgress(0.5, detail = "Running PAM50.parker method...")
                res <- do.call(BreastSubtypeR::BS_parker, args)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "cIHC") {
                incProgress(0.5, detail = "Running cIHC method...")
                res <- BreastSubtypeR::BS_cIHC(se_obj = reactive_files$data_input$se_NC, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "cIHC.itr") {
                incProgress(0.5, detail = "Running cIHC.itr method...")
                res <- BreastSubtypeR::BS_cIHC.itr(se_obj = reactive_files$data_input$se_NC, iteration = input$iteration, ratio = input$ratio, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "PCAPAM50") {
                incProgress(0.5, detail = "Running PCAPAM50 method...")
                res <- BreastSubtypeR::BS_PCAPAM50(se_obj = reactive_files$data_input$se_NC, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "ssBC") {
                incProgress(0.5, detail = "Running ssBC method...")
                res <- BreastSubtypeR::BS_ssBC(se_obj = reactive_files$data_input$se_NC, s = input$s, hasClinical = input$hasClinical)

                output_res <- res$score.ROR
            }

            if (input$BSmethod == "AIMS") {
                incProgress(0.5, detail = "Running AIMS method...")

                data("BreastSubtypeRobj", package = "BreastSubtypeR")
                res <- BreastSubtypeR::BS_AIMS(se_obj = reactive_files$data_input$se_SSP)

                res$BS.all <- data.frame(
                    PatientID = rownames(res$cl),
                    BS = res$cl[, 1]
                )

                output_res <- res$BS.all
            }

            if (input$BSmethod == "sspbc") {
                incProgress(0.5, detail = "Running sspbc method...")

                res_sspbc <- BreastSubtypeR::BS_sspbc(reactive_files$data_input$se_SSP, ssp.name = "ssp.pam50")

                BS.all <- data.frame(PatientID = rownames(res_sspbc), BS = res_sspbc[, 1], row.names = rownames(res_sspbc))

                if (input$Subtype == "TRUE") {
                    res_sspbc.Subtype <- BreastSubtypeR::BS_sspbc(reactive_files$data_input$se_SSP, ssp.name = "ssp.subtype")

                    BS.all$BS.Subtype <- res_sspbc.Subtype[, 1]
                }
                res <- list()
                res$BS.all <- BS.all
                output_res <- BS.all
            }


            incProgress(0.9, detail = "Finalizing the results...")
            ## Allow downloading the result
            reactive_files$output_res <- output_res
            output$download <- downloadHandler(
                filename = function() {
                    paste("results-", input$BSmethod, ".txt", sep = "")
                },
                content = function(file) {
                    write.table(reactive_files$output_res, file, row.names = F, sep = "\t", quote = F)
                }
            )

            ## Allow visualization
            out <- data.frame(
                PatientID = res$BS.all$PatientID,
                Subtype = if (input$BSmethod == "sspbc" && input$Subtype == "TRUE") {
                    res$BS.all$BS.Subtype
                } else {
                    res$BS.all$BS
                }
            )
            output$pie1 <- renderPlot(BreastSubtypeR::Vis_pie(out))

            matrix <- if (input$BSmethod == "sspbc" || input$BSmethod == "AIMS") {
                log2(SummarizedExperiment::assay(reactive_files$data_input$se_SSP) + 1) %>% as.matrix()
            } else {
                SummarizedExperiment::assay(reactive_files$data_input$se_NC) %>% as.matrix()
            }

            if (input$BSmethod == "sspbc" || input$BSmethod == "AIMS") {
                rownames(matrix) <- NULL
            }

            output$heat2 <- renderPlot(BreastSubtypeR::Vis_heatmap(matrix, out = out))

            output$plotSection <- renderUI({
                req(input$run) # Wait until the run button is clicked

                # Once 'Run' is clicked, show the plot layout
                bslib::layout_columns(
                    col_width = 2,
                    bslib::card(plotOutput("pie1")),
                    bslib::card(plotOutput("heat2"))
                )
            })

            incProgress(1, detail = "Analysis is complete.")
        })

        # Show a notification when the analysis is complete
        # showNotification(  paste( "Analysis complete" ,"You can exit or continue to run another subtyping method.", sep = "\n") , type = "message", duration = 5)
        showNotification(HTML(paste(c(
            "Analysis is complete.",
            "Please review your results."
        ), collapse = "<br>")), type = "message", duration = 5)
    })
}
