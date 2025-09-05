# Define server logic
server <- function(input, output, session) {
  
  # Optional logo (kept from your version)
  output$logo <- renderImage({ list(src = "logo.svg", height = "100%", inline = FALSE) }, deleteFile = FALSE)
  
  # --- helpers ---
  raw_mode <- reactive({ identical(input$is_raw_counts, "raw") })
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  has_clinical <- reactive({ isTRUE((input$hasClinical %||% FALSE)) })
  
  .validate_has_clinical <- function(se, want) {
    out <- list(use = FALSE, se = se, msg = NULL)
    if (!want || is.null(se)) return(out)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    missing <- setdiff(c("TSIZE","NODE"), names(cd))
    if (length(missing) > 0) { out$msg <- paste("missing column(s):", paste(missing, collapse = ", ")); return(out) }
    # coerce to numeric & re-attach
    for (nm in c("TSIZE","NODE")) {
      if (!is.numeric(cd[[nm]])) cd[[nm]] <- suppressWarnings(as.numeric(cd[[nm]]))
    }
    SummarizedExperiment::colData(se)$TSIZE <- cd$TSIZE
    SummarizedExperiment::colData(se)$NODE  <- cd$NODE
    if (any(!is.finite(cd$TSIZE)) || any(!is.finite(cd$NODE))) {
      out$msg <- "non-numeric or NA in TSIZE/NODE"; return(out)
    }
    out$use <- TRUE; out$se <- se; out
  }
  
  # --- reactive storage ---
  reactive_files <- shiny::reactiveValues(
    GEX = NULL, clinic = NULL, anno = NULL,
    RawCounts = NULL, data_input = NULL, output_res = NULL
  )
  
  # --- dynamic help: GEX ---
  output$gex_help <- renderUI({
    if (raw_mode()) {
      tagList(
        helpText(
          strong("Raw counts mode:"),
          tags$ul(
            tags$li(tags$b("Input:"), " genes × samples ", tags$em("integer counts"), "."),
            tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
            tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
            tags$li(tags$b("Processing:"), " NC → log2-CPM (upper-quartile); SSP → FPKM (linear).")
          )
        )
      )
    } else {
      tagList(
        helpText(
          strong("Normalized data mode:"),
          tags$ul(
            tags$li(tags$b("Input:"), " genes × samples ", tags$em("log2-normalized expression"), "."),
            tags$li(tags$b("Examples:"), " ", tags$code("log2(FPKM+1)"), ", microarray/nCounter log2."),
            tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
            tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
            tags$li(tags$b("Processing:"), " NC → uses log2 values as-is; SSP → back-transforms to linear scale.")
          )
        )
      )
    }
  })
  
  # --- dynamic help: Clinical ---
  output$clin_help <- renderUI({
    use_clin <- has_clinical()
    req(TRUE)
    if (use_clin) {
      tagList(
        helpText(
          strong("Clinical metadata:"),
          tags$ul(
            tags$li(tags$b("Required:"),
                    " ",
                    tags$code("PatientID"),
                    " (must match gene-expression column names), and ",
                    tags$code("ER"),
                    " coded as ",
                    tags$code("ER+"),
                    " or ",
                    tags$code("ER-"),
                    ".")
          ),
          tags$div(strong("Additional fields (used when ‘Use clinical variables (ROR)’ is checked):")),
          tags$ul(
            tags$li(tags$code("TSIZE"), " — tumor size (",
                    tags$code("0"), " = \u2264 2 cm; ",
                    tags$code("1"), " = > 2 cm)"),
            tags$li(tags$code("NODE"), " — lymph node status (numeric: ",
                    tags$code("0"), " = negative; ",
                    tags$code("\u2265 1"), " = positive)")
          )
        )
      )
    } else {
      tagList(
        helpText(
          strong("Clinical metadata:"),
          tags$ul(
            tags$li(tags$b("Required:"),
                    " ",
                    tags$code("PatientID"),
                    " (must match gene-expression column names), and ",
                    tags$code("ER"),
                    " coded as ",
                    tags$code("ER+"),
                    " or ",
                    tags$code("ER-"),
                    ".")
          ),
          tags$div(strong("Optional (recommended if you plan to compute ROR):")),
          tags$ul(
            tags$li(tags$code("TSIZE"), " and ", tags$code("NODE"), "."),
            tags$li("Other columns are preserved in ", tags$code("colData()"), ".")
          )
        )
      )
    }
  })
  
  # --- dynamic help: Annotation ---
  output$anno_help <- renderUI({
    if (raw_mode()) {
      tagList(
        helpText(
          strong("Raw counts mode:"),
          tags$ul(
            tags$li("Annotation must include ", tags$code("probe"), " (matches GEX row names)."),
            tags$li("Annotation must include ", tags$code("ENTREZID"), " (mandatory)."),
            tags$li("Annotation must include ", tags$code("Length"), " (gene length in bp, mandatory).")
          )
        )
      )
    } else {
      tagList(
        helpText(
          strong("Normalized data mode:"),
          tags$ul(
            tags$li("Annotation must include ", tags$code("probe"), " (matches GEX row names)."),
            tags$li(tags$code("ENTREZID"), " is mandatory.")
          )
        )
      )
    }
  })
  
  # --- read files + Mapping ---
  shiny::observeEvent(input$map, {
    req(input$GEX, input$clinic, input$anno)
    
    withProgress(message = "Running Mapping...", value = 0, {
      incProgress(0.2, detail = "Loading gene expression data...")
      
      .read_table <- function(inFile) {
        ext <- tools::file_ext(inFile$name)
        sep <- if (tolower(ext) == "csv") "," else "\t"
        utils::read.table(inFile$datapath, header = TRUE, sep = sep,
                          row.names = NULL, check.names = FALSE,
                          quote = "", comment.char = "", stringsAsFactors = FALSE, fill = TRUE)
      }
      
      # GEX
      gex_df <- .read_table(input$GEX)
      if (!any(grepl("^probe$", names(gex_df), ignore.case = TRUE))) {
        if (!is.numeric(gex_df[[1]])) {
          rn <- gex_df[[1]]; gex_df[[1]] <- NULL; rownames(gex_df) <- rn
        }
      } else {
        rownames(gex_df) <- gex_df$probe; gex_df$probe <- NULL
      }
      gex_mat <- as.matrix(sapply(gex_df, function(x) suppressWarnings(as.numeric(x))))
      rownames(gex_mat) <- rownames(gex_df)
      if (anyNA(gex_mat)) {
        showNotification("Warning: NAs introduced while coercing GEX to numeric.", type = "warning")
      }
      reactive_files$GEX <- gex_mat
      
      # RawCounts flag
      reactive_files$RawCounts <- raw_mode()
      showNotification(sprintf("Data type: %s", if (reactive_files$RawCounts) "Raw counts" else "Normalized (log2)"),
                       type = "message")
      
      # Clinical
      incProgress(0.45, detail = "Loading clinical data...")
      clinic_df <- .read_table(input$clinic)
      if (!"PatientID" %in% names(clinic_df)) stop("Clinical file must contain a 'PatientID' column.")
      rownames(clinic_df) <- clinic_df$PatientID
      reactive_files$clinic <- clinic_df
      
      # Annotation
      incProgress(0.65, detail = "Loading feature annotation...")
      anno_df <- .read_table(input$anno)
      if (!all(c("probe","ENTREZID") %in% names(anno_df))) stop("Annotation file must contain 'probe' and 'ENTREZID' columns.")
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
        updateRadioButtons(session, "is_raw_counts", selected = "norm")
        removeModal()
        showNotification("Switched to normalized data mode", type = "warning")
      })
      
      # Align & validate
      incProgress(0.8, detail = "Aligning samples and features...")
      samples <- intersect(colnames(reactive_files$GEX), reactive_files$clinic$PatientID)
      if (length(samples) == 0) stop("No overlapping sample IDs between GEX columns and clinical 'PatientID'. Check headers.")
      probeID <- intersect(rownames(reactive_files$GEX), reactive_files$anno$probe)
      if (length(probeID) == 0) stop("No overlapping probes between GEX rownames and annotation 'probe'.")
      
      validate(need(!raw_mode() || "Length" %in% names(reactive_files$anno),
                    "Raw counts selected but annotation lacks ‘Length’ (bp)."))
      
      se_obj <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(counts = reactive_files$GEX[probeID, samples, drop = FALSE]),
        rowData = S4Vectors::DataFrame(reactive_files$anno[probeID, , drop = FALSE]),
        colData = S4Vectors::DataFrame(reactive_files$clinic[samples, , drop = FALSE])
      )
      
      # Run Mapping
      tryCatch({
        output_text <- capture.output({
          message("Mapping() -> RawCounts = ", reactive_files$RawCounts)
          data_input <- Mapping(se_obj = se_obj,
                                RawCounts = reactive_files$RawCounts,
                                impute = TRUE, verbose = TRUE)
          cat("Step 1 is complete.\nYou may now proceed to Step 2.\n")
        })
        reactive_files$data_input <- data_input
        showNotification(HTML(paste(output_text, collapse = "<br>")),
                         type = "message", duration = NULL)
        
        # quick peek at what Mapping returned
        nn <- names(reactive_files$data_input)
        showNotification(paste("Mapping returned:", paste(nn, collapse = ", ")),
                         type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("Error in Mapping:", e$message), type = "error", duration = NULL)
      })
      
      incProgress(1, detail = "Completed.")
    })
  })
  
  # --- run subtyping ---
  shiny::observeEvent(input$run, {
    req(reactive_files$data_input)
    
    # Normalize Mapping outputs: prefer x_*; fallback to se_* if present
    se_NC  <- reactive_files$data_input$x_NC  %||% reactive_files$data_input$se_NC
    se_SSP <- reactive_files$data_input$x_SSP %||% reactive_files$data_input$se_SSP
    
    is_nc_method <- input$BSmethod %in% c("PAM50.parker","cIHC","cIHC.itr","PCAPAM50","ssBC")
    want_clin <- isTRUE(input$hasClinical)
    chk <- .validate_has_clinical(if (is_nc_method) se_NC else NULL, want_clin)
    use_clin <- want_clin && chk$use
    se_NC    <- chk$se
    
    if (want_clin && !use_clin && is_nc_method) {
      showNotification(paste("Clinical variables disabled for this run —", chk$msg),
                       type = "warning", duration = 5)
    }
    
    showNotification(
      sprintf("Running %s | hasClinical=%s",
              input$BSmethod,
              if (use_clin) "TRUE" else "FALSE"),
      type = "message", duration = 3
    )
    
    withProgress(message = "Performing analysis...", value = 0, {
      incProgress(0.2, detail = "Initializing the method...")
      
      res <- NULL
      output_res <- NULL
      
      if (input$BSmethod == "PAM50.parker") {
        req(se_NC)
        args <- list(
          se_obj      = se_NC,
          calibration = input$calibration,
          hasClinical = use_clin
        )
        if (input$calibration == "Internal") {
          args$internal <- input$internal
        } else if (input$calibration == "External" & input$external == "Given.mdns") {
          req(input$medians)
          inFile <- input$medians
          ext <- tolower(tools::file_ext(inFile$name))
          medians <- if (ext == "csv") {
            read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",")
          } else {
            read.table(inFile$datapath, header = TRUE, row.names = NULL, sep = "\t",
                       fill = TRUE, stringsAsFactors = FALSE)
          }
          args$external <- input$external
          args$medians  <- medians
        } else if (input$calibration == "External") {
          args$external <- input$external
        }
        incProgress(0.5, detail = "Running PAM50.parker...")
        res <- do.call(BreastSubtypeR::BS_parker, args)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "cIHC") {
        req(se_NC)
        incProgress(0.5, detail = "Running cIHC...")
        res <- BreastSubtypeR::BS_cIHC(se_obj = se_NC, hasClinical = use_clin)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "cIHC.itr") {
        req(se_NC)
        incProgress(0.5, detail = "Running cIHC.itr...")
        res <- BreastSubtypeR::BS_cIHC.itr(
          se_obj = se_NC, iteration = input$iteration,
          ratio = input$ratio, hasClinical = use_clin
        )
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "PCAPAM50") {
        req(se_NC)
        incProgress(0.5, detail = "Running PCAPAM50...")
        res <- BreastSubtypeR::BS_PCAPAM50(se_obj = se_NC, hasClinical = use_clin)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "ssBC") {
        req(se_NC)
        incProgress(0.5, detail = "Running ssBC...")
        res <- BreastSubtypeR::BS_ssBC(se_obj = se_NC, s = input$s, hasClinical = use_clin)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "AIMS") {
        req(se_SSP)
        incProgress(0.5, detail = "Running AIMS...")
        data("BreastSubtypeRobj", package = "BreastSubtypeR")
        res <- BreastSubtypeR::BS_AIMS(se_obj = se_SSP)
        if (!is.null(res$cl) && is.matrix(res$cl) && ncol(res$cl) >= 1) {
          res$BS.all <- data.frame(
            PatientID = rownames(res$cl),
            BS = res$cl[, 1, drop = TRUE],
            row.names = rownames(res$cl),
            check.names = FALSE
          )
          output_res <- res$BS.all
        } else {
          output_res <- data.frame()
        }
      }
      
      if (input$BSmethod == "sspbc") {
        req(se_SSP)
        incProgress(0.5, detail = "Running sspbc...")
        res_sspbc <- BreastSubtypeR::BS_sspbc(se_SSP, ssp.name = "ssp.pam50")
        BS.all <- data.frame(
          PatientID = rownames(res_sspbc),
          BS = res_sspbc[, 1, drop = TRUE],
          row.names = rownames(res_sspbc),
          check.names = FALSE
        )
        if (identical(input$Subtype, "TRUE")) {
          res_sspbc.Subtype <- BreastSubtypeR::BS_sspbc(se_SSP, ssp.name = "ssp.subtype")
          BS.all$BS.Subtype <- res_sspbc.Subtype[, 1, drop = TRUE]
        }
        res <- list(BS.all = BS.all)
        output_res <- BS.all
      }
      
      incProgress(0.9, detail = "Finalizing the results...")
      
      # Save results for download
      reactive_files$output_res <- output_res
      output$download <- downloadHandler(
        filename = function() paste0("results-", input$BSmethod, ".txt"),
        content = function(file) {
          write.table(reactive_files$output_res, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      )
      
      # Visualization: only if we have a table with PatientID + label
      out <- NULL
      if (!is.null(res) && !is.null(res$BS.all) && nrow(res$BS.all) > 0) {
        labcol <- if ("BS" %in% names(res$BS.all)) "BS" else if ("Subtype" %in% names(res$BS.all)) "Subtype" else names(res$BS.all)[2]
        out <- data.frame(PatientID = res$BS.all$PatientID, Subtype = res$BS.all[[labcol]], check.names = FALSE)
      } else if (!is.null(res) && !is.null(res$cl) && is.matrix(res$cl) && ncol(res$cl) >= 1) {
        out <- data.frame(PatientID = rownames(res$cl), Subtype = res$cl[, 1, drop = TRUE], check.names = FALSE)
      }
      
      if (is.null(out) || nrow(out) == 0) {
        showNotification("No per-sample subtype table returned; skipping plots.", type = "warning", duration = 5)
        output$plotSection <- renderUI(NULL)
      } else {
        mat <- if (input$BSmethod %in% c("sspbc", "AIMS")) {
          req(se_SSP); log2(SummarizedExperiment::assay(se_SSP) + 1)
        } else {
          req(se_NC);  SummarizedExperiment::assay(se_NC)
        }
        if (input$BSmethod %in% c("sspbc", "AIMS")) rownames(mat) <- NULL
        
        output$pie1  <- renderPlot(BreastSubtypeR::Vis_pie(out))
        output$heat2 <- renderPlot(BreastSubtypeR::Vis_heatmap(as.matrix(mat), out = out))
        
        output$plotSection <- renderUI({
          bslib::layout_columns(
            col_width = 2,
            bslib::card(plotOutput("pie1")),
            bslib::card(plotOutput("heat2"))
          )
        })
      }
      
      incProgress(1, detail = "Analysis is complete.")
    })
    
    showNotification(HTML(paste(c("Analysis is complete.","Please review your results."), collapse = "<br>")),
                     type = "message", duration = 5)
  })
  
  # Tiny log when toggling the checkbox
  observeEvent(input$hasClinical, ignoreInit = TRUE, {
    showNotification(sprintf("Use clinical variables (ROR): %s",
                             if (isTRUE(input$hasClinical)) "ON" else "OFF"),
                     type = "message", duration = 2)
  })
}
