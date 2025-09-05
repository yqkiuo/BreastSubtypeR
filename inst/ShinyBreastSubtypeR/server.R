# Define server logic
server <- function(input, output, session) {
  # --- helpers ---
  raw_mode     <- reactive({ identical(input$is_raw_counts, "raw") })
  `%||%`       <- function(a, b) if (!is.null(a)) a else b
  has_clinical <- reactive({ isTRUE((input$hasClinical %||% FALSE)) })
  k_is_4class  <- reactive({ identical(input$k_subtypes, "4") })
  
  # Dynamically restrict subtype choices for AIMS (AIMS = 5-class only)
  observeEvent(input$BSmethod, {
    if (identical(input$BSmethod, "AIMS")) {
      updateRadioButtons(
        session, "k_subtypes",
        choices = c("5 classes (includes Normal-like)" = "5"),
        selected = "5"
      )
    } else {
      # restore both choices, keep previous selection if valid
      prev <- isolate(input$k_subtypes) %||% "5"
      if (!prev %in% c("4","5")) prev <- "5"
      updateRadioButtons(
        session, "k_subtypes",
        choices = c("5 classes (includes Normal-like)" = "5",
                    "4 classes (excludes Normal-like)" = "4"),
        selected = prev
      )
    }
  }, ignoreInit = FALSE)
  
  .validate_has_clinical <- function(se, want) {
    out <- list(use = FALSE, se = se, msg = NULL)
    if (!want || is.null(se)) return(out)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    missing <- setdiff(c("TSIZE","NODE"), names(cd))
    if (length(missing) > 0) { out$msg <- paste("missing column(s):", paste(missing, collapse = ", ")); return(out) }
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
  
  # --- dynamic help panels (unchanged from your version) ---
  output$gex_help <- renderUI({
    badge <- function(txt) tags$span(
      txt, style = "font-size:12px; padding:2px 6px; background:#f1f3f5; border-radius:12px; margin-left:6px;"
    )
    if (raw_mode()) {
      tagList(
        div(class = "method-help",
            tags$b("Input requirements"), badge("Selected: Raw counts"),
            tags$br(),
            tags$ul(
              tags$li(tags$b("Format:"), " genes × samples ", tags$em("integer raw counts (RNA-seq)"), "."),
              tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
              tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
              tags$li(tags$b("Processing:"), " NC → log2-CPM (upper-quartile); SSP → FPKM (linear)."),
              tags$li(tags$small("Note: Annotation must include ", tags$code("ENTREZID"),
                                 " and ", tags$code("Length"), " (bp). See Feature annotation."))
            )
        )
      )
    } else {
      tagList(
        div(class = "method-help",
            tags$b("Input requirements"), badge("Selected: Normalized (log2)"),
            tags$br(),
            tags$ul(
              tags$li(tags$b("Format:"), " genes × samples ", tags$em("log2-normalized expression"), "."),
              tags$li(tags$b("Examples:"), " ", tags$code("log2(FPKM+1)"), ", microarray/nCounter log2."),
              tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
              tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
              tags$li(tags$b("Processing:"), " NC → uses log2 values as-is; SSP → back-transforms to linear scale.")
            )
        )
      )
    }
  })
  
  output$clin_help <- renderUI({
    use_clin <- has_clinical()
    if (use_clin) {
      tagList(
        div(class = "method-help",
            tags$b("Clinical data requirements"),
            tags$br(),
            HTML("<ul>
              <li><b>Required:</b> <code>PatientID</code> (matches GEX columns) and <code>ER</code> coded as <code>ER+</code> or <code>ER-</code>.</li>
              <li><b>Additional (for ROR):</b> <code>TSIZE</code> (0 = ≤ 2 cm; 1 = > 2 cm) and <code>NODE</code> (0 = negative; ≥ 1 = positive).</li>
             </ul>")
        )
      )
    } else {
      tagList(
        div(class = "method-help",
            tags$b("Clinical data requirements"),
            tags$br(),
            HTML("<ul>
              <li><b>Required:</b> <code>PatientID</code> and <code>ER</code> (ER+/ER-).</li>
              <li><b>Optional:</b> include <code>TSIZE</code> and <code>NODE</code> if you plan to compute ROR (NC methods).</li>
             </ul>")
        )
      )
    }
  })
  
  output$anno_help <- renderUI({
    if (raw_mode()) {
      tagList(
        div(class = "method-help",
            tags$b("Annotation requirements"),
            tags$br(),
            HTML("<ul>
              <li><code>probe</code> — matches GEX row names.</li>
              <li><code>ENTREZID</code> — mandatory.</li>
              <li><code>Length</code> — gene length (bp), mandatory for raw counts.</li>
             </ul>")
        )
      )
    } else {
      tagList(
        div(class = "method-help",
            tags$b("Annotation requirements"),
            tags$br(),
            HTML("<ul>
              <li><code>probe</code> — matches GEX row names.</li>
              <li><code>ENTREZID</code> — mandatory.</li>
             </ul>")
        )
      )
    }
  })
  
  # --- method descriptions (your restored citations) ---
  output$method_help <- renderUI({
    m <- input$BSmethod
    if (is.null(m)) return(NULL)
    .p <- function(title, items) {
      lis <- Map(function(k, v) tags$li(HTML(paste0("<b>", k, ":</b> ", v))), names(items), items)
      div(class = "method-help", tags$b(title), tags$br(), tags$ul(lis))
    }
    switch(m,
           "AUTO Mode" = .p(
             "AUTO Mode (cohort-aware selection)",
             list(
               "Category" = "NC-/SSP-based",
               "IHC Input Requirement" = "IHC ER and/or HER2, or TN status",
               "Description" = "Evaluates cohort diagnostics (e.g., receptor-status distribution, subtype purity, subgroup sizes) and programmatically disables classifiers whose assumptions are likely violated — reducing misclassification in skewed or small cohorts.",
               "Notes" = "AUTO may include PAM50 variants, cIHC / cIHC.itr, PCAPAM50, ssBC/ssBC.v2, AIMS, sspbc as appropriate.",
               "Key refs" = "Yang et al., 2025 NAR Genom Bioinform."
             )
           ),
           "PAM50" = .p(
             "PAM50 family (NC-based)",
             list(
               "Variants" = "parker.original | genefu.scale | genefu.robust",
               "IHC Input Requirement" = "Not required",
               "Description" = "Nearest-centroid subtyping using the PAM50 gene set; genefu variants provide alternative internal scaling/centering.",
               "Key refs" = "Parker et al., 2009 JCO; Gendoo et al., 2016 Bioinformatics."
             )
           ),
           "cIHC" = .p(
             "cIHC",
             list(
               "Category" = "NC-based",
               "IHC Input Requirement" = "IHC ER status",
               "Description" = "Conventional ER-balancing using IHC prior to PAM50-style classification.",
               "Key refs" = "Ciriello et al., 2015 Cell; Raj-Kumar et al., 2019 Sci Rep."
             )
           ),
           "cIHC.itr" = .p(
             "cIHC.itr",
             list(
               "Category" = "NC-based",
               "IHC Input Requirement" = "IHC ER status",
               "Description" = "Iterative ER-balancing variant.",
               "Key refs" = "Curtis et al., 2012 Nature."
             )
           ),
           "PCAPAM50" = .p(
             "PCAPAM50",
             list(
               "Category" = "NC-based",
               "IHC Input Requirement" = "IHC ER status → ESR1 axis",
               "Description" = "PCA on ESR1 to achieve ER-aware centering before PAM50; improves consistency with IHC.",
               "Key ref" = "Raj-Kumar et al., 2019 Sci Rep."
             )
           ),
           "ssBC" = .p(
             "ssBC",
             list(
               "Category" = "NC-based",
               "IHC Input Requirement" = "IHC ER and/or HER2, or TN status",
               "Description" = "Subgroup-specific gene-centering PAM50.",
               "Key ref" = "Zhao et al., 2015 Breast Cancer Res."
             )
           ),
           "AIMS" = .p(
             "AIMS",
             list(
               "Category" = "SSP-based",
               "IHC Input Requirement" = "Not required",
               "Description" = "Absolute Intrinsic Molecular Subtyping using pairwise gene rules.",
               "Key ref" = "Paquet & Hallett, 2015 JNCI."
             )
           ),
           "sspbc" = .p(
             "sspbc",
             list(
               "Category" = "SSP-based",
               "IHC Input Requirement" = "Not required",
               "Description" = "RNA-seq single-sample predictors for subtypes and ROR; 5-class vs 4-class via model choice.",
               "Key ref" = "Staaf et al., 2022 npj Breast Cancer."
             )
           ),
           NULL
    )
  })
  
  # Tiny toasts
  observeEvent(input$k_subtypes, ignoreInit = TRUE, {
    showNotification(sprintf("Subtype classes: %s",
                             if (k_is_4class()) "4 (Normal-like excluded)" else "5 (includes Normal-like)"),
                     type = "message", duration = 2
    )
  })
  observeEvent(input$hasClinical, ignoreInit = TRUE, {
    showNotification(sprintf("Use clinical variables (ROR): %s",
                             if (isTRUE(input$hasClinical)) "ON" else "OFF"),
                     type = "message", duration = 2
    )
  })
  
  # --- data upload + Mapping ---
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
      showNotification(sprintf("Data type: %s",
                               if (reactive_files$RawCounts) "Raw counts" else "Normalized (log2)"),
                       type = "message"
      )
      
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
      
      # Run Mapping with proper RawCounts
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
      }, error = function(e) {
        showNotification(paste("Error in Mapping:", e$message), type = "error", duration = NULL)
      })
      
      incProgress(1, detail = "Completed.")
    })
  })
  
  # --- run subtyping (AUTO + single methods) ---
  shiny::observeEvent(input$run, {
    req(reactive_files$data_input)
    
    se_NC  <- reactive_files$data_input$x_NC  %||% reactive_files$data_input$se_NC
    se_SSP <- reactive_files$data_input$x_SSP %||% reactive_files$data_input$se_SSP
    want_4 <- k_is_4class()
    
    # ---------- AUTO Mode ----------
    if (identical(input$BSmethod, "AUTO Mode")) {
      want_clin <- isTRUE(input$hasClinical)
      chk <- .validate_has_clinical(se_NC, want_clin)
      se_NC    <- chk$se
      use_clin <- want_clin && chk$use
      
      withProgress(message = "Performing analysis (AUTO Mode)...", value = 0, {
        incProgress(0.5, detail = "BS_Multi() running...")
        result_auto <- BreastSubtypeR::BS_Multi(
          data_input  = reactive_files$data_input,
          methods     = "AUTO",
          Subtype     = want_4,     # 4 vs 5 passed here
          hasClinical = use_clin
        )
        reactive_files$output_res <- result_auto$res_subtypes
        
        output$download <- downloadHandler(
          filename = function() "results-AUTO.txt",
          content  = function(file) {
            write.table(reactive_files$output_res, file, row.names = FALSE, sep = "\t", quote = FALSE)
          }
        )
        
        output$plotSection <- renderUI({
          bslib::layout_columns(col_width = 2, bslib::card(plotOutput("multi_plot", height = 480)))
        })
        output$multi_plot <- renderPlot({ BreastSubtypeR::Vis_Multi(result_auto$res_subtypes) })
        
        incProgress(1, detail = "AUTO Mode complete.")
      })
      return(invisible(NULL))
    }
    
    # ---------- Single methods ----------
    is_nc_method <- input$BSmethod %in% c("PAM50","cIHC","cIHC.itr","PCAPAM50","ssBC")
    want_clin <- isTRUE(input$hasClinical)
    chk <- .validate_has_clinical(if (is_nc_method) se_NC else NULL, want_clin)
    se_NC    <- chk$se
    use_clin <- want_clin && chk$use
    if (want_clin && !use_clin && is_nc_method) {
      showNotification(paste("Clinical variables disabled —", chk$msg), type = "warning", duration = 5)
    }
    
    showNotification(sprintf("Running %s | Subtype=%s | ROR=%s",
                             input$BSmethod,
                             if (want_4) "4-class" else "5-class",
                             if (use_clin) "ON" else "OFF"),
                     type = "message", duration = 3)
    
    withProgress(message = "Performing analysis...", value = 0, {
      incProgress(0.25, detail = "Initializing...")
      
      res <- NULL
      output_res <- NULL
      
      if (input$BSmethod == "PAM50") {
        req(se_NC)
        cal <- input$calibration %||% "Internal"
        args <- list(se_obj = se_NC, hasClinical = use_clin, Subtype = want_4)
        if (identical(cal, "Internal")) {
          args$calibration <- "Internal"; args$internal <- input$internal
        } else if (identical(cal, "External")) {
          args$calibration <- "External"
          if (identical(input$external, "Given.mdns")) {
            req(input$medians)
            inFile <- input$medians
            ext <- tolower(tools::file_ext(inFile$name))
            medians <- if (ext == "csv") {
              read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",")
            } else {
              read.table(inFile$datapath, header = TRUE, row.names = NULL, sep = "\t",
                         fill = TRUE, stringsAsFactors = FALSE)
            }
            args$external <- "Given.mdns"; args$medians  <- medians
          } else {
            args$external <- input$external
          }
        }
        incProgress(0.55, detail = "BS_parker...")
        res <- do.call(BreastSubtypeR::BS_parker, args)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "cIHC") {
        req(se_NC)
        incProgress(0.55, detail = "BS_cIHC...")
        res <- BreastSubtypeR::BS_cIHC(se_obj = se_NC, Subtype = want_4, hasClinical = use_clin)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "cIHC.itr") {
        req(se_NC)
        incProgress(0.55, detail = "BS_cIHC.itr...")
        res <- BreastSubtypeR::BS_cIHC.itr(
          se_obj = se_NC, iteration = input$iteration,
          ratio = input$ratio, Subtype = want_4, hasClinical = use_clin
        )
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "PCAPAM50") {
        req(se_NC)
        incProgress(0.55, detail = "BS_PCAPAM50...")
        res <- BreastSubtypeR::BS_PCAPAM50(se_obj = se_NC, Subtype = want_4, hasClinical = use_clin)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "ssBC") {
        req(se_NC)
        incProgress(0.55, detail = "BS_ssBC...")
        res <- BreastSubtypeR::BS_ssBC(se_obj = se_NC, s = input$s, Subtype = want_4, hasClinical = use_clin)
        output_res <- res$score.ROR
      }
      
      if (input$BSmethod == "AIMS") {
        # AIMS always 5-class, the UI already forces "5"
        req(se_SSP)
        incProgress(0.55, detail = "BS_AIMS...")
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
        incProgress(0.55, detail = "BS_sspbc...")
        model <- if (want_4) "ssp.subtype" else "ssp.pam50"  # 4 vs 5
        res_sspbc <- BreastSubtypeR::BS_sspbc(se_obj = se_SSP, ssp.name = model)
        BS.all <- data.frame(
          PatientID = rownames(res_sspbc),
          BS = res_sspbc[, 1, drop = TRUE],
          row.names = rownames(res_sspbc),
          check.names = FALSE
        )
        res <- list(BS.all = BS.all)
        output_res <- BS.all
      }
      
      incProgress(0.9, detail = "Finalizing...")
      
      # ---------- Build a single Subtype column for plotting ----------
      out <- NULL
      if (!is.null(res) && !is.null(res$BS.all) && nrow(res$BS.all) > 0) {
        # Prefer BS.Subtype (4-class) when requested and available; otherwise use BS or Subtype
        labcol <- if (want_4 && "BS.Subtype" %in% names(res$BS.all)) {
          "BS.Subtype"
        } else if ("BS" %in% names(res$BS.all)) {
          "BS"
        } else if ("Subtype" %in% names(res$BS.all)) {
          "Subtype"
        } else {
          names(res$BS.all)[2]
        }
        out <- data.frame(
          PatientID = res$BS.all$PatientID,
          Subtype   = res$BS.all[[labcol]],
          check.names = FALSE
        )
      } else if (!is.null(res) && !is.null(res$cl) && is.matrix(res$cl) && ncol(res$cl) >= 1) {
        # AIMS path (cl matrix)
        out <- data.frame(
          PatientID = rownames(res$cl),
          Subtype   = res$cl[, 1, drop = TRUE],
          check.names = FALSE
        )
      }
      
      # ---------- Save & download (export both labels when present) ----------
      if (!is.null(out) && nrow(out) > 0) {
        # Start with selected label
        download_tbl <- data.frame(PatientID = out$PatientID, Subtype = out$Subtype, check.names = FALSE)
        
        # If the model returned both 5-class and 4-class, include them as well
        if (!is.null(res$BS.all)) {
          if ("BS" %in% names(res$BS.all))        download_tbl$BS_5class <- res$BS.all$BS
          if ("BS.Subtype" %in% names(res$BS.all)) download_tbl$BS_4class <- res$BS.all$BS.Subtype
        }
        
        reactive_files$output_res <- download_tbl
      } else {
        reactive_files$output_res <- data.frame()
      }
      
      output$download <- downloadHandler(
        filename = function() paste0("results-", input$BSmethod, ".txt"),
        content  = function(file) {
          write.table(reactive_files$output_res, file, row.names = FALSE, sep = "\t", quote = FALSE)
        }
      )
      
      # ---------- Visualization ----------
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
          bslib::layout_columns(col_width = 2,
                                bslib::card(plotOutput("pie1")),
                                bslib::card(plotOutput("heat2")))
        })
      }
      
      incProgress(1, detail = "Analysis complete.")
    })
  })
}
