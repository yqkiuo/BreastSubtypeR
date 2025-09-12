# --- Safe CSS unit helper (optional) ---
safeCssUnit <- function(x, fallback = "auto") {
  if (is.null(x) || (length(x) == 0)) return(fallback)
  tryCatch(shiny::validateCssUnit(x),
           error = function(e) if (is.numeric(x) && is.finite(x)) paste0(x, "px") else fallback)
}

# --- bslib::callout shim ---
bs_callout <- function(title, body, color = c("info","success","warning","danger","secondary")) {
  color <- match.arg(color)
  if (requireNamespace("bslib", quietly = TRUE)) {
    ns <- asNamespace("bslib")
    if (exists("callout", envir = ns, inherits = FALSE)) {
      return(get("callout", envir = ns)(title = title, body, color = color))
    }
  }
  cl <- switch(color, "success"="alert alert-success","warning"="alert alert-warning",
               "danger"="alert alert-danger","secondary"="alert alert-secondary","info"="alert alert-info")
  shiny::tags$div(class = cl, role = "alert", shiny::tags$div(shiny::tags$b(title)), body)
}

server <- function(input, output, session) {
  
  # --- small helpers ---
  raw_mode     <- reactive({ identical(input$is_raw_counts, "raw") })
  `%||%`       <- function(a, b) if (!is.null(a)) a else b
  has_clinical <- reactive({ isTRUE((input$hasClinical %||% FALSE)) })
  k_is_4class  <- reactive({ identical(input$k_subtypes, "4") })

  # --- DOI / version info for exports ---
  PKG_VER      <- as.character(utils::packageVersion("BreastSubtypeR"))
  CITATION_DOI <- "10.1093/nargab/lqaf131"
  CITATION_NOTE <- sprintf(
    "BreastSubtypeR v%s — NAR Genomics and Bioinformatics (2025) DOI: %s",
    PKG_VER, CITATION_DOI
  )
  
  # Write a TSV with one or more header lines (each prefixed with "# ")
  .write_tsv_with_header <- function(path, df, header_lines = CITATION_NOTE) {
    con <- file(path, open = "wt", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    if (length(header_lines)) writeLines(paste0("# ", header_lines), con)
    utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  # Map PAM50 UI choices to a clear export header
  .pam50_variant_from_internal <- function(internal) {
    switch(internal,
           "medianCtr" = "parker.original",
           "meanCtr"   = "genefu.scale",
           "qCtr"      = "genefu.robust",
           "parker.original"
    )
  }
  
  .method_header_line <- function(input, want_4, use_clin) {
    method_label <- input$BSmethod %||% "Unknown"
    
    if (identical(method_label, "PAM50")) {
      cal <- input$calibration %||% "None"
      if (identical(cal, "Internal")) {
        internal <- input$internal %||% "medianCtr"
        variant  <- .pam50_variant_from_internal(internal)
        method_label <- sprintf("PAM50 (variant: %s, calibration: Internal (%s))", variant, internal)
      } else if (identical(cal, "External")) {
        ext <- input$external %||% "Unknown"
        ext_label <- if (identical(ext, "Given.mdns")) "Custom (Given.mdns)" else ext
        method_label <- sprintf("PAM50 (variant: parker.original, calibration: External (%s))", ext_label)
      } else {
        method_label <- "PAM50 (variant: parker.original, calibration: None)"
      }
    }
    
    sprintf("Method: %s | Classes: %s | ROR: %s",
            method_label,
            if (want_4) "4" else "5",
            if (use_clin) "ON" else "OFF")
  }
  
  # Build a safe data.frame with PatientID prepended (avoid cbind(..., check.names=...))
  .build_out_df <- function(ids, tab) {
    tab <- as.data.frame(tab, stringsAsFactors = FALSE, check.names = FALSE)
    data.frame(
      PatientID = ids,
      tab,
      row.names = NULL,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  # Build Calls-only and Full metrics for single-method results
  .build_single_exports <- function(res, want_4) {
    stopifnot(!is.null(res$BS.all), "PatientID" %in% names(res$BS.all))
    calls_col <- if (want_4 && "BS.Subtype" %in% names(res$BS.all)) "BS.Subtype"
    else if ("BS" %in% names(res$BS.all)) "BS"
    else setdiff(names(res$BS.all), "PatientID")[1]
    nm <- if (want_4) "Call_4class" else "Call_5class"
    
    calls_tbl <- data.frame(PatientID = res$BS.all$PatientID,
                            stringsAsFactors = FALSE, check.names = FALSE)
    calls_tbl[[nm]] <- res$BS.all[[calls_col]]
    
    full_tbl <- as.data.frame(res$BS.all, stringsAsFactors = FALSE, check.names = FALSE)
    if ("BS" %in% names(full_tbl))        { full_tbl$Call_5class <- full_tbl$BS;         full_tbl$BS <- NULL }
    if ("BS.Subtype" %in% names(full_tbl)){ full_tbl$Call_4class <- full_tbl$BS.Subtype; full_tbl$BS.Subtype <- NULL }
    
    if (!is.null(res$score.ROR) && nrow(res$score.ROR) > 0) {
      sr <- as.data.frame(res$score.ROR, stringsAsFactors = FALSE, check.names = FALSE)
      if (!"PatientID" %in% names(sr)) {
        if ("sample" %in% names(sr)) names(sr)[names(sr) == "sample"] <- "PatientID" else sr$PatientID <- rownames(sr)
      }
      keep <- setdiff(names(sr), names(full_tbl))
      if (length(keep)) {
        full_tbl <- merge(full_tbl, sr[, c("PatientID", keep), drop = FALSE],
                          by = "PatientID", all.x = TRUE, sort = FALSE)
      }
    }
    list(calls = calls_tbl, full = full_tbl)
  }
  
  # Named-list badges
  badge <- function(txt) tags$span(
    txt, style = "font-size:12px; padding:2px 6px; background:#f1f3f5; border-radius:12px; margin-left:6px;"
  )
  
  .validate_has_clinical <- function(se, want) {
    out <- list(use = FALSE, se = se, msg = NULL)
    if (!want || is.null(se)) return(out)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    missing <- setdiff(c("TSIZE","NODE"), names(cd))
    if (length(missing) > 0) { out$msg <- paste("missing column(s):", paste(missing, collapse = ", ")); return(out) }
    for (nm in c("TSIZE","NODE")) if (!is.numeric(cd[[nm]])) cd[[nm]] <- suppressWarnings(as.numeric(cd[[nm]]))
    SummarizedExperiment::colData(se)$TSIZE <- cd$TSIZE
    SummarizedExperiment::colData(se)$NODE  <- cd$NODE
    if (any(!is.finite(cd$TSIZE)) || any(!is.finite(cd$NODE))) { out$msg <- "non-numeric or NA in TSIZE/NODE"; return(out) }
    out$use <- TRUE; out$se <- se; out
  }
  
  # --- reactive storage ---
  reactive_files <- shiny::reactiveValues(
    GEX = NULL, clinic = NULL, anno = NULL,
    RawCounts = NULL, data_input = NULL,
    download_calls = NULL, download_full = NULL
  )
  
  # Minimal DOI helper for compact links
  doi_tag <- function(doi) {
    sprintf("<a href='https://doi.org/%s' target='_blank' rel='noopener noreferrer'>doi</a>", doi)
  }
  
  # --- dynamic help panels ---
  output$gex_help <- renderUI({
    if (raw_mode()) {
      tagList(
        div(class = "method-help",
            tags$b("Expression matrix requirements"), badge("Selected: Raw counts"),
            tags$br(),
            tags$ul(
              tags$li(tags$b("Format:"), " genes × samples ", tags$em("integer raw counts (RNA-seq)"), "."),
              tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
              tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
              tags$li(tags$b("Processing:"), " NC → log2-CPM (upper-quartile); SSP → FPKM (linear)."),
              tags$li(tags$small("Note: Annotation must include ", tags$code("ENTREZID"),
                                 " and ", tags$code("Length"), " (bp). See Feature annotation."))
            ))
      )
    } else {
      tagList(
        div(class = "method-help",
            tags$b("Expression matrix requirements"), badge("Selected: Normalized (log2)"),
            tags$br(),
            tags$ul(
              tags$li(tags$b("Format:"), " genes × samples ", tags$em("log2-normalized expression"), "."),
              tags$li(tags$b("Examples:"), " ", tags$code("log2(FPKM+1)"), ", microarray/nCounter log2."),
              tags$li(tags$b("Rows:"), " gene IDs matching ", tags$code("annotation$probe"), "."),
              tags$li(tags$b("Columns:"), " sample IDs matching ", tags$code("clinical$PatientID"), "."),
              tags$li(tags$b("Processing:"), " NC → uses log2 values as-is; SSP → back-transforms to linear scale.")
            ))
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
              <li><b>Minimum:</b> <code>PatientID</code> (matches GEX columns).</li>
              <li><b>Method-specific:</b> 
                cIHC / cIHC.itr / PCAPAM50 need <code>ER</code> (values <code>ER+</code>/<code>ER-</code>). 
                ssBC requires <code>ER</code> (and <code>HER2</code> for <code>ER.v2</code>, values <code>HER2+</code>/<code>HER2-</code>) or <code>TN</code> (values <code>TN</code>/<code>nonTN</code>) depending on subgroup.
              </li>
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
              <li><b>Minimum:</b> <code>PatientID</code>.</li>
              <li><b>Method-specific:</b> 
                cIHC / cIHC.itr / PCAPAM50 need <code>ER</code> (values <code>ER+</code>/<code>ER-</code>). 
                ssBC requires <code>ER</code> (and <code>HER2</code> for <code>ER.v2</code>) or <code>TN</code> depending on subgroup.
              </li>
              <li><b>Optional for ROR (NC methods):</b> include <code>TSIZE</code> and <code>NODE</code>.</li>
             </ul>")
        )
      )
    }
  })
  
  output$anno_help <- renderUI({
    if (raw_mode()) {
      tagList(
        div(class = "method-help",
            tags$b("Feature annotation requirements"),
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
            tags$b("Feature annotation requirements"),
            tags$br(),
            HTML("<ul>
              <li><code>probe</code> — matches GEX row names.</li>
              <li><code>ENTREZID</code> — mandatory.</li>
             </ul>")
        )
      )
    }
  })
  
  output$calib_help <- renderUI({
    if (!identical(input$BSmethod, "PAM50")) return(NULL)
    
    tagList(
      div(class = "method-help",
          tags$b("Calibration notes (PAM50)"),
          tags$ul(
            tags$li(HTML("<b>None</b>: no centering (use only if your expression scale already matches the training data).")),
            
            tags$li(
              HTML("<b>Internal</b>: center your cohort; choose one of:"),
              tags$ul(
                tags$li(HTML("<code>medianCtr</code> — gene-wise <i>median</i> centering (Parker original).")),
                tags$li(HTML("<code>meanCtr</code> — gene-wise z-score (mean&nbsp;0, sd&nbsp;1).")),
                tags$li(HTML("<code>qCtr</code> — robust quantile re-centering (mq&nbsp;&approx;&nbsp;0.05)."))
              )
            ),
            
            tags$li(
              HTML("<b>External</b>: subtract <i>reference medians</i> from a training cohort/platform."),
              tags$ul(
                tags$li(HTML("<b>Reference medians</b> selector: choose a <i>Built-in</i> preset (e.g., <code>RNAseq.V2</code>, <code>nCounter</code>, <code>Agilent_244K</code>) or <b>Custom (upload file…)</b>.")),
                tags$li(HTML("<b>Custom upload</b>: CSV/TXT with two columns — <code>X</code> (PAM50 gene symbol) and <code>Given.mdns</code> (median log2 expression). Exactly 50 unique symbols; extra rows are ignored; gene symbols are case-sensitive."))
              )
            ),
            
            tags$li(HTML("<b>Gene overlap</b>: all steps operate on the intersection of your genes and PAM50; missing genes are ignored."))
          )
      )
    )
  })
  
  # --- method descriptions ---
  output$method_help <- renderUI({
    m <- input$BSmethod; if (is.null(m)) return(NULL)
    .p <- function(title, items) {
      lis <- Map(function(k, v) tags$li(HTML(paste0("<b>", k, ":</b> ", v))), names(items), items)
      div(class = "method-help", tags$b(title), tags$br(), tags$ul(lis))
    }
    switch(m,
           "AUTO Mode" = .p("AUTO Mode (cohort-aware selection)", list(
             "Category" = "NC-/SSP-based",
             "IHC Input Requirement" = "ER status and/or HER2, or TN status",
             "Description" = "Evaluates cohort diagnostics and disables classifiers with violated assumptions.",
             "Notes" = "ROR is disabled in AUTO.",
             "Key refs" = paste0(
               "Yang et al., 2025 NAR Genom Bioinform. (", doi_tag("10.1093/nargab/lqaf131"), ")"
             )
           )),
           "PAM50" = .p("PAM50 family (NC-based)", list(
             "Variants" = "parker.original | genefu.scale | genefu.robust",
             "IHC Input Requirement" = "Not required",
             "Description" = "Nearest-centroid subtyping using the PAM50 gene set; multiple centering options.",
             "Key refs" = paste0(
               "Parker et al., 2009 JCO (", doi_tag("10.1200/JCO.2008.18.1370"), "); ",
               "Gendoo et al., 2016 Bioinformatics (", doi_tag("10.1093/bioinformatics/btv693"), ")"
             )
           )),
           "cIHC" = .p("cIHC", list(
             "Category" = "NC-based",
             "IHC Input Requirement" = "ER status",
             "Description" = "Conventional ER-balancing with PAM50-style classification.",
             "Key refs" = paste0(
               "Ciriello et al., 2015 Cell (", doi_tag("10.1016/j.cell.2015.09.033"), ") "
             )
           )),
           "cIHC.itr" = .p("cIHC.itr", list(
             "Category" = "NC-based",
             "IHC Input Requirement" = "ER status",
             "Description" = "Iterative ER-balancing variant.",
             "Key refs" = paste0(
               "Curtis et al., 2012 Nature (", doi_tag("10.1038/nature10983"), ")"
             )
           )),
           "PCAPAM50" = .p("PCAPAM50", list(
             "Category" = "NC-based",
             "IHC Input Requirement" = "ER status (IHC) → ESR1 axis",
             "Description" = "PCA on ESR1 to achieve ER-aware centering before PAM50.",
             "Key refs" = paste0(
               "Raj-Kumar et al., 2019 Sci Rep (", doi_tag("10.1038/s41598-019-44339-4"), ")"
             )
           )),
           "ssBC" = .p("ssBC", list(
             "Category" = "NC-based",
             "IHC Input Requirement" = "ER (and HER2 for ER.v2) or TN depending on subgroup.",
             "Description" = "Subgroup-specific gene-centering PAM50.",
             "Key refs" = paste0(
               "Zhao et al., 2015 BCR (", doi_tag("10.1186/s13058-015-0520-4"), "); ",
               "Fernandez-Martinez et al., 2020 JCO (", doi_tag("10.1200/JCO.20.01276"), ")"
             )
           )),
           "AIMS" = .p("AIMS", list(
             "Category" = "SSP-based",
             "IHC Input Requirement" = "Not required",
             "Description" = "Absolute Intrinsic Molecular Subtyping (pairwise rules).",
             "Key refs" = paste0(
               "Paquet & Hallett, 2015 JNCI (", doi_tag("10.1093/jnci/dju357"), ")"
             )
           )),
           "sspbc" = .p("sspbc", list(
             "Category" = "SSP-based",
             "IHC Input Requirement" = "Not required",
             "Description" = "SCAN-B SSP models (4- or 5-class via model).",
             "Key refs" = paste0(
               "Staaf et al., 2022 NPJ Breast Cancer (", doi_tag("10.1038/s41523-022-00465-3"), ")"
             )
           )),
           NULL)
    
  })
  
  # --- Smart AUTO chip (dynamic) ---
  output$auto_chip <- renderUI({
    # fall back icon if FA6 is missing
    wand_ok <- tryCatch({ htmltools::tagList(icon("wand-magic-sparkles")); TRUE }, error = function(...) FALSE)
    icn <- if (wand_ok) icon("wand-magic-sparkles") else icon("magic")
    
    chk <- auto_requirements()  # uses existing reactive
    cls <- if (isTRUE(chk$ready)) "chip auto ready" else "chip auto blocked"
    state <- if (isTRUE(chk$ready)) "ready" else "action needed"
    note  <- if (isTRUE(chk$show)) (chk$msg %||% "OK") else "Run Step 1 first"
    
    # little “Why AUTO?” link opens a modal with cohort-aware details
    htmltools::tags$span(
      class = cls,
      icn,
      htmltools::HTML("<b>AUTO:</b> assumption-aware selection"),
      htmltools::span(class = "chip-note", paste("—", state)),
      actionLink("show_auto_modal", label = "Why AUTO?", icon = icon("circle-question"), style = "margin-left:6px;")
    )
  })
  
  # Lock 4/5-class choice for AIMS (AIMS = 5-class only)
  observeEvent(input$BSmethod, {
    if (identical(input$BSmethod, "AIMS")) {
      updateRadioButtons(
        session, "k_subtypes",
        choices  = list("5 classes (includes Normal-like)" = "5"),
        selected = "5"
      )
    } else {
      prev <- isolate(input$k_subtypes)
      if (!prev %in% c("4","5")) prev <- "5"
      updateRadioButtons(
        session, "k_subtypes",
        choices  = list("5 classes (includes Normal-like)" = "5",
                        "4 classes (excludes Normal-like)" = "4"),
        selected = prev
      )
    }
  }, ignoreInit = FALSE)
  
  # Methods where "Full metrics" should be hidden
  block_full <- reactive({
    input$BSmethod %in% c("AUTO Mode", "AIMS", "sspbc")
  })
  
  # Dynamic export selector (shows only "Calls only" for AUTO/AIMS/SSPBC)
  output$export_selector <- renderUI({
    show_only_calls <- isTRUE(block_full())
    choices <- if (show_only_calls) {
      list("Calls only" = "calls")
    } else {
      list("Calls only" = "calls", "Full metrics (incl. ROR when available)" = "full")
    }
    # default to "calls" if option not present or not yet set
    sel <- if (show_only_calls) "calls" else (isolate(input$export_kind) %||% "calls")
    radioButtons("export_kind", "Export content", choices = choices, inline = TRUE, selected = sel)
  })
  
  
  # --- Cohort-aware AUTO explainer modal ---
  observeEvent(input$show_auto_modal, {
    di <- reactive_files$data_input
    have_data <- !is.null(di) && !is.null(di$se_NC)
    details_ui <- if (!have_data) {
      # general explanation (no data yet)
      tagList(
        tags$p("AUTO analyzes your cohort to choose compatible classifiers and avoid misclassification in skewed or small cohorts."),
        tags$ul(
          tags$li("Checks ER/HER2 presence and coding (ER+/ER-, HER2+/HER2-)."),
          tags$li("Screens for subgroup imbalance (e.g., ER+ ≫ ER-), subtype purity and subgroup size."),
          tags$li("Disables methods whose assumptions are violated; runs robust alternatives (e.g., ssBC, AIMS/SSPBC)."),
          tags$li("Keeps outputs per method.")
        ),
        tags$p(tags$em("Tip: Run Step 1 so AUTO can show a cohort snapshot here."))
      )
    } else {
      se <- di$se_NC
      ph <- as.data.frame(SummarizedExperiment::colData(se))
      n  <- nrow(ph)
      ER_tab   <- if ("ER"   %in% names(ph)) table(ph$ER, useNA = "ifany")   else NULL
      HER2_tab <- if ("HER2" %in% names(ph)) table(ph$HER2, useNA = "ifany") else NULL
      
      # Try to preview which methods AUTO would keep (best-effort)
      kept <- tryCatch({
        gm <- get_methods(ph)  # from package
        gm$methods
      }, error = function(e) NULL)
      
      tagList(
        tags$p("AUTO inspects your cohort and selects only the classifiers whose assumptions match your data. Quick snapshot:"),
        tags$ul(
          tags$li(HTML(sprintf("<b>Samples:</b> %s", n))),
          if (!is.null(ER_tab))
            tags$li(HTML(sprintf("<b>ER distribution:</b> %s",
                                 paste(sprintf("%s=%s", names(ER_tab), as.integer(ER_tab)), collapse = ", ")))),
          if (!is.null(HER2_tab))
            tags$li(HTML(sprintf("<b>HER2 distribution:</b> %s",
                                 paste(sprintf("%s=%s", names(HER2_tab), as.integer(HER2_tab)), collapse = ", "))))
        ),
        if (!is.null(kept))
          tags$p(HTML(sprintf("<b>Selected methods for this cohort:</b> %s", paste(kept, collapse = " · ")))),
        tags$hr(),
        tags$p("How AUTO decides"),
        tags$ul(
          tags$li("Checks ER/HER2 imbalance, subtype purity, and subgroup size."),
          tags$li("Disables classic NC-based PAM50 variants when ER/HER2 appear markedly imbalanced."),
          tags$li("Prefers subgroup-specific gene centring (ssBC / ssBC.v2) when a subtype-specific cohort is detected (e.g., ER+/HER2-, HER2+)."),
          tags$li("Retains cohort-agnostic SSPs (AIMS / SSPBC) as robust baselines."),
          tags$li("Returns per-method calls; an entropy score summarises agreement across methods."),
          tags$li("ROR is not computed in AUTO.")
        )
      )
      
    }
    
    showModal(modalDialog(
      title = "AUTO: assumption-aware method selection",
      size  = "l",
      easyClose = TRUE,
      details_ui,
      footer = tagList(modalButton("Close"))
    ))
  })
  
  
  # Tiny toasts
  observeEvent(input$k_subtypes, ignoreInit = TRUE, {
    showNotification(sprintf("Subtype classes: %s",
                             if (k_is_4class()) "4 (Normal-like excluded)" else "5 (includes Normal-like)"),
                     type = "message", duration = 2)
  })
  observeEvent(input$hasClinical, ignoreInit = TRUE, {
    showNotification(sprintf("Use clinical variables (ROR): %s",
                             if (isTRUE(input$hasClinical)) "ON" else "OFF"),
                     type = "message", duration = 2)
  })
  
  observeEvent(input$about, {
    showModal(modalDialog(
      title = "About & Citation",
      easyClose = TRUE, footer = NULL,
      HTML(paste0(
        "<p><b>BreastSubtypeR</b> — intrinsic molecular subtyping for breast cancer (R/Bioconductor).</p>",
        "<p><b>Please cite:</b><br/>",
        "Yang Q., Hartman J., Sifakis E.G. ",
        "<em>BreastSubtypeR: A Unified R/Bioconductor Package for Intrinsic Molecular Subtyping in Breast Cancer Research</em>. ",
        "<em>NAR Genomics and Bioinformatics</em> (2025). ",
        "<a href='https://doi.org/10.1093/nargab/lqaf131' target='_blank' rel='noopener noreferrer'>https://doi.org/10.1093/nargab/lqaf131</a></p>"
      ))
    ))
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
      if (any(grepl("^probe$", names(gex_df), ignore.case = TRUE))) {
        probe_col <- grep("^probe$", names(gex_df), ignore.case = TRUE)[1]
        rownames(gex_df) <- gex_df[[probe_col]]
        gex_df[[probe_col]] <- NULL
      } else {
        rn <- gex_df[[1]]; gex_df[[1]] <- NULL; rownames(gex_df) <- as.character(rn)
      }
      gex_mat <- as.matrix(suppressWarnings(sapply(gex_df, function(x) as.numeric(x))))
      rownames(gex_mat) <- rownames(gex_df)
      if (anyNA(gex_mat)) showNotification("Warning: NAs introduced while coercing GEX to numeric.", type = "warning")
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
      if (length(samples) == 0) stop("No overlapping sample IDs between GEX columns and clinical 'PatientID'.")
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
          data_input <- Mapping(se_obj = se_obj, RawCounts = reactive_files$RawCounts, impute = TRUE, verbose = TRUE)
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
  
  # --- persistent AUTO preflight status (panel in UI) ---
  auto_requirements <- reactive({
    if (!identical(input$BSmethod, "AUTO Mode")) return(list(show = FALSE))
    di <- reactive_files$data_input
    if (is.null(di) || is.null(di$se_NC)) return(list(show = TRUE, ready = FALSE, msg = "Run Step 1 (Preprocess & map) first."))
    pheno <- as.data.frame(SummarizedExperiment::colData(di$se_NC))
    msgs <- character(0)
    miss <- setdiff(c("ER","HER2"), names(pheno))
    if (length(miss)) msgs <- c(msgs, paste0("Missing columns: ", paste(miss, collapse = ", ")))
    if ("ER" %in% names(pheno)) {
      badER <- setdiff(na.omit(unique(pheno$ER)), c("ER+","ER-"))
      if (length(badER)) msgs <- c(msgs, paste0("Invalid ER values: ", paste(badER, collapse = ", ")))
    }
    if ("HER2" %in% names(pheno)) {
      badH <- setdiff(na.omit(unique(pheno$HER2)), c("HER2+","HER2-"))
      if (length(badH)) msgs <- c(msgs, paste0("Invalid HER2 values: ", paste(badH, collapse = ", ")))
    }
    list(show = TRUE, ready = length(msgs) == 0, msg = if (length(msgs)) paste(msgs, collapse = "; ") else "OK")
  })
  
  output$auto_preflight <- renderUI({
    chk <- auto_requirements()
    if (!isTRUE(chk$show)) return(NULL)
    if (isTRUE(chk$ready)) {
      bs_callout("AUTO preflight: ready",
                 "ER and HER2 present with valid coding (ER+/ER-, HER2+/HER2-). You can run AUTO.",
                 "success")
    } else {
      bs_callout("AUTO preflight: action needed",
                 HTML(sprintf("<p><b>Fix these before running AUTO:</b></p><ul><li>%s</li></ul><p>Or pick a single method from the dropdown.</p>", chk$msg)),
                 "danger")
    }
  })
  
  # --- run subtyping (AUTO + single methods) ---
  shiny::observeEvent(input$run, {
    req(reactive_files$data_input)
    
    se_NC  <- reactive_files$data_input$x_NC  %||% reactive_files$data_input$se_NC
    se_SSP <- reactive_files$data_input$x_SSP %||% reactive_files$data_input$se_SSP
    want_4 <- k_is_4class()
    if (identical(input$BSmethod, "AIMS")) {
      want_4 <- FALSE  # enforce 5-class for AIMS
    }
    
    # ---------- AUTO Mode ----------
    if (identical(input$BSmethod, "AUTO Mode")) {
      # Gate AUTO with preflight
      chk <- auto_requirements()
      if (!isTRUE(chk$ready)) {
        showNotification(HTML(sprintf(
          "AUTO requires <code>ER</code> and <code>HER2</code> coded as ER+/ER- and HER2+/HER2-. Issue: %s.", chk$msg)),
          type = "error", duration = 10)
        showModal(modalDialog(
          title = "AUTO Mode requirements not met",
          div(HTML(sprintf(
            "<p><b>AUTO</b> requires both <code>ER</code> and <code>HER2</code> in the clinical data.</p>
             <ul><li>ER must be coded <code>ER+</code>/<code>ER-</code>.</li>
                 <li>HER2 must be coded <code>HER2+</code>/<code>HER2-</code>.</li></ul>
             <p><b>Issue:</b> %s</p>
             <p>You can either fix the clinical file and re-run, or choose a single method from the dropdown.</p>", chk$msg))),
          easyClose = TRUE, footer = tagList(modalButton("OK")))
        )
        return(invisible(NULL))
      }
      
      # ROR is disabled in AUTO
      use_clin <- FALSE
      
      withProgress(message = "Performing analysis (AUTO Mode)...", value = 0, {
        incProgress(0.5, detail = "BS_Multi() running...")
        result_auto <- BreastSubtypeR::BS_Multi(
          data_input  = reactive_files$data_input,
          methods     = "AUTO",
          Subtype     = want_4,   # TRUE -> 4-class; FALSE -> 5-class
          hasClinical = use_clin
        )
        
        # Table for plotting (chosen k)
        tab <- if (isTRUE(want_4)) result_auto$res_subtypes.Subtype else result_auto$res_subtypes
        
        # ---- Build exports for AUTO ----
        # Calls-only = chosen k table (+ entropy)
        calls_tbl_auto <- .build_out_df(rownames(tab), tab)
        
        # Full metrics = chosen k table; if the other k exists, append it with suffixes
        other_tab <- if (isTRUE(want_4)) (result_auto$res_subtypes %||% NULL) else (result_auto$res_subtypes.Subtype %||% NULL)
        full_tbl_auto <- calls_tbl_auto
        if (!is.null(other_tab)) {
          a <- calls_tbl_auto
          b <- .build_out_df(rownames(other_tab), other_tab)
          # add suffixes to b's method columns + entropy
          nonid <- setdiff(names(b), "PatientID")
          names(b)[names(b) %in% nonid] <- paste0(names(b)[names(b) %in% nonid],
                                                  if (want_4) "_5class" else "_4class")
          # also suffix a
          nonid.a <- setdiff(names(a), "PatientID")
          names(a)[names(a) %in% nonid.a] <- paste0(names(a)[names(a) %in% nonid.a],
                                                    if (want_4) "_4class" else "_5class")
          full_tbl_auto <- merge(a, b, by = "PatientID", all = TRUE, sort = FALSE)
        }
        
        # Save export tables
        reactive_files$download_calls <- calls_tbl_auto
        reactive_files$download_full  <- NULL
        
        # Downloads
        output$download <- downloadHandler(
          filename = function() {
            kind <- if (!isTRUE(block_full()) && identical(input$export_kind, "full")) "full" else "calls"
            sprintf("results-AUTO-%s-%s.txt", if (want_4) "4class" else "5class", kind)
          },
          content = function(file) {
            want_full <- identical(input$export_kind, "full") && !isTRUE(block_full())
            obj <- if (want_full && !is.null(reactive_files$download_full) && nrow(reactive_files$download_full) > 0) {
              reactive_files$download_full
            } else {
              reactive_files$download_calls
            }
            hdr <- c(
              CITATION_NOTE,
              sprintf("Method: AUTO | Classes: %s | ROR: OFF", if (want_4) "4" else "5"),
              sprintf("Export: %s | Date: %s | R: %s",
                      if (want_full) "full" else "calls",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                      getRversion())
            )
            .write_tsv_with_header(file, obj, header_lines = hdr)
          }
        )
        
        # Visualization
        output$plotSection <- renderUI({
          bslib::layout_columns(col_width = 2, bslib::card(plotOutput("multi_plot", height = 480)))
        })
        output$multi_plot <- renderPlot({
          BreastSubtypeR::Vis_Multi(tab)   # unchanged; expects the per-method calls matrix
        })
        
        incProgress(1, detail = "AUTO Mode complete.")
      })
      return(invisible(NULL))
    }
    
    # ---------- Single methods ----------
    is_nc_method <- input$BSmethod %in% c("PAM50","cIHC","cIHC.itr","PCAPAM50","ssBC")
    want_clin <- isTRUE(input$hasClinical)
    chkclin <- .validate_has_clinical(if (is_nc_method) se_NC else NULL, want_clin)
    se_NC    <- chkclin$se
    use_clin <- want_clin && chkclin$use
    if (want_clin && !use_clin && is_nc_method) {
      showNotification(paste("Clinical variables disabled —", chkclin$msg), type = "warning", duration = 5)
    }
    
    showNotification(sprintf("Running %s | Subtype=%s | ROR=%s",
                             input$BSmethod,
                             if (want_4) "4-class" else "5-class",
                             if (use_clin) "ON" else "OFF"),
                     type = "message", duration = 3)
    
    withProgress(message = "Performing analysis...", value = 0, {
      incProgress(0.25, detail = "Initializing...")
      
      res <- NULL
      
      if (input$BSmethod == "PAM50") {
        req(se_NC)
        cal  <- input$calibration %||% "Internal"
        args <- list(se_obj = se_NC, hasClinical = use_clin, Subtype = want_4)
        
        if (identical(cal, "Internal")) {
          args$calibration <- "Internal"
          args$internal    <- input$internal
          
        } else if (identical(cal, "External")) {
          args$calibration <- "External"
          
          if (identical(input$external, "Given.mdns")) {
            req(input$medians)
            inFile <- input$medians
            ext <- tolower(tools::file_ext(inFile$name))
            medians <- if (ext == "csv") {
              read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",",
                       check.names = FALSE)
            } else {
              read.table(inFile$datapath, header = TRUE, row.names = NULL, sep = "\t",
                         fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE,
                         quote = "", comment.char = "")
            }
            
            # --- VALIDATE custom medians ---
            if (ncol(medians) < 2) {
              showNotification("Custom medians file must have 2 columns: gene symbol + value.",
                               type = "error", duration = 7)
              return(invisible(NULL))
            }
            # Count unique gene symbols in first column
            uniq_genes <- length(unique(medians[[1]]))
            if (uniq_genes < 40) {
              showNotification(
                paste("Custom medians: only", uniq_genes,
                      "unique gene symbols found. Expect ~50 PAM50 genes."),
                type = "warning", duration = 7
              )
            }
            # Gently warn if column names aren’t X/Given.mdns
            if (!all(c("X", "Given.mdns") %in% names(medians))) {
              showNotification(
                "Tip: Columns will be treated as (X, Given.mdns). Rename in file to avoid ambiguity.",
                type = "message", duration = 6
              )
            }
            # -----------------------------------------
            
            args$external <- "Given.mdns"
            args$medians  <- medians
            
          } else {
            args$external <- input$external
          }
          
        } else {
          args$calibration <- "None"
        }
        
        incProgress(0.55, detail = "BS_parker...")
        res <- do.call(BreastSubtypeR::BS_parker, args)
      }
      
      
      if (input$BSmethod == "cIHC") {
        req(se_NC)
        incProgress(0.55, detail = "BS_cIHC...")
        res <- BreastSubtypeR::BS_cIHC(se_obj = se_NC, Subtype = want_4, hasClinical = use_clin)
      }
      
      if (input$BSmethod == "cIHC.itr") {
        req(se_NC)
        incProgress(0.55, detail = "BS_cIHC.itr...")
        res <- BreastSubtypeR::BS_cIHC.itr(
          se_obj = se_NC, iteration = input$iteration, ratio = input$ratio,
          Subtype = want_4, hasClinical = use_clin
        )
      }
      
      if (input$BSmethod == "PCAPAM50") {
        req(se_NC)
        incProgress(0.55, detail = "BS_PCAPAM50...")
        res <- tryCatch(
          BreastSubtypeR::BS_PCAPAM50(se_obj = se_NC, Subtype = want_4, hasClinical = use_clin),
          error = function(e) { showNotification("PCAPAM50 failed on this dataset; skipping.", type = "warning", duration = 6); NULL }
        )
      }
      
      if (input$BSmethod == "ssBC") {
        req(se_NC)
        cd <- as.data.frame(SummarizedExperiment::colData(se_NC))
        need_cols <- switch(input$s, "ER"=c("ER"), "ER.v2"=c("ER","HER2"), "TN"=c("TN"), "TN.v2"=c("TN"))
        miss <- setdiff(need_cols, names(cd))
        if (length(miss)) { showNotification(paste("ssBC requires:", paste(need_cols, collapse=", "),
                                                   "| missing:", paste(miss, collapse=", ")), type="error", duration=7)
          return(invisible(NULL)) }
        if (input$s %in% c("ER","ER.v2")) {
          badER <- setdiff(na.omit(unique(cd$ER)), c("ER+","ER-"))
          if (length(badER)) showNotification(sprintf("ER values should be 'ER+' or 'ER-'. Found: %s", paste(badER, collapse=", ")),
                                              type="warning", duration=7)
        }
        if (input$s == "ER.v2") {
          badH <- setdiff(na.omit(unique(cd$HER2)), c("HER2+","HER2-"))
          if (length(badH)) showNotification(sprintf("HER2 values should be 'HER2+' or 'HER2-'. Found: %s", paste(badH, collapse=", ")),
                                             type="warning", duration=7)
        }
        if (input$s %in% c("TN","TN.v2")) {
          badTN <- setdiff(na.omit(unique(cd$TN)), c("TN","nonTN"))
          if (length(badTN)) showNotification(sprintf("TN values should be 'TN' or 'nonTN'. Found: %s", paste(badTN, collapse=", ")),
                                              type="warning", duration=7)
        }
        incProgress(0.55, detail = "BS_ssBC...")
        res <- BreastSubtypeR::BS_ssBC(se_obj = se_NC, s = input$s, Subtype = want_4, hasClinical = use_clin)
      }
      
      if (input$BSmethod == "AIMS") {
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
        }
      }
      
      if (input$BSmethod == "sspbc") {
        req(se_SSP)
        incProgress(0.55, detail = "BS_sspbc...")
        model <- if (want_4) "ssp.subtype" else "ssp.pam50"
        res_sspbc <- BreastSubtypeR::BS_sspbc(se_obj = se_SSP, ssp.name = model)
        res <- list(BS.all = data.frame(PatientID = rownames(res_sspbc),
                                        BS = res_sspbc[, 1, drop = TRUE],
                                        row.names = rownames(res_sspbc),
                                        check.names = FALSE))
      }
      
      incProgress(0.9, detail = "Finalizing...")
      
      # ---------- Visualization (unchanged API expectations) ----------
      out <- NULL
      if (!is.null(res) && !is.null(res$BS.all) && nrow(res$BS.all) > 0) {
        labcol <- if (want_4 && "BS.Subtype" %in% names(res$BS.all)) "BS.Subtype"
        else if ("BS" %in% names(res$BS.all)) "BS"
        else setdiff(names(res$BS.all), "PatientID")[1]
        out <- data.frame(PatientID = res$BS.all$PatientID,
                          Subtype   = res$BS.all[[labcol]],
                          check.names = FALSE)
      } else if (!is.null(res) && !is.null(res$cl) && is.matrix(res$cl) && ncol(res$cl) >= 1) {
        out <- data.frame(PatientID = rownames(res$cl),
                          Subtype   = res$cl[, 1, drop = TRUE],
                          check.names = FALSE)
      }
      
      # ---------- Exports (Calls-only / Full) ----------
      if (!is.null(res) && !is.null(res$BS.all) && nrow(res$BS.all) > 0) {
        ex <- .build_single_exports(res, want_4)
        reactive_files$download_calls <- ex$calls
        reactive_files$download_full  <- ex$full
        
        # Hide "full" for AIMS/SSPBC
        if (input$BSmethod %in% c("AIMS", "sspbc")) {
          reactive_files$download_full <- NULL
        }
      } else {
        reactive_files$download_calls <- data.frame()
        reactive_files$download_full  <- NULL
      }
      
      output$download <- downloadHandler(
        filename = function() {
          kind <- if (!isTRUE(block_full()) && identical(input$export_kind, "full")) "full" else "calls"
          sprintf("results-%s-%s-%s.txt", input$BSmethod,
                  if (want_4) "4class" else "5class", kind)
        },
        content = function(file) {
          want_full <- identical(input$export_kind, "full") && !isTRUE(block_full())
          obj <- if (want_full && !is.null(reactive_files$download_full) && nrow(reactive_files$download_full) > 0) {
            reactive_files$download_full
          } else {
            reactive_files$download_calls
          }
          
          hdr <- c(
            CITATION_NOTE,
            .method_header_line(input, want_4, use_clin),
            sprintf("Export: %s | Date: %s | R: %s",
                    if (want_full) "full" else "calls",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                    getRversion())
          )
          
          .write_tsv_with_header(file, obj, header_lines = hdr)
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
