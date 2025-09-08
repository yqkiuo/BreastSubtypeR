# Define UI for iBreastSubtypeR
app_theme <- bslib::bs_theme(
  version = 5,
  base_font    = bslib::font_google("Inter"),
  heading_font = bslib::font_google("Poppins"),
  primary = "#FF69B4" # matches download button
)

ui <- bslib::page_fluid(
  theme = app_theme,
  
  # --- Global CSS fixes & minor polish ---
  tags$head(tags$style(HTML("
    /* Keep dropdowns visible above cards, but below modals */
    .bslib-card, .card, .card-body { overflow: visible !important; }
    .selectize-control, .selectize-dropdown, .dropdown-menu { z-index: 1040 !important; } /* < modal */
    .modal-backdrop { z-index: 1050 !important; }
    .modal { z-index: 1060 !important; }

    /* Center the top heading */
    .app-title { text-align: center; margin: 16px 0 8px; }
    .app-title h2 { margin: 0; }

    /* Styled info box for help panels */
    .method-help {
      margin: 6px 0 12px; padding: 10px 12px;
      border-left: 4px solid #e9ecef; background: #fafbfc; border-radius: 6px;
    }
    .method-help ul { margin-bottom: 0; }
  "))),
  
  # --- Centered heading ---
  tags$div(class = "app-title",
           tags$h2("Interactive Breast Cancer Intrinsic Molecular Subtyping")
  ),
  
  # --- Welcome / hero card ---
  bslib::card(
    tags$div(
      style = "display: flex; align-items: center;",
      tags$div(
        style = "flex: 1; text-align: left; margin-right: 20px;",
        bslib::card_image("logo.svg", height = "180px")
      ),
      tags$div(
        style = "flex: 3; text-align: center;",
        tags$div(style = "margin-bottom: 10px;",
                 bslib::card_header("Welcome to iBreastSubtypeR!")
        ),
        "iBreastSubtypeR provides an interactive interface for intrinsic molecular subtyping of breast cancer using nearest-centroid (NC) and single-sample predictor (SSP) methods. It standardizes I/O handling and presents results compatible with R/Bioconductor workflows.",
        tags$div(style = "margin-top: 10px;",
                 bslib::card_footer("Enjoy your subtyping journey!")
        )
      )
    )
  ),
  
  #### Step 1
  h3("Step 1 · Upload your data"),
  
  bslib::layout_column_wrap(
    col_width = 3,
    
    # 1) Gene expression
    bslib::card(
      bslib::card_header("1) Gene expression (GEX)"),
      fileInput(
        "GEX",
        "Upload expression matrix",
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")
      ),
      radioButtons(
        "is_raw_counts", "Data type",
        choices  = c("Normalized (log2)" = "norm", "Raw counts (RNA-seq)" = "raw"),
        selected = "norm", inline = TRUE
      ),
      uiOutput("gex_help")     # "Requirements" panel (mode-aware)
    ),
    
    # 2) Clinical data
    bslib::card(
      bslib::card_header("2) Clinical data"),
      fileInput(
        "clinic",
        "Upload clinical table",
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")
      ),
      uiOutput("clin_help")    # "Clinical data requirements"
    ),
    
    # 3) Feature annotation
    bslib::card(
      bslib::card_header("3) Feature annotation"),
      fileInput(
        "anno",
        "Upload annotation table",
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")
      ),
      uiOutput("anno_help")    # "Feature annotation requirements"
    )
  ),
  
  # Preprocess / Map
  bslib::card(
    actionButton("map", "Preprocess & map", icon = icon("sliders")),
    helpText("Runs Mapping() to align IDs and normalize (as needed). When complete, continue to Step 2.")
  ),
  
  #### Step 2
  h3("Step 2 · Choose method & parameters"),
  
  bslib::card(
    style = "overflow: visible;",
    
    # Subtyping method + AUTO
    selectizeInput(
      "BSmethod", "Subtyping method",
      choices = c(
        "Run all (AUTO Mode)" = "AUTO Mode",
        "PAM50 (parker.original | genefu.scale | genefu.robust)" = "PAM50",
        "cIHC"     = "cIHC",
        "cIHC.itr" = "cIHC.itr",
        "PCAPAM50" = "PCAPAM50",
        "ssBC"     = "ssBC",
        "AIMS"     = "AIMS",
        "sspbc"    = "sspbc"
      ),
      selected = "AUTO Mode",
      width = "100%",
      options = list(openOnFocus = TRUE, dropdownParent = "body")
    ),
    
    # --- NEW: persistent AUTO preflight panel ---
    uiOutput("auto_preflight"),  # appears only when AUTO Mode is selected
    
    # Global 4 vs 5 classes (AIMS always 5)
    radioButtons(
      "k_subtypes", "Subtype classes",
      choices = c("5 classes (includes Normal-like)" = "5",
                  "4 classes (excludes Normal-like)" = "4"),
      selected = "5", inline = TRUE
    ),
    
    # Live per-method help
    uiOutput("method_help"),
    
    # ROR checkbox (NC methods + AUTO)
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50' || input.BSmethod == 'cIHC' || input.BSmethod == 'cIHC.itr' || input.BSmethod == 'PCAPAM50' || input.BSmethod == 'ssBC'",
      checkboxInput("hasClinical", "Use clinical variables (ROR)", value = FALSE)
    ),
    
    # PAM50 calibration row with side-by-side controls
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50'",
      bslib::card(
        style = "overflow: visible;",
        bslib::layout_columns(
          col_widths = c(4, 8),
          
          selectizeInput(
            "calibration", "Calibration strategy",
            choices = c("None" = "None", "External" = "External", "Internal" = "Internal"),
            selected = "Internal", width = "100%",
            options = list(openOnFocus = TRUE, dropdownParent = "body")
          ),
          
          div(
            style = "position: relative; overflow: visible;",
            
            # External
            conditionalPanel(
              condition = "input.calibration == 'External'",
              tagList(
                selectizeInput(
                  "external", "External calibration method",
                  choices = c(
                    "Given.mdns" = "Given.mdns",
                    "nCounter" = "nCounter",
                    "RNAseq.Freeze.20120907" = "RNAseq.Freeze.20120907",
                    "totalRNA.FFPE.20151111" = "totalRNA.FFPE.20151111",
                    "RNAseq.V2"  = "RNAseq.V2",
                    "RNAseq.V1"  = "RNAseq.V1",
                    "GC.4x44Kcustom" = "GC.4x44Kcustom",
                    "Agilent_244K"   = "Agilent_244K",
                    "commercial_1x44k_postMeanCollapse_WashU"    = "commercial_1x44k_postMeanCollapse_WashU",
                    "commercial_4x44k_postMeanCollapse_WashU_v2" = "commercial_4x44k_postMeanCollapse_WashU_v2",
                    "htp1.5_WU_update" = "htp1.5_WU_update",
                    "arrayTrain_postMeanCollapse" = "arrayTrain_postMeanCollapse"
                  ),
                  selected = "RNAseq.V2", width = "100%",
                  options = list(placeholder = "Choose a reference set…", openOnFocus = TRUE, dropdownParent = "body")
                ),
                conditionalPanel(
                  condition = "input.external == 'Given.mdns'",
                  fileInput(
                    "medians", "Upload Given.mdns file",
                    accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")
                  )
                )
              )
            ),
            
            # Internal
            conditionalPanel(
              condition = "input.calibration == 'Internal'",
              selectizeInput(
                "internal", "Internal calibration method",
                choices = c(
                  "medianCtr (parker.original)"   = "medianCtr",
                  "meanCtr (genefu.scale)"        = "meanCtr",
                  "qCtr (genefu.robust)"   = "qCtr"
                ),
                selected = "medianCtr", width = "100%",
                options = list(placeholder = "Pick a calibration…", openOnFocus = TRUE, dropdownParent = "body")
              )
            )
          )
        )
      )
    ),

    uiOutput("calib_help"),
    
    
    # cIHC / PCAPAM50 / AIMS (no extra UI blocks)
    conditionalPanel(condition = "input.BSmethod == 'cIHC'", div()),
    conditionalPanel(
      condition = "input.BSmethod == 'cIHC.itr'",
      bslib::layout_columns(
        col_widths = c(6, 6),
        numericInput("iteration", label = "Iterations", value = 100, min = 10, step = 10),
        numericInput("ratio", label = "ER+ training ratio", value = 54/64, min = 0, max = 1, step = 0.01)
      )
    ),
    conditionalPanel(condition = "input.BSmethod == 'PCAPAM50'", div()),
    conditionalPanel(
      condition = "input.BSmethod == 'ssBC'",
      bslib::layout_column_wrap(
        selectInput(
          "s", "Subgroup",
          choices = list("ER" = "ER", "ER.v2" = "ER.v2", "TN" = "TN", "TN.v2" = "TN.v2"),
          selected = "ER.v2"
        )
      )
    ),
    conditionalPanel(condition = "input.BSmethod == 'AIMS'", div()),
    conditionalPanel(condition = "input.BSmethod == 'sspbc'", div()),
    
    bslib::card(
      actionButton("run", "Run subtyping", icon = icon("play-circle"))
    )
  ),
  
  # Visualization (rendered after run)
  uiOutput("plotSection"),
  
  # Footer
  fluidRow(
    column(
      12, align = "center",
      div(
        style = "margin-top: 20px; margin-bottom: 20px;",
        downloadButton(
          "download",
          HTML("Download results<br>(.tsv)"),
          style = "
            width: 220px;
            height: auto;
            padding: 10px 16px;
            white-space: normal;
            line-height: 1.2;
            background-color: #FF69B4;
            color: white;
            border: none;
          "
        )
      )
    )
  )
)
