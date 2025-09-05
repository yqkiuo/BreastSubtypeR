# Define UI for iBreastSubtypeR
ui <- bslib::page_fluid(
  titlePanel("Interactive Breast Cancer Intrinsic Molecular Subtyping"),
  
  bslib::card(
    tags$div(
      style = "display: flex; align-items: center;",
      tags$div(
        style = "flex: 1; text-align: left; margin-right: 20px;",
        bslib::card_image("logo.svg", height = "180px")
      ),
      tags$div(
        style = "flex: 3; text-align: center;",
        tags$div(style = "margin-bottom: 10px;", bslib::card_header("Welcome to iBreastSubtypeR!")),
        "iBreastSubtypeR provides an interactive interface for intrinsic molecular subtyping of breast cancer using nearest-centroid (NC) and single-sample predictor (SSP) methods. It standardizes I/O handling and presents results compatible with R/Bioconductor workflows.",
        tags$div(style = "margin-top: 10px;", bslib::card_footer("Enjoy your subtyping journey!"))
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
      uiOutput("gex_help")
    ),
    
    # 2) Clinical
    bslib::card(
      bslib::card_header("2) Clinical metadata"),
      fileInput(
        "clinic",
        "Upload clinical table",
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")
      ),
      uiOutput("clin_help")
    ),
    
    # 3) Feature annotation
    bslib::card(
      bslib::card_header("3) Feature annotation"),
      fileInput(
        "anno",
        "Upload annotation table",
        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")
      ),
      uiOutput("anno_help")
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
    selectInput(
      "BSmethod", "Subtyping method",
      choices = list(
        "PAM50.parker (NC)" = "PAM50.parker",
        "cIHC (NC)"         = "cIHC",
        "cIHC.itr (NC)"     = "cIHC.itr",
        "PCAPAM50 (NC)"     = "PCAPAM50",
        "ssBC (NC)"         = "ssBC",
        "AIMS (SSP)"        = "AIMS",
        "sspbc (SSP)"       = "sspbc"
      ),
      selected = "PAM50.parker"
    ),
    
    # One checkbox for ROR-enabled methods
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50.parker' || input.BSmethod == 'cIHC' || input.BSmethod == 'cIHC.itr' || input.BSmethod == 'PCAPAM50' || input.BSmethod == 'ssBC'",
      checkboxInput("hasClinical", "Use clinical variables (ROR)", value = FALSE)
    ),
    
    # PAM50.parker options
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50.parker'",
      bslib::layout_column_wrap(
        selectInput(
          "calibration", "Calibration strategy",
          choices = list("None" = "None", "External" = "External", "Internal" = "Internal"),
          selected = "Internal"
        )
      ),
      bslib::card(
        conditionalPanel(
          condition = "input.calibration == 'External'",
          bslib::layout_column_wrap(
            selectInput(
              "external", "External reference set",
              choices = list(
                "Given.mdns" = "Given.mdns",
                "RNAseq.V2" = "RNAseq.V2",
                "RNAseq.V1" = "RNAseq.V1",
                "GC.4x44Kcustom" = "GC.4x44Kcustom",
                "Agilent_244K" = "Agilent_244K",
                "commercial_1x44k_postMeanCollapse_WashU" = "commercial_1x44k_postMeanCollapse_WashU",
                "commercial_4x44k_postMeanCollapse_WashU_v2" = "commercial_4x44k_postMeanCollapse_WashU_v2",
                "htp1.5_WU_update" = "htp1.5_WU_update",
                "arrayTrain_postMeanCollapse" = "arrayTrain_postMeanCollapse"
              ),
              selected = "RNAseq.V2"
            )
          )
        ),
        conditionalPanel(
          condition = "input.calibration == 'External' && input.external == 'Given.mdns'",
          fileInput(
            "medians", "Upload Given.mdns file",
            accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")
          )
        ),
        conditionalPanel(
          condition = "input.calibration == 'Internal'",
          selectInput(
            "internal", "Internal calibration method",
            choices = list(
              "Parker (median-centering)"            = "medianCtr",
              "genefu (mean scaling)"                = "meanCtr",
              "genefu (quantile/robust centering)"   = "qCtr"
            ),
            selected = "medianCtr"
          )
        )
      )
    ),
    
    # cIHC (no extra controls)
    conditionalPanel(
      condition = "input.BSmethod == 'cIHC'",
      div()
    ),
    
    # cIHC.itr controls
    conditionalPanel(
      condition = "input.BSmethod == 'cIHC.itr'",
      bslib::layout_column_wrap(
        numericInput("iteration", label = "Iterations", value = 100, min = 10, step = 10),
        numericInput("ratio", label = "ER+ training ratio", value = 54/64, min = 0, max = 1, step = 0.01)
      )
    ),
    
    # PCAPAM50 (no extra controls)
    conditionalPanel(
      condition = "input.BSmethod == 'PCAPAM50'",
      div()
    ),
    
    # ssBC controls
    conditionalPanel(
      condition = "input.BSmethod == 'ssBC'",
      bslib::layout_column_wrap(
        selectInput(
          "s", "Subgroup",
          choices = list("ER" = "ER", "ER.v2" = "ER.v2", "TN" = "TN"),
          selected = "ER.v2"
        )
      )
    ),
    
    # AIMS (no extra controls)
    conditionalPanel(
      condition = "input.BSmethod == 'AIMS'",
      div()
    ),
    
    # sspbc controls
    conditionalPanel(
      condition = "input.BSmethod == 'sspbc'",
      bslib::layout_column_wrap(
        selectInput(
          "Subtype",
          "Subtype method (TRUE = 4 subtypes; FALSE = 5 subtypes)",
          choices = list("5 subtypes" = "FALSE", "4 subtypes" = "TRUE"),
          selected = "FALSE"
        )
      )
    ),
    
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
