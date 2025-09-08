# Define UI for iBreastSubtypeR
app_theme <- bslib::bs_theme(
  version = 5,
  base_font    = bslib::font_google("Inter"),
  heading_font = bslib::font_google("Poppins"),
  primary = "#FF69B4"
)

ui <- bslib::page_fluid(
  theme = app_theme,
  
  # --- Global CSS fixes ---
  tags$head(tags$style(HTML("
    .bslib-card, .card, .card-body { overflow: visible !important; }
    .selectize-control, .selectize-dropdown, .dropdown-menu { z-index: 1040 !important; }
    .modal-backdrop { z-index: 1050 !important; }
    .modal { z-index: 1060 !important; }
    .app-title { text-align: center; margin: 16px 0 8px; }
    .app-title h2 { margin: 0; }
    .method-help { margin: 6px 0 12px; padding: 10px 12px;
      border-left: 4px solid #e9ecef; background: #fafbfc; border-radius: 6px; }
    .method-help ul { margin-bottom: 0; }
  "))),
  
  tags$head(tags$style(HTML("
  /* HERO */
  .hero-card {
    background: linear-gradient(135deg, #ffe3f0 0%, #ffffff 60%);
    border: 1px solid #f1f3f5; border-radius: 12px;
  }
  .hero-inner { display: grid; grid-template-columns: 180px 1fr; gap: 16px; align-items: center; }
  @media (max-width: 768px) { .hero-inner { grid-template-columns: 1fr; text-align:center; } }

  .hero-title { margin: 4px 0 6px; font-weight: 600; }
  .hero-subtitle { margin: 0 0 10px; color: #4a4f55; }

  .chips { display: flex; flex-wrap: wrap; gap: 6px; margin: 8px 0 6px; }
  .chip {
    font-size: 12px; padding: 4px 8px; background: #f8f9fa;
    border: 1px solid #e9ecef; border-radius: 999px;
    display: inline-flex; align-items: center; gap: 6px;
  }
  .chip i { opacity: 0.8; }

  .feature-list { margin: 8px 0 0 0; padding-left: 18px; }
  .feature-list li { margin: 3px 0; }

  .hero-ctas { display:flex; gap:10px; margin-top:12px; flex-wrap:wrap; }
  .btn-pink {
    background-color:#FF69B4; border:none; color:white;
  }
  .btn-outline { border:1px solid #dee2e6; background:#fff; }
  
  .chip.auto { border:1px solid; padding:6px 10px; border-radius:999px; display:inline-flex; gap:8px; align-items:center; }
.chip.auto.ready   { background:#e6f7ed; border-color:#2fb170; }
.chip.auto.blocked { background:#fff5f5; border-color:#e03131; }
.chip .chip-note   { font-size:11px; opacity:0.75; }

"))),
  
  tags$head(
    tags$script(HTML("
    (function () {
      var lastY = 0;
      var restoreOnIdle = false;

      function saveScroll() {
        lastY = window.scrollY || window.pageYOffset || 0;
      }
      function restoreScroll() {
        window.scrollTo(0, lastY || 0);
      }
      function scrollToMap() {
        var el = document.getElementById('map');
        if (el && el.scrollIntoView) {
          el.scrollIntoView({ behavior: 'smooth', block: 'center' });
        }
      }

      // ---- File inputs: prevent jump on open/close; after anno selection, jump to Map
      ['GEX','clinic','anno'].forEach(function(id) {
        $(document).on('mousedown click focusin', 'label[for=\"' + id + '\"], #' + id, function() {
          saveScroll();
          setTimeout(restoreScroll, 0);
        });
        $(document).on('change', '#' + id, function() {
          setTimeout(function() {
            restoreScroll();
            if (id === 'anno') scrollToMap();
          }, 50);
        });
        $(document).on('focusout', '#' + id, function() {
          setTimeout(restoreScroll, 0);
        });
      });

      // ---- Map button: keep position through withProgress
      $(document).on('mousedown click', '#map', function() {
        saveScroll();
        restoreOnIdle = true;
      });

      // ---- When Shiny finishes the server work, restore where we were
      $(document).on('shiny:idle', function() {
        if (restoreOnIdle) {
          setTimeout(function() {
            restoreScroll();
            scrollToMap();
            restoreOnIdle = false;
          }, 50);
        }
      });
    })();
  "))
  ),
  
  
  # --- Heading ---
  tags$div(class = "app-title", tags$h2("Interactive Breast Cancer Intrinsic Molecular Subtyping")),
  
  # --- Welcome / hero card ---
  bslib::card(
    class = "hero-card",
    bslib::card_body(
      div(class = "hero-inner",
          # left: logo
          div(bslib::card_image("logo.svg", height = "160px")),
          # right: text
          div(
            h3(class = "hero-title", "Welcome to iBreastSubtypeR"),
            div(class = "hero-subtitle",
                "An interactive companion to the BreastSubtypeR Bioconductor package for intrinsic molecular subtyping."
            ),
            # chips row (icons & labels)
            div(class = "chips",
                span(
                  class = "chip", icon("bullseye"),
                  HTML("<b>NC-based:</b> PAM50 (parker.original | genefu.scale | genefu.robust), cIHC / cIHC.itr, PCAPAM50, ssBC / ssBC.v2")
                ),
                span(
                  class = "chip", icon("vial"),
                  HTML("<b>SSP-based:</b> AIMS, SSPBC")
                ),
                uiOutput("auto_chip"),  # smart AUTO chip (dynamic)
                span(
                  class = "chip", icon("chart-line"),
                  HTML("<b>ROR:</b> research-use Risk of Recurrence (uses TSIZE, NODE; NC methods only)")
                )
            ),
            # short bullets
            tags$ul(class = "feature-list",
                    tags$li("Assumption-aware AUTO reduces bias in ER/HER2-skewed or subtype-specific cohorts."),
                    tags$li(HTML("Standardised mapping & normalisation across platforms (<i>log<sub>2</sub></i>-CPM for NC; FPKM for SSP).")),
                    tags$li("Choose 5-class (incl. Normal-like) or 4-class; AIMS is 5-class only."),
                    tags$li("All computation runs locally; exports are Bioconductor-ready.")
            ),
            # CTAs
            div(class = "hero-ctas",
                tags$a(class = "btn btn-pink", href = "#step1", icon("upload"), "Go to Step 1"),
                tags$a(class = "btn btn-outline", href = "#step2", icon("play-circle"), "Go to Step 2")
            )
          )
      )
    )
  ),
  
  #### Step 1
  h3(id = "step1", "Step 1 · Upload your data"),
  
  bslib::layout_column_wrap(
    col_width = 3,
    
    # 1) Gene expression
    bslib::card(
      bslib::card_header("1) Gene expression (GEX)"),
      fileInput("GEX", "Upload expression matrix",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")),
      radioButtons("is_raw_counts", "Data type",
                   choices  = list("Normalized (log2)" = "norm", "Raw counts (RNA-seq)" = "raw"),
                   selected = "norm", inline = TRUE),
      uiOutput("gex_help")
    ),
    
    # 2) Clinical data
    bslib::card(
      bslib::card_header("2) Clinical data"),
      fileInput("clinic", "Upload clinical table",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")),
      uiOutput("clin_help")
    ),
    
    # 3) Feature annotation
    bslib::card(
      bslib::card_header("3) Feature annotation"),
      fileInput("anno", "Upload annotation table",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt")),
      uiOutput("anno_help")
    )
  ),
  
  bslib::card(
    actionButton("map", "Preprocess & map", icon = icon("sliders")),
    helpText("Runs Mapping() to align IDs and normalize (as needed). When complete, continue to Step 2.")
  ),
  
  #### Step 2
  h3(id = "step2", "Step 2 · Choose method & parameters"),
  
  bslib::card(
    style = "overflow: visible;",
    
    # Subtyping method + AUTO
    selectizeInput(
      "BSmethod", "Subtyping method",
      choices = list(
        "Run all (AUTO Mode)" = "AUTO Mode",
        "PAM50 (parker.original | genefu.scale | genefu.robust)" = "PAM50",
        "cIHC" = "cIHC",
        "cIHC.itr" = "cIHC.itr",
        "PCAPAM50" = "PCAPAM50",
        "ssBC" = "ssBC",
        "AIMS" = "AIMS",
        "sspbc" = "sspbc"
      ),
      selected = "AUTO Mode", width = "100%",
      options = list(openOnFocus = TRUE, dropdownParent = "body")
    ),
    
    uiOutput("auto_preflight"),
    
    # Global 4 vs 5 classes (AIMS always 5)
    radioButtons(
      "k_subtypes", "Subtype classes",
      choices = list("5 classes (includes Normal-like)" = "5",
                     "4 classes (excludes Normal-like)" = "4"),
      selected = "5", inline = TRUE
    ),
    
    # Live per-method help
    uiOutput("method_help"),
    
    # ROR checkbox (NC methods + AUTO -> disabled for AUTO)
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50' || input.BSmethod == 'cIHC' || input.BSmethod == 'cIHC.itr' || input.BSmethod == 'PCAPAM50' || input.BSmethod == 'ssBC'",
      checkboxInput("hasClinical", "Use clinical variables (ROR)", value = FALSE)
    ),
    
    # PAM50 calibration row
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50'",
      bslib::card(
        style = "overflow: visible;",
        bslib::layout_columns(
          col_widths = c(4, 8),
          
          selectizeInput(
            "calibration", "Calibration strategy",
            choices = list("None" = "None", "External" = "External", "Internal" = "Internal"),
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
                  choices = list(
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
                  options = list(placeholder = "Choose a reference set…",
                                 openOnFocus = TRUE, dropdownParent = "body")
                ),
                conditionalPanel(
                  condition = "input.external == 'Given.mdns'",
                  fileInput("medians", "Upload Given.mdns file",
                            accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt"))
                )
              )
            ),
            
            # Internal
            conditionalPanel(
              condition = "input.calibration == 'Internal'",
              selectizeInput(
                "internal", "Internal calibration method",
                choices = list(
                  "medianCtr (parker.original)" = "medianCtr",
                  "meanCtr (genefu.scale)"      = "meanCtr",
                  "qCtr (genefu.robust)"        = "qCtr"
                ),
                selected = "medianCtr", width = "100%",
                options = list(openOnFocus = TRUE, dropdownParent = "body")
              )
            )
          )
        )
      )
    ),
    
    uiOutput("calib_help"),
    
    # cIHC.itr params
    conditionalPanel(
      condition = "input.BSmethod == 'cIHC.itr'",
      bslib::layout_columns(
        col_widths = c(6, 6),
        numericInput("iteration", label = "Iterations", value = 100, min = 10, step = 10),
        numericInput("ratio", label = "ER+ training ratio", value = 54/64, min = 0, max = 1, step = 0.01)
      )
    ),
    
    # ssBC subgroup
    conditionalPanel(
      condition = "input.BSmethod == 'ssBC'",
      bslib::layout_column_wrap(
        selectInput("s", "Subgroup",
                    choices = list("ER" = "ER", "ER.v2" = "ER.v2", "TN" = "TN", "TN.v2" = "TN.v2"),
                    selected = "ER.v2")
      )
    ),
    
    bslib::card(
      actionButton("run", "Run subtyping", icon = icon("play-circle"))
    )
  ),
  
  # Visualization (rendered after run)
  uiOutput("plotSection"),
  
  # --- Export controls + Download ---
  bslib::card(
    bslib::layout_columns(
      col_widths = c(8,4),
      uiOutput("export_selector"),
      div(style = "text-align:center;",
          downloadButton(
            "download",
            HTML("Download results<br>(.tsv)"),
            style = "
            width: 220px; height: auto; padding: 10px 16px;
            white-space: normal; line-height: 1.2;
            background-color: #FF69B4; color: white; border: none;"
          )
      )
    )
  )
  
)
