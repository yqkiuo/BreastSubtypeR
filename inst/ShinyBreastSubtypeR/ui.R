# Define UI for iBreastSubtypeR

ui = page_fluid(

  titlePanel("interactive Breast Cancer Intrinsic Subtype (iBreastSubtypeR)"),

  card(
    card_header("Welcom to iBreastSubtypeR"),
    "BreastSubtypeR integrates intrinsic molecular subtyping methods for breast cancer, 
    including nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches. 
    It employs standardized input and output formats, 
    offering a unified framework that is highly compatible with other R packages in the gene expression profiling field.",
    card_image("logo.svg", height = "150px"),
    card_footer("Enjoy this subtyping journey.")
  ),
  
  #### input your data ####
  # Main page section with a title
  h3("Please input your data"),  # Add a title here

  ## main page, input gene expression etc.
  layout_column_wrap(

    col_width = 3,
    ## gene expression
    card(
      card_header("Gene expression table"),
      fileInput("GEX", "choose file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".txt")
                )
      #checkboxInput("Log", "Log2", value = FALSE)
    ),
    
    ## clinic information
    card(
      card_header("Clinical table"),
      fileInput("clinic", "choose file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".txt")
                ), ## csv/text
    ),
    
    ## feature information
    card(
      card_header("Feature table"),
      fileInput("anno", "choose file",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  ".txt")
      ), ## csv/text
    )
  ),
  
  #### Map button ####
  card(
    actionButton("map", "Step 1: Map" )
  )
  , ## layout 3

  #### Please choose your method ####
  # Main page section with a title
  h3("Please choose your method"),  # Add a title here
  
  ## options
  card(
    card_header( "Subtyping method setting"),

    selectInput(
      "BSmethod",
      "Select subtyping method",
      choices = list("PAM50.parker" = "PAM50.parker",
                     "cIHC" = "cIHC",
                     "cIHC.itr" = "cIHC.itr",
                     "PCAPAM50" = "PCAPAM50",
                     "ssBC" ="ssBC" ,
                     "AIMS" = "AIMS",
                     "sspbc" ="sspbc"),
      selected = "PAM50.parker"
    ),

    # Conditional panels for each choice
    conditionalPanel(
      condition = "input.BSmethod == 'PAM50.parker' ", ## parker

      layout_column_wrap(
        ##calibration
        selectInput(
          "hasClinical",
          "Has Clinical", choices =  list( "FALSE" = "FALSE", "TRUE" = "TRUE" ),
          selected = "FALSE"),
        selectInput(
          "calibration",
          "calibration",
          choices = list("None" = "None","External" = "External", "Internal" = "Internal"),
          selected = "Internal")
      ),
      
      
      card(
        conditionalPanel(
          condition = "input.calibration == 'External'",
          layout_column_wrap(
            selectInput(
              "external",
              "External",
              choices = list("Given.mdns" = "Given.mdns",
                             "RNAseq.V2" = "RNAseq.V2","RNAseq.V1"  = "RNAseq.V1" ,
                             "GC.4x44Kcustom" = "GC.4x44Kcustom", "Agilent_244K" = "Agilent_244K",
                             "commercial_1x44k_postMeanCollapse_WashU"  = "commercial_1x44k_postMeanCollapse_WashU" , "commercial_4x44k_postMeanCollapse_WashU_v2" = "commercial_4x44k_postMeanCollapse_WashU_v2",
                             "htp1.5_WU_update"  =  "htp1.5_WU_update", "arrayTrain_postMeanCollapse" = "arrayTrain_postMeanCollapse"
                             ),
              selected = "RNAseq.V2")
          )
          ),
        
        # Show a new file input field when calibration is "External" and external is "Given.mdns"
        conditionalPanel(
          condition = "input.calibration == 'External' && input.external == 'Given.mdns'",
          fileInput("medians", "Upload Given.mdns File", 
                    accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".txt"))
        ),

        conditionalPanel(
          condition = "input.calibration == 'Internal'",
        selectInput(
          "internal",
          "Internal",
          choices = list("median" = "medianCtr", "mean" = "meanCtr", "quantile" = "qCtr"),
          selected = "median")
      )
      
      )
      
      ),

    conditionalPanel(
      condition = "input.BSmethod == 'cIHC'",
      layout_column_wrap(
        selectInput(
          "hasClinical",
          "Has Clinical",
          choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
          selected = "FALSE"
        )
      )
    ),

    conditionalPanel(
      condition = "input.BSmethod == 'cIHC.itr'",
      layout_column_wrap(
        numericInput("iteration", label = "Iteration", value = 100),
        numericInput("ratio", label = "Ratio", value = 54/64),
        selectInput(
          "hasClinical",
          "Has Clinical",
          choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
          selected = "FALSE"
        )
      )
    ),

    conditionalPanel(
      condition = "input.BSmethod == 'PCAPAM50'",
      layout_column_wrap(
        selectInput(
          "hasClinical",
          "Has Clinical",
          choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
          selected = "FALSE"
        )
      )
    ),

    conditionalPanel(
      condition = "input.BSmethod == 'ssBC'",
      layout_column_wrap(
        selectInput(
          "s",
          "subgroup",
          choices = list("ER" = "ER", "TN" = "TN", "ER_JAMA" = "ER_JAMA", "HER2+" = "HER2+", "TNBC" = "TNBC"),
          selected = "ER_JAMA"
        ),
        selectInput(
          "hasClinical",
          "Has Clinical",
          choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
          selected = "FALSE"
        )
      )
    ),

    conditionalPanel(
      condition = "input.BSmethod == 'AIMS'",
    ),

    conditionalPanel(
      condition = "input.BSmethod == 'sspbc'",
      layout_column_wrap(
        selectInput(
          "Subtype",
          "Subtype",
          choices = list("FALSE" = "FALSE", "TRUE" = "TRUE"),
          selected = "FALSE"
        )
      )
    ),

    card(
    actionButton("run", "Step 2: Run" )
    )
  ), ## card option
  
  # Visualization layout placed separately and conditionally rendered
  uiOutput("plotSection"),
  
  # Footer with download and reset buttons
  fluidRow(
    column(12, align = "center",
           div(style = "margin-top: 20px; margin-bottom: 20px;",
               downloadButton("download", "", style = "width: 150px; height: 40px;"))
    )
  )
  
)

