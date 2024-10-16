library(shiny)
library(shinyjs) ## js action !!!
library(bslib)

# Define UI for app that draws a histogram ----

ui = page_fluid(


  titlePanel("Breast Cancer Intrinsic Subtype"),


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
  ), ## layout 3

  
  ## Map button
  card(
    actionButton("map", "Map" )
  ),

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
          "hasClinical", choices =  list( "FALSE" = "FALSE", "TRUE" = "TRUE" ),
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
              choices = list("Given.mdns" = "Given.mdns"),
              selected = "Given.mdns")
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
    ),

    card(
    actionButton("run", "Run" )
    )
  ), ## card option
  
  
  # Footer with download and reset buttons
  fluidRow(
    column(12, align = "center",
           div(style = "margin-top: 20px; margin-bottom: 20px;",
               downloadButton("download", "", style = "width: 150px; height: 40px;"),
               actionButton("reset", "Reset", style = "width: 150px; height: 40px; margin-left: 10px;")
           )
    )
  ),

  # Visualization layout placed separately and conditionally rendered
  uiOutput("plotSection")
  
  
)

