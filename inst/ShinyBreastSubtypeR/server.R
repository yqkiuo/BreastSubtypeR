# Define server logic required to draw a histogram ----

library(BreastSubtypeR)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)

library(dplyr)
library(magrittr)
library(stringr)


server <- function(input, output) {
  
  output$logo = renderImage( 
    { 
      list(src = "logo.svg", height = "100%",inline = FALSE) 
    }, 
    deleteFile = FALSE 
  ) 
  
  
  #### upload data ####

  # Reactive values to store uploaded data and results
  reactive_files = reactiveValues(GEX = NULL, clinic = NULL, anno = NULL, data_input = NULL, output_res = NULL )

  observeEvent( input$map, {
    req(input$GEX)
    inFile = input$GEX
    fileExt = tools::file_ext(inFile$name)
    
    # Read the file based on its extension
    reactive_files$GEX = if (fileExt == "csv") {
      read.csv(inFile$datapath, header = TRUE, row.names = 1, sep = ",")
    } else if (fileExt == "txt") {
      read.table(inFile$datapath, header = TRUE,row.names = 1, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file type")
    }
    
    
    # Read the file based on its extension
    req(input$clinic)
    inFile = input$clinic
    fileExt = tools::file_ext(inFile$name)
    
    reactive_files$clinic = if (fileExt == "csv") {
      read.csv(inFile$datapath, header = TRUE, row.names = NULL,  sep = ",")
    } else if (fileExt == "txt") {
      read.table(inFile$datapath, header = TRUE,row.names = NULL, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file type")
    }
    rownames(reactive_files$clinic) = reactive_files$clinic$patientID
    
    # Read the file based on its extension
    req(input$anno)
    inFile = input$anno
    fileExt = tools::file_ext(inFile$name)
    
    reactive_files$anno = if (fileExt == "csv") {
      read.csv(inFile$datapath, header = TRUE, row.names = NULL,   sep = ",")
    } else if (fileExt == "txt") {
      read.table(inFile$datapath, header = TRUE, row.names = NULL,  sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported file type")
    }
    
    ## Run Mapping function
    data_input = Mapping(gene_expression_matrix = reactive_files$GEX, featuredata = reactive_files$anno, impute = TRUE, verbose = TRUE )
    reactive_files$data_input = data_input
  })
  
  #### perform analysis ####
   observeEvent( input$run, {

 
     if ( input$BSmethod == "PAM50.parker" ){
       
       if(input$calibration == "External" ){
         
         req(input$medians)
         inFile = input$medians
         fileExt = tools::file_ext(inFile$name)
         
         # Read the file based on its extension
         medians = if (fileExt == "csv") {
           read.csv(inFile$datapath, header = TRUE, row.names = NULL, sep = ",")
         } else if (fileExt == "txt") {
           read.table(inFile$datapath, header = TRUE,row.names = NULL,sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
         } else {
           stop("Unsupported file type")
         }
         
       }
       
       if ( input$calibration == "Internal" ) {
         res = BreastSubtypeR::BS_parker(gene_expression_matrix = reactive_files$data_input$x_NC.log, phenodata= reactive_files$clinic, calibration = input$calibration, 
                                         internal = input$internal,  
                                         hasClinical = input$hasClinical
         )
       }
       if ( input$calibration == "External" ) {
         res = BreastSubtypeR::BS_parker(gene_expression_matrix = reactive_files$data_input$x_NC.log, phenodata= reactive_files$clinic, calibration = input$calibration, 
                                         external = input$external,  
                                         hasClinical = input$hasClinical, medians = medians
         )
       }
       
       output_res = res$score.ROR
     }
     
     if (input$BSmethod == "cIHC"){

       res = BreastSubtypeR::BS_cIHC(gene_expression_matrix = reactive_files$data_input$x_NC.log, phenodata= reactive_files$clinic,  hasClinical = input$hasClinical )
       
       output_res = res$score.ROR
     }
     
     if (input$BSmethod == "cIHC.itr"){

       res= BreastSubtypeR::BS_cIHC.itr(gene_expression_matrix = reactive_files$data_input$x_NC.log, phenodata= reactive_files$clinic, iteration = input$iteration, ratio = input$ratio, hasClinical = input$hasClinical)
       
       output_res = res$score.ROR
     }
       
     if( input$BSmethod == "PCAPAM50"){

       
       res= BreastSubtypeR::BS_PCAPAM50(gene_expression_matrix = reactive_files$data_input$x_NC.log, phenodata= reactive_files$clinic,  hasClinical = input$hasClinical )
       
       output_res = res$score.ROR
     }
       
     if (input$BSmethod == "ssBC"){

       res = BreastSubtypeR::BS_ssBC(gene_expression_matrix = reactive_files$data_input$x_NC.log, phenodata= reactive_files$clinic, s = input$s, hasClinical = input$hasClinical)
       
       output_res = res$score.ROR
     }
       
     if ( input$BSmethod == "AIMS") {
       
       data("genes.signature")
       genes = as.character( genes.signature$EntrezGene.ID[which( genes.signature$AIMS == "Yes" )])

       res = BreastSubtypeR::BS_AIMS( reactive_files$data_input$x_SSP[ colnames(reactive_files$data_input$x_SSP) %in% genes,], genes)
       
       res$BS.all = data.frame( PatientID = rownames(res$cl) ,
                                     BS = res$cl[,1])
       
       if( Prosigna){
         res$BS.all$BS.prosigna = res$BS.all$BS
       }
       
       output_res = res$BS.all
     }
     
     if (input$BSmethod == "sspbc"){

       res = BreastSubtypeR::BS_sspbc( as.matrix(reactive_files$data_input$x_SSP), ssp.name = "ssp.pam50")
       BS.all = data.frame( PatinetID = rownames(res),
                            BS = res,
                            row.names = rownames(res)
                            )

       
       if(Prosigna) {
         res_sspbc.prosigna = BS_sspbc( gene_expression_matrix = as.matrix(reactive_files$data_input$x_SSP), ssp.name= "ssp.subtype"  )
         BS.all$BS.prosigna = res_sspbc.prosigna[,1]
       } 
       
       output_res = BS.all
     }


       ## Allow downloading the result
     reactive_files$output_res = output_res
     output$download = downloadHandler(
       filename = function() {
         paste("results-", input$BSmethod, ".txt", sep = "")
         },
       content = function(file) {
         write.table(reactive_files$output_res, file, row.names = F, sep = "\t", quote = F)
         } )
  
       ## Allow visualization
       out= data.frame(PatientID = res$BS.all$PatientID,
                         Subtype = res$BS.all$BS )
       output$pie1 = renderPlot( BreastSubtypeR::Vis_pie(out) )
       output$pie2 = renderPlot( BreastSubtypeR::Vis_pie(out) )
      
      # Example function to generate plots conditionally based on the run button
      output$plotSection = renderUI({
        req(input$run)  # Wait until the run button is clicked
        
        # Once 'Run' is clicked, show the plot layout
        layout_columns(
          col_width = 2,
          card(plotOutput("pie1")),
          card(plotOutput("pie2"))
        )
      })
   

})

  # 
  # #### Reset functionality ####
  # observeEvent(input$reset, {
  #   # Clear the reactive values
  #   data_input$GEX = NULL
  #   data_input$clinic = NULL
  #   data_input$anno = NULL
  #   data_input$output_res = NULL
  #   
  #   # Reset all input fields
  #   updateFileInput(session, "GEX", label = "choose file", value = NULL)
  #   updateFileInput(session, "clinic", label = "choose file", value = NULL)
  #   updateFileInput(session, "anno", label = "choose file", value = NULL)
  #   updateSelectInput(session, "BSmethod", selected = "PAM50.parker")
  #   
  #   # Reset plot outputs
  #   output$pie1 = renderPlot(NULL)
  #   output$pie2 = renderPlot(NULL)
  #   
  #   # Clear any additional fields or reset other inputs as necessary
  # })
  
}