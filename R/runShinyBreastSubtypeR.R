
#' @title iBreastSubtypeR 
#' 
#' @description Starts an interactive BreastSubtypeR shiny web app.
#' 
#' BreastSubtypeR integrates intrinsic molecular subtyping methods for breast cancer, 
#' including nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches. 
#' It employs standardized input and output formats, 
#' providing a unified framework that is highly compatible with other R packages in the gene expression profiling field.   
#' 
#' The `iBreastSubtypeR()` function launches an interactive Shiny web application. This app enables users 
#' to configure the arguments of subtyping functions and
#' execute subtyping on their local computer. For detailed descriptions of the arguments,
#' including their default and alternative values, please refer to the manual pages of the respective functions.
#' 
#' Step 1:   
#' 
#' The input data may be loaded from the users workspace or by selecting a CSV/text
#' file for the expression data, a CSV/text file for clinical information and 
#' a CSV/text file for feature annotations.    
#' 
#' After loading the necessary files, users can click the "Map" button once and wait for the notification. 
#' If the Mapping() function runs successfully, a message will appear stating, "You can proceed to Step 2."
#' 
#' 
#' 
#' Step 2:  
#' 
#' Users can select any subtyping method and adjust the relevant parameters to conduct their analysis. 
#' Once the analysis is complete, you will see a message indicating, "Analysis complete." 
#' Two visualizations will be displayed, and you will have the option to download the results as a text file. 
#' If you wish to continue your analysis, you can directly run another method without needing to repeat Step 1.
#' 
#'    
#' 
#' @usage iBreastSubtypeR()
#' 
#' @return A table with subtyping and ROR score 
#' 
#'
#' @aliases iBreastSubtypeR
#' 
#' @name iBreastSubtypeR
#' 
#' @rdname iBreastSubtypeR
#' 
#' @keywords BreastSubtypeR shiny
#' @examples
#' 
#' iBreastSubtypeR() ## this will open your browser with the BreastSubtypeR shiny web app

#'
#' @export
iBreastSubtypeR = function() {
  
  shinydeps = c("shiny","bslib")
  maskshinydeps = shinydeps %in% installed.packages()
  if (any(!maskshinydeps))
    stop(sprintf("Please install the following packages to use the BreastSubtypeR WebApp:\n\n  %s\n",
                 paste(shinydeps[!maskshinydeps], collapse=", ")))
  
  appDir = system.file("ShinyBreastSubtypeR", package="BreastSubtypeR")
  if (appDir == "")
    stop("The BreastSubtypeR Shiny app cannot be found within the package.")
  
  
  suppressMessages(suppressWarnings({
    require(bslib, quietly = T)
    require(shiny, quietly = T)
    require(BreastSubtypeR, quietly = T)
  }))

  ## increase file limit
  options(shiny.maxRequestSize=1000*1024^2) 
  
  
  shiny::runApp(appDir)
}
