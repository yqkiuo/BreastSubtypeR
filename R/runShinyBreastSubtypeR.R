
#' @title iBreastSubtypeR 
#' 
#' @description Starts an interactive BreastSubtypeR shiny web app.
#' 
#' BreastSubtypeR integrates intrinsic molecular subtyping methods for breast cancer, 
#' including nearest-centroid (NC-based) and single-sample predictor (SSP-based) approaches. 
#' It employs standardized input and output formats, 
#' offering a unified framework that is highly compatible with other R packages in the gene expression profiling field.   
#' 
#' The `iBreastSubtypeR()` function starts an interactive shiny web app that allows
#' the user to configure the arguments of subtyping functions and
#' runs it on the local computer. Please see the manual page of the functions for a description of the arguments and their
#' default and alternative values.   
#' 
#' The input data may be loaded from the users workspace or by selecting a CSV/text
#' file for the expression data, a CSV/text file for clinical information and 
#' a CSV/text file for feature annotations.    
#' 
#' After loading necessary files, users only need to run/click "Mapping" button once and waiting for the notification.
#' If Mapping() function runs well, it will say "You can go for Step 2.". 
#' Then, users can select any subtyping method and change relevant parameters to perform analysis. You can see "Analysis complete." when the analysis is finished. 
#' Two plots would be visualized as well. Users can download the result as text file.
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
  
  
  require(bslib, quietly = T)
  require(shiny, quietly = T)
  require(BreastSubtypeR, quietly = T)
  
  ## increase file limit
  options(shiny.maxRequestSize=1000*1024^2) 
  
  
  shiny::runApp(appDir)
}
