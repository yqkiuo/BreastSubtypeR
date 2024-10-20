
#' @title iBreastSubtypeR 
#' 
#' @description Starts an interactive BreastSubtypeR shiny web app.
#' 
#' BreastSubtypeR collected breast cancer subtyping methods, 
#' near centrioids based (NC-based) and single sample prediction based (SSP-based). 
#' It employs standardized input and output and provides a unified framework that 
#' is highly compatible with other R packages in the gene expression profiling field.   
#' 
#' The `iBreastSubtypeR()` function starts an interactive shiny web app that allows
#' the user to configure the arguments of subtyping functions and
#' runs it on the computer. Please see the manual page of the functions for a description of the arguments and their
#' default and alternative values.   
#' 
#' The input data may be loaded from the users workspace or by selecting a CSV/text
#' file for the expression data, a CSV/text file for clinical information and 
#' a CSV/text file for feature annotations. 
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
#' \dontrun{
#' res <- iBreastSubtypeR() ## this will open your browser with the GSVA shiny web app
#' }
#'
#' @export
iBreastSubtypeR = function() {
  
  shinydeps = c("shiny", "data.table","magrittr")
  maskshinydeps = shinydeps %in% installed.packages()
  if (any(!maskshinydeps))
    stop(sprintf("Please install the following packages to use the GSVA WebApp:\n\n  %s\n",
                 paste(shinydeps[!maskshinydeps], collapse=", ")))
  
  appDir = system.file("ShinyBreastSubtypeR", package="BreastSubtypeR")
  if (appDir == "")
    stop("The BreastSubtypeR Shiny app cannot be found within the package.")
  
  
  shiny::runApp(appDir)
}
