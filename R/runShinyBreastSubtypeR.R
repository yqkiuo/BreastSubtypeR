#' @title iBreastSubtypeR
#'
#' @description Launches the **iBreastSubtypeR** Shiny application, a graphical interface for the BreastSubtypeR package.
#'
#' The app allows users to run intrinsic molecular subtyping of breast cancer
#' using both nearest-centroid (NC-based) and single-sample predictor (SSP-based)
#' methods. It standardizes input and output handling, providing a reproducible
#' framework compatible with downstream R/Bioconductor tools.
#'
#' ## Workflow
#'
#' **Step 1 – Data mapping**
#' 
#' - Load input data from the current R session or upload expression, clinical,
#'   and feature annotation files (CSV/text format).
#' - Run the `Mapping()` step by clicking *Map Now*. Once complete, a message
#'   confirms: *"You may now proceed to Step 2."*
#'
#' **Step 2 – Subtyping analysis**
#'
#' - Select one or more subtyping methods (including `AUTO`) and configure
#'   parameters.
#' - Run the analysis. Upon completion, a message confirms: *"Analysis is complete."*
#' - Two diagnostic visualizations are displayed. Results can be downloaded as a
#'   text file. Users may re-run other methods directly without repeating Step 1.
#'
#' ## Notes
#'
#' - See `?BS_Multi` and `?Mapping` for detailed parameter descriptions.
#' - iBreastSubtypeR is intended for local use; no data are uploaded externally.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#'
#' @usage iBreastSubtypeR()
#'
#' @return
#' Launches the interactive Shiny application. The app provides subtyping calls
#' and ROR scores, downloadable as text files.
#'
#' @aliases iBreastSubtypeR
#' @name iBreastSubtypeR
#' @rdname iBreastSubtypeR
#'
#' @keywords BreastSubtypeR Shiny
#' 
#' @examples
#' \donttest{
#'   library(BreastSubtypeR)
#'   # Opens the BreastSubtypeR Shiny app in your browser
#'   iBreastSubtypeR()
#' }
#'
#' @export
iBreastSubtypeR <- function() {
    shinydeps <- c("shiny", "bslib")
    maskshinydeps <- shinydeps %in% installed.packages()
    if (any(!maskshinydeps)) {
        stop(
            sprintf(
                "Please install the following packages before running iBreastSubtypeR:\n\n  %s\n",
                paste(shinydeps[!maskshinydeps], collapse = ", ")
            )
        )
    }

    appDir <- system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR")
    if (appDir == "") {
        stop("The Shiny app directory 'ShinyBreastSubtypeR' was not found in the package.")
    }


    ## increase file upload limit (1 GB)
    options(shiny.maxRequestSize = 1000 * 1024^2)


    shiny::runApp(appDir, display.mode = "normal")
}
