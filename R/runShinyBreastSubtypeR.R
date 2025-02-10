#' @title iBreastSubtypeR
#'
#' @description Starts an interactive BreastSubtypeR shiny web app.
#'
#'   BreastSubtypeR integrates intrinsic molecular subtyping methods for breast
#'   cancer, including nearest-centroid (NC-based) and single-sample predictor
#'   (SSP-based) approaches. It employs standardized input and output formats,
#'   providing a unified framework that is highly compatible with other R
#'   packages in the gene expression profiling field.
#'
#'   The `iBreastSubtypeR()` function launches an interactive Shiny web
#'   application. This app enables users to configure the arguments of subtyping
#'   functions and execute subtyping on their local computer. For detailed
#'   descriptions of the arguments, including their default and alternative
#'   values, please refer to the manual pages of the respective functions.
#'
#'   Step 1:
#'
#'   The input data can be loaded from the user's workspace or by selecting a
#'   CSV/text file for the expression data, a CSV/text file for clinical
#'   information, and a CSV/text file for feature annotations.
#'
#'   After loading the necessary files, users can click the "Map Now" button
#'   once and wait for the notification. If the Mapping() function runs
#'   successfully, a message will appear stating, "You may now proceed to Step
#'   2."
#'
#'
#'
#'   Step 2:
#'
#'   Users can select the desired subtyping method and adjust the relevant
#'   parameters to conduct their analysis. Once the analysis is complete, a
#'   message will indicate, "Analysis is complete." Two visualizations will be
#'   displayed, and you will have the option to download the results as a text
#'   file. If you wish to continue your analysis, you can directly run another
#'   method without needing to repeat Step 1.
#'
#' @import bslib
#' @import shiny
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
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
#' @keywords BreastSubtypeR Shiny
#' @examples
#'
#' library(BreastSubtypeR)
#'
#' # This will open your browser with the BreastSubtypeR shiny web app
#' \donttest{
#' iBreastSubtypeR()
#' }
#'
#' @export
iBreastSubtypeR <- function() {
    shinydeps <- c("shiny", "bslib")
    maskshinydeps <- shinydeps %in% installed.packages()
    if (any(!maskshinydeps)) {
        stop(
            sprintf(
                "Please install the following packages :\n\n  %s\n",
                paste(shinydeps[!maskshinydeps], collapse = ", ")
            )
        )
    }

    appDir <- system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR")
    if (appDir == "") {
        stop("The iBreastSubtypeR cannot be found within the package.")
    }


    ## increase file limit
    options(shiny.maxRequestSize = 1000 * 1024^2)


    shiny::runApp(appDir)
}
