#' Launch the iBreastSubtypeR Shiny app
#'
#' @description Starts the Shiny UI bundled with the BreastSubtypeR package.
#' The launcher can (optionally) attach Shiny/Bslib so UI/server can use
#' unqualified functions like `tags`, `icon`, `fileInput`, etc.
#'
#' @param attach Character vector of packages to attach before launch.
#'        Defaults to c("shiny","bslib"). Set to character(0) to skip attaching.
#' @param attach_tidyverse Logical; if TRUE and tidyverse is installed, it will
#'        be attached quietly for the session (purely optional).
#' @param max_upload_mb Numeric; Shiny upload size limit (in MB). Default 1000.
#' @return The value returned by `shiny::runApp()` (usually `invisible(NULL)`).
#'
#' @examples
#' if (interactive()) {
#'     iBreastSubtypeR()
#'     iBreastSubtypeR(attach = character(0))
#' }
#'
#' @export
iBreastSubtypeR <- function(
        attach = c("shiny", "bslib"),
        attach_tidyverse = FALSE,
        max_upload_mb = 1000) {
    # Load Shiny/bslib (and optionally tidyverse) namespaces for this session.
    # ui.R/server.R use qualified calls and shiny::runApp() attaches shiny
    # itself, so loading the namespaces is sufficient.
    .load_app_dependencies(attach)
    if (isTRUE(attach_tidyverse) && "tidyverse" %in% rownames(utils::installed.packages())) {
        .load_app_dependencies("tidyverse")
    }

    # App directory shipped inside the package
    appDir <- system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR")
    if (identical(appDir, "") || !dir.exists(appDir)) {
        stop("Shiny app directory 'inst/ShinyBreastSubtypeR' not found in the package.", call. = FALSE)
    }

    # Increase upload limit (MB -> bytes)
    options(shiny.maxRequestSize = max_upload_mb * 1024^2)

    shiny::runApp(appDir, display.mode = "normal")
}

#' Check that app dependencies are installed and load their namespaces
#'
#' @param pkgs Character vector of package names.
#' @return `invisible(NULL)`. Stops with the list of missing packages, or with
#'   the name of a package that is installed but cannot be loaded.
#' @keywords internal
#' @noRd
.load_app_dependencies <- function(pkgs) {
    pkgs <- unique(pkgs)
    if (!length(pkgs)) {
        return(invisible(NULL))
    }
    installed <- rownames(utils::installed.packages())
    missing <- setdiff(pkgs, installed)
    if (length(missing)) {
        stop(
            sprintf(
                "Please install required package(s) before launching the app: %s",
                paste(missing, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    for (p in pkgs) {
        loaded <- suppressPackageStartupMessages(
            requireNamespace(p, quietly = TRUE)
        )
        if (!isTRUE(loaded)) {
            stop(
                sprintf("Package '%s' is installed but could not be loaded.", p),
                call. = FALSE
            )
        }
    }
    invisible(NULL)
}

#' (Deprecated) Run iBreastSubtypeR
#'
#' Internal wrapper kept for back-compat. Use [iBreastSubtypeR()] directly.
#' @keywords internal
#' @noRd
runShinyBreastSubtypeR <- function(...) {
    iBreastSubtypeR(...)
}
