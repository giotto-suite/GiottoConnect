# Internal Python environment management for GiottoConnect
# All Python module handles are cached in .gc_python_env to avoid repeated
# reticulate::import() calls across function invocations in the same session.

.gc_python_env <- new.env(parent = emptyenv())

# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

#' Set up the conda environment and cache module handles
#'
#' @param envname Character. Name of the conda environment.
#' @param conda Character. Path to conda binary, or "auto".
#' @noRd
.setup_python_env <- function(
    envname = "giotto-connect",
    conda   = "auto"
) {
    required_pkgs <- c("anndata", "zarr", "vitessce", "spatialdata", "numpy")

    # Create env + install packages if the env doesn't already exist,
    # otherwise just activate it.
    existing_envs <- tryCatch(
        reticulate::conda_list(conda = conda)$name,
        error = function(e) character(0)
    )

    if (!envname %in% existing_envs) {
        message("[GiottoConnect] Creating conda env '", envname, "' ...")
        reticulate::conda_install(
            envname  = envname,
            packages = required_pkgs,
            conda    = conda,
            pip      = TRUE   # vitessce is pip-installable; spatialdata too
        )
    }

    reticulate::use_condaenv(envname, required = TRUE)

    # Cache module handles
    .gc_python_env$anndata     <- reticulate::import("anndata",     convert = FALSE)
    .gc_python_env$zarr        <- reticulate::import("zarr",        convert = FALSE)
    .gc_python_env$spatialdata <- reticulate::import("spatialdata", convert = FALSE)
    .gc_python_env$vitessce    <- reticulate::import("vitessce",    convert = FALSE)
    .gc_python_env$numpy       <- reticulate::import("numpy",       convert = FALSE)
    .gc_python_env$pandas      <- reticulate::import("pandas",      convert = FALSE)

    invisible(NULL)
}

#' Ensure Python env is initialized; set up lazily on first call.
#' @noRd
.ensure_python <- function() {
    if (is.null(.gc_python_env$zarr)) {
        .setup_python_env()
    }
}

# -------------------------------------------------------------------------
# Exported user-facing setup function
# -------------------------------------------------------------------------

#' Set up or verify the Python environment required by GiottoConnect
#'
#' Call this once per session before using any GiottoConnect export or
#' visualization functions. It creates a conda environment named
#' `"giotto-connect"` (if it does not already exist) and installs the
#' required Python packages: `anndata`, `zarr`, `vitessce`, `spatialdata`,
#' `numpy`.
#'
#' @param envname Character. Name of the conda environment.
#'   Default: `"giotto-connect"`.
#' @param conda Character. Path to the conda binary, or `"auto"` to let
#'   reticulate find it automatically. Default: `"auto"`.
#'
#' @return Invisibly `NULL`. Called for side effects.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_giotto_connect_env()
#' }
setup_giotto_connect_env <- function(
    envname = "giotto-connect",
    conda   = "auto"
) {
    .setup_python_env(envname = envname, conda = conda)
    message("[GiottoConnect] Python environment '", envname, "' is ready.")
    invisible(NULL)
}
