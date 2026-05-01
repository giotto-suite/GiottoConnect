# Stage 1: Export a Giotto object to SpatialData-Zarr
#
# Uses GiottoClass::giottoToSpatialData() as the primary conversion path.
# This keeps Giotto object introspection logic in GiottoClass and lets
# GiottoConnect focus on the downstream export/deployment steps.
#
# TODO(giotto-upstream): giottoToSpatialData() is the canonical conversion
# path. Any missing slot exports (UMAP in obsm, polygon shapes, cell metadata)
# should be added to giottoToSpatialData() in GiottoClass rather than
# worked around here in GiottoConnect.

#' Export a Giotto object to SpatialData-Zarr
#'
#' Converts a Giotto object to a SpatialData store and writes it as a Zarr
#' archive that Vitessce can consume. The output directory will contain a
#' `<dataset_name>.zarr/` folder with a consolidated `.zmetadata` file.
#'
#' @param gobject A Giotto object (`giotto` class from GiottoClass).
#' @param output_dir Character. Path to the directory where the zarr store
#'   will be written.
#' @param dataset_name Character. Base name for the output zarr store
#'   (without `.zarr` extension). Default: `"dataset"`.
#' @param spat_unit Character. Spatial unit to export (e.g., `"cell"`).
#'   Passed to `giottoToSpatialData()`. Default: `"cell"`.
#' @param feat_type Character. Feature type to export (e.g., `"rna"`).
#'   Passed to `giottoToSpatialData()`. Default: `"rna"`.
#' @param include_expression Logical. If `TRUE`, include the expression
#'   matrix in the export. Not yet implemented in v1; reserved for future
#'   use. Default: `FALSE`.
#' @param include_images Logical. If `TRUE`, include spatial images (as
#'   OME-Zarr). Not yet implemented in v1; reserved for future use.
#'   Default: `FALSE`.
#' @param overwrite Logical. If `TRUE`, overwrite an existing zarr store at
#'   the target path. Default: `FALSE`.
#' @param python_env Character. Name of the conda environment to use.
#'   Default: `"giotto-connect"`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return Invisibly returns the path to the written zarr store.
#' @export
#'
#' @examples
#' \dontrun{
#' library(GiottoData)
#' gobject <- GiottoData::loadGiottoMini("merfish")
#' zarr_path <- export_giotto_to_zarr(
#'     gobject    = gobject,
#'     output_dir = tempdir(),
#'     dataset_name = "mini_merfish"
#' )
#' }
export_giotto_to_zarr <- function(
    gobject,
    output_dir,
    dataset_name       = "dataset",
    spat_unit          = "cell",
    feat_type          = "rna",
    include_expression = FALSE,
    include_images     = FALSE,
    overwrite          = FALSE,
    python_env         = "giotto-connect",
    verbose            = TRUE
) {
    # -- Input validation -----------------------------------------------------
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(output_dir)
    checkmate::assert_string(dataset_name, min.chars = 1L)
    checkmate::assert_string(spat_unit)
    checkmate::assert_string(feat_type)
    checkmate::assert_flag(include_expression)
    checkmate::assert_flag(include_images)
    checkmate::assert_flag(overwrite)
    checkmate::assert_flag(verbose)

    # -- Prepare output path --------------------------------------------------
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    zarr_path <- file.path(output_dir, paste0(dataset_name, ".zarr"))

    if (dir.exists(zarr_path)) {
        if (!overwrite) {
            stop(
                "[GiottoConnect] Zarr store already exists at: ", zarr_path,
                "\n  Use overwrite = TRUE to replace it."
            )
        }
        unlink(zarr_path, recursive = TRUE)
    }

    # -- Ensure Python env is ready -------------------------------------------
    .ensure_python()

    # -- Convert Giotto → SpatialData -----------------------------------------
    if (verbose) message("[GiottoConnect] Stage 1: Converting Giotto object to SpatialData ...")

    # giottoToSpatialData() returns a Python SpatialData object via reticulate.
    # If this function needs to be updated to include UMAP embeddings or
    # polygon shapes, those changes belong in GiottoClass.
    sdata <- GiottoClass::giottoToSpatialData(
        gobject   = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # -- Write SpatialData to zarr --------------------------------------------
    if (verbose) message("[GiottoConnect]   Writing zarr to: ", zarr_path)

    sdata$write(zarr_path)

    # -- Consolidate metadata (required by Vitessce) --------------------------
    if (verbose) message("[GiottoConnect]   Consolidating zarr metadata ...")

    zarr_mod <- .gc_python_env$zarr
    zarr_mod$consolidate_metadata(zarr_path)

    if (verbose) message("[GiottoConnect] Stage 1 complete: ", zarr_path)

    invisible(zarr_path)
}
