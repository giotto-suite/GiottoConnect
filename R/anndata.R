#' @import GiottoClass
#' @import GiottoUtils
NULL




## anndata object ####


#' @title Convert an AnnData .h5ad file to a Giotto object
#' @name anndataToGiotto
#' @description Convert an AnnData (.h5ad) file with the mapping that fits
#' its structure. The file is read directly (hdf5r, no python): it is
#' profiled, the matching bundled mapping is selected (or the one given) and
#' the conversion runs through [convertWithMapping()].
#' @param anndata_path path to the .h5ad file
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id.
#' By default the bundled mapping that applies to the file is used.
#' @param source `NULL` for an in-memory giotto object, or a directory for a
#' GiottoDisk-backed object
#' @param verbose logical
#' @returns giotto object; see [mappingProvenance()] for what was converted
#' and which losses apply
#' @seealso [findMapping()], [convertWithMapping()], [readbackCheck()]
#' @export
anndataToGiotto <- function(anndata_path, mapping = NULL, source = NULL, verbose = FALSE) {
    profile <- adapterProfile(.adapter_for("anndata"), anndata_path)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "anndata")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this AnnData file:",
                paste(found$reasons, collapse = "; ")), call. = FALSE)
        }
        if (found$status == "repair_needed") {
            warning(wrap_txt(sprintf(
                "Structure differs from mapping '%s' (fingerprint %s, expected %s); converting with re-bound params",
                found$mapping$mapping$id, found$fingerprint,
                found$mapping$mapping$fingerprint %null% "unset"
            )), call. = FALSE)
        }
        mapping <- found$mapping
    }
    convertWithMapping(anndata_path, mapping, profile = profile, source = source,
        verbose = verbose)
}
