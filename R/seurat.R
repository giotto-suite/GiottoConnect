#' @import GiottoClass
#' @import GiottoUtils
NULL




## Seurat object ####



#' @title Convert a Giotto object to a Seurat object
#' @name giottoToSeurat
#' @description Convert a giotto object (in memory or GiottoDisk-backed) to a
#' Seurat v5 object with the bundled mapping that fits its structure (see
#' [findMapping()]). The active spatial unit and feature type are converted:
#' expression (raw counts and normalized values), cell and feature metadata,
#' spatial locations and polygons as an FOV, dimension reductions and
#' nearest-neighbour networks.
#' @param gobject giotto object
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id
#' @param verbose logical
#' @returns Seurat object; see [mappingProvenance()] for what was converted and
#' which losses apply
#' @seealso [seuratToGiotto()], [convertWithMapping()], [readbackCheck()]
#' @export
giottoToSeurat <- function(gobject, mapping = NULL, verbose = FALSE) {
    profile <- profileObject(gobject)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "giotto")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this giotto object:",
                paste(found$reasons, collapse = "; ")), call. = FALSE)
        }
        mapping <- found$mapping
    }
    convertWithMapping(gobject, mapping, profile = profile, verbose = verbose)
}




#' @title Convert a Seurat object to a Giotto object
#' @name seuratToGiotto
#' @description Convert a Seurat object with the mapping that fits its
#' structure. The object is profiled, the matching mapping is selected from
#' those bundled with GiottoConnect (or the one given), and the conversion runs
#' through [convertWithMapping()]. Object versions are handled by choosing a
#' mapping, not by separate functions.
#' @param sobject Seurat object
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id.
#' By default the bundled mapping that applies to `sobject` is used.
#' @param source `NULL` for an in-memory giotto object, or a directory for a
#' GiottoDisk-backed object
#' @param verbose logical
#' @returns giotto object; see [mappingProvenance()] for what was converted
#' and which losses apply
#' @seealso [findMapping()], [convertWithMapping()], [validateMapping()]
#' @export
seuratToGiotto <- function(sobject, mapping = NULL, source = NULL, verbose = FALSE) {
    profile <- profileObject(sobject)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "seurat")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this Seurat object:",
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
    convertWithMapping(sobject, mapping, profile = profile, source = source,
        verbose = verbose)
}
