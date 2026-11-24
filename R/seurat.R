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
#' @param spat_unit,feat_type spatial unit and feature type to export (default:
#' the active ones)
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id
#' @param verbose logical
#' @returns Seurat object; see [mappingProvenance()] for what was converted and
#' which losses apply
#' @seealso [seuratToGiotto()], [convertWithMapping()], [readbackCheck()]
#' @export
giottoToSeurat <- function(gobject, spat_unit = NULL, feat_type = NULL, mapping = NULL,
    verbose = FALSE) {
    profile <- adapterProfile(.adapter_for("giotto"), gobject)
    pinned <- .pinned_pair(spat_unit, feat_type)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "giotto", target = "seurat")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this giotto object:",
                paste(found$reasons, collapse = "; ")), call. = FALSE)
        }
        mapping <- found$mapping
    }
    convertWithMapping(gobject, mapping, profile = profile, verbose = verbose, params = pinned)
}

# the spat_unit / feat_type a giotto export is pinned to (NULL: the active ones)
.pinned_pair <- function(spat_unit = NULL, feat_type = NULL) {
    out <- list(spat_unit = spat_unit, feat_type = feat_type)
    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out)) out else NULL
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
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- Matrix::rsparsematrix(20, 30, density = 0.3,
#'         rand.x = function(n) stats::rpois(n, 5) + 1)
#'     dimnames(counts) <- list(paste0("gene", 1:20), paste0("cell", 1:30))
#'     s <- SeuratObject::CreateSeuratObject(counts)
#'     # coordinates kept as metadata columns are found as well
#'     s$x <- stats::runif(30, 0, 500)
#'     s$y <- stats::runif(30, 0, 500)
#'     g <- seuratToGiotto(s)
#'     g
#'     mappingProvenance(g)$losses
#' }
#' @export
seuratToGiotto <- function(sobject, mapping = NULL, source = NULL, verbose = FALSE) {
    profile <- adapterProfile(.adapter_for("seurat"), sobject)
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
