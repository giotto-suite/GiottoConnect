#' @title Convert a SpatialExperiment to a Giotto object
#' @name spatialExperimentToGiotto
#' @description Convert a SpatialExperiment or SpatialFeatureExperiment with the
#' mapping that fits its structure (the layouts SpatialExperimentIO, VisiumIO
#' and SpatialFeatureExperiment's readers produce). The object is profiled, the
#' matching bundled mapping is selected (or the one given) and the conversion
#' runs through [convertWithMapping()].
#' @param x SpatialExperiment or SpatialFeatureExperiment
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id.
#' By default the bundled mapping that applies to `x` is used.
#' @param source `NULL` for an in-memory giotto object, or a directory for a
#' GiottoDisk-backed object
#' @param verbose logical
#' @returns giotto object; see [mappingProvenance()] for what was converted
#' and which losses apply
#' @seealso [findMapping()], [convertWithMapping()], [readbackCheck()]
#' @examples
#' if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- Matrix::rsparsematrix(20, 30, density = 0.3,
#'         rand.x = function(n) stats::rpois(n, 5) + 1)
#'     dimnames(counts) <- list(paste0("gene", 1:20), paste0("cell", 1:30))
#'     spe <- SpatialExperiment::SpatialExperiment(
#'         assays = list(counts = counts),
#'         spatialCoords = cbind(x = stats::runif(30, 0, 500), y = stats::runif(30, 0, 500))
#'     )
#'     g <- spatialExperimentToGiotto(spe)
#'     g
#'     mappingProvenance(g)$mapping
#' }
#' @export
spatialExperimentToGiotto <- function(x, mapping = NULL, source = NULL, verbose = FALSE) {
    profile <- adapterProfile(.adapter_for("spatialexperiment"), x)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "spatialexperiment")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this SpatialExperiment:",
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
    convertWithMapping(x, mapping, profile = profile, source = source, verbose = verbose)
}

#' @title Convert a Giotto object to a SpatialExperiment
#' @name giottoToSpatialExperiment
#' @description Convert one spatial unit and feature type of a giotto object
#' (in memory or GiottoDisk-backed) to a SpatialExperiment, or to a
#' SpatialFeatureExperiment that also holds the polygons, transcripts and
#' spatial networks. Raw counts become the `counts` assay, normalized values
#' `logcounts`; reductions use scater's names (`PCA`, `UMAP`, `TSNE`) with
#' PCA loadings and variances as attributes; nearest-neighbour networks
#' become colPairs that keep their Giotto name and type, so
#' [spatialExperimentToGiotto()] brings them back. Coordinates are in image
#' space (y down).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type to convert
#' (default: the active ones)
#' @param sfe `"auto"` builds a SpatialFeatureExperiment when the spatial unit
#' has polygons, the feature type has transcripts or the spatial unit has
#' spatial networks (and SpatialFeatureExperiment is installed); `"yes"` or
#' `"no"` decide it
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id
#' (overrides `sfe`)
#' @param verbose logical
#' @returns SpatialExperiment or SpatialFeatureExperiment; see
#' [mappingProvenance()] for what was converted
#' @seealso [spatialExperimentToGiotto()], [readbackCheck()]
#' @export
giottoToSpatialExperiment <- function(gobject, spat_unit = NULL, feat_type = NULL,
    sfe = c("auto", "yes", "no"), mapping = NULL, verbose = FALSE) {
    sfe <- match.arg(sfe)
    profile <- adapterProfile(.adapter_for("giotto"), gobject)
    if (is.null(mapping)) {
        if (sfe == "auto") {
            sfe <- if (requireNamespace("SpatialFeatureExperiment", quietly = TRUE) &&
                .has_sfe_parts(gobject, spat_unit, feat_type)) "yes" else "no"
        }
        mapping <- if (sfe == "yes") "giotto-sfe" else "giotto-spe"
    }
    convertWithMapping(gobject, mapping, profile = profile, verbose = verbose,
        params = .pinned_pair(spat_unit, feat_type))
}

# polygons, transcripts or spatial networks for the pair being converted
.has_sfe_parts <- function(gobject, spat_unit = NULL, feat_type = NULL) {
    spat_unit <- spat_unit %null% GiottoClass::activeSpatUnit(gobject)
    feat_type <- feat_type %null% GiottoClass::activeFeatType(gobject)
    poly <- spat_unit %in% names(methods::slot(gobject, "spatial_info"))
    pts <- feat_type %in% names(methods::slot(gobject, "feat_info"))
    nets <- length(methods::slot(gobject, "spatial_network")[[spat_unit]]) > 0L
    poly || pts || nets
}
