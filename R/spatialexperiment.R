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
