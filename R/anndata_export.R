#' @title Convert a Giotto object to an AnnData .h5ad file
#' @name giottoToAnnData
#' @description Write one spatial unit and feature type of a giotto object
#' (in memory or GiottoDisk-backed) as an AnnData `.h5ad` file in scanpy's
#' layout, without python (the file is written with anndataR). `X` holds the
#' normalized values (raw counts when there are none) and `layers/counts` the
#' raw counts; other expression values become layers. Reductions go to
#' `obsm/X_<name>` (PCA loadings to `varm/PCs`), networks to
#' `obsp/<name>_distances` / `_connectivities` with a `uns/<name>` record, so
#' [anndataToGiotto()] brings them back under the same names.
#' @param gobject giotto object
#' @param file path of the `.h5ad` file to write
#' @param spat_unit,feat_type spatial unit and feature type to export (default:
#' the active ones)
#' @param mapping optional `giottoMapping`, mapping file or bundled mapping id
#' @param verbose logical
#' @returns the path of the written file (invisibly); see [mappingProvenance()]
#' for what was converted
#' @seealso [anndataToGiotto()], [readbackCheck()]
#' @export
giottoToAnnData <- function(gobject, file, spat_unit = NULL, feat_type = NULL, mapping = NULL,
    verbose = FALSE) {
    profile <- adapterProfile(.adapter_for("giotto"), gobject)
    pinned <- .pinned_pair(spat_unit, feat_type)
    if (is.null(mapping)) {
        found <- findMapping(profile, ecosystem = "giotto", target = "anndata")
        if (found$status == "none") {
            stop(wrap_txt("No bundled mapping applies to this giotto object:",
                paste(found$reasons, collapse = "; ")), call. = FALSE)
        }
        mapping <- found$mapping
    }
    invisible(convertWithMapping(gobject, mapping, profile = profile, verbose = verbose,
        params = pinned, output = file))
}
