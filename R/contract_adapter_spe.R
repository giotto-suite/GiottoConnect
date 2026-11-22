# SpatialExperiment adapter. Covers SpatialFeatureExperiment (a subclass) too.
# Objects are profiled structurally like any S4 object; accessors go through
# the SummarizedExperiment / SingleCellExperiment / SpatialExperiment /
# SpatialFeatureExperiment getters and return plain R types (data.frames,
# matrices), whatever DataFrame or DelayedArray class holds the data.

setClass("speAdapter", contains = "giottoConnectAdapter",
    prototype = list(ecosystem = "spatialexperiment", package = "SpatialExperiment"))

setMethod("adapterListPaths", "speAdapter", function(adapter) c(
    Assays = "assays/data/listData",
    Reductions = "int_colData/listData/reducedDims/listData",
    AltExps = "int_colData/listData/altExps/listData",
    ColPairs = "int_colData/listData/colPairs/listData",
    # SpatialFeatureExperiment
    ColGeometries = "int_colData/listData/colGeometries/listData",
    RowGeometries = "int_elementMetadata/listData/rowGeometries/listData",
    # column graphs of one sample
    SpatialGraphs = "int_metadata/spatialGraphs/listData/{on}/col",
    children = "{on}", columns = "{on}"
))
setMethod("adapterAccessors", "speAdapter", function(adapter) names(.SPE_ACCESSORS))
setMethod("adapterObjectLists", "speAdapter", function(adapter) {
    c("Assays", "Reductions", "AltExps", "ColGeometries", "RowGeometries")
})
setMethod("adapterProfile", "speAdapter", function(adapter, x, ...) {
    profileObject(x, n_cells = ncol(x), n_feat = nrow(x), ...)
})
setMethod("adapterRead", "speAdapter", function(adapter, x, node, params) {
    args <- lapply(node$get$args %null% list(), function(a) {
        if (is.character(a)) .pf_bind(a, params) else a
    })
    if (identical(node$get$fn, "path")) return(.resolve_path(x, .pf_bind(node$get$on, params)))
    do.call(.SPE_ACCESSORS[[node$get$fn]], c(list(x), args))
})

# accessors ####

# a matrix R can use: dense or dgCMatrix (DelayedArray/HDF5 assays are read)
.spe_matrix <- function(m) {
    if (methods::is(m, "Matrix") || is.matrix(m)) return(m)
    methods::as(m, "dgCMatrix")
}

.spe_frame <- function(df, ids) {
    out <- as.data.frame(df, optional = TRUE)
    rownames(out) <- ids
    out
}

.spe_ns <- function(pkg) {
    package_check(pkg)
    asNamespace(pkg)
}

.SPE_ACCESSORS <- list(
    speAssay = function(x, i) .spe_matrix(.spe_ns("SummarizedExperiment")$assay(x, i)),
    speColData = function(x) .spe_frame(.spe_ns("SummarizedExperiment")$colData(x), colnames(x)),
    speRowData = function(x) .spe_frame(.spe_ns("SummarizedExperiment")$rowData(x), rownames(x)),
    speSpatialCoords = function(x) {
        m <- .spe_ns("SpatialExperiment")$spatialCoords(x)
        rownames(m) <- colnames(x)
        m
    },
    # keeps the attributes scater stores (rotation, varExplained, percentVar)
    speReducedDim = function(x, type) {
        m <- .spe_ns("SingleCellExperiment")$reducedDim(x, type)
        rownames(m) <- colnames(x)
        m
    },
    speAltExpAssay = function(x, e, i = "counts") {
        ae <- .spe_ns("SingleCellExperiment")$altExp(x, e)
        .spe_matrix(.spe_ns("SummarizedExperiment")$assay(ae, i))
    },
    speAltExpRowData = function(x, e) {
        ae <- .spe_ns("SingleCellExperiment")$altExp(x, e)
        .spe_frame(.spe_ns("SummarizedExperiment")$rowData(ae), rownames(ae))
    },
    # sf data.frame of one colGeometry (cell or nucleus segmentation, centroids)
    sfeColGeometry = function(x, type) {
        .spe_ns("SpatialFeatureExperiment")$colGeometry(x, type)
    },
    # transcripts as x, y, feat_ID (a MULTIPOINT per feature in rowGeometries)
    sfeTxSpots = function(x, type = "txSpots") {
        g <- .spe_ns("SpatialFeatureExperiment")$rowGeometry(x, type)
        xy <- sf::st_coordinates(g)
        ids <- rownames(x)[match(rownames(g), rownames(x))]
        if (anyNA(ids)) ids <- rownames(g)
        grp <- xy[, ncol(xy)]
        data.frame(x = xy[, "X"], y = xy[, "Y"], feat_ID = ids[grp], stringsAsFactors = FALSE)
    },
    sfeSpatialGraph = function(x, type, sample_id = 1L, MARGIN = 2L) {
        SFE <- .spe_ns("SpatialFeatureExperiment")
        SFE$spatialGraph(x, type = type, MARGIN = MARGIN, sample_id = sample_id)
    }
)
