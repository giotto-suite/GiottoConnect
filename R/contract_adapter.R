# Ecosystem adapters and conversion targets. The engine (params, for_each,
# ops, hooks, validation, readback, repair) is ecosystem-neutral; everything
# that depends on how one ecosystem stores and reads its data sits behind an
# adapter (source side) or a target (output side). A mapping names its
# adapter with `mapping.source.ecosystem` and its target with
# `mapping.target.contract`.
#
# Adapter  — how to profile an input, list names, which accessors a mapping
#            may call and how a node is read.
# Target   — how prepared parts become the output object, which target keys
#            a part needs, and how the output is checked and read back.

# adapters ####

setClass("giottoConnectAdapter", representation("VIRTUAL",
    ecosystem = "character", package = "character"))

# Seurat objects, read through exported SeuratObject accessors
setClass("seuratAdapter", contains = "giottoConnectAdapter",
    prototype = list(ecosystem = "seurat", package = "SeuratObject"))

# giotto objects (in memory or GiottoDisk-backed), read through GiottoClass getters
setClass("giottoAdapter", contains = "giottoConnectAdapter",
    prototype = list(ecosystem = "giotto", package = "GiottoClass"))

setGeneric("adapterProfile", function(adapter, x, ...) standardGeneric("adapterProfile"))
setGeneric("adapterListPaths", function(adapter) standardGeneric("adapterListPaths"))
setGeneric("adapterAccessors", function(adapter) standardGeneric("adapterAccessors"))
setGeneric("adapterObjectLists", function(adapter) standardGeneric("adapterObjectLists"))
setGeneric("adapterRead", function(adapter, x, node, params) standardGeneric("adapterRead"))
# a GiottoDisk streaming input for a node, or NULL to read it into R
setGeneric("adapterStream", function(adapter, x, node, params) standardGeneric("adapterStream"))
# paths a node's accessor reads besides the node itself (e.g. the feature
# names of a matrix); they count as mapped
setGeneric("adapterReads", function(adapter, node, params) standardGeneric("adapterReads"))

# defaults: in-R objects profiled structurally, nodes read by an exported
# function of the adapter's package, by `[[`, or structurally ("path", for
# repair patches only)
setMethod("adapterProfile", "giottoConnectAdapter", function(adapter, x, ...) {
    profileObject(x, ...)
})

setMethod("adapterRead", "giottoConnectAdapter", function(adapter, x, node, params) {
    get <- node$get
    if (identical(get$fn, "path")) return(.resolve_path(x, .pf_bind(get$on, params)))
    fn <- if (identical(get$fn, "[[")) base::`[[` else getExportedValue(adapter@package, get$fn)
    target <- if (is.null(get$on)) x else .resolve_path(x, .pf_bind(get$on, params))
    args <- lapply(get$args %null% list(), function(a) {
        if (is.character(a)) .pf_bind(a, params) else a
    })
    do.call(fn, c(list(target), args))
})

setMethod("adapterStream", "giottoConnectAdapter", function(adapter, x, node, params) NULL)
setMethod("adapterReads", "giottoConnectAdapter", function(adapter, node, params) character())

setMethod("adapterListPaths", "seuratAdapter", function(adapter) c(
    Assays = "assays", Layers = "assays/{on}/layers",
    # Assay5 layers or v3 Assay slots
    AssayLayers = "assays/{on}/layers|assays/{on}",
    Reductions = "reductions", Graphs = "graphs", Neighbors = "neighbors",
    Images = "images", Boundaries = "images/{on}/boundaries", Molecules = "images/{on}/molecules",
    children = "{on}", columns = "{on}"
))
setMethod("adapterAccessors", "seuratAdapter", function(adapter) c(
    "LayerData", "Embeddings", "Loadings", "Stdev", "GetTissueCoordinates", "Misc", "Key",
    "[[", "path"
))
# a reduction's embeddings come with its loadings and stdev (read by the
# seurat_reduction_misc hook)
setMethod("adapterReads", "seuratAdapter", function(adapter, node, params) {
    p <- .pf_bind(node$path, params)
    if (!identical(node$get$fn, "Embeddings") || !grepl("/cell.embeddings$", p)) return(character())
    paste0(sub("cell.embeddings$", "", p), c("feature.loadings", "stdev"))
})
setMethod("adapterObjectLists", "seuratAdapter", function(adapter) {
    c("Assays", "Images", "Reductions", "Graphs", "Neighbors")
})

setMethod("adapterListPaths", "giottoAdapter", function(adapter) c(
    SpatUnits = "expression", FeatTypes = "expression/{on}", Expressions = "expression/{on}",
    SpatLocs = "spatial_locs/{on}", Polygons = "spatial_info",
    DimMethods = "dimension_reduction/cells/{on}", NNNames = "nn_network/{on}",
    # every cells-side reduction, whatever its method (name and method may differ)
    DimNames = "dimension_reduction/cells/{on}/*",
    SpatNets = "spatial_network/{on}", FeatInfo = "feat_info",
    # every nn network by name; its type (kNN, sNN) is the parent
    NNAll = "nn_network/{on}/*",
    Images = "images",
    children = "{on}", columns = "{on}"
))
setMethod("adapterAccessors", "giottoAdapter", function(adapter) c(
    "getExpression", "getCellMetadata", "getFeatureMetadata", "getSpatialLocations",
    "getPolygonInfo", "getDimReduction", "getNearestNetwork", "getSpatialNetwork",
    "getFeatureInfo", "getGiottoImage", "[[", "path"
))
setMethod("adapterObjectLists", "giottoAdapter", function(adapter) {
    c("SpatUnits", "FeatTypes", "Polygons", "DimMethods", "NNNames")
})
# a giotto object knows its cell and feature counts; no inference needed
setMethod("adapterProfile", "giottoAdapter", function(adapter, x, ...) {
    n <- tryCatch(list(c = length(GiottoClass::spatIDs(x)), f = length(GiottoClass::featIDs(x))),
        error = function(e) list(c = NULL, f = NULL))
    profileObject(x, n_cells = n$c, n_feat = n$f, ...)
})


# targets ####

setClass("giottoConnectTarget", representation("VIRTUAL", contract = "character"))
setClass("giottoTarget", contains = "giottoConnectTarget", prototype = list(contract = "giotto"))
setClass("seuratTarget", contains = "giottoConnectTarget", prototype = list(contract = "seurat"))

# `source`: GiottoDisk project of a giotto result; `output`: file a file
# target (e.g. .h5ad) writes
setGeneric("targetBuild", function(target, prepared, mapping, prov, source = NULL, output = NULL) {
    standardGeneric("targetBuild")
})
setGeneric("targetCheckKeys", function(target, mapping, out) standardGeneric("targetCheckKeys"))
setGeneric("targetSemantic", function(target, mapping, profile, result, out) {
    standardGeneric("targetSemantic")
})
setGeneric("targetComparators", function(target) standardGeneric("targetComparators"))

setMethod("targetBuild", "giottoTarget", function(target, prepared, mapping, prov, source = NULL,
    output = NULL) {
    contract <- giottoContract()
    src <- if (is.null(source)) NULL else {
        package_check("GiottoDisk")
        GiottoDisk::sourceCreate(source, type = "gDirSource")
    }
    parts <- lapply(prepared, function(pp) {
        data <- pp$data
        # streamed parts go straight from the input into the project
        if (methods::is(data, "exprInput")) {
            data <- GiottoDisk::sourceWrite(src, data, store_type = "parquetExpr")
        }
        .build_part(pp$spec, data, contract, mapping)
    })
    g <- .assemble(parts, source = src)
    .record_conversion(g, mapping_id = prov$mapping, mapping_version = prov$version, prov = prov)
})
setMethod("targetCheckKeys", "giottoTarget", function(target, mapping, out) {
    .vm_targets_giotto(mapping, out)
})
setMethod("targetSemantic", "giottoTarget", function(target, mapping, profile, result, out) {
    .vm_semantic(mapping, profile, result, out)
})
setMethod("targetComparators", "giottoTarget", function(target) .RB_COMPARE)

setMethod("targetBuild", "seuratTarget", function(target, prepared, mapping, prov, source = NULL,
    output = NULL) {
    if (!is.null(source)) stop(wrap_txt("`source` applies to Giotto targets only"), call. = FALSE)
    .record_seurat_conversion(.assemble_seurat(prepared), prov)
})
setMethod("targetCheckKeys", "seuratTarget", function(target, mapping, out) {
    .vm_targets_seurat(mapping, out)
})
setMethod("targetSemantic", "seuratTarget", function(target, mapping, profile, result, out) {
    .vm_semantic_seurat(mapping, profile, result, out)
})
setMethod("targetComparators", "seuratTarget", function(target) .RB_COMPARE_SEURAT)


# registries ####

.ADAPTERS <- list(seurat = "seuratAdapter", giotto = "giottoAdapter", anndata = "anndataAdapter",
    spatialexperiment = "speAdapter")
.TARGETS <- list(giotto = "giottoTarget", seurat = "seuratTarget", anndata = "anndataTarget",
    spatialexperiment = "speTarget")

.adapter_for <- function(x) {
    eco <- if (is.list(x)) x$mapping$source$ecosystem else x
    cls <- .ADAPTERS[[eco %null% ""]]
    if (is.null(cls)) {
        stop(wrap_txt(sprintf("No adapter for ecosystem '%s' (available: %s)", eco,
            paste(names(.ADAPTERS), collapse = ", "))), call. = FALSE)
    }
    methods::new(cls)
}

.target_for <- function(x) {
    con <- if (is.list(x)) x$mapping$target$contract else x
    cls <- .TARGETS[[con %null% ""]]
    if (is.null(cls)) {
        stop(wrap_txt(sprintf("No target for contract '%s' (available: %s)", con,
            paste(names(.TARGETS), collapse = ", "))), call. = FALSE)
    }
    methods::new(cls)
}

#' @title Ecosystem adapters and targets
#' @name adapterRegistry
#' @description List the ecosystems the mapping engine can read (adapters),
#' the accessors each lets a mapping call, and the outputs it can build
#' (targets).
#' @returns data.frame with columns `role`, `name`, `class`, `package`,
#' `accessors`
#' @examples
#' adapterRegistry()
#' @export
adapterRegistry <- function() {
    a <- lapply(names(.ADAPTERS), function(nm) {
        ad <- .adapter_for(nm)
        data.frame(role = "adapter", name = nm, class = .ADAPTERS[[nm]], package = ad@package,
            accessors = paste(adapterAccessors(ad), collapse = ", "), stringsAsFactors = FALSE)
    })
    t <- lapply(names(.TARGETS), function(nm) {
        data.frame(role = "target", name = nm, class = .TARGETS[[nm]], package = NA_character_,
            accessors = NA_character_, stringsAsFactors = FALSE)
    })
    do.call(rbind, c(a, t))
}
