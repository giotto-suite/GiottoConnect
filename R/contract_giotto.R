# The Giotto contract: for each part a mapping can produce, the subobject
# class, constructor, setter and nesting keys. GiottoClass keeps this
# implicit (S4 dispatch, list_* tables); .contract_build() declares the
# part -> function names once and derives everything else from GiottoClass
# itself by building a tiny object.

.CONTRACT_SCHEMA_VERSION <- "0.1.0"

.CONTRACT_PARTS <- list(
    expression = list(slot = "expression", constructor = "createExprObj", data_arg = "expression_data",
        setter = "setExpression", list_fn = "list_expression"),
    cell_metadata = list(slot = "cell_metadata", constructor = "createCellMetaObj", data_arg = "metadata",
        setter = "setCellMetadata", list_fn = "list_cell_metadata"),
    feat_metadata = list(slot = "feat_metadata", constructor = "createFeatMetaObj", data_arg = "metadata",
        setter = "setFeatureMetadata", list_fn = "list_feat_metadata"),
    spatial_locs = list(slot = "spatial_locs", constructor = "createSpatLocsObj", data_arg = "coordinates",
        setter = "setSpatialLocations", list_fn = "list_spatial_locations"),
    spatial_info = list(slot = "spatial_info", constructor = "createGiottoPolygonsFromDfr", data_arg = "segmdfr",
        setter = "setPolygonInfo", list_fn = "list_spatial_info"),
    dim_reduction = list(slot = "dimension_reduction", constructor = "createDimObj", data_arg = "coordinates",
        setter = "setDimReduction", list_fn = "list_dim_reductions"),
    nn_network = list(slot = "nn_network", constructor = "createNearestNetObj", data_arg = "network",
        setter = "setNearestNetwork", list_fn = "list_nearest_networks"),
    spatial_network = list(slot = "spatial_network", constructor = "createSpatNetObj", data_arg = "network",
        setter = "setSpatialNetwork", list_fn = "list_spatial_networks"),
    # transcripts; `method` names the S4 method whose arguments apply
    feature_info = list(slot = "feat_info", constructor = "createGiottoPoints", method = "data.frame",
        data_arg = "x", setter = "setFeatureInfo", list_fn = "list_feature_info"),
    # a SpatRaster placed in spatial coordinates (its extent)
    images = list(slot = "images", constructor = "createGiottoLargeImage", data_arg = "raster_object",
        setter = "setGiottoImage", list_fn = "list_images")
)

# columns of list_* tables that are nesting keys (the rest are descriptive)
.CONTRACT_KEY_COLS <- c("spat_unit", "feat_type", "name", "type", "nn_type",
    "method", "spat_info", "feat_info")

#' @title The Giotto contract
#' @name giottoContract
#' @description Machine-readable description of the Giotto parts a mapping
#' can target: for each part the subobject class, constructor (and its
#' arguments), setter, nesting keys, matrix orientation and whether the part
#' is written to disk when the giotto object has a GiottoDisk source.
#' @param source `"auto"` (bundled contract for the installed GiottoClass
#' version, otherwise generated live), `"bundled"` or `"live"`
#' @returns a `giottoContract` list
#' @examples
#' ct <- giottoContract("bundled")
#' names(ct$parts)
#' @export
giottoContract <- function(source = c("auto", "bundled", "live")) {
    source <- match.arg(source)
    files <- list.files(system.file("contract", package = "GiottoConnect"),
        pattern = "^giotto-contract-.*\\.json$", full.names = TRUE
    )
    vers <- sub("^giotto-contract-(.*)\\.json$", "\\1", basename(files))
    installed <- if (requireNamespace("GiottoClass", quietly = TRUE)) {
        as.character(utils::packageVersion("GiottoClass"))
    } else {
        NA_character_
    }
    latest_bundled <- function() {
        if (!length(files)) {
            stop(wrap_txt("No bundled Giotto contract found"), call. = FALSE)
        }
        files[order(numeric_version(vers), decreasing = TRUE)][1L]
    }
    ct <- switch(source,
        live = .contract_build(),
        bundled = .contract_read(latest_bundled()),
        auto = {
            if (!is.na(installed) && installed %in% vers) {
                .contract_read(files[vers == installed])
            } else if (!is.na(installed)) {
                warning(wrap_txt(sprintf(
                    "No bundled contract for GiottoClass %s (bundled: %s); generated live",
                    installed, paste(vers, collapse = ", ")
                )), call. = FALSE)
                .contract_build()
            } else {
                f <- latest_bundled()
                warning(wrap_txt(sprintf(
                    "GiottoClass is not installed; using bundled contract %s",
                    basename(f)
                )), call. = FALSE)
                .contract_read(f)
            }
        }
    )
    structure(ct, class = "giottoContract")
}

#' @export
print.giottoContract <- function(x, ...) {
    cat(sprintf("<giottoContract> GiottoClass %s (schema %s, generated %s)\n",
        x$giottoclass_version, x$schema_version, x$generated_at))
    for (nm in names(x$parts)) {
        p <- x$parts[[nm]]
        cat(sprintf("  %-14s %-15s %-28s keys: %-28s%s\n", nm, p$class,
            p$constructor, paste(p$nesting, collapse = "/"),
            if (isTRUE(p$disk_backable)) " [disk]" else ""))
    }
    invisible(x)
}

.contract_read <- function(file) {
    jsonlite::read_json(file, simplifyVector = TRUE)
}

.contract_write <- function(ct, dir) {
    f <- file.path(dir, sprintf("giotto-contract-%s.json", ct$giottoclass_version))
    jsonlite::write_json(unclass(ct), f, auto_unbox = TRUE, pretty = TRUE)
    f
}

# Build the contract from the installed GiottoClass.
.contract_build <- function() {
    package_check("GiottoClass")
    ns <- asNamespace("GiottoClass")
    tiny <- .contract_tiny_parts()
    g <- GiottoClass::setGiotto(GiottoClass::giotto(), unname(tiny),
        verbose = FALSE
    )
    parts <- lapply(names(.CONTRACT_PARTS), function(nm) {
        d <- .CONTRACT_PARTS[[nm]]
        ctor <- get(d$constructor, envir = ns)
        fm <- if (is.null(d$method)) formals(ctor) else .contract_method_formals(d$constructor, d$method, ns)
        no_default <- vapply(fm, function(a) {
            is.symbol(a) && !nzchar(as.character(a))
        }, logical(1))
        if (!d$data_arg %in% names(fm)) {
            stop(sprintf("%s has no argument '%s'", d$constructor, d$data_arg), call. = FALSE)
        }
        listed <- get(d$list_fn, envir = ns)(g)
        list(
            slot = d$slot,
            class = class(tiny[[nm]])[1L],
            constructor = d$constructor,
            data_arg = d$data_arg,
            args = setdiff(names(fm), "..."),
            required_args = setdiff(names(fm)[no_default], "..."),
            setter = d$setter,
            list_fn = d$list_fn,
            nesting = intersect(names(listed), .CONTRACT_KEY_COLS),
            disk_backable = .contract_setter_writes_source(d$setter, ns)
        )
    })
    names(parts) <- names(.CONTRACT_PARTS)

    e <- tiny$expression
    orient <- if (identical(GiottoClass::spatIDs(e), colnames(e[]))) {
        "features_x_cells"
    } else {
        "cells_x_features"
    }
    list(
        contract = "giotto",
        schema_version = .CONTRACT_SCHEMA_VERSION,
        giottoclass_version = as.character(utils::packageVersion("GiottoClass")),
        generated_at = format(Sys.Date()),
        manifest_schema = .contract_manifest_schema(),
        orientation = list(expression = orient),
        parts = parts
    )
}

# A setter is disk-backable when one of its methods routes through the
# giotto object's GiottoDisk source.
# formals of an S4 method, including those of its `.local` function
.contract_method_formals <- function(generic, signature, ns) {
    m <- methods::getMethod(generic, signature, where = ns)
    fm <- formals(m)
    b <- body(m)
    for (e in as.list(b)) {
        if (is.call(e) && identical(e[[1L]], as.name("<-")) && identical(e[[2L]], as.name(".local"))) {
            return(formals(eval(e[[3L]])))
        }
    }
    fm
}

.contract_setter_writes_source <- function(setter, ns) {
    src <- if (methods::isGeneric(setter, where = ns)) {
        m <- methods::findMethods(setter, where = ns)
        unlist(lapply(m, function(mm) deparse(body(mm))))
    } else {
        deparse(body(get(setter, envir = ns)))
    }
    any(grepl("@source", src, fixed = TRUE))
}

.contract_manifest_schema <- function() {
    f <- list.files(system.file("schema", package = "GiottoClass"),
        pattern = "^giotto-manifest-.*\\.json$"
    )
    if (length(f)) sub("\\.json$", "", sort(f, decreasing = TRUE)[1L]) else NA_character_
}

# one minimal instance per part, built only with GiottoClass constructors
.contract_tiny_parts <- function() {
    cells <- c("c1", "c2")
    feats <- c("f1", "f2", "f3")
    m <- matrix(c(1, 0, 2, 0, 3, 1), nrow = 3L, dimnames = list(feats, cells))
    dt <- data.table::data.table
    sq <- dt(
        poly_ID = rep(cells, each = 4L),
        x = c(0, 1, 1, 0, 2, 3, 3, 2), y = c(0, 0, 1, 1, 0, 0, 1, 1)
    )
    emb <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2L,
        dimnames = list(cells, c("Dim.1", "Dim.2"))
    )
    suppressMessages(list(
        expression = GiottoClass::createExprObj(m, name = "raw",
            spat_unit = "cell", feat_type = "rna"),
        cell_metadata = GiottoClass::createCellMetaObj(dt(cell_ID = cells, a = 1:2),
            spat_unit = "cell", feat_type = "rna", verbose = FALSE),
        feat_metadata = GiottoClass::createFeatMetaObj(dt(feat_ID = feats, b = 1:3),
            spat_unit = "cell", feat_type = "rna", verbose = FALSE),
        spatial_locs = GiottoClass::createSpatLocsObj(
            dt(sdimx = c(0.5, 2.5), sdimy = c(0.5, 0.5), cell_ID = cells),
            name = "raw", spat_unit = "cell", verbose = FALSE),
        spatial_info = GiottoClass::createGiottoPolygonsFromDfr(sq, name = "cell",
            calc_centroids = FALSE, verbose = FALSE),
        dim_reduction = GiottoClass::createDimObj(emb, name = "pca", method = "pca",
            spat_unit = "cell", feat_type = "rna"),
        nn_network = GiottoClass::createNearestNetObj(name = "kNN.pca",
            network = dt(from = "c1", to = "c2", weight = 1, distance = 0),
            nn_type = "kNN", spat_unit = "cell", feat_type = "rna"),
        spatial_network = GiottoClass::createSpatNetObj(
            network = .edges_to_igraph(dt(from = "c1", to = "c2", distance = 2, weight = 0.5), cells),
            name = "Delaunay_network", spat_unit = "cell"),
        feature_info = GiottoClass::createGiottoPoints(
            data.frame(x = c(0.2, 2.2, 2.7), y = c(0.5, 0.4, 0.6), feat_ID = feats),
            feat_type = "rna", verbose = FALSE),
        images = GiottoClass::createGiottoLargeImage(
            terra::rast(array(1:12, c(2L, 2L, 3L)), extent = terra::ext(0, 3, -1, 1)),
            name = "image", use_rast_ext = TRUE, verbose = FALSE)
    ))
}
