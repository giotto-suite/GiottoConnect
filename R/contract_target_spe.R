# SpatialExperiment as a mapping target. Prepared parts are assembled into a
# SpatialExperiment, or a SpatialFeatureExperiment when the mapping's target
# class says so. Target keys per part:
#   expression      assay          ("counts", "logcounts", ...)
#   cell_metadata   -              colData (table with a cell_ID column)
#   feat_metadata   -              rowData (table with a feat_ID column)
#   spatial_locs    -              spatialCoords (cell_ID, sdimx, sdimy)
#   dim_reduction   reducedDim     PCA loadings / eigenvalues as scater's
#                                  "rotation" / "varExplained" attributes
#   nn_network      colPair, nn_type   SelfHits with distance and weight
#   spatial_network colPair            (SPE) or spatialGraph (SFE, listw)
#                   | spatialGraph
#   spatial_info    colGeometry    (SFE) sf polygons named by cell
#   feature_info    rowGeometry    (SFE) a MULTIPOINT per feature
#   images          imgData        image_id; list(image, scale) as a loaded
#                                  SpatialImage with its scaleFactor
# Edge networks keep their Giotto name and type in metadata(<hits>)$giotto.
# The conversion record is stored in metadata(x)$giottoConnect.

.SPE_TARGET <- list(
    expression = "assay",
    cell_metadata = character(),
    feat_metadata = character(),
    spatial_locs = character(),
    dim_reduction = "reducedDim",
    nn_network = c("colPair", "nn_type"),
    spatial_network = list("colPair", "spatialGraph"),
    spatial_info = "colGeometry",
    feature_info = "rowGeometry",
    images = "imgData"
)
# parts only a SpatialFeatureExperiment holds
.SFE_ONLY_KEYS <- c("spatialGraph", "colGeometry", "rowGeometry")

setClass("speTarget", contains = "giottoConnectTarget", prototype = list(contract = "spatialexperiment"))

setMethod("targetBuild", "speTarget", function(target, prepared, mapping, prov, source = NULL,
    output = NULL) {
    if (!is.null(source)) stop(wrap_txt("`source` applies to Giotto targets only"), call. = FALSE)
    sfe <- identical(mapping$mapping$target$class, "SpatialFeatureExperiment")
    x <- .assemble_spe(prepared, sfe = sfe)
    md <- S4Vectors::metadata(x)
    md$giottoConnect <- prov
    S4Vectors::`metadata<-`(x, value = md)
})
setMethod("targetCheckKeys", "speTarget", function(target, mapping, out) {
    .vm_targets_spe(mapping, out)
})
setMethod("targetSemantic", "speTarget", function(target, mapping, profile, result, out) {
    .vm_semantic_spe(mapping, profile, result, out)
})
setMethod("targetComparators", "speTarget", function(target) .RB_COMPARE_SPE)

# edge table -> SelfHits over `cells` with distance and weight columns
.spe_edges_to_hits <- function(e, cells, giotto) {
    e <- as.data.frame(e)
    i <- match(e$from, cells)
    j <- match(e$to, cells)
    ok <- !is.na(i) & !is.na(j)
    num <- function(v) if (is.null(v)) rep(NA_real_, sum(ok)) else as.numeric(v[ok])
    h <- S4Vectors::SelfHits(i[ok], j[ok], nnode = length(cells), distance = num(e$distance),
        weight = num(e$weight))
    S4Vectors::metadata(h)$giotto <- giotto
    h
}

# edge table -> spdep listw (binary neighbours, Giotto weights as glist)
.spe_edges_to_listw <- function(e, cells, giotto) {
    package_check("spdep")
    e <- as.data.frame(e)
    i <- match(e$from, cells)
    j <- match(e$to, cells)
    ok <- !is.na(i) & !is.na(j)
    w <- if (is.null(e$weight)) rep(1, sum(ok)) else as.numeric(e$weight[ok])
    # both directions, as spdep neighbour lists are
    from <- c(i[ok], j[ok])
    to <- c(j[ok], i[ok])
    w <- c(w, w)
    nb <- lapply(seq_along(cells), function(k) {
        hit <- which(from == k)
        o <- order(to[hit])
        list(to = as.integer(to[hit][o]), w = w[hit][o])
    })
    neigh <- lapply(nb, function(v) if (length(v$to)) v$to else 0L)
    glist <- lapply(nb, function(v) if (length(v$to)) v$w else NULL)
    class(neigh) <- "nb"
    attr(neigh, "region.id") <- cells
    attr(neigh, "sym") <- TRUE
    lw <- spdep::nb2listw(neigh, glist = glist, style = "B", zero.policy = TRUE)
    attr(lw, "giotto") <- giotto
    lw
}

# a MULTIPOINT per feature (empty for features without transcripts)
.spe_points_to_sf <- function(d, feats) {
    d <- as.data.frame(d)
    by <- split(seq_len(nrow(d)), factor(d$feat_ID, levels = feats))
    geo <- lapply(by, function(k) sf::st_multipoint(cbind(d$x[k], d$y[k])))
    out <- sf::st_sf(geometry = sf::st_sfc(geo))
    rownames(out) <- feats
    out
}

# list(image, scale) -> a SpatialExperiment image placed by its scale factor,
# or for a SpatialFeatureExperiment a SpatRasterImage placed by its extent in
# spatialCoords units (y down: the raster's first row is the image's bottom
# row, as terra reads a PNG)
.spe_image <- function(img, sfe = FALSE) {
    if (!sfe) return(SpatialExperiment::SpatialImage(grDevices::as.raster(img$image)))
    a <- img$image
    r <- terra::rast(a[rev(seq_len(dim(a)[1L])), , , drop = FALSE] * 255)
    terra::ext(r) <- c(0, dim(a)[2L] / img$scale, 0, dim(a)[1L] / img$scale)
    SpatialFeatureExperiment::SpatRasterImage(r)
}

# an imgData image back as an H x W x 3 array of 0..1 values (top row first)
.spe_image_array <- function(im) {
    if (methods::is(im, "SpatRasterImage")) {
        a <- terra::as.array(methods::as(im, "SpatRaster"))
        return(a[rev(seq_len(dim(a)[1L])), , seq_len(min(3L, dim(a)[3L])), drop = FALSE] / 255)
    }
    m <- as.matrix(SpatialExperiment::imgRaster(im))
    rgb <- grDevices::col2rgb(m) / 255
    array(t(rgb), c(nrow(m), ncol(m), 3L))
}

.assemble_spe <- function(prepared, sfe = FALSE) {
    package_check("SpatialExperiment", repository = "Bioc:SpatialExperiment")
    if (sfe) package_check("SpatialFeatureExperiment", repository = "Bioc:SpatialFeatureExperiment")
    by_part <- function(type) Filter(function(pp) pp$spec$part == type, prepared)
    ex <- by_part("expression")
    if (!length(ex)) stop(wrap_txt("A SpatialExperiment needs at least one expression part"), call. = FALSE)
    mat <- function(m) {
        if (is.matrix(m)) return(m)
        if (methods::is(m, "Matrix")) return(methods::as(m, "CsparseMatrix"))
        methods::as(m, "dgCMatrix")
    }
    assays <- lapply(ex, function(pp) mat(pp$data))
    names(assays) <- vapply(ex, function(pp) pp$spec$to$assay, character(1))
    feats <- rownames(assays[[1L]])
    cells <- colnames(assays[[1L]])
    assays <- lapply(assays, function(m) m[feats, cells, drop = FALSE])

    table_rows <- function(d, id, ids) {
        d <- as.data.frame(d)
        rn <- as.character(d[[id]])
        d[[id]] <- NULL
        d <- d[match(ids, rn), , drop = FALSE]
        rownames(d) <- ids
        S4Vectors::DataFrame(d, check.names = FALSE)
    }
    cd <- S4Vectors::DataFrame(row.names = cells)
    rd <- S4Vectors::DataFrame(row.names = feats)
    for (pp in by_part("cell_metadata")) cd <- table_rows(pp$data, "cell_ID", cells)
    for (pp in by_part("feat_metadata")) rd <- table_rows(pp$data, "feat_ID", feats)

    xy <- matrix(0, length(cells), 2L, dimnames = list(cells, c("x", "y")))
    for (pp in by_part("spatial_locs")) {
        d <- as.data.frame(pp$data)
        d <- d[match(cells, d$cell_ID), , drop = FALSE]
        xy <- cbind(x = d$sdimx, y = d$sdimy)
        rownames(xy) <- cells
    }
    geoms <- list()
    for (pp in by_part("spatial_info")) {
        g <- pp$data[match(cells, rownames(pp$data)), , drop = FALSE]
        rownames(g) <- cells
        geoms[[pp$spec$to$colGeometry]] <- g
    }
    x <- SpatialExperiment::SpatialExperiment(assays = assays, colData = cd, rowData = rd, spatialCoords = xy)
    if (sfe) {
        # geometries are added afterwards: given to the constructor, they
        # would replace spatialCoords by the polygon centroids
        x <- SpatialFeatureExperiment::toSpatialFeatureExperiment(x)
        for (nm in names(geoms)) {
            x <- SpatialFeatureExperiment::`colGeometry<-`(x, nm, value = geoms[[nm]])
        }
    }

    for (pp in by_part("dim_reduction")) {
        d <- pp$data
        m <- if (methods::is(d, "dimObj")) d[] else as.matrix(d)
        m <- m[cells, , drop = FALSE]
        misc <- if (methods::is(d, "dimObj")) d@misc else NULL
        if (length(misc$loadings)) attr(m, "rotation") <- as.matrix(misc$loadings)
        if (length(misc$eigenvalues)) attr(m, "varExplained") <- as.numeric(misc$eigenvalues)
        SingleCellExperiment::reducedDim(x, pp$spec$to$reducedDim) <- m
    }
    for (pp in by_part("nn_network")) {
        to <- pp$spec$to
        SingleCellExperiment::colPair(x, to$colPair) <- .spe_edges_to_hits(pp$data, cells,
            list(name = to$name %null% to$colPair, type = to$nn_type))
    }
    for (pp in by_part("spatial_network")) {
        to <- pp$spec$to
        giotto <- list(name = to$name %null% to$colPair %null% to$spatialGraph, type = "spatial")
        if (!is.null(to$spatialGraph)) {
            x <- SpatialFeatureExperiment::`colGraph<-`(x, to$spatialGraph,
                value = .spe_edges_to_listw(pp$data, cells, giotto))
        } else {
            SingleCellExperiment::colPair(x, to$colPair) <- .spe_edges_to_hits(pp$data, cells, giotto)
        }
    }
    imgs <- by_part("images")
    if (length(imgs)) {
        SpatialExperiment::imgData(x) <- S4Vectors::DataFrame(
            sample_id = rep(x$sample_id[1L], length(imgs)),
            image_id = unname(vapply(imgs, function(pp) pp$spec$to$imgData, character(1))),
            data = I(unname(lapply(imgs, function(pp) .spe_image(pp$data, sfe)))),
            scaleFactor = unname(vapply(imgs, function(pp) pp$data$scale, numeric(1))))
    }
    for (pp in by_part("feature_info")) {
        x <- SpatialFeatureExperiment::`rowGeometry<-`(x, pp$spec$to$rowGeometry,
            value = .spe_points_to_sf(pp$data, feats))
    }
    x
}


# validation ####

# every part carries the target keys its SpatialExperiment place needs; SFE
# places need a SpatialFeatureExperiment target class
.vm_targets_spe <- function(mapping, out) {
    sfe <- identical(mapping$mapping$target$class, "SpatialFeatureExperiment")
    for (p in mapping$parts) {
        opts <- .SPE_TARGET[[p$part]]
        if (is.null(opts)) {
            out$add("target", p$id, "FAIL", sprintf("SpatialExperiment has no place for %s", p$part))
            next
        }
        if (!is.list(opts)) opts <- list(opts)
        keys <- setdiff(names(p$to), "name")
        fits <- vapply(opts, function(need) setequal(keys, need), logical(1))
        lvl <- if (any(fits)) "OK" else "FAIL"
        msg <- if (lvl == "OK") sprintf("%s -> SpatialExperiment", p$part) else
            sprintf("SpatialExperiment %s needs keys %s (has: %s)", p$part,
                paste(vapply(opts, paste, character(1), collapse = ", "), collapse = " or "),
                paste(keys, collapse = ", "))
        if (lvl == "OK" && !sfe && any(keys %in% .SFE_ONLY_KEYS)) {
            lvl <- "FAIL"
            msg <- sprintf("%s needs target class SpatialFeatureExperiment", paste(keys, collapse = ", "))
        }
        out$add("target", p$id, lvl, msg)
    }
}

.vm_semantic_spe <- function(mapping, profile, x, out) {
    cells <- colnames(x)
    if (!is.na(profile$n_cells)) {
        out$add("dims", NA, if (length(cells) == profile$n_cells) "OK" else "WARN",
            sprintf("%d cells in the SpatialExperiment, %d in source", length(cells), profile$n_cells))
    }
    an <- SummarizedExperiment::assayNames(x)
    if ("counts" %in% an) {
        m <- SummarizedExperiment::assay(x, "counts")[, utils::head(cells, 99L), drop = FALSE]
        v <- if (methods::is(m, "sparseMatrix")) m@x else as.vector(m)
        ok <- all(v == round(v)) && all(v >= 0)
        out$add("raw_counts", "counts", if (ok) "OK" else "FAIL",
            sprintf("counts are %sinteger-valued and non-negative", if (ok) "" else "NOT "))
    }
    xy <- SpatialExperiment::spatialCoords(x)
    rng <- max(diff(range(xy[, 1])), diff(range(xy[, 2])))
    out$add("coord_units", "spatialCoords", if (rng < 50 || rng > 1e6) "WARN" else "OK",
        sprintf("coordinate range %.3g", rng))
    for (r in SingleCellExperiment::reducedDimNames(x)) {
        n <- nrow(SingleCellExperiment::reducedDim(x, r))
        out$add("dim_reduction", r, if (n == length(cells)) "OK" else "FAIL",
            sprintf("%d rows for %d cells", n, length(cells)))
    }
    out$add("validity", class(x)[1L], if (isTRUE(methods::validObject(x, test = TRUE))) "OK" else "FAIL",
        "validObject()")
}


# readback: compare the source parts with the built object ####

.rb_spe_edges <- function(h, cells) {
    paste(cells[S4Vectors::from(h)], cells[S4Vectors::to(h)])
}

.RB_COMPARE_SPE <- list(
    expression = function(spec, data, x, n_cells) {
        data <- .rb_as_matrix(data)
        got <- SummarizedExperiment::assay(x, spec$to$assay)
        use <- .rb_sample(colnames(data), n_cells)
        list(
            .rb_row("placement", TRUE, paste0("assay ", spec$to$assay)),
            .rb_row("ids", setequal(colnames(got), colnames(data)) && setequal(rownames(got), rownames(data)),
                sprintf("%d cells, %d features", ncol(got), nrow(got))),
            .rb_row("values", .rb_same(as.numeric(as.matrix(got[rownames(data), use, drop = FALSE])),
                as.numeric(as.matrix(data[, use, drop = FALSE]))), sprintf("%d sampled cells", length(use)))
        )
    },
    cell_metadata = function(spec, data, x, n_cells) {
        cd <- as.data.frame(SummarizedExperiment::colData(x), optional = TRUE)
        .rb_table(data, cbind(cell_ID = rownames(cd), cd), "cell_ID")
    },
    feat_metadata = function(spec, data, x, n_cells) {
        rd <- as.data.frame(SummarizedExperiment::rowData(x), optional = TRUE)
        .rb_table(data, cbind(feat_ID = rownames(rd), rd), "feat_ID")
    },
    spatial_locs = function(spec, data, x, n_cells) {
        xy <- SpatialExperiment::spatialCoords(x)
        d <- as.data.frame(data)
        d <- d[match(colnames(x), d$cell_ID), , drop = FALSE]
        list(
            .rb_row("placement", TRUE, "spatialCoords"),
            .rb_row("values", .rb_same(c(xy[, 1], xy[, 2]), c(d$sdimx, d$sdimy)), sprintf("%d cells", nrow(xy)))
        )
    },
    dim_reduction = function(spec, data, x, n_cells) {
        m <- SingleCellExperiment::reducedDim(x, spec$to$reducedDim)
        misc <- if (methods::is(data, "dimObj")) data@misc else NULL
        coords <- if (methods::is(data, "dimObj")) data[] else as.matrix(data)
        list(
            .rb_row("placement", TRUE, paste0("reducedDim ", spec$to$reducedDim)),
            .rb_row("values", .rb_same(as.numeric(m), as.numeric(coords[colnames(x), , drop = FALSE])), ""),
            if (length(misc$loadings)) {
                .rb_row("loadings", .rb_same(as.numeric(attr(m, "rotation")), as.numeric(misc$loadings)),
                    sprintf("%d features", nrow(misc$loadings)))
            }
        )
    },
    nn_network = function(spec, data, x, n_cells) {
        h <- SingleCellExperiment::colPair(x, spec$to$colPair)
        got <- .rb_spe_edges(h, colnames(x))
        d <- as.data.frame(data)
        exp <- paste(d$from, d$to)
        k <- match(exp, got)
        list(
            .rb_row("placement", TRUE, paste0("colPair ", spec$to$colPair)),
            .rb_row("edges", setequal(got, exp), sprintf("%d edges", length(exp))),
            .rb_row("weight", .rb_same(S4Vectors::mcols(h)$weight[k], as.numeric(d$weight)), "")
        )
    },
    spatial_network = function(spec, data, x, n_cells) {
        d <- as.data.frame(data)
        if (!is.null(spec$to$colPair)) {
            h <- SingleCellExperiment::colPair(x, spec$to$colPair)
            got <- .rb_spe_edges(h, colnames(x))
            where <- paste0("colPair ", spec$to$colPair)
        } else {
            lw <- SpatialFeatureExperiment::colGraph(x, spec$to$spatialGraph)
            ids <- colnames(x)
            from <- rep(seq_along(lw$neighbours), lengths(lw$neighbours))
            to <- unlist(lw$neighbours)
            keep <- to > 0L & from < to
            got <- paste(ids[from[keep]], ids[to[keep]])
            where <- paste0("spatialGraph ", spec$to$spatialGraph)
            # listw neighbours are symmetric: each pair once, in column order
            i <- match(d$from, ids)
            j <- match(d$to, ids)
            d <- data.frame(from = ids[pmin(i, j)], to = ids[pmax(i, j)])
        }
        exp <- unique(paste(d$from, d$to))
        list(
            .rb_row("placement", TRUE, where),
            .rb_row("edges", setequal(got, exp), sprintf("%d edges", length(exp)))
        )
    },
    spatial_info = function(spec, data, x, n_cells) {
        g <- SpatialFeatureExperiment::colGeometry(x, spec$to$colGeometry)
        a <- sf::st_bbox(data)
        b <- sf::st_bbox(g)
        list(
            .rb_row("placement", TRUE, paste0("colGeometry ", spec$to$colGeometry)),
            .rb_row("ids", setequal(rownames(g), rownames(data)), sprintf("%d polygons", nrow(g))),
            .rb_row("geometry", .rb_same(as.numeric(a), as.numeric(b)), "bounding box")
        )
    },
    images = function(spec, data, x, n_cells) {
        d <- SpatialExperiment::imgData(x)
        i <- match(spec$to$imgData, d$image_id)
        got <- if (!is.na(i)) .spe_image_array(d$data[[i]])
        list(
            .rb_row("placement", !is.na(i), paste0("imgData ", spec$to$imgData)),
            # an R raster holds 8-bit colours
            .rb_row("values", !is.na(i) && .rb_same(got, data$image, tol = 1 / 255),
                sprintf("%s px", paste(dim(data$image)[2:1], collapse = " x "))),
            .rb_row("scale", !is.na(i) && .rb_same(d$scaleFactor[i], data$scale), sprintf("scaleFactor %.4g", data$scale))
        )
    },
    feature_info = function(spec, data, x, n_cells) {
        g <- SpatialFeatureExperiment::rowGeometry(x, spec$to$rowGeometry)
        got <- vapply(sf::st_geometry(g), function(p) nrow(as.matrix(p)), integer(1))
        names(got) <- rownames(g)
        a <- table(as.character(data$feat_ID))
        a <- a[names(a) %in% rownames(g)]
        list(
            .rb_row("placement", TRUE, paste0("rowGeometry ", spec$to$rowGeometry)),
            .rb_row("ids", all(as.integer(a) == got[names(a)]) && sum(got) == sum(a),
                sprintf("%d transcripts of %d features", sum(a), length(a)))
        )
    }
)
