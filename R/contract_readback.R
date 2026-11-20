#' @title Read a converted object back and compare it with its source
#' @name readbackCheck
#' @description For every part a mapping converted, read the Giotto content
#' back through Giotto getters and compare it with what the mapping prepared
#' from the source (the accessor output after hooks and ops): placement under
#' the declared nesting keys, IDs, shape, columns and values. Disk-backed
#' parts are compared on a sample of cells read through the store.
#'
#' Each comparison is `MATCH`, `DECLARED` (the difference is covered by a loss
#' declared in the mapping or recorded during conversion) or `UNDECLARED`.
#' @param x the source object
#' @param result object converted from `x` (giotto, or Seurat for a
#' Giotto-to-Seurat mapping)
#' @param mapping mapping used for the conversion. Defaults to the one
#' recorded in `result` (see [mappingProvenance()]).
#' @param n_cells number of cells sampled for expression values
#' @returns a `giottoReadback` data.frame with columns `part_id`, `part`,
#' `aspect`, `status`, `detail`; attribute `pass` is `TRUE` when nothing is
#' `UNDECLARED`
#' @export
readbackCheck <- function(x, result, mapping = NULL, n_cells = 99L) {
    prov <- mappingProvenance(result)
    if (is.null(mapping)) {
        if (is.null(prov)) {
            stop(wrap_txt("`result` has no conversion record; pass `mapping`"), call. = FALSE)
        }
        mapping <- prov$mapping
    }
    mapping <- if (inherits(mapping, "giottoMapping")) mapping else readMapping(mapping)
    if (!is.null(prov) && !identical(prov$version, mapping$mapping$version)) {
        warning(wrap_txt(sprintf("result was converted with %s v%s, checking against v%s",
            prov$mapping, prov$version, mapping$mapping$version)), call. = FALSE)
    }
    package_check(.adapter_for(mapping)@package)
    profile <- adapterProfile(.adapter_for(mapping), x)
    res <- .resolve_parts(mapping, profile)
    mapping <- res$mapping
    params <- res$params
    mapping <- .expand_for_each(mapping, profile, params)
    ctx <- .new_ctx(mapping, params, profile)
    prepared <- .prepare_parts(x, mapping, params, ctx)
    losses <- .collect_losses(mapping, ctx, params$values)

    compare <- targetComparators(.target_for(mapping))
    rows <- list()
    for (pp in prepared) {
        spec <- pp$spec
        res <- tryCatch(
            compare[[spec$part]](spec, pp$data, result, n_cells),
            error = function(e) list(.rb_row("placement", FALSE, conditionMessage(e), "placement"))
        )
        for (r in Filter(Negate(is.null), res)) {
            status <- if (r$ok) {
                "MATCH"
            } else if (.rb_declared(losses, spec$id, r$field)) {
                "DECLARED"
            } else {
                "UNDECLARED"
            }
            rows[[length(rows) + 1L]] <- data.frame(part_id = spec$id, part = spec$part,
                aspect = r$aspect, status = status, detail = r$detail, stringsAsFactors = FALSE)
        }
    }
    out <- do.call(rbind, rows)
    structure(out, class = c("giottoReadback", "data.frame"),
        pass = !any(out$status == "UNDECLARED"))
}

#' @export
print.giottoReadback <- function(x, ...) {
    n <- table(factor(x$status, levels = c("MATCH", "DECLARED", "UNDECLARED")))
    cat(sprintf("<giottoReadback> %s: %d MATCH, %d DECLARED, %d UNDECLARED\n",
        if (isTRUE(attr(x, "pass"))) "PASS" else "FAIL",
        n[["MATCH"]], n[["DECLARED"]], n[["UNDECLARED"]]))
    show <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
    print(show[, c("status", "part_id", "aspect", "detail")], row.names = FALSE, right = FALSE)
    invisible(x)
}

.rb_row <- function(aspect, ok, detail, field = aspect) {
    list(aspect = aspect, ok = isTRUE(ok), detail = detail, field = field)
}

# a difference is declared when a loss names the part (or the whole mapping)
# and the affected field
.rb_declared <- function(losses, part_id, field) {
    if (!nrow(losses)) return(FALSE)
    hit <- (is.na(losses$part_id) | losses$part_id == part_id) & !is.na(losses$field) &
        losses$field == field
    any(hit)
}

.rb_sample <- function(ids, n) {
    if (length(ids) <= n) return(ids)
    ids[as.integer(seq(1L, length(ids), length.out = n))]
}

.rb_as_matrix <- function(m) {
    if (methods::is(m, "Matrix") || is.matrix(m)) m else methods::as(m, "dgCMatrix")
}

.rb_same <- function(a, b, tol = 1e-8) isTRUE(all.equal(a, b, tolerance = tol, check.attributes = FALSE))


# comparators: list of .rb_row() per part ####

.RB_COMPARE <- list(
    expression = function(spec, data, g, n_cells) {
        to <- spec$to
        e <- GiottoClass::getExpression(g, values = to$name, spat_unit = to$spat_unit,
            feat_type = to$feat_type, output = "exprObj")
        cells <- GiottoClass::spatIDs(e)
        feats <- GiottoClass::featIDs(e)
        use <- .rb_sample(colnames(data), n_cells)
        got <- .rb_as_matrix(e[][, use, drop = FALSE])
        list(
            .rb_row("placement", TRUE, sprintf("%s/%s/%s", to$spat_unit, to$feat_type, to$name)),
            .rb_row("ids", setequal(cells, colnames(data)) && setequal(feats, rownames(data)),
                sprintf("%d cells, %d features", length(cells), length(feats))),
            .rb_row("values", .rb_same(as.numeric(got[rownames(data), use]),
                as.numeric(data[, use])), sprintf("%d sampled cells", length(use)))
        )
    },
    cell_metadata = function(spec, data, g, n_cells) {
        dt <- GiottoClass::pDataDT(g, spat_unit = spec$to$spat_unit, feat_type = spec$to$feat_type)
        .rb_table(data, dt, "cell_ID")
    },
    feat_metadata = function(spec, data, g, n_cells) {
        dt <- GiottoClass::fDataDT(g, spat_unit = spec$to$spat_unit, feat_type = spec$to$feat_type)
        .rb_table(data, dt, "feat_ID")
    },
    spatial_locs = function(spec, data, g, n_cells) {
        sl <- GiottoClass::getSpatialLocations(g, spat_unit = spec$to$spat_unit,
            name = spec$to$name, output = "data.table")
        .rb_table(data, sl, "cell_ID")
    },
    spatial_info = function(spec, data, g, n_cells) {
        p <- GiottoClass::getPolygonInfo(g, name = spec$to$name, return_giottoPolygon = TRUE)
        sv <- if (inherits(p[], "SpatVector")) p[] else GiottoClass::as.terra(p[])
        ids <- sv$poly_ID
        geo <- data.table::as.data.table(terra::geom(sv))
        geo$poly_ID <- ids[geo$geom]
        exp <- as.data.frame(data)
        bb <- function(d) {
            r <- lapply(split(d[, c("x", "y")], d$poly_ID), function(v) {
                c(xmin = min(v$x), xmax = max(v$x), ymin = min(v$y), ymax = max(v$y))
            })
            data.frame(poly_ID = names(r), do.call(rbind, r), row.names = NULL)
        }
        a <- bb(exp)
        b <- bb(as.data.frame(geo))
        m <- merge(a, b, by = "poly_ID")
        moved <- pmax(abs(m$xmin.x - m$xmin.y), abs(m$xmax.x - m$xmax.y),
            abs(m$ymin.x - m$ymin.y), abs(m$ymax.x - m$ymax.y)) > 1e-8
        list(
            .rb_row("placement", TRUE, sprintf("spatial_info/%s", spec$to$name)),
            .rb_row("ids", setequal(ids, unique(exp$poly_ID)), sprintf("%d polygons", length(ids))),
            .rb_row("geometry", !any(moved), sprintf("%d of %d polygons with a different extent",
                sum(moved), nrow(m)))
        )
    },
    feature_info = function(spec, data, g, n_cells) {
        have <- GiottoClass::list_feature_info(g)$feat_info
        fts <- intersect(unlist(spec$to$feat_type), have)
        pts <- do.call(rbind, lapply(fts, function(ft) {
            p <- GiottoClass::getFeatureInfo(g, feat_type = ft, return_giottoPoints = TRUE)
            sv <- if (inherits(p[], "SpatVector")) p[] else GiottoClass::as.terra(p[])
            xy <- terra::crds(sv)
            data.frame(feat_ID = sv$feat_ID, x = xy[, 1], y = xy[, 2])
        }))
        # transcripts per feature and their mean position
        per <- function(d) {
            f <- as.character(d$feat_ID)
            ids <- sort(unique(f))
            list(feat_ID = ids, n = as.integer(table(f)[ids]),
                mx = as.numeric(tapply(d$x, f, mean)[ids]), my = as.numeric(tapply(d$y, f, mean)[ids]))
        }
        a <- per(data)
        b <- per(pts)
        list(
            .rb_row("placement", length(fts) > 0L, paste(fts, collapse = ", ")),
            .rb_row("ids", identical(a$feat_ID, b$feat_ID) && identical(a$n, b$n),
                sprintf("%d transcripts of %d features", nrow(data), length(a$feat_ID))),
            .rb_row("values", identical(a$feat_ID, b$feat_ID) && .rb_same(c(a$mx, a$my), c(b$mx, b$my)),
                "mean x, y per feature")
        )
    },
    dim_reduction = function(spec, data, g, n_cells) {
        to <- spec$to
        m <- GiottoClass::getDimReduction(g, spat_unit = to$spat_unit, feat_type = to$feat_type,
            reduction_method = to$method, name = to$name, output = "matrix")
        list(
            .rb_row("placement", TRUE, sprintf("%s/%s/%s/%s", to$spat_unit, to$feat_type,
                to$method, to$name)),
            .rb_row("ids", setequal(rownames(m), rownames(data)), sprintf("%d cells", nrow(m))),
            .rb_row("shape", identical(dim(m), dim(data)), paste(dim(m), collapse = "x")),
            .rb_row("values", .rb_same(as.numeric(m[rownames(data), ]), as.numeric(data)), ""),
            if (!is.null(to$misc$loadings) || !is.null(to$misc$eigenvalues)) {
                d <- GiottoClass::getDimReduction(g, spat_unit = to$spat_unit, feat_type = to$feat_type,
                    reduction_method = to$method, name = to$name, output = "dimObj")
                .rb_row("loadings", identical(d@misc$loadings, to$misc$loadings) &&
                    identical(d@misc$eigenvalues, to$misc$eigenvalues), "loadings and eigenvalues")
            }
        )
    },
    spatial_network = function(spec, data, g, n_cells) {
        to <- spec$to
        e <- GiottoClass::getSpatialNetwork(g, spat_unit = to$spat_unit, name = to$name,
            output = "networkDT")
        und <- function(f, t) ifelse(f < t, paste(f, t, sep = "\r"), paste(t, f, sep = "\r"))
        ka <- und(data$from, data$to)
        kb <- und(e$from, e$to)
        list(
            .rb_row("placement", TRUE, sprintf("%s/%s", to$spat_unit, to$name)),
            .rb_row("edges", setequal(ka, kb), sprintf("%d edges", nrow(e))),
            .rb_row("distance", .rb_same(data$distance[order(ka)], e$distance[match(sort(ka), kb)]),
                "", field = "distance")
        )
    },
    nn_network = function(spec, data, g, n_cells) {
        to <- spec$to
        e <- GiottoClass::getNearestNetwork(g, spat_unit = to$spat_unit, feat_type = to$feat_type,
            nn_type = to$nn_type, name = to$name, output = "data.table")
        key <- function(f, t) paste(f, t, sep = "\r")
        und <- function(f, t) ifelse(f < t, key(f, t), key(t, f))
        same_dir <- setequal(key(e$from, e$to), key(data$from, data$to))
        same_und <- setequal(und(e$from, e$to), und(data$from, data$to))
        m <- merge(data.table::data.table(k = und(data$from, data$to), w = data$weight,
            d = data$distance), data.table::data.table(k = und(e$from, e$to), w = e$weight,
            d = e$distance), by = "k", allow.cartesian = TRUE)
        list(
            .rb_row("placement", TRUE, sprintf("%s/%s/%s/%s", to$spat_unit, to$feat_type,
                to$nn_type, to$name)),
            .rb_row("edges", same_und, sprintf("%d edges, %s", nrow(e),
                if (same_dir) "same direction" else "same as undirected pairs")),
            .rb_row("weight", .rb_same(m$w.x, m$w.y), ""),
            .rb_row("distance", identical(is.na(m$d.x), is.na(m$d.y)) &&
                .rb_same(m$d.x[!is.na(m$d.x)], m$d.y[!is.na(m$d.y)]),
                sprintf("%d NA in source, %d NA in giotto", sum(is.na(m$d.x)), sum(is.na(m$d.y))),
                field = "distance")
        )
    }
)

# table parts: ids, columns and per-column values aligned on the id column
.rb_table <- function(exp, got, id) {
    exp <- as.data.frame(exp)
    got <- as.data.frame(got)
    ids_ok <- setequal(exp[[id]], got[[id]])
    miss <- setdiff(names(exp), names(got))
    extra <- setdiff(names(got), names(exp))
    got <- got[match(exp[[id]], got[[id]]), , drop = FALSE]
    shared <- setdiff(intersect(names(exp), names(got)), id)
    bad <- shared[!vapply(shared, function(cl) {
        a <- exp[[cl]]
        b <- got[[cl]]
        if (is.numeric(a) && is.numeric(b)) .rb_same(a, b) else identical(as.character(a), as.character(b))
    }, logical(1))]
    list(
        .rb_row("ids", ids_ok, sprintf("%d rows", nrow(exp))),
        .rb_row("columns", !length(miss), if (length(miss) || length(extra)) {
            lst <- function(v) if (length(v)) paste(v, collapse = ",") else "-"
            sprintf("missing: %s; extra: %s", lst(miss), lst(extra))
        } else sprintf("%d columns", length(shared))),
        .rb_row("values", !length(bad), if (length(bad)) paste("differ:", paste(bad, collapse = ",")) else "")
    )
}
