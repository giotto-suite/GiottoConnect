#' @title Validate a mapping against an object
#' @name validateMapping
#' @description Check a mapping against a profiled object (structural checks)
#' and, when `result` is given, check the converted giotto object (semantic
#' checks). No conversion is run here.
#'
#' Structural: schema, parameter binding, source paths, target keys under the
#' Giotto contract, hooks, and data in the source that is neither mapped nor
#' declared lossy. Semantic: raw expression is integer-valued and
#' non-negative, cell IDs agree across parts, dimensions match the profile,
#' coordinate ranges are plausible for microns, reductions and networks index
#' the cells.
#' @param mapping a `giottoMapping`, mapping file or bundled mapping id
#' @param profile `giottoProfile` of the source object
#' @param params optional bound params (as returned internally by the engine)
#' @param result optional giotto object produced with this mapping
#' @param checks `"all"`, `"structural"` or `"semantic"` (the latter requires
#' `result`)
#' @returns a `giottoValidation` data.frame with columns `check`, `part_id`,
#' `level` (OK, WARN, FAIL) and `message`
#' @export
validateMapping <- function(mapping, profile, params = NULL, result = NULL,
    checks = c("all", "structural", "semantic")) {
    checks <- match.arg(checks)
    mapping <- if (inherits(mapping, "giottoMapping")) mapping else readMapping(mapping)
    checkmate::assert_class(profile, "giottoProfile")
    params <- params %null% .bind_params(mapping, profile)
    out <- .vm_collector()

    if (checks != "semantic") {
        for (pr in .mapping_schema_check(unclass(mapping))) out$add("schema", NA, "FAIL", pr)
        .vm_params(mapping, params, out)
        .vm_nodes(mapping, profile, params, out)
        .vm_targets(mapping, out)
        .vm_hooks(mapping, out)
        .vm_unmapped(mapping, profile, params, out)
    }
    if (checks != "structural" && !is.null(result)) .vm_semantic(mapping, profile, result, out)

    res <- out$table()
    class(res) <- c("giottoValidation", "data.frame")
    res
}

#' @export
print.giottoValidation <- function(x, ...) {
    n <- table(factor(x$level, levels = c("OK", "WARN", "FAIL")))
    cat(sprintf("<giottoValidation> %d OK, %d WARN, %d FAIL\n", n[["OK"]], n[["WARN"]], n[["FAIL"]]))
    show <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
    show$part_id[is.na(show$part_id)] <- ""
    print(show[, c("level", "check", "part_id", "message")], row.names = FALSE, right = FALSE)
    invisible(x)
}

.vm_collector <- function() {
    rows <- list()
    list(
        add = function(check, part_id, level, message) {
            rows[[length(rows) + 1L]] <<- data.frame(check = check,
                part_id = part_id %null% NA_character_, level = level,
                message = message, stringsAsFactors = FALSE)
        },
        table = function() {
            if (!length(rows)) {
                return(data.frame(check = character(), part_id = character(),
                    level = character(), message = character()))
            }
            do.call(rbind, rows)
        }
    )
}

# which params each part depends on (through its source node)
.vm_part_params <- function(mapping, part) {
    node <- mapping$source_nodes[[part$from]]
    if (is.null(node)) return(character())
    txt <- paste(node$path, node$get$on %null% "", paste(unlist(node$get$args), collapse = " "))
    m <- regmatches(txt, gregexpr("\\{[^}]+\\}", txt))[[1L]]
    unique(gsub("[{}]", "", m))
}

.vm_params <- function(mapping, params, out) {
    for (nm in names(mapping$params)) {
        pr <- mapping$params[[nm]]
        if (nm %in% names(params$values)) {
            out$add("param", NA, "OK", sprintf("%s = '%s'", nm, params$values[[nm]]))
        } else {
            needed <- any(vapply(mapping$parts, function(p) {
                !isTRUE(p$optional) && nm %in% .vm_part_params(mapping, p)
            }, logical(1)))
            lvl <- if (needed) "FAIL" else if (isTRUE(pr$optional)) "WARN" else "FAIL"
            out$add("param", NA, lvl, sprintf("%s could not be bound (%s %s)", nm,
                pr$list, pr$rule))
        }
    }
}

.vm_nodes <- function(mapping, profile, params, out) {
    for (p in mapping$parts) {
        node <- mapping$source_nodes[[p$from]]
        if (is.null(node)) {
            out$add("source", p$id, "FAIL", sprintf("unknown source node '%s'", p$from))
            next
        }
        path <- .pf_bind(node$path, params$values)
        lvl_missing <- if (isTRUE(p$optional)) "WARN" else "FAIL"
        if (grepl("{", path, fixed = TRUE)) {
            out$add("source", p$id, lvl_missing, sprintf("'%s' has unbound params", node$path))
        } else if (!path %in% profile$nodes$path) {
            out$add("source", p$id, lvl_missing, sprintf("path '%s' not found", path))
        } else {
            out$add("source", p$id, "OK", path)
        }
    }
}

.vm_targets <- function(mapping, out) {
    ct <- giottoContract()
    for (p in mapping$parts) {
        cp <- ct$parts[[p$part]]
        keys <- names(p$to)
        unknown <- setdiff(keys, c(cp$nesting, cp$args))
        need <- setdiff(intersect(cp$nesting, c("spat_unit", "feat_type", "name", "nn_type")),
            c(keys, if (p$part == "spatial_info") c("spat_unit", "feat_type", "spat_info")))
        if (p$part == "spatial_info" && is.null(p$to$name)) need <- c(need, "name")
        if (length(unknown)) {
            out$add("target", p$id, "FAIL", sprintf("'%s' is not a key or argument of %s",
                paste(unknown, collapse = "', '"), cp$constructor))
        } else if (length(need)) {
            out$add("target", p$id, "FAIL", sprintf("missing nesting key(s): %s",
                paste(need, collapse = ", ")))
        } else {
            out$add("target", p$id, "OK", sprintf("%s via %s", p$part, cp$constructor))
        }
    }
}

.vm_hooks <- function(mapping, out) {
    reg <- hookRegistry()
    for (p in mapping$parts) {
        if (is.null(p$hook)) next
        h <- reg[reg$hook == p$hook, , drop = FALSE]
        if (!nrow(h)) {
            out$add("hook", p$id, "FAIL", sprintf("hook '%s' is not registered", p$hook))
        } else if (!h$available) {
            out$add("hook", p$id, "FAIL", sprintf("hook '%s' needs %s", p$hook, h$packages))
        } else {
            out$add("hook", p$id, "OK", p$hook)
        }
    }
}

# numeric matrices and tables that no part reads and no loss declares
.vm_unmapped <- function(mapping, profile, params, out) {
    nd <- profile$nodes
    empty <- grepl("(^|x)0(x|$)", nd$shape)
    data_nodes <- nd$path[!empty & (
        (nd$kind %in% c("sparse_matrix", "dense_matrix") & nd$dtype %in% c("double", "integer")) |
            nd$kind == "data.frame")
    ]
    covered <- c(
        vapply(mapping$source_nodes, function(n) .pf_bind(n$path, params$values), character(1)),
        unlist(lapply(mapping$parts, function(p) lapply(p$lossy, function(l) l$path))),
        unlist(lapply(mapping$lossy, function(l) .pf_bind(l$path %null% "", params$values)))
    )
    covered <- covered[!grepl("{", covered, fixed = TRUE) & nzchar(covered)]
    inside <- function(p) any(p == covered | startsWith(p, paste0(covered, "/")) |
        startsWith(covered, paste0(p, "/")))
    # report subtree roots only
    left <- data_nodes[!vapply(data_nodes, inside, logical(1))]
    left <- left[!vapply(left, function(p) any(startsWith(p, paste0(setdiff(left, p), "/"))), logical(1))]
    for (p in left) {
        r <- nd[nd$path == p, ]
        out$add("unmapped", NA, "WARN", sprintf("%s (%s %s) is neither mapped nor declared lossy",
            p, r$class, r$shape))
    }
}


# semantic checks on the converted object ####

.vm_semantic <- function(mapping, profile, g, out) {
    G <- asNamespace("GiottoClass")
    cells_expr <- GiottoClass::spatIDs(g, spat_unit = "cell")
    feats <- GiottoClass::featIDs(g, feat_type = "rna")

    raw <- tryCatch(GiottoClass::getExpression(g, values = "raw", spat_unit = "cell",
        feat_type = "rna", output = "exprObj"), error = function(e) NULL)
    if (!is.null(raw)) {
        v <- .vm_expr_sample(raw[], cells_expr)
        ok_int <- all(v == round(v))
        ok_pos <- all(v >= 0)
        out$add("raw_counts", NA, if (ok_int && ok_pos) "OK" else "FAIL",
            sprintf("raw expression %s integer-valued, %s non-negative (sampled %d values)",
                if (ok_int) "is" else "is NOT", if (ok_pos) "is" else "is NOT", length(v)))
    }

    if (!is.na(profile$n_cells)) {
        lvl <- if (length(cells_expr) == profile$n_cells) "OK" else "WARN"
        out$add("dims", NA, lvl, sprintf("%d cells in giotto, %d in source", length(cells_expr),
            profile$n_cells))
    }
    if (!is.na(profile$n_feat)) {
        lvl <- if (length(feats) == profile$n_feat) "OK" else "WARN"
        out$add("dims", NA, lvl, sprintf("%d features in giotto, %d in source", length(feats),
            profile$n_feat))
    }

    ids <- list(
        metadata = tryCatch(GiottoClass::pDataDT(g, spat_unit = "cell", feat_type = "rna")$cell_ID,
            error = function(e) NULL),
        spatial_locs = tryCatch(GiottoClass::getSpatialLocations(g, spat_unit = "cell",
            output = "data.table")$cell_ID, error = function(e) NULL),
        polygons = tryCatch(GiottoClass::spatIDs(GiottoClass::getPolygonInfo(g,
            name = "cell", return_giottoPolygon = TRUE)), error = function(e) NULL)
    )
    for (nm in names(ids)) {
        if (is.null(ids[[nm]])) next
        ov <- mean(cells_expr %in% ids[[nm]])
        out$add("cell_ids", nm, if (ov >= 0.99) "OK" else "FAIL",
            sprintf("%.1f%% of expression cells found in %s", 100 * ov, nm))
    }

    sl <- tryCatch(GiottoClass::getSpatialLocations(g, spat_unit = "cell", output = "data.table"),
        error = function(e) NULL)
    if (!is.null(sl) && nrow(sl)) {
        rng <- max(diff(range(sl$sdimx)), diff(range(sl$sdimy)))
        lvl <- if (rng < 50 || rng > 1e6) "WARN" else "OK"
        msg <- if (rng < 50) "looks like mm, not micron" else if (rng > 1e6) "implausibly large for micron" else "plausible for micron"
        out$add("coord_units", "locs", lvl, sprintf("coordinate range %.3g: %s", rng, msg))
    }

    dr <- GiottoClass::list_dim_reductions(g)
    for (i in seq_len(NROW(dr))) {
        m <- GiottoClass::getDimReduction(g, spat_unit = dr$spat_unit[i], feat_type = dr$feat_type[i],
            reduction_method = dr$dim_type[i], name = dr$name[i], output = "matrix")
        lvl <- if (nrow(m) == length(cells_expr) && all(rownames(m) %in% cells_expr)) "OK" else "FAIL"
        out$add("dim_reduction", dr$name[i], lvl, sprintf("%d rows for %d cells", nrow(m),
            length(cells_expr)))
    }
    nn <- GiottoClass::list_nearest_networks(g)
    for (i in seq_len(NROW(nn))) {
        e <- GiottoClass::getNearestNetwork(g, spat_unit = nn$spat_unit[i], feat_type = nn$feat_type[i],
            nn_type = nn$nn_type[i], name = nn$name[i], output = "data.table")
        nodes <- unique(c(e$from, e$to))
        lvl <- if (all(nodes %in% cells_expr)) "OK" else "FAIL"
        out$add("nn_network", nn$name[i], lvl, sprintf("%d edges over %d nodes, all cells: %s",
            nrow(e), length(nodes), lvl == "OK"))
    }
    invisible(G)
}

# values from the first cells, for in-memory or disk matrices; GiottoDisk
# materialises at most 100 along one axis without an explicit override
.vm_expr_sample <- function(m, cells, n_cells = 99L) {
    use <- utils::head(cells, n_cells)
    sub <- m[, use, drop = FALSE]
    if (!methods::is(sub, "Matrix") && !is.matrix(sub)) sub <- methods::as(sub, "dgCMatrix")
    v <- if (methods::is(sub, "sparseMatrix")) sub@x else as.vector(sub)
    v
}
