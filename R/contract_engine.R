# The mapping engine: bind params on a profile, read source nodes through the
# source package's accessors, apply ops/hooks, build Giotto subobjects with
# the GiottoClass constructors and assemble them with setGiotto(). Disk
# backing comes from the giotto object's GiottoDisk source, so the same parts
# list builds an in-memory or a disk-backed object.

# listing verbs -> profile paths, per source ecosystem
.LIST_PATHS <- list(
    seurat = c(
        Assays = "assays", Layers = "assays/{on}/layers",
        Reductions = "reductions", Graphs = "graphs", Neighbors = "neighbors",
        Images = "images", Boundaries = "images/{on}/boundaries",
        children = "{on}", columns = "{on}"
    ),
    giotto = c(
        SpatUnits = "expression", FeatTypes = "expression/{on}", Expressions = "expression/{on}",
        SpatLocs = "spatial_locs/{on}", Polygons = "spatial_info",
        DimMethods = "dimension_reduction/cells/{on}", NNNames = "nn_network/{on}",
        children = "{on}", columns = "{on}"
    )
)

# assembly order: data before the parts that index into it
.PART_ORDER <- c("expression", "cell_metadata", "feat_metadata", "spatial_locs",
    "spatial_info", "dim_reduction", "nn_network")

#' @title Convert an object with a mapping
#' @name convertWithMapping
#' @description Convert an external object into a giotto object using only a
#' mapping file and the engine: parameters are bound on the object's profile,
#' each source node is read through its declared accessor, ops or a hook are
#' applied and the result is built with the GiottoClass constructor for the
#' target part. With `source`, the giotto object is backed by a GiottoDisk
#' project and disk-backable parts (expression, polygons, networks) are
#' written there.
#' @param x object to convert
#' @param mapping a `giottoMapping`, mapping file or bundled mapping id
#' @param profile optional `giottoProfile` of `x`
#' @param source `NULL` for an in-memory giotto object, or a directory for a
#' GiottoDisk project (requires GiottoDisk)
#' @param validate logical. Run [validateMapping()] before (structural) and
#' after (semantic) conversion; FAIL stops, WARN warns.
#' @param verbose logical
#' @returns giotto object. The conversion (mapping id, version, fingerprint,
#' params and declared losses) is recorded in the object history.
#' @export
convertWithMapping <- function(x, mapping, profile = NULL, source = NULL,
    validate = TRUE, verbose = FALSE) {
    mapping <- if (inherits(mapping, "giottoMapping")) mapping else readMapping(mapping)
    profile <- profile %null% profileObject(x)
    pkg <- mapping$mapping$source$package
    package_check(pkg)

    params <- .bind_params(mapping, profile)
    # fingerprints describe the mapping as written, not its expansion
    fp <- .mapping_fingerprint(mapping, profile, params$values)
    mapping <- .expand_for_each(mapping, profile, params)
    if (isTRUE(validate)) {
        pre <- validateMapping(mapping, profile, params = params, checks = "structural")
        .stop_on_fail(pre, "structural validation")
    }

    ctx <- .new_ctx(mapping, params, profile)
    prepared <- .prepare_parts(x, mapping, params, ctx, verbose = verbose)
    prov <- list(
        mapping = mapping$mapping$id,
        version = mapping$mapping$version,
        status = mapping$mapping$status,
        fingerprint = fp,
        params = params$values,
        parts = names(prepared),
        targets = lapply(prepared, function(pp) pp$spec$to),
        losses = .collect_losses(mapping, ctx, params$values)
    )
    if (identical(mapping$mapping$target$contract, "seurat")) {
        if (!is.null(source)) stop(wrap_txt("`source` applies to Giotto targets only"), call. = FALSE)
        g <- .record_seurat_conversion(.assemble_seurat(prepared), prov)
    } else {
        contract <- giottoContract()
        parts <- lapply(prepared, function(pp) .build_part(pp$spec, pp$data, contract, mapping))
        g <- .assemble(parts, source = source)
        g <- .record_conversion(g, mapping_id = prov$mapping, mapping_version = prov$version,
            prov = prov)
    }

    if (isTRUE(validate)) {
        post <- validateMapping(mapping, profile, params = params, result = g,
            checks = "semantic")
        .stop_on_fail(post, "semantic validation")
    }
    g
}

#' @title Conversion record of a mapped giotto object
#' @name mappingProvenance
#' @description Retrieve what [convertWithMapping()] recorded in the giotto
#' object history: mapping id, version, status, fingerprint, bound params,
#' converted parts and losses.
#' @param gobject giotto object, or a Seurat object produced by a mapping
#' @returns list, or `NULL` if the object was not produced by a mapping
#' @export
mappingProvenance <- function(gobject) {
    if (!inherits(gobject, "giotto")) {
        misc <- tryCatch(methods::slot(gobject, "misc"), error = function(e) NULL)
        return(misc$giottoConnect)
    }
    hist <- gobject@parameters
    for (h in rev(hist)) {
        a <- attr(h, "attachments")
        if (!is.null(a$giottoConnect)) return(a$giottoConnect)
    }
    NULL
}

# recorded through GiottoClass so it shows up in objHistory(); the recorded
# call is this function's (small) argument list
.record_conversion <- function(gobject, mapping_id, mapping_version, prov) {
    GiottoClass::update_giotto_params(gobject,
        description = "_convertWithMapping", toplevel = 1L,
        attachments = list(giottoConnect = prov)
    )
}


# Read, hook and transform every bound part. Returns a list (in assembly
# order) of list(spec, data): spec is the part with its final, bound `to`
# keys; data is exactly what the constructor receives.
.prepare_parts <- function(x, mapping, params, ctx, verbose = FALSE) {
    pkg <- mapping$mapping$source$package
    ctx$x <- x
    paths <- ctx$profile$nodes$path
    out <- list()
    for (p in mapping$parts) {
        node <- mapping$source_nodes[[p$from]]
        present <- .part_bound(p, node, params) && .pf_bind(node$path, params$values) %in% paths
        if (!present) {
            if (isTRUE(p$optional)) next
            stop(wrap_txt(sprintf("Part '%s': source '%s' is not bound or not present", p$id,
                .pf_bind(node$path, params$values))), call. = FALSE)
        }
        ctx$part <- p$id
        vmsg(.v = verbose, sprintf("part %s <- %s", p$id, .pf_bind(node$path, params$values)))
        p$to <- lapply(p$to, function(v) if (is.character(v)) .pf_bind(v, params$values) else v)
        ctx$part_spec <- p
        data <- .read_node(x, node, params$values, pkg)
        if (!is.null(p$hook)) {
            data <- .hook_get(p$hook)(data, params$values, ctx)
            if (inherits(data, "giottoHookResult")) {
                p$to <- utils::modifyList(p$to, data$to)
                data <- data$data
            }
        }
        data <- .apply_ops(data, p$ops, params$values, p$id)
        for (l in p$lossy) ctx$add_loss(l$kind, l$field %null% l$path, l$note, declared = TRUE)
        out[[p$id]] <- list(spec = p, data = data)
    }
    types <- vapply(out, function(pp) pp$spec$part, character(1))
    out[order(match(types, .PART_ORDER))]
}


# for_each ####

# Expand template parts (`for_each`) into one concrete, optional part per
# matching child of the listing path, e.g. one dim_reduction part per
# reduction. Each copy gets its own source node with the iteration value
# bound, so the expanded mapping is an ordinary mapping.
.expand_for_each <- function(mapping, profile, params) {
    if (!any(vapply(mapping$parts, function(p) !is.null(p$for_each), logical(1)))) {
        return(mapping)
    }
    lists <- .LIST_PATHS[[mapping$mapping$source$ecosystem]]
    nodes <- profile$nodes
    out <- list()
    for (p in mapping$parts) {
        fe <- p$for_each
        if (is.null(fe)) {
            out[[length(out) + 1L]] <- p
            next
        }
        on <- .pf_bind(fe$on %null% "", params$values)
        if (grepl("{", on, fixed = TRUE)) next
        parent <- gsub("{on}", on, lists[[fe$list]], fixed = TRUE)
        kids <- nodes[.pf_parent(nodes$path) == parent & nodes$path != parent &
            !endsWith(nodes$path, "/*"), , drop = FALSE]
        if (!nrow(kids)) next
        name <- sub(".*/", "", kids$path)
        keep <- if (is.null(fe$where)) rep(TRUE, nrow(kids)) else
            .param_where(fe$where, kids, name, nodes)
        vals <- name[keep]
        if (!is.null(fe$exclude)) vals <- vals[!grepl(fe$exclude, vals)]
        for (v in vals) {
            sub_v <- function(x) gsub(sprintf("{%s}", fe$as), v, x, fixed = TRUE)
            q <- p
            q$for_each <- NULL
            q <- rapply(q, sub_v, classes = "character", how = "replace")
            q$optional <- TRUE
            node_id <- sprintf("%s__%s", p$from, v)
            mapping$source_nodes[[node_id]] <- rapply(mapping$source_nodes[[p$from]], sub_v,
                classes = "character", how = "replace")
            q$from <- node_id
            out[[length(out) + 1L]] <- q
        }
    }
    mapping$parts <- out
    mapping
}


# params ####

# Bind each param on the profile, in declaration order. Returns list(values,
# unbound) with values a named list of strings.
.bind_params <- function(mapping, profile) {
    eco <- mapping$mapping$source$ecosystem
    lists <- .LIST_PATHS[[eco]]
    if (is.null(lists)) stop(wrap_txt(sprintf("No listing paths for ecosystem '%s'", eco)))
    nodes <- profile$nodes
    values <- list()
    unbound <- character()
    for (nm in names(mapping$params)) {
        pr <- mapping$params[[nm]]
        on <- if (is.null(pr$on)) NULL else .pf_bind(pr$on, values)
        if (!is.null(pr$value)) pr$value <- .pf_bind(pr$value, values)
        if (!is.null(on) && grepl("{", on, fixed = TRUE)) {
            unbound <- c(unbound, nm)
            next
        }
        parent <- gsub("{on}", on %null% "", lists[[pr$list]], fixed = TRUE)
        kids <- nodes[.pf_parent(nodes$path) == parent & nodes$path != parent &
            !endsWith(nodes$path, "/*"), , drop = FALSE]
        val <- .param_pick(pr, kids, nodes)
        if (is.na(val)) unbound <- c(unbound, nm) else values[[nm]] <- val
    }
    list(values = values, unbound = unbound)
}

.param_pick <- function(pr, kids, nodes) {
    if (!nrow(kids)) return(NA_character_)
    name <- sub(".*/", "", kids$path)
    ok <- switch(pr$rule,
        first = seq_len(nrow(kids)) == 1L,
        name_is = name == pr$value,
        value_of = {
            v <- nodes$value[nodes$path == .pf_bind(pr$path, list())]
            hit <- length(v) == 1L && !is.na(v) & name == v
            if (!is.null(pr$where)) hit & .param_where(pr$where, kids, name, nodes) else hit
        },
        first_of_class = kids$class == pr$class,
        first_where = .param_where(pr$where, kids, name, nodes)
    )
    if (!any(ok)) NA_character_ else name[which(ok)[1L]]
}

.param_where <- function(w, kids, name, nodes) {
    has <- function(h, roles) {
        hs <- strsplit(ifelse(is.na(h), "", h), ",", fixed = TRUE)
        vapply(hs, function(x) all(roles %in% x), logical(1))
    }
    ok <- rep(TRUE, nrow(kids))
    if (!is.null(w$class)) ok <- ok & kids$class == w$class
    if (!is.null(w$kind)) ok <- ok & kids$kind == w$kind
    if (!is.null(w$role)) ok <- ok & has(kids$role_hints, unlist(w$role))
    if (!is.null(w$has_child_role)) {
        # data held below the node: matrices only, not scalar slots
        ok <- ok & vapply(kids$path, function(p) {
            sub <- nodes[startsWith(nodes$path, paste0(p, "/")) &
                nodes$kind %in% c("sparse_matrix", "dense_matrix"), , drop = FALSE]
            any(has(sub$role_hints, unlist(w$has_child_role)))
        }, logical(1))
    }
    if (!is.null(w$name_matches) || !is.null(w$or_role)) {
        by_name <- if (is.null(w$name_matches)) FALSE else grepl(w$name_matches, name)
        by_role <- if (is.null(w$or_role)) FALSE else has(kids$role_hints, unlist(w$or_role))
        ok <- ok & (by_name | by_role)
    }
    ok
}

# a part is bound when its source path, accessor and op arguments have no
# unbound placeholders left
.part_bound <- function(p, node, params) {
    txt <- c(node$path, node$get$on, unlist(node$get$args), unlist(p$ops),
        names(unlist(lapply(p$ops, `[[`, "map"))))
    txt <- txt[is.character(txt)]
    !any(grepl("{", .pf_bind(txt, params$values), fixed = TRUE))
}


# source nodes ####

# walk a concrete "/" path through S4 slots, list elements and environments
.resolve_path <- function(x, path) {
    if (!nzchar(path)) return(x)
    for (seg in strsplit(path, "/", fixed = TRUE)[[1L]]) {
        seg <- gsub("%2F", "/", seg, fixed = TRUE)
        x <- if (isS4(x)) {
            attr(x, seg, exact = TRUE)
        } else if (is.environment(x)) {
            get0(seg, envir = x, inherits = FALSE)
        } else {
            x[[seg]]
        }
        if (is.null(x)) {
            stop(wrap_txt(sprintf("Path '%s' not found (at '%s')", path, seg)), call. = FALSE)
        }
    }
    x
}

.read_node <- function(x, node, params, pkg) {
    get <- node$get
    # "path": the node at `on`, read structurally. Only patches use it; a
    # verified mapping reads through the source package's accessors.
    if (identical(get$fn, "path")) return(.resolve_path(x, .pf_bind(get$on, params)))
    fn <- if (identical(get$fn, "[[")) {
        base::`[[`
    } else {
        getExportedValue(pkg, get$fn)
    }
    target <- if (is.null(get$on)) x else .resolve_path(x, .pf_bind(get$on, params))
    args <- lapply(get$args %null% list(), function(a) {
        if (is.character(a)) .pf_bind(a, params) else a
    })
    do.call(fn, c(list(target), args))
}


# ops ####

.apply_ops <- function(data, ops, params, part_id) {
    for (o in ops) {
        o <- rapply(o, function(v) .pf_bind(v, params), classes = "character", how = "replace")
        if (!is.null(o$map)) names(o$map) <- .pf_bind(names(o$map), params)
        data <- tryCatch(.apply_op(data, o), error = function(e) {
            stop(sprintf("Part '%s', op '%s': %s", part_id, o$op, conditionMessage(e)),
                call. = FALSE)
        })
    }
    data
}

.apply_op <- function(data, o) {
    tab <- is.data.frame(data)
    switch(o$op,
        select_cols = {
            miss <- setdiff(o$cols, colnames(data))
            if (length(miss)) stop("missing column(s): ", paste(miss, collapse = ", "))
            data[, unlist(o$cols), drop = FALSE]
        },
        exclude_cols = data[, setdiff(colnames(data), unlist(o$cols)), drop = FALSE],
        rename = {
            old <- names(o$map)
            miss <- setdiff(old, colnames(data))
            if (length(miss)) stop("missing column(s): ", paste(miss, collapse = ", "))
            cn <- colnames(data)
            cn[match(old, cn)] <- unlist(o$map)
            colnames(data) <- cn
            data
        },
        transpose = t(data),
        scale = {
            if (!is.null(o$cols) && tab) {
                for (cl in unlist(o$cols)) data[[cl]] <- data[[cl]] * o$factor
                data
            } else {
                data * o$factor
            }
        },
        cast = .cast(data, o$to),
        set_colnames = {
            if (length(o$cols) != ncol(data)) stop("set_colnames needs one name per column")
            colnames(data) <- unlist(o$cols)
            data
        },
        rownames_to_col = {
            if (!tab) stop("rownames_to_col needs a table")
            rn <- rownames(data)
            # an existing column of that name is reused when it holds the row
            # names, otherwise kept as "<col>.source"
            if (o$col %in% names(data)) {
                if (identical(as.character(data[[o$col]]), rn)) {
                    data <- data[, setdiff(names(data), o$col), drop = FALSE]
                } else {
                    names(data)[names(data) == o$col] <- paste0(o$col, ".source")
                }
            }
            out <- cbind(stats::setNames(data.frame(rn, stringsAsFactors = FALSE), o$col),
                data)
            rownames(out) <- NULL
            out
        },
        set = {
            if (!tab) stop("set needs a table")
            data[[o$field]] <- o$value
            data
        },
        stop("unknown op")
    )
}

.cast <- function(data, to) {
    switch(to,
        table = as.data.frame(as.matrix(data)),
        sparse = if (methods::is(data, "Matrix") || is.matrix(data)) {
            methods::as(data, "CsparseMatrix")
        } else {
            methods::as(data, "dgCMatrix")
        },
        integer = {
            v <- if (methods::is(data, "sparseMatrix")) data@x else as.vector(data)
            if (any(v != round(v), na.rm = TRUE)) stop("values are not integer-valued")
            if (methods::is(data, "sparseMatrix")) data else {
                storage.mode(data) <- "integer"
                data
            }
        },
        double = {
            if (methods::is(data, "sparseMatrix")) data else {
                storage.mode(data) <- "double"
                data
            }
        },
        character = {
            storage.mode(data) <- "character"
            data
        }
    )
}


# build and assemble ####

.build_part <- function(p, data, contract, mapping) {
    cp <- contract$parts[[p$part]]
    if (is.data.frame(data)) data <- data.table::as.data.table(data)
    if (p$part == "dim_reduction") data <- as.matrix(data)
    ctor <- getExportedValue("GiottoClass", cp$constructor)
    args <- stats::setNames(list(data), cp$data_arg)
    to <- p$to
    for (k in intersect(names(to), cp$args)) args[[k]] <- to[[k]]
    if ("misc" %in% cp$args) {
        args$misc <- list(giottoConnect = list(mapping = mapping$mapping$id,
            version = mapping$mapping$version, part = p$id))
    }
    if ("verbose" %in% cp$args) args$verbose <- FALSE
    if (p$part == "spatial_info") args$calc_centroids <- FALSE
    suppressMessages(do.call(ctor, args))
}

.assemble <- function(parts, source = NULL) {
    g <- if (is.null(source)) {
        GiottoClass::giotto()
    } else {
        package_check("GiottoDisk")
        src <- GiottoDisk::sourceCreate(source, type = "gDirSource")
        GiottoClass::giotto(source = src)
    }
    g <- GiottoClass::setGiotto(g, unname(parts), initialize = FALSE, verbose = FALSE)
    methods::initialize(g)
}


# context, losses, fingerprints ####

.new_ctx <- function(mapping, params, profile) {
    ctx <- new.env(parent = emptyenv())
    ctx$losses <- list()
    ctx$x <- NULL
    ctx$profile <- profile
    ctx$part <- NA_character_
    ctx$add_loss <- function(kind, field, note, declared = FALSE) {
        ctx$losses[[length(ctx$losses) + 1L]] <- data.frame(
            part_id = ctx$part, kind = kind, field = field %null% NA_character_,
            note = note, source = if (declared) "declared" else "recorded",
            stringsAsFactors = FALSE
        )
    }
    ctx
}

.collect_losses <- function(mapping, ctx, params) {
    top <- lapply(mapping$lossy, function(l) data.frame(
        part_id = NA_character_, kind = l$kind,
        field = .pf_bind(l$path %null% l$field %null% NA_character_, params),
        note = l$note, source = "declared", stringsAsFactors = FALSE
    ))
    out <- do.call(rbind, c(top, ctx$losses))
    if (is.null(out)) {
        out <- data.frame(part_id = character(), kind = character(), field = character(),
            note = character(), source = character())
    }
    out
}

# templated paths a mapping depends on: the source nodes its parts read (not
# those only alternatives use) and the applies_to predicates
.mapping_paths <- function(mapping) {
    fixed <- Filter(function(p) is.null(p$for_each), mapping$parts)
    used <- unique(vapply(fixed, `[[`, character(1), "from"))
    unique(c(
        vapply(mapping$source_nodes[used], `[[`, character(1), "path"),
        vapply(mapping$mapping$source$applies_to$predicates %null% list(),
            `[[`, character(1), "path")
    ))
}

.mapping_fingerprint <- function(mapping, profile, params) {
    fingerprint(profile, paths = .mapping_paths(mapping), params = params)
}

.stop_on_fail <- function(v, stage) {
    fails <- v[v$level == "FAIL", , drop = FALSE]
    warns <- v[v$level == "WARN", , drop = FALSE]
    if (nrow(warns)) {
        warning(sprintf("%s:\n%s", stage, paste0("  - ", .fmt_checks(warns), collapse = "\n")),
            call. = FALSE)
    }
    if (nrow(fails)) {
        stop(sprintf("%s failed:\n%s", stage, paste0("  - ", .fmt_checks(fails), collapse = "\n")),
            call. = FALSE)
    }
    invisible(v)
}

.fmt_checks <- function(v) {
    sprintf("[%s]%s %s", v$check, ifelse(is.na(v$part_id), "", paste0(" ", v$part_id)), v$message)
}
