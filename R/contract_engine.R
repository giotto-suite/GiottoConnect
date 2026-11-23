# The mapping engine: bind params on a profile, read source nodes through the
# source package's accessors, apply ops/hooks, build Giotto subobjects with
# the GiottoClass constructors and assemble them with setGiotto(). Disk
# backing comes from the giotto object's GiottoDisk source, so the same parts
# list builds an in-memory or a disk-backed object.

# assembly order: data before the parts that index into it
.PART_ORDER <- c("expression", "cell_metadata", "feat_metadata", "spatial_locs",
    "spatial_info", "feature_info", "dim_reduction", "nn_network", "spatial_network")

#' @title Convert an object with a mapping
#' @name convertWithMapping
#' @description Convert an external object into a giotto object using only a
#' mapping file and the engine: parameters are bound on the object's profile,
#' each source node is read through its declared accessor, ops or a hook are
#' applied and the result is built with the GiottoClass constructor for the
#' target part. A part whose source is absent switches to the first of its
#' declared alternatives that is present; otherwise an optional part is
#' skipped. With `source`, the giotto object is backed by a GiottoDisk
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
#' @param params optional named list of param values that replace the
#' mapping's rules, e.g. `list(spat_unit = "nucleus", feat_type = "protein")`
#' to export a chosen spatial unit and feature type; recorded in the
#' provenance and reused by [readbackCheck()]
#' @param output file to write, for file targets (e.g. `.h5ad`)
#' @returns giotto object. The conversion (mapping id, version, fingerprint,
#' params, losses, alternatives used and absent optional parts) is recorded
#' in the object history.
#' @export
convertWithMapping <- function(x, mapping, profile = NULL, source = NULL,
    validate = TRUE, verbose = FALSE, params = NULL, output = NULL) {
    pinned <- params
    mapping <- if (inherits(mapping, "giottoMapping")) mapping else readMapping(mapping)
    adapter <- .adapter_for(mapping)
    target <- .target_for(mapping)
    package_check(adapter@package)
    profile <- profile %null% adapterProfile(adapter, x)

    params <- .bind_params(mapping, profile, pinned)
    # fingerprints describe the mapping as written, not its expansion
    fp <- .mapping_fingerprint(mapping, profile, params$values)
    res <- .resolve_parts(mapping, profile, pinned)
    mapping <- res$mapping
    params <- res$params
    mapping <- .expand_for_each(mapping, profile, params)
    if (isTRUE(validate)) {
        pre <- validateMapping(mapping, profile, params = params, checks = "structural")
        .stop_on_fail(pre, "structural validation")
    }

    ctx <- .new_ctx(mapping, params, profile)
    prepared <- .prepare_parts(x, mapping, params, ctx, verbose = verbose,
        stream = !is.null(source))
    prov <- list(
        mapping = mapping$mapping$id,
        version = mapping$mapping$version,
        status = mapping$mapping$status,
        fingerprint = fp,
        params = params$values,
        parts = names(prepared),
        targets = lapply(prepared, function(pp) pp$spec$to),
        losses = .collect_losses(mapping, ctx, params$values),
        alternatives = res$switched,
        absent = ctx$absent,
        pinned = pinned
    )
    g <- targetBuild(target, prepared, mapping, prov, source = source, output = output)

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
#' @param gobject giotto object, a Seurat object produced by a mapping, or the
#' path of an `.h5ad` file written by [giottoToAnnData()]
#' @returns list, or `NULL` if the object was not produced by a mapping
#' @export
mappingProvenance <- function(gobject) {
    if (.is_h5ad_path(gobject)) return(.anndata_provenance(gobject))
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
# Bind params and, for a part whose source is absent, switch to the first
# declared alternative that binds. repairMapping() can stamp the same switch
# into a new mapping version.
.resolve_parts <- function(mapping, profile, forced = NULL) {
    params <- .bind_params(mapping, profile, forced)
    switched <- character()
    for (i in seq_along(mapping$parts)) {
        p <- mapping$parts[[i]]
        if (!length(p$alternatives) || !is.null(p$for_each)) next
        node <- mapping$source_nodes[[p$from]]
        if (.part_bound(p, node, params) &&
            .pf_bind(node$path, params$values) %in% profile$nodes$path) next
        alt <- .rule_alternative(mapping, p, profile)
        if (is.null(alt)) next
        mapping <- .apply_alternative(mapping, i, alt)
        params <- .bind_params(mapping, profile, forced)
        switched <- c(switched, sprintf("%s: %s", p$id, alt$detail))
    }
    list(mapping = mapping, params = params, switched = switched)
}

# stream = TRUE: expression parts without ops or hook may come back as a
# GiottoDisk input (see adapterStream) instead of an R matrix
.prepare_parts <- function(x, mapping, params, ctx, verbose = FALSE, stream = FALSE) {
    adapter <- .adapter_for(mapping)
    ctx$x <- x
    paths <- ctx$profile$nodes$path
    out <- list()
    for (p in mapping$parts) {
        node <- .drop_optional_args(mapping$source_nodes[[p$from]], mapping, params)
        present <- .part_bound(p, node, params) && .pf_bind(node$path, params$values) %in% paths
        if (!present) {
            if (isTRUE(p$optional)) {
                ctx$absent <- c(ctx$absent, p$id)
                next
            }
            stop(wrap_txt(sprintf("Part '%s': source '%s' is not bound or not present", p$id,
                .pf_bind(node$path, params$values))), call. = FALSE)
        }
        ctx$part <- p$id
        vmsg(.v = verbose, sprintf("part %s <- %s", p$id, .pf_bind(node$path, params$values)))
        p$to <- lapply(p$to, function(v) if (is.character(v)) .pf_bind(v, params$values) else v)
        ctx$part_spec <- p
        streamable <- stream && p$part == "expression" && is.null(p$hook) && !length(p$ops)
        data <- if (streamable) adapterStream(adapter, x, node, params$values)
        if (!is.null(data)) {
            for (l in p$lossy) ctx$add_loss(l$kind, l$field %null% l$path, l$note, declared = TRUE)
            out[[p$id]] <- list(spec = p, data = data)
            next
        }
        data <- adapterRead(adapter, x, node, params$values)
        if (!is.null(p$hook)) {
            data <- .hook_get(p$hook)(data, params$values, ctx)
            if (inherits(data, "giottoHookResult")) {
                p$to <- utils::modifyList(p$to, data$to)
                data <- data$data
            }
        }
        data <- .apply_ops(data, p$ops, params$values, p$id)
        # an optional part with nothing left after its ops (e.g. no control
        # features to select) is absent
        if (isTRUE(p$optional) && !is.null(dim(data)) && nrow(data) == 0L) {
            ctx$absent <- c(ctx$absent, p$id)
            next
        }
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
    lists <- adapterListPaths(.adapter_for(mapping))
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
        kids <- .list_kids(nodes, lists, fe$list, on)
        if (!nrow(kids)) next
        name <- sub(".*/", "", kids$path)
        keep <- if (is.null(fe$where)) rep(TRUE, nrow(kids)) else
            .param_where(fe$where, kids, name, nodes)
        if (!is.null(fe$exclude)) keep <- keep & !grepl(.pf_bind(fe$exclude, params$values), name)
        for (j in which(keep)) {
            v <- name[j]
            kid <- kids$path[j]
            label <- fe$rename[[v]] %null% (if (is.null(fe$strip)) v else sub(fe$strip, "", v))
            sub_v <- function(x) {
                # {<as>_path}: the item's profile path; {<as>_parent}: its parent's name
                x <- gsub(sprintf("{%s_path}", fe$as), kid, x, fixed = TRUE)
                x <- gsub(sprintf("{%s_parent}", fe$as), sub(".*/", "", .pf_parent(kid)), x, fixed = TRUE)
                x <- gsub(sprintf("{%s_label}", fe$as), label, x, fixed = TRUE)
                gsub(sprintf("{%s}", fe$as), v, x, fixed = TRUE)
            }
            q <- p
            q$for_each <- NULL
            q <- rapply(q, sub_v, classes = "character", how = "replace")
            q$optional <- TRUE
            # a name can repeat under different parents (a `*` listing level)
            node_id <- sprintf("%s__%s", p$from, if (sum(name[keep] == v) > 1L) gsub("/", "__", kid) else v)
            mapping$source_nodes[[node_id]] <- rapply(mapping$source_nodes[[p$from]], sub_v,
                classes = "character", how = "replace")
            q$from <- node_id
            out[[length(out) + 1L]] <- q
        }
    }
    mapping$parts <- out
    mapping
}


# children of a listing verb's path(s); "a|b" lists the children of both
.list_kids <- function(nodes, lists, verb, on) {
    parents <- strsplit(gsub("{on}", on %null% "", lists[[verb]], fixed = TRUE), "|", fixed = TRUE)[[1L]]
    # a `*` segment stands for every child at that level
    glob <- grepl("(^|/)\\*(/|$)", parents)
    if (any(glob)) {
        pat <- paste0("^", gsub("\\*", "[^/]+", gsub("([.+?^${}()|\\[\\]\\\\])", "\\\\\\1", parents[glob])), "$")
        hit <- unique(unlist(lapply(pat, function(p) grep(p, nodes$path, value = TRUE))))
        parents <- c(parents[!glob], hit)
    }
    nodes[.pf_parent(nodes$path) %in% parents & !nodes$path %in% parents &
        !endsWith(nodes$path, "/*"), , drop = FALSE]
}


# params ####

# Bind each param on the profile, in declaration order. Returns list(values,
# unbound) with values a named list of strings.
# `forced`: named values that replace the rule of those params (e.g. a
# spat_unit / feat_type chosen by the caller)
.bind_params <- function(mapping, profile, forced = NULL) {
    lists <- adapterListPaths(.adapter_for(mapping))
    nodes <- profile$nodes
    values <- list()
    unbound <- character()
    for (nm in names(mapping$params)) {
        if (!is.null(forced[[nm]])) {
            values[[nm]] <- as.character(forced[[nm]])
            next
        }
        pr <- mapping$params[[nm]]
        on <- if (is.null(pr$on)) NULL else .pf_bind(pr$on, values)
        if (!is.null(pr$value)) pr$value <- .pf_bind(pr$value, values)
        if (!is.null(on) && grepl("{", on, fixed = TRUE)) {
            unbound <- c(unbound, nm)
            next
        }
        kids <- .list_kids(nodes, lists, pr$list, on)
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
    if (!any(ok)) return(NA_character_)
    hit <- which(ok)[1L]
    if (identical(pr$as, "path")) kids$path[hit] else name[hit]
}

.param_where <- function(w, kids, name, nodes) {
    has <- function(h, roles) {
        hs <- strsplit(ifelse(is.na(h), "", h), ",", fixed = TRUE)
        vapply(hs, function(x) all(roles %in% x), logical(1))
    }
    ok <- rep(TRUE, nrow(kids))
    if (!is.null(w$class)) ok <- ok & kids$class %in% unlist(w$class)
    if (!is.null(w$kind)) ok <- ok & kids$kind %in% unlist(w$kind)
    if (!is.null(w$role)) ok <- ok & has(kids$role_hints, unlist(w$role))
    if (!is.null(w$not_role)) {
        hs <- strsplit(ifelse(is.na(kids$role_hints), "", kids$role_hints), ",", fixed = TRUE)
        ok <- ok & !vapply(hs, function(x) any(unlist(w$not_role) %in% x), logical(1))
    }
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
# accessor args that only refer to optional params which did not bind are
# left out, so the accessor uses its default
.drop_optional_args <- function(node, mapping, params) {
    args <- node$get$args
    if (!length(args)) return(node)
    opt <- names(Filter(function(pr) isTRUE(pr$optional), mapping$params))
    gone <- intersect(params$unbound, opt)
    keep <- vapply(args, function(a) {
        if (!is.character(a)) return(TRUE)
        ph <- gsub("[{}]", "", regmatches(a, gregexpr("\\{[^}]+\\}", a))[[1L]])
        !length(ph) || !all(ph %in% gone)
    }, logical(1))
    node$get$args <- args[keep]
    node
}

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
        # keep rows whose names match `pattern` (regex), or the others with invert
        select_rows = {
            hit <- grepl(o$pattern, rownames(data))
            if (isTRUE(o$invert)) hit <- !hit
            data[hit, , drop = FALSE]
        },
        replace = {
            for (cl in unlist(o$cols)) data[[cl]] <- gsub(o$pattern, o$replacement, data[[cl]], fixed = TRUE)
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
    # Giotto keeps spatial networks as undirected igraphs over all cells
    if (p$part == "spatial_network" && is.data.frame(data)) {
        data <- .edges_to_igraph(data, attr(data, "nodes"))
    }
    ctor <- getExportedValue("GiottoClass", cp$constructor)
    args <- stats::setNames(list(data), cp$data_arg)
    to <- p$to
    for (k in intersect(names(to), cp$args)) args[[k]] <- to[[k]]
    if ("misc" %in% cp$args) {
        # kept next to misc a hook set (e.g. PCA loadings)
        args$misc <- c(args$misc, list(giottoConnect = list(mapping = mapping$mapping$id,
            version = mapping$mapping$version, part = p$id)))
    }
    if ("verbose" %in% cp$args) args$verbose <- FALSE
    if (p$part == "spatial_info") args$calc_centroids <- FALSE
    if (p$part == "feature_info" && !is.null(args$split_keyword)) {
        # split only on keywords that occur; the first feat_type takes the rest
        kw <- unlist(args$split_keyword)
        hit <- vapply(kw, function(k) any(grepl(k, data$feat_ID, fixed = TRUE)), logical(1))
        args$feat_type <- c(unlist(args$feat_type)[1L], unlist(args$feat_type)[-1L][hit])
        args$split_keyword <- if (any(hit)) as.list(kw[hit])
        if (!any(hit)) args$feat_type <- args$feat_type[1L]
    }
    suppressMessages(do.call(ctor, args))
}

# edge table (from, to, distance, weight) -> undirected igraph; `nodes`
# defaults to the cells on an edge
.edges_to_igraph <- function(edges, nodes = NULL) {
    e <- as.data.frame(edges)[, intersect(c("from", "to", "distance", "weight"), names(edges))]
    v <- data.frame(name = nodes %null% unique(c(e$from, e$to)))
    igraph::graph_from_data_frame(e, directed = FALSE, vertices = v)
}

.assemble <- function(parts, source = NULL) {
    g <- if (is.null(source)) GiottoClass::giotto() else GiottoClass::giotto(source = source)
    # a constructor may return several subobjects (points split by feat_type)
    parts <- do.call(c, lapply(unname(parts), function(p) if (is.list(p) && !isS4(p)) p else list(p)))
    g <- GiottoClass::setGiotto(g, parts, initialize = FALSE, verbose = FALSE)
    methods::initialize(g)
}


# context, losses, fingerprints ####

.new_ctx <- function(mapping, params, profile) {
    ctx <- new.env(parent = emptyenv())
    ctx$losses <- list()
    ctx$absent <- character()
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
    # a declared drop only applies when the object holds data at that path
    # (an empty slot or list holds nothing)
    nodes <- ctx$profile$nodes
    paths <- nodes$path
    data_kinds <- c("sparse_matrix", "dense_matrix", "data.frame", "vector", "factor")
    held <- Filter(function(l) {
        if (is.null(l$path)) return(TRUE)
        p <- .pf_bind(l$path, params)
        below <- startsWith(paths, paste0(p, "/"))
        any(below) || any(paths == p & nodes$kind %in% data_kinds)
    }, mapping$lossy)
    top <- lapply(held, function(l) data.frame(
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
# source paths of the fixed parts. For fingerprints: required parts only,
# since optional parts may be absent without the structure having changed
# (their absence is recorded in the provenance instead).
.mapping_paths <- function(mapping, optional = FALSE) {
    fixed <- Filter(function(p) is.null(p$for_each) && (optional || !isTRUE(p$optional)),
        mapping$parts)
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
