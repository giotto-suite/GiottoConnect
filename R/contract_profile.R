# Structural profiler for external objects. Walks any R object through
# attributes() so the owning package (e.g. SeuratObject) is never loaded.
# On S4 nodes only isS4/typeof/class/attr/attributes/unclass are used:
# inherits(), length(), dim() and is.matrix() dispatch and load the
# class's namespace.

.PROFILER_VERSION <- "0.1.0"

# packages whose objects are recorded but not descended into
.PROFILE_STOP_PKGS <- c(
    "BPCells", "DelayedArray", "HDF5Array", "TileDBArray", "terra", "sp", "sf"
)

#' @title Profile the structure of an object
#' @name profileObject
#' @description Describe the structure of any R object (S4, lists,
#' environments, matrices, tables) as a flat node table, without loading the
#' package that defines it. Profiles are used to choose a mapping, to compute
#' structural fingerprints and to detect drift between object versions.
#' @param x object to profile
#' @param max_depth maximum depth walked
#' @param max_children maximum number of named children walked per node;
#' unnamed lists are represented by their first element only
#' @param sample_n maximum number of values inspected per matrix or vector
#' for role hints
#' @param n_cells,n_feat optional known numbers of cells and features. When
#' `NULL` they are inferred from the largest count-like matrix.
#' @returns a `giottoProfile`: a list with `source`, `profiled_at`,
#' `profiler_version`, `n_cells`, `n_feat`, `size_from`, `notes` and `nodes`
#' (a data.frame with columns `path`, `kind`, `class`, `dtype`, `shape`,
#' `n_children`, `truncated`, `role_hints`, `sample`, `value`; `value` holds
#' scalar vectors and is not part of the fingerprint)
#' @examples
#' p <- profileObject(list(counts = matrix(1:6, 2), meta = data.frame(a = 1:3)))
#' p
#' @export
profileObject <- function(x,
    max_depth = 8L,
    max_children = 50L,
    sample_n = 10000L,
    n_cells = NULL,
    n_feat = NULL) {
    env <- new.env(parent = emptyenv())
    env$rows <- list()
    env$notes <- character()
    .pf_walk(x, "", 0L, env,
        max_depth = max_depth, max_children = max_children,
        sample_n = sample_n
    )
    nodes <- as.data.frame(data.table::rbindlist(env$rows))

    size <- .pf_infer_size(nodes, n_cells = n_cells, n_feat = n_feat)
    nodes$role_hints <- .pf_size_hints(nodes, size$n_cells, size$n_feat)
    nodes$.flags <- NULL

    structure(list(
        source = .pf_source(x),
        profiled_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        profiler_version = .PROFILER_VERSION,
        n_cells = size$n_cells,
        n_feat = size$n_feat,
        size_from = size$from,
        notes = unique(env$notes),
        nodes = nodes
    ), class = "giottoProfile")
}

#' @export
print.giottoProfile <- function(x, ...) {
    src <- x$source
    ver <- if (is.na(src$version)) "" else paste0(" ", src$version)
    pkg <- if (is.na(src$package)) "" else paste0(" (", src$package, ver, ")")
    cat(sprintf("<giottoProfile> %s%s\n", src$class, pkg))
    if (!is.na(x$n_cells)) {
        cat(sprintf("  %s cells x %s features (from %s)\n",
            x$n_cells, x$n_feat, x$size_from))
    }
    nd <- x$nodes
    cat(sprintf("  %d nodes, depth %d\n", nrow(nd),
        max(.pf_depth(nd$path))))
    top <- nd[.pf_depth(nd$path) == 1L, , drop = FALSE]
    for (i in seq_len(nrow(top))) {
        p <- top$path[i]
        line <- .pf_describe_children(nd, p, top[i, ])
        if (!is.null(line)) cat(sprintf("  %-12s: %s\n", p, line))
    }
    if (length(x$notes)) cat("  notes:", paste(x$notes, collapse = "; "), "\n")
    invisible(x)
}

#' @title Structural fingerprint of a profile
#' @name fingerprint
#' @description Hash the structure of a profiled object: the sorted set of
#' (templated path, kind, dtype, ndim). Sizes and values are excluded, so two
#' datasets with the same structure share a fingerprint. Parameter values in
#' paths are replaced by their `{param}` placeholder first.
#' @param profile a `giottoProfile`
#' @param paths templated paths to restrict the fingerprint to (plus their
#' ancestors). Missing paths are part of the hash. Ignored when
#' `scope = "full"`.
#' @param params named list of bound parameter values, e.g.
#' `list(assay = "Xenium")`. For `scope = "paths"` they fill the `{param}`
#' placeholders of `paths`. For `scope = "full"` their values are replaced by
#' placeholders in the profile, as whole path segments or as `_`/`.`-delimited
#' tokens (`Xenium_nn`, `FindNeighbors.Xenium.pca`); pass only object names
#' there, since a layer name like `data` would also match `meta.data`.
#' @param scope `"paths"` (default when `paths` is given) or `"full"`
#' @returns character, an xxhash64 digest
#' @examples
#' p <- profileObject(list(a = list(counts = matrix(1:4, 2))))
#' fingerprint(p, scope = "full")
#' fingerprint(p, paths = "a/counts")
#' @export
fingerprint <- function(profile,
    paths = NULL,
    params = NULL,
    scope = if (is.null(paths)) "full" else "paths") {
    checkmate::assert_class(profile, "giottoProfile")
    scope <- match.arg(scope, c("paths", "full"))
    nd <- profile$nodes
    if (scope == "full") {
        tpath <- .pf_template(nd$path, params)
        sig <- paste(tpath, nd$kind, nd$dtype, .pf_ndim(nd$shape), sep = "|")
    } else {
        if (is.null(paths)) stop(wrap_txt("`paths` is required when scope = 'paths'"))
        # bind templated paths forward and look the concrete paths up
        want <- unique(unlist(lapply(paths, .pf_with_ancestors)))
        found <- match(.pf_bind(want, params), nd$path)
        hit <- !is.na(found)
        sig <- c(
            paste(want[hit], nd$kind[found[hit]], nd$dtype[found[hit]],
                .pf_ndim(nd$shape[found[hit]]), sep = "|"),
            paste(want[!hit], "(missing)", sep = "|")
        )
    }
    # radix sort is C-locale ordering, so the hash does not depend on LC_COLLATE
    digest::digest(paste(sort(unique(sig), method = "radix"), collapse = "\n"),
        algo = "xxhash64", serialize = FALSE
    )
}

#' @title Compare two profiles
#' @name diffProfiles
#' @description Classify structural changes between two profiles. Removed
#' and added subtrees are reported once, at their root. A removed subtree that
#' reappears elsewhere with the same kind, dtype, dimensionality and size roles
#' is reported as `renamed_key` (same parent) or `moved` (different parent).
#' @param old,new `giottoProfile` objects
#' @param params_old,params_new optional bound object names (e.g. assay
#' names), templated as in [fingerprint()] before comparing
#' @returns a `giottoProfileDiff` data.frame with columns `change` (`added`,
#' `removed`, `moved`, `renamed_key`, `kind_changed`, `dtype_changed`),
#' `path`, `new_path`, `kind`, `class`, `detail`
#' @examples
#' a <- profileObject(list(x = list(counts = matrix(1:4, 2))))
#' b <- profileObject(list(x = list(raw = matrix(1:4, 2))))
#' diffProfiles(a, b)
#' @export
diffProfiles <- function(old, new, params_old = NULL, params_new = NULL) {
    checkmate::assert_class(old, "giottoProfile")
    checkmate::assert_class(new, "giottoProfile")
    a <- old$nodes
    b <- new$nodes
    a$tpath <- .pf_template(a$path, params_old)
    b$tpath <- .pf_template(b$path, params_new)

    out <- list()
    add <- function(change, path, new_path, kind, class, detail) {
        out[[length(out) + 1L]] <<- data.frame(
            change = change, path = path, new_path = new_path, kind = kind,
            class = class, detail = detail, stringsAsFactors = FALSE
        )
    }

    common <- intersect(a$tpath, b$tpath)
    ia <- match(common, a$tpath)
    ib <- match(common, b$tpath)
    for (k in seq_along(common)) {
        ra <- a[ia[k], ]
        rb <- b[ib[k], ]
        if (ra$kind != rb$kind) {
            add("kind_changed", ra$tpath, rb$tpath, rb$kind, rb$class,
                paste(ra$kind, "->", rb$kind))
        } else if (!identical(ra$dtype, rb$dtype)) {
            add("dtype_changed", ra$tpath, rb$tpath, rb$kind, rb$class,
                paste(ra$dtype, "->", rb$dtype))
        }
    }

    removed <- .pf_subtree_roots(setdiff(a$tpath, b$tpath))
    added <- .pf_subtree_roots(setdiff(b$tpath, a$tpath))
    sig_a <- vapply(removed, .pf_subtree_sig, character(1), nodes = a)
    sig_b <- vapply(added, .pf_subtree_sig, character(1), nodes = b)
    matched_b <- character()
    for (i in seq_along(removed)) {
        cand <- added[sig_b == sig_a[i] & !added %in% matched_b]
        if (length(cand) == 1L) {
            ra <- a[match(removed[i], a$tpath), ]
            same_parent <- .pf_parent(removed[i]) == .pf_parent(cand)
            add(if (same_parent) "renamed_key" else "moved", removed[i], cand,
                ra$kind, ra$class, "same structure")
            matched_b <- c(matched_b, cand)
        } else {
            ra <- a[match(removed[i], a$tpath), ]
            add("removed", removed[i], NA_character_, ra$kind, ra$class,
                .pf_subtree_detail(removed[i], a))
        }
    }
    for (p in setdiff(added, matched_b)) {
        rb <- b[match(p, b$tpath), ]
        add("added", NA_character_, p, rb$kind, rb$class,
            .pf_subtree_detail(p, b))
    }

    res <- if (length(out)) {
        do.call(rbind, out)
    } else {
        data.frame(
            change = character(), path = character(), new_path = character(),
            kind = character(), class = character(), detail = character()
        )
    }
    class(res) <- c("giottoProfileDiff", "data.frame")
    res
}

#' @export
print.giottoProfileDiff <- function(x, ...) {
    if (!nrow(x)) {
        cat("<giottoProfileDiff> no structural changes\n")
        return(invisible(x))
    }
    cat(sprintf("<giottoProfileDiff> %d change(s)\n", nrow(x)))
    print(as.data.frame(unclass(x)[c("change", "path", "new_path", "kind", "detail")],
        stringsAsFactors = FALSE), row.names = FALSE)
    invisible(x)
}


# walker ####

.pf_walk <- function(x, path, depth, env, max_depth, max_children, sample_n) {
    info <- .pf_node(x, sample_n)
    children <- NULL
    truncated <- FALSE
    if (info$descend && depth < max_depth) {
        children <- .pf_children(x, info$kind)
        if (!is.null(children) && length(children) > max_children) {
            truncated <- TRUE
            children <- children[seq_len(max_children)]
        }
    } else if (info$descend) {
        truncated <- TRUE
    }
    if (!is.null(info$note)) env$notes <- c(env$notes, info$note)

    env$rows[[length(env$rows) + 1L]] <- list(
        path = path, kind = info$kind, class = info$class, dtype = info$dtype,
        shape = info$shape, n_children = as.integer(info$n_children),
        truncated = truncated, role_hints = NA_character_,
        sample = info$sample, value = info$value, .flags = info$flags
    )
    for (nm in names(children)) {
        .pf_walk(children[[nm]], .pf_join(path, nm), depth + 1L, env,
            max_depth = max_depth, max_children = max_children,
            sample_n = sample_n
        )
    }
}

# describe a single node using only dispatch-free primitives
.pf_node <- function(x, sample_n) {
    cls <- if (is.null(x)) "NULL" else class(x)[1L]
    pkg <- attr(class(x), "package") %null% NA_character_
    s4 <- isS4(x)
    a <- attributes(x)
    out <- list(
        kind = "other", class = cls, dtype = NA_character_,
        shape = NA_character_, n_children = 0L, descend = FALSE,
        sample = NA_character_, value = NA_character_, flags = "", note = NULL
    )
    if (is.null(x)) {
        out$kind <- "null"
        return(out)
    }
    if (!is.na(pkg) && pkg %in% .PROFILE_STOP_PKGS) {
        out$kind <- if (typeof(x) == "externalptr") "external_pointer" else "S4"
        out$note <- sprintf("profiling %s requires package %s", cls, pkg)
        return(out)
    }
    tp <- typeof(x)
    if (tp == "externalptr") {
        out$kind <- "external_pointer"
        return(out)
    }
    if (tp %in% c("closure", "builtin", "special")) {
        out$kind <- "function"
        return(out)
    }
    if (tp == "environment") {
        out$kind <- "environment"
        out$n_children <- length(ls(x, all.names = TRUE))
        out$descend <- TRUE
        return(out)
    }
    # sparse matrix by structure: Dim + p/i/j + x (dgCMatrix, Graph, ...)
    if (s4 && !is.null(a$Dim) && (!is.null(a$i) || !is.null(a$j))) {
        out$kind <- "sparse_matrix"
        out$shape <- paste(a$Dim, collapse = "x")
        vals <- a$x
        out$dtype <- if (is.null(vals)) "logical" else typeof(vals)
        out$flags <- .pf_value_flags(vals, sample_n)
        out$sample <- .pf_sample_names(a$Dimnames)
        if (!is.null(a$Dimnames) && any(lengths(a$Dimnames) > 0L)) {
            out$flags <- .pf_add_flag(out$flags, "has_dimnames")
        }
        return(out)
    }
    if (!is.null(a$dim) && length(a$dim) == 2L &&
        tp %in% c("logical", "integer", "double", "complex", "character")) {
        out$kind <- "dense_matrix"
        out$shape <- paste(a$dim, collapse = "x")
        out$dtype <- tp
        vals <- unclass(x)
        attributes(vals) <- NULL
        out$flags <- .pf_value_flags(vals, sample_n)
        out$sample <- .pf_sample_names(a$dimnames)
        if (!is.null(a$dimnames) && any(lengths(a$dimnames) > 0L)) {
            out$flags <- .pf_add_flag(out$flags, "has_dimnames")
        }
        return(out)
    }
    if (!s4 && is.data.frame(x)) {
        out$kind <- "data.frame"
        out$shape <- paste(nrow(x), ncol(x), sep = "x")
        out$n_children <- ncol(x)
        out$descend <- TRUE
        out$sample <- .pf_short(names(x))
        num <- vapply(x, is.numeric, logical(1))
        if (sum(num) >= 2L) out$flags <- "two_numeric_cols"
        return(out)
    }
    if (s4) {
        out$kind <- "S4"
        slots <- setdiff(names(a), "class")
        out$n_children <- length(slots)
        out$descend <- TRUE
        out$sample <- .pf_short(slots)
        return(out)
    }
    if (inherits(x, "package_version")) {
        out$kind <- "vector"
        out$dtype <- "character"
        out$shape <- "1"
        out$sample <- as.character(x)
        return(out)
    }
    if (is.factor(x)) {
        out$kind <- "factor"
        out$dtype <- "factor"
        out$shape <- as.character(length(x))
        out$sample <- .pf_short(levels(x))
        return(out)
    }
    if (tp == "list") {
        out$kind <- "list"
        out$n_children <- length(x)
        out$descend <- TRUE
        out$sample <- .pf_short(names(x))
        return(out)
    }
    if (is.atomic(x)) {
        out$kind <- "vector"
        out$dtype <- tp
        out$shape <- as.character(length(x))
        # scalar values can bind params (e.g. Seurat's active.assay)
        if (length(x) == 1L && !is.na(x)) out$value <- as.character(x)
        vals <- x
        attributes(vals) <- NULL
        out$flags <- .pf_value_flags(vals, sample_n)
        out$sample <- .pf_short(if (!is.null(names(x))) names(x) else vals)
        return(out)
    }
    out
}

# children as a named list; unnamed lists collapse to a representative "*"
.pf_children <- function(x, kind) {
    switch(kind,
        S4 = {
            a <- attributes(x)
            a[setdiff(names(a), "class")]
        },
        environment = mget(sort(ls(x, all.names = TRUE)), envir = x),
        data.frame = as.list(x),
        list = {
            nms <- names(x)
            if (is.null(nms) || all(!nzchar(nms))) {
                if (length(x)) list(`*` = x[[1L]]) else NULL
            } else {
                nms[!nzchar(nms)] <- "*"
                stats::setNames(as.list(x), nms)
            }
        },
        NULL
    )
}

# value flags from a strided sample of at most n values
.pf_value_flags <- function(v, n) {
    if (is.null(v) || !length(v) || !is.numeric(v)) return("")
    if (length(v) > n) v <- v[as.integer(seq(1L, length(v), length.out = n))]
    v <- v[!is.na(v)]
    if (!length(v)) return("")
    flags <- character()
    if (all(v == round(v))) flags <- c(flags, "integer_valued")
    if (all(v >= 0)) flags <- c(flags, "nonneg")
    if (all(v %in% c(0, 1))) flags <- c(flags, "binary_valued")
    paste(flags, collapse = ",")
}


# size inference and size-based role hints ####

.pf_infer_size <- function(nodes, n_cells = NULL, n_feat = NULL) {
    from <- NA_character_
    if (is.null(n_cells) || is.null(n_feat)) {
        mats <- nodes[nodes$kind %in% c("sparse_matrix", "dense_matrix") &
            grepl("integer_valued", nodes$.flags) &
            grepl("nonneg", nodes$.flags) &
            !grepl("binary_valued", nodes$.flags), , drop = FALSE]
        if (nrow(mats)) {
            dims <- do.call(rbind, lapply(strsplit(mats$shape, "x"), as.numeric))
            dims <- matrix(dims, ncol = 2L)
            best <- which.max(dims[, 1] * dims[, 2])
            d <- dims[best, ]
            from <- mats$path[best]
            # the shallowest table whose row count matches one dimension marks
            # the cells (object-level cell metadata sits above per-assay
            # feature metadata), then a table named for cells; otherwise the
            # larger dimension
            tabs <- nodes[nodes$kind == "data.frame", , drop = FALSE]
            trows <- as.numeric(sub("x.*", "", tabs$shape))
            hit <- trows %in% d
            cells_dim <- if (any(hit)) {
                tdepth <- .pf_depth(tabs$path)
                top <- hit & tdepth == min(tdepth[hit])
                rows_top <- unique(trows[top])
                # tie between cell and feature tables at the same depth:
                # decide by their names
                named <- top & grepl("cell", tabs$path, ignore.case = TRUE) &
                    !grepl("feat|gene|var", tabs$path, ignore.case = TRUE)
                if (length(rows_top) == 1L && sum(d == rows_top) == 1L) {
                    match(rows_top, d)
                } else if (length(unique(trows[named])) == 1L && sum(d == unique(trows[named])) == 1L) {
                    match(unique(trows[named]), d)
                } else {
                    which.max(d)
                }
            } else {
                which.max(d)
            }
            n_cells <- n_cells %null% d[cells_dim]
            n_feat <- n_feat %null% d[3L - cells_dim]
        }
    } else {
        from <- "user"
    }
    list(
        n_cells = if (is.null(n_cells)) NA_integer_ else as.integer(n_cells),
        n_feat = if (is.null(n_feat)) NA_integer_ else as.integer(n_feat),
        from = from
    )
}

.pf_size_hints <- function(nodes, n_cells, n_feat) {
    hints <- strsplit(nodes$.flags, ",", fixed = TRUE)
    for (i in seq_len(nrow(nodes))) {
        sh <- nodes$shape[i]
        if (is.na(sh)) next
        d <- as.numeric(strsplit(sh, "x", fixed = TRUE)[[1L]])
        h <- hints[[i]]
        if (!is.na(n_cells)) {
            if (d[1] == n_cells) h <- c(h, "rows_eq_ncells")
            if (length(d) == 2L && d[2] == n_cells) h <- c(h, "cols_eq_ncells")
        }
        if (!is.na(n_feat)) {
            if (d[1] == n_feat) h <- c(h, "rows_eq_nfeat")
            if (length(d) == 2L && d[2] == n_feat) h <- c(h, "cols_eq_nfeat")
        }
        hints[[i]] <- h
    }
    vapply(hints, function(h) {
        h <- h[nzchar(h)]
        if (length(h)) paste(h, collapse = ",") else NA_character_
    }, character(1))
}


# paths ####

.pf_join <- function(parent, name) {
    name <- gsub("/", "%2F", name, fixed = TRUE)
    if (!nzchar(parent)) name else paste(parent, name, sep = "/")
}

.pf_depth <- function(path) {
    ifelse(nzchar(path), lengths(strsplit(path, "/", fixed = TRUE)), 0L)
}

.pf_parent <- function(path) {
    ifelse(grepl("/", path, fixed = TRUE), sub("/[^/]*$", "", path), "")
}

.pf_with_ancestors <- function(path) {
    parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
    vapply(seq_along(parts), function(i) {
        paste(parts[seq_len(i)], collapse = "/")
    }, character(1))
}

# replace bound param values by {param}, as whole segments or as tokens
# delimited by "_" or "." within a segment ("Xenium_nn", "nCount_Xenium",
# "FindNeighbors.Xenium.pca"). Pass only object names (assays, images,
# reductions): a layer value like "data" would also match "meta.data".
.pf_template <- function(paths, params) {
    if (!length(params)) return(paths)
    vals <- unlist(params, use.names = TRUE)
    vals <- vals[nzchar(vals)]
    vals <- vals[order(nchar(vals), decreasing = TRUE)]
    for (nm in names(vals)) {
        pat <- sprintf("(?<=^|/|_|\\.)%s(?=_|/|$|\\.)",
            .pf_regex_escape(vals[[nm]]))
        paths <- gsub(pat, sprintf("{%s}", nm), paths, perl = TRUE)
    }
    paths
}

.pf_regex_escape <- function(x) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

# substitute {param} placeholders by bound values
.pf_bind <- function(path, params) {
    for (nm in names(params)) {
        path <- gsub(sprintf("{%s}", nm), params[[nm]], path, fixed = TRUE)
    }
    path
}

.pf_ndim <- function(shape) {
    ifelse(is.na(shape), 0L, lengths(strsplit(shape, "x", fixed = TRUE)))
}

# keep only paths whose parent is not itself in the set
.pf_subtree_roots <- function(paths) {
    if (!length(paths)) return(character())
    paths[!.pf_parent(paths) %in% paths]
}

# structure of a subtree relative to its root, sizes as roles only
.pf_subtree_sig <- function(root, nodes) {
    inside <- nodes$tpath == root |
        startsWith(nodes$tpath, paste0(root, "/"))
    sub <- nodes[inside, , drop = FALSE]
    rel <- substring(sub$tpath, nchar(root) + 1L)
    roles <- vapply(strsplit(ifelse(is.na(sub$role_hints), "", sub$role_hints), ","),
        function(h) paste(sort(grep("_eq_", h, value = TRUE)), collapse = "+"),
        character(1))
    sig <- paste(rel, sub$kind, sub$dtype, .pf_ndim(sub$shape), roles, sep = "|")
    paste(sort(sig), collapse = "\n")
}

.pf_subtree_detail <- function(root, nodes) {
    n <- sum(startsWith(nodes$tpath, paste0(root, "/")))
    if (n) sprintf("subtree of %d node(s)", n + 1L) else ""
}


# source and printing helpers ####

.pf_source <- function(x) {
    cls <- class(x)[1L]
    pkg <- attr(class(x), "package") %null% NA_character_
    ver <- NA_character_
    a <- attributes(x)
    # Seurat-style objects record their package version in a slot
    if (!is.null(a$version) && inherits(a$version, "package_version")) {
        ver <- as.character(a$version)
    } else if (!is.na(pkg) && requireNamespace("utils", quietly = TRUE)) {
        ver <- tryCatch(as.character(utils::packageVersion(pkg)),
            error = function(e) NA_character_
        )
    }
    list(class = cls, package = pkg, version = ver)
}

.pf_short <- function(v, n = 3L) {
    if (is.null(v) || !length(v)) return(NA_character_)
    v <- as.character(utils::head(v, n))
    s <- paste(ifelse(nchar(v) > 20L, paste0(substr(v, 1L, 17L), "..."), v),
        collapse = ", ")
    s
}

.pf_sample_names <- function(dn) {
    if (is.null(dn)) return(NA_character_)
    parts <- vapply(dn, function(d) {
        if (is.null(d) || !length(d)) "" else .pf_short(d, 2L)
    }, character(1))
    if (all(!nzchar(parts))) NA_character_ else paste(parts, collapse = " | ")
}

.pf_add_flag <- function(flags, f) {
    if (!nzchar(flags)) f else paste(flags, f, sep = ",")
}

.pf_describe_children <- function(nodes, path, row) {
    if (row$kind %in% c("data.frame", "sparse_matrix", "dense_matrix")) {
        return(sprintf("%s %s: %s", row$kind, row$shape,
            row$sample %na% ""))
    }
    kids <- nodes[.pf_parent(nodes$path) == path &
        .pf_depth(nodes$path) == 2L, , drop = FALSE]
    if (!nrow(kids)) return(NULL)
    lab <- sprintf("%s [%s]", sub(".*/", "", kids$path), kids$class)
    if (length(lab) > 6L) lab <- c(lab[1:6], sprintf("... (%d)", length(lab)))
    paste(lab, collapse = ", ")
}
