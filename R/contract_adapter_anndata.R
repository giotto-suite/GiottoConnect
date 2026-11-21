# AnnData adapter. The input is the path to an .h5ad file; nothing is loaded
# up front. The profile walks the HDF5 groups and reads only attributes,
# shapes and small value samples. Accessors decode one AnnData element each.
#
# Orientation: AnnData stores obs x var, row-major. Profile shapes are given
# as stored (obs x var). `h5adMatrix` returns var x obs, the orientation of R
# packages: for CSR and dense storage that costs no transpose (hdf5r reports
# dims reversed). obsm and obsp elements keep obs as rows.

setClass("anndataAdapter", contains = "giottoConnectAdapter",
    prototype = list(ecosystem = "anndata", package = "hdf5r"))

setMethod("adapterListPaths", "anndataAdapter", function(adapter) c(
    # X, raw/X and every layer; bind with `as: path`
    Matrices = "|layers|raw", Layers = "layers",
    Obsm = "obsm", Obsp = "obsp", Varm = "varm", Uns = "uns",
    # categories of var/feature_types (10x: Gene Expression, Antibody Capture, ...)
    FeatureTypes = "var/feature_types",
    children = "{on}", columns = "{on}"
))
setMethod("adapterAccessors", "anndataAdapter", function(adapter) c(
    "h5adMatrix", "h5adFrame", "h5adVar", "h5adObsm", "h5adObsp", "h5adValue", "path"
))
setMethod("adapterObjectLists", "anndataAdapter", function(adapter) {
    c("Layers", "Obsm", "Obsp")
})
setMethod("adapterProfile", "anndataAdapter", function(adapter, x, ...) {
    .h5ad_profile(.h5ad_path(x), ...)
})
# every accessor reads the element at `get.on`, else at the node's path
setMethod("adapterRead", "anndataAdapter", function(adapter, x, node, params) {
    get <- node$get
    el <- .pf_bind(get$on %null% node$path, params)
    args <- lapply(get$args %null% list(), function(a) {
        if (is.character(a)) .pf_bind(a, params) else a
    })
    h5 <- .h5ad_open(.h5ad_path(x))
    on.exit(h5$close_all())
    fn <- .H5AD_ACCESSORS[[get$fn]]
    do.call(fn, c(list(h5, el), args))
})

# matrices are named by their var table; feature metadata for raw/X counts
# reads raw/var
setMethod("adapterReads", "anndataAdapter", function(adapter, node, params) {
    fn <- node$get$fn
    if (identical(fn, "h5adMatrix")) return(.h5ad_var_group(.pf_bind(node$path, params)))
    # X_pca carries its loadings (read by the anndata_pca_misc hook)
    if (identical(fn, "h5adObsm") && identical(.pf_bind(node$path, params), "obsm/X_pca")) {
        return("varm/PCs")
    }
    # a squidpy spatial graph is read with its <key>distances partner
    if (identical(fn, "h5adObsp")) {
        p <- .pf_bind(node$path, params)
        if (grepl("connectivities$", p)) return(sub("connectivities$", "distances", p))
    }
    m <- node$get$args$matrix
    if (identical(fn, "h5adVar") && is.character(m) && startsWith(.pf_bind(m, params), "raw/")) {
        return("raw/var")
    }
    character()
})

# X, raw/X and layers stream into a GiottoDisk project with h5adInput()
setMethod("adapterStream", "anndataAdapter", function(adapter, x, node, params) {
    # a feature-type subset is read into R
    if (!is.null(node$get$args$feature_type)) return(NULL)
    if (!identical(node$get$fn, "h5adMatrix") || !requireNamespace("GiottoDisk", quietly = TRUE) ||
        !exists("h5adInput", envir = asNamespace("GiottoDisk"))) {
        return(NULL)
    }
    el <- .pf_bind(node$get$on %null% node$path, params)
    layer <- if (el == "X") NULL else if (el == "raw/X") "raw" else sub("^layers/", "", el)
    suppressMessages(GiottoDisk::h5adInput(.h5ad_path(x), layer = layer))
})

.h5ad_path <- function(x) {
    if (!is.character(x) || length(x) != 1L || !file.exists(x)) {
        stop(wrap_txt("AnnData input must be the path to an existing .h5ad file"), call. = FALSE)
    }
    normalizePath(x)
}


# element readers ####

.h5ad_open <- function(path) {
    package_check("hdf5r")
    hdf5r::H5File$new(path, mode = "r")
}

# walk one level at a time: H5Lexists errors on a missing intermediate group
# or one that is a dataset
.h5ad_exists <- function(h5, p) {
    parts <- strsplit(p, "/", fixed = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    cur <- character()
    for (j in seq_along(parts)) {
        cur <- paste(c(cur, parts[[j]]), collapse = "/")
        if (!h5$exists(cur)) return(FALSE)
        # only groups have children
        if (j < length(parts) && !.h5ad_is_group(h5[[cur]])) return(FALSE)
    }
    TRUE
}

.h5ad_attr <- function(node, name) {
    if (!node$attr_exists(name)) return(NULL)
    hdf5r::h5attr(node, name)
}

.h5ad_encoding <- function(node) {
    enc <- .h5ad_attr(node, "encoding-type")
    if (!is.null(enc)) return(enc)
    if (.h5ad_is_group(node)) "dict" else "array"
}

.h5ad_is_group <- function(node) inherits(node, c("H5Group", "H5File"))

# obs / var index: the group's `_index` attribute names the index dataset
.h5ad_index <- function(h5, group) {
    g <- h5[[group]]
    as.character(g[[.h5ad_attr(g, "_index") %null% "_index"]][])
}

.h5ad_feat_ids <- function(h5, group) .unique_feat_ids(.h5ad_index(h5, group))

# duplicated feature names get "--1", "--2", ... (as GiottoDisk readers)
.unique_feat_ids <- function(ids) {
    dup <- ids %in% ids[duplicated(ids)]
    if (!any(dup)) return(ids)
    ids[dup] <- paste0(ids[dup], "--", stats::ave(seq_along(ids)[dup], ids[dup], FUN = seq_along))
    ids
}

.h5ad_read_array <- function(node) {
    x <- node$read()
    if (length(node$dims) == 2L) x <- t(x)
    # bool enum on older hdf5r
    if (is.factor(x) && all(levels(x) %in% c("FALSE", "TRUE"))) {
        x <- as.logical(as.character(x))
    }
    x
}

# anndata < 0.8: codes in the column, categories in a sibling `__categories`
# group under the same name
.h5ad_legacy_categorical <- function(h5, node, path) {
    codes <- as.integer(node$read())
    cat_path <- file.path(dirname(path), "__categories", basename(path))
    if (!.h5ad_exists(h5, cat_path)) {
        warning(wrap_txt(sprintf("'%s' is a pre-0.8 categorical without categories; returning its codes",
            path)), call. = FALSE)
        return(codes)
    }
    cats <- as.character(h5[[cat_path]][])
    codes[codes < 0L] <- NA_integer_
    factor(cats[codes + 1L], levels = cats)
}

# one column-like element (array, string-array, categorical, nullable-*)
.h5ad_read_elem <- function(h5, path) {
    node <- h5[[path]]
    enc <- .h5ad_encoding(node)
    switch(enc,
        # scalars (anndata >= 0.8 uns)
        "string" = , "numeric-scalar" = , "numeric" = .h5ad_read_array(node),
        "array" = , "string-array" = {
            if (node$attr_exists("categories")) return(.h5ad_legacy_categorical(h5, node, path))
            .h5ad_read_array(node)
        },
        "categorical" = {
            codes <- as.integer(node[["codes"]][])
            cats <- as.character(node[["categories"]][])
            codes[codes < 0L] <- NA_integer_
            factor(cats[codes + 1L], levels = cats)
        },
        "nullable-integer" = , "nullable-boolean" = , "nullable-string-array" = {
            v <- node[["values"]][]
            v[as.logical(node[["mask"]][])] <- NA
            v
        },
        stop(wrap_txt(sprintf("'%s' has unsupported encoding '%s'", path, enc)), call. = FALSE)
    )
}

# obs, var or raw/var -> data.frame with the index as rownames
.h5ad_frame <- function(h5, path) {
    g <- h5[[path]]
    if (!identical(.h5ad_encoding(g), "dataframe")) {
        stop(wrap_txt(sprintf("'%s' is not an AnnData dataframe", path)), call. = FALSE)
    }
    ids <- .h5ad_index(h5, path)
    if (basename(path) == "var") ids <- .unique_feat_ids(ids)
    cols <- as.character(.h5ad_attr(g, "column-order") %null% character())
    out <- lapply(cols, function(cn) .h5ad_read_elem(h5, paste(path, cn, sep = "/")))
    names(out) <- cols
    out <- if (length(out)) as.data.frame(out, optional = TRUE, stringsAsFactors = FALSE) else
        data.frame(row.names = seq_along(ids))
    rownames(out) <- ids
    out
}

# the feature table that belongs to a matrix
.h5ad_var_group <- function(path) if (startsWith(path, "raw/")) "raw/var" else "var"

# Feature metadata for the features of `matrix`: var, or for raw/X the
# raw/var table (which may hold more features) with the var columns added by
# feature ID; a column in both takes the var value where var has the feature.
.h5ad_var <- function(h5, path, matrix = NULL, feature_type = NULL) {
    out <- .h5ad_var_table(h5, path, matrix)
    if (is.null(feature_type) || !"feature_types" %in% names(out)) return(out)
    out[as.character(out$feature_types) %in% feature_type, , drop = FALSE]
}

.h5ad_var_table <- function(h5, path, matrix = NULL) {
    v <- .h5ad_frame(h5, path)
    if (is.null(matrix) || !startsWith(matrix, "raw/") || !.h5ad_exists(h5, "raw/var")) return(v)
    r <- .h5ad_frame(h5, "raw/var")
    ids <- union(rownames(r), rownames(v))
    ir <- match(ids, rownames(r))
    iv <- match(ids, rownames(v))
    cols <- union(names(r), names(v))
    out <- lapply(cols, function(cn) {
        a <- if (cn %in% names(r)) r[[cn]][ir]
        b <- if (cn %in% names(v)) v[[cn]][iv]
        if (is.null(a)) return(b)
        if (is.null(b)) return(a)
        if (is.factor(a) || is.factor(b)) {
            a <- as.character(a)
            b <- as.character(b)
        }
        a[!is.na(iv)] <- b[!is.na(iv)]
        a
    })
    names(out) <- cols
    out <- as.data.frame(out, optional = TRUE, stringsAsFactors = FALSE)
    rownames(out) <- ids
    out
}

# a csr/csc group as a dgCMatrix with the stored (rows x cols) layout
.h5ad_sparse <- function(node, transposed = FALSE) {
    enc <- .h5ad_encoding(node)
    shape <- as.integer(.h5ad_attr(node, "shape"))
    i <- as.integer(node[["indices"]][])
    p <- as.integer(node[["indptr"]][])
    x <- as.double(node[["data"]][])
    # CSR rows are the compressed axis: read as CSC it is the transpose
    by_row <- identical(enc, "csr_matrix")
    dims <- if (by_row) rev(shape) else shape
    m <- Matrix::sparseMatrix(i = i, p = p, x = x, dims = dims, index1 = FALSE, repr = "C")
    # `m` is t(stored) for CSR, stored for CSC
    if (by_row == transposed) m else Matrix::t(m)
}

# the feature_types of a var table, or NULL
.h5ad_feature_types <- function(h5, group) {
    p <- paste(group, "feature_types", sep = "/")
    if (!.h5ad_exists(h5, p)) return(NULL)
    as.character(.h5ad_read_elem(h5, p))
}

# rows of one feature type (all rows when `feature_type` is NULL)
.h5ad_ft_rows <- function(h5, group, feature_type) {
    if (is.null(feature_type)) return(NULL)
    ft <- .h5ad_feature_types(h5, group)
    if (is.null(ft)) return(NULL)
    which(ft == feature_type)
}

# X, raw/X or layers/<l> -> var x obs, named; `feature_type` keeps the rows
# of one var/feature_types category
.h5ad_matrix <- function(h5, path, feature_type = NULL) {
    node <- h5[[path]]
    enc <- .h5ad_encoding(node)
    m <- if (enc %in% c("csr_matrix", "csc_matrix")) {
        .h5ad_sparse(node, transposed = TRUE)
    } else if (identical(enc, "array")) {
        # dims come back reversed: var x obs as read
        node$read()
    } else {
        stop(wrap_txt(sprintf("'%s' has unsupported matrix encoding '%s'", path, enc)), call. = FALSE)
    }
    dimnames(m) <- list(.h5ad_feat_ids(h5, .h5ad_var_group(path)), .h5ad_index(h5, "obs"))
    rows <- .h5ad_ft_rows(h5, .h5ad_var_group(path), feature_type)
    if (is.null(rows)) m else m[rows, , drop = FALSE]
}

# obsm/<k> -> obs x k matrix (or the data.frame it is stored as)
.h5ad_obsm <- function(h5, path) {
    node <- h5[[path]]
    enc <- .h5ad_encoding(node)
    if (identical(enc, "dataframe")) return(.h5ad_frame(h5, path))
    m <- if (enc %in% c("csr_matrix", "csc_matrix")) as.matrix(.h5ad_sparse(node)) else .h5ad_read_array(node)
    if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
    dimnames(m) <- list(.h5ad_index(h5, "obs"), paste0("Dim.", seq_len(ncol(m))))
    m
}

# obsp/<k> -> obs x obs dgCMatrix, rows = from
.h5ad_obsp <- function(h5, path) {
    node <- h5[[path]]
    enc <- .h5ad_encoding(node)
    m <- if (enc %in% c("csr_matrix", "csc_matrix")) .h5ad_sparse(node) else
        methods::as(.h5ad_read_array(node), "CsparseMatrix")
    ids <- .h5ad_index(h5, "obs")
    dimnames(m) <- list(ids, ids)
    m
}

# uns values: datasets as vectors, groups as named lists
.h5ad_value <- function(h5, path) {
    node <- h5[[path]]
    if (!.h5ad_is_group(node)) return(.h5ad_read_elem(h5, path))
    enc <- .h5ad_encoding(node)
    if (!enc %in% c("dict", "raw", "anndata")) return(.h5ad_any(h5, path))
    kids <- names(node)
    out <- lapply(kids, function(k) .h5ad_value(h5, paste(path, k, sep = "/")))
    names(out) <- kids
    out
}

# structural read: decode by encoding and location
.h5ad_any <- function(h5, path) {
    node <- h5[[path]]
    enc <- .h5ad_encoding(node)
    top <- sub("/.*", "", path)
    if (identical(enc, "dataframe")) return(.h5ad_frame(h5, path))
    if (path == "X" || top %in% c("layers", "raw") && enc %in% c("csr_matrix", "csc_matrix", "array") &&
        !(.h5ad_is_group(node) && enc == "dict")) {
        return(.h5ad_matrix(h5, path))
    }
    if (top == "obsm") return(.h5ad_obsm(h5, path))
    if (top == "obsp") return(.h5ad_obsp(h5, path))
    .h5ad_value(h5, path)
}

.H5AD_ACCESSORS <- list(
    h5adMatrix = .h5ad_matrix, h5adFrame = .h5ad_frame, h5adVar = .h5ad_var, h5adObsm = .h5ad_obsm,
    h5adObsp = .h5ad_obsp, h5adValue = .h5ad_value, path = .h5ad_any
)


# profile ####

# HDF5 groups and datasets -> profile node rows. Only attributes, dims and
# value samples are read.
.h5ad_profile <- function(path, max_depth = 8L, max_children = 50L, sample_n = 10000L, ...) {
    h5 <- .h5ad_open(path)
    on.exit(h5$close_all())
    env <- new.env(parent = emptyenv())
    env$rows <- list()
    env$notes <- character()
    .h5ad_walk(h5, h5, "", 0L, env, max_depth = max_depth, max_children = max_children,
        sample_n = sample_n)
    # categories of var/feature_types as child nodes, for the FeatureTypes listing
    ft <- .h5ad_feature_types(h5, "var")
    for (cat in unique(ft)) {
        env$rows[[length(env$rows) + 1L]] <- list(
            path = .pf_join("var/feature_types", cat), kind = "category", class = "category",
            dtype = NA_character_, shape = as.character(sum(ft == cat)), n_children = 0L,
            truncated = FALSE, role_hints = NA_character_, sample = NA_character_,
            value = NA_character_, .flags = "")
    }
    n_cells <- if (.h5ad_exists(h5, "obs")) length(.h5ad_index(h5, "obs")) else NULL
    n_feat <- if (.h5ad_exists(h5, "var")) length(.h5ad_index(h5, "var")) else NULL
    enc_ver <- .h5ad_attr(h5, "encoding-version")
    source <- list(class = "AnnData", package = "anndata",
        version = if (is.null(enc_ver)) NA_character_ else as.character(enc_ver))
    .pf_finish(env$rows, env$notes, source, n_cells = n_cells, n_feat = n_feat)
}

.h5ad_walk <- function(h5, node, path, depth, env, max_depth, max_children, sample_n) {
    info <- .h5ad_node(h5, node, path, sample_n)
    children <- character()
    truncated <- FALSE
    if (info$descend) {
        children <- info$children %null% setdiff(names(node), "__categories")
        if (depth >= max_depth) {
            truncated <- length(children) > 0L
            children <- character()
        } else if (length(children) > max_children) {
            truncated <- TRUE
            children <- children[seq_len(max_children)]
        }
    }
    env$rows[[length(env$rows) + 1L]] <- list(
        path = path, kind = info$kind, class = info$class, dtype = info$dtype,
        shape = info$shape, n_children = as.integer(info$n_children),
        truncated = truncated, role_hints = NA_character_,
        sample = info$sample, value = info$value, .flags = info$flags
    )
    for (nm in children) {
        child <- if (nzchar(path)) paste(path, nm, sep = "/") else nm
        .h5ad_walk(h5, h5[[child]], .pf_join(path, nm), depth + 1L, env,
            max_depth = max_depth, max_children = max_children, sample_n = sample_n)
    }
}

.h5ad_node <- function(h5, node, path, sample_n) {
    out <- list(kind = "other", class = NA_character_, dtype = NA_character_,
        shape = NA_character_, n_children = 0L, descend = FALSE, children = NULL,
        sample = NA_character_, value = NA_character_, flags = "")
    enc <- .h5ad_encoding(node)
    out$class <- enc
    named <- sub("/.*", "", path) %in% c("X", "layers", "raw", "obsm", "obsp")
    if (.h5ad_is_group(node)) {
        if (identical(enc, "dataframe")) {
            idx <- .h5ad_attr(node, "_index") %null% "_index"
            cols <- as.character(.h5ad_attr(node, "column-order") %null% character())
            n <- if (node$exists(idx)) node[[idx]]$dims else 0L
            out$kind <- "data.frame"
            out$shape <- paste(n, length(cols), sep = "x")
            out$n_children <- length(cols)
            out$descend <- TRUE
            out$children <- cols
            out$sample <- .pf_short(cols)
            num <- vapply(cols, function(cn) {
                cn <- node[[cn]]
                !.h5ad_is_group(cn) && is.null(.h5ad_attr(cn, "categories")) &&
                    .h5ad_dtype(cn) %in% c("integer", "double")
            }, logical(1))
            if (sum(num) >= 2L) out$flags <- "two_numeric_cols"
        } else if (enc %in% c("csr_matrix", "csc_matrix")) {
            out$kind <- "sparse_matrix"
            out$shape <- paste(.h5ad_attr(node, "shape"), collapse = "x")
            out$dtype <- .h5ad_dtype(node[["data"]])
            out$flags <- .pf_value_flags(.h5ad_sample(node[["data"]], sample_n), sample_n)
            if (named) out$flags <- .pf_add_flag(out$flags, "has_dimnames")
        } else if (identical(enc, "categorical")) {
            cats <- node[["categories"]]
            out$kind <- "factor"
            out$dtype <- "factor"
            out$shape <- as.character(node[["codes"]]$dims)
            out$sample <- .pf_short(cats[seq_len(min(3L, cats$dims))])
        } else if (startsWith(enc, "nullable-")) {
            v <- node[["values"]]
            out$kind <- "vector"
            out$dtype <- .h5ad_dtype(v)
            out$shape <- as.character(v$dims)
        } else if (enc %in% c("dict", "raw", "anndata")) {
            out$kind <- "list"
            out$n_children <- length(setdiff(names(node), "__categories"))
            out$descend <- TRUE
            out$sample <- .pf_short(names(node))
        }
        return(out)
    }
    # datasets
    dims <- node$dims
    out$dtype <- .h5ad_dtype(node)
    if (!is.null(.h5ad_attr(node, "categories"))) {
        out$kind <- "factor"
        out$dtype <- "factor"
        out$shape <- as.character(dims)
        return(out)
    }
    if (length(dims) == 2L) {
        out$kind <- "dense_matrix"
        out$shape <- paste(rev(dims), collapse = "x")
        if (out$dtype %in% c("integer", "double", "logical")) {
            out$flags <- .pf_value_flags(.h5ad_sample(node, sample_n), sample_n)
        }
        if (named) out$flags <- .pf_add_flag(out$flags, "has_dimnames")
        return(out)
    }
    out$kind <- "vector"
    n <- if (length(dims)) prod(dims) else 1L
    out$shape <- as.character(n)
    # null dataspaces (e.g. None in uns) cannot be read
    head <- tryCatch(if (!length(dims)) node$read() else if (n) node[seq_len(min(n, 3L))],
        error = function(e) NULL)
    if (is.null(head)) {
        out$kind <- "null"
        return(out)
    }
    if (is.factor(head)) head <- as.logical(as.character(head))
    if (n == 1L && length(head) == 1L && !is.na(head)) out$value <- as.character(head)
    out$sample <- .pf_short(head)
    if (out$dtype %in% c("integer", "double", "logical")) {
        out$flags <- .pf_value_flags(.h5ad_sample(node, sample_n), sample_n)
    }
    out
}

.h5ad_dtype <- function(node) {
    switch(class(node$get_type())[1L],
        H5T_INTEGER = "integer", H5T_FLOAT = "double", H5T_STRING = "character",
        H5T_LOGICAL = "logical", H5T_ENUM = "logical", H5T_COMPOUND = "compound",
        "other")
}

# about n values from ten evenly spaced contiguous runs (one hyperslab each)
.h5ad_sample <- function(node, n) {
    dims <- node$dims
    total <- prod(dims)
    if (!total) return(NULL)
    if (length(dims) == 2L) {
        # dims reversed: columns are obs; read whole obs rows
        per <- max(1L, as.integer(n %/% (10 * dims[1L])))
        starts <- unique(as.integer(seq(1, max(1, dims[2L] - per + 1), length.out = 10L)))
        v <- unlist(lapply(starts, function(s) node[, s:min(dims[2L], s + per - 1L)]))
    } else {
        per <- max(1L, as.integer(n %/% 10))
        if (total <= n) return(as.vector(node$read()))
        starts <- unique(as.integer(seq(1, total - per + 1, length.out = 10L)))
        v <- unlist(lapply(starts, function(s) node[s:(s + per - 1L)]))
    }
    if (is.factor(v)) v <- as.logical(as.character(v))
    as.numeric(v)
}
