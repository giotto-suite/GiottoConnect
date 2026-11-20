# Mapping files are validated against inst/schema/mapping.schema.json by a
# small hand-written checker. It supports only the JSON Schema constructs the
# mapping schema uses: type, required, enum, properties, items,
# additionalProperties (false or a schema) and local $ref.

# arguments each op requires / accepts
.OP_ARGS <- list(
    select_cols = list(required = "cols"),
    exclude_cols = list(required = "cols"),
    rename = list(required = "map"),
    transpose = list(required = character()),
    scale = list(required = "factor", optional = "cols"),
    cast = list(required = "to"),
    rownames_to_col = list(required = "col"),
    set_colnames = list(required = "cols"),
    set = list(required = c("field", "value")),
    replace = list(required = c("cols", "pattern", "replacement"))
)

.mapping_schema <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            f <- system.file("schema", "mapping.schema.json",
                package = "GiottoConnect", mustWork = TRUE
            )
            cache <<- jsonlite::read_json(f, simplifyVector = FALSE)
        }
        cache
    }
})

# Validate a parsed mapping. Returns a character vector of problems, each
# prefixed with the YAML path it concerns. Empty when valid.
.mapping_schema_check <- function(m, schema = .mapping_schema()) {
    probs <- .schema_check(m, schema, "", schema)
    c(probs, .mapping_ops_check(m), .mapping_adapter_check(m))
}

# the source ecosystem and target must be registered, and every accessor a
# source node calls must be one its adapter allows
.mapping_adapter_check <- function(m) {
    eco <- m$mapping$source$ecosystem
    con <- m$mapping$target$contract
    probs <- character()
    if (is.character(con) && length(con) == 1L && is.null(.TARGETS[[con]])) {
        probs <- c(probs, sprintf("mapping.target.contract: no target '%s' (available: %s)", con,
            paste(names(.TARGETS), collapse = ", ")))
    }
    if (!is.character(eco) || length(eco) != 1L) return(probs)
    if (is.null(.ADAPTERS[[eco]])) {
        return(c(probs, sprintf("mapping.source.ecosystem: no adapter for '%s' (available: %s)", eco,
            paste(names(.ADAPTERS), collapse = ", "))))
    }
    adapter <- .adapter_for(eco)
    verbs <- names(adapterListPaths(adapter))
    # malformed entries are reported by the schema check
    objs <- function(x) Filter(is.list, x)
    lists <- c(
        lapply(objs(m$params), function(pr) list(at = "params", v = pr$list)),
        unlist(lapply(objs(m$parts), function(p) c(
            if (is.list(p$for_each)) list(list(at = sprintf("parts.%s.for_each", p$id), v = p$for_each$list)),
            unlist(lapply(objs(p$alternatives), function(a) lapply(objs(a$params), function(pr) {
                list(at = sprintf("parts.%s.alternatives.%s.params", p$id, a$id), v = pr$list)
            })), recursive = FALSE)
        )), recursive = FALSE)
    )
    for (l in lists) {
        if (is.character(l$v) && length(l$v) == 1L && !l$v %in% verbs) {
            probs <- c(probs, sprintf("%s: list '%s' is not a listing of the %s adapter (allowed: %s)",
                l$at, l$v, eco, paste(verbs, collapse = ", ")))
        }
    }
    allowed <- adapterAccessors(adapter)
    for (nm in names(m$source_nodes)) {
        fn <- m$source_nodes[[nm]]$get$fn
        if (is.character(fn) && length(fn) == 1L && !fn %in% allowed) {
            probs <- c(probs, sprintf("source_nodes.%s.get.fn: '%s' is not an accessor of the %s adapter (allowed: %s)",
                nm, fn, eco, paste(allowed, collapse = ", ")))
        }
    }
    probs
}

.schema_check <- function(x, s, path, root) {
    if (!is.null(s[["$ref"]])) {
        s <- .schema_resolve(s[["$ref"]], root)
    }
    if (!length(s)) return(character())
    p <- if (nzchar(path)) path else "<root>"

    if (!is.null(s$type)) {
        types <- unlist(s$type)
        if (!any(vapply(types, .schema_is_type, logical(1), x = x))) {
            return(sprintf("%s: expected %s, got %s", p,
                paste(types, collapse = " or "), .schema_describe(x)))
        }
    }
    probs <- character()
    if (!is.null(s$enum)) {
        allowed <- unlist(s$enum)
        bad <- setdiff(as.character(unlist(x)), allowed)
        if (length(bad)) {
            probs <- c(probs, sprintf("%s: '%s' is not one of: %s", p,
                bad[1L], paste(allowed, collapse = ", ")))
        }
    }
    if (.schema_is_type("object", x) && length(x) && !is.null(names(x))) {
        for (req in unlist(s$required)) {
            if (is.null(x[[req]])) {
                probs <- c(probs, sprintf("%s: missing required field '%s'", p, req))
            }
        }
        props <- s$properties %null% list()
        extra <- setdiff(names(x), names(props))
        ap <- s$additionalProperties
        if (length(extra) && isFALSE(ap)) {
            probs <- c(probs, sprintf("%s: unknown field '%s' (allowed: %s)", p,
                extra, paste(names(props), collapse = ", ")))
        }
        for (nm in names(x)) {
            sub <- if (nm %in% names(props)) props[[nm]] else if (is.list(ap)) ap else NULL
            if (!is.null(sub)) {
                probs <- c(probs, .schema_check(x[[nm]], sub, .schema_path(path, nm), root))
            }
        }
    } else if (.schema_is_type("object", x) && length(unlist(s$required))) {
        probs <- c(probs, sprintf("%s: missing required field '%s'", p,
            unlist(s$required)[1L]))
    }
    if (!is.null(s$items) && .schema_is_type("array", x)) {
        items <- if (is.list(x)) x else as.list(x)
        for (i in seq_along(items)) {
            probs <- c(probs, .schema_check(items[[i]], s$items,
                .schema_index(path, i, items[[i]]), root))
        }
    }
    probs
}

.schema_resolve <- function(ref, root) {
    if (!startsWith(ref, "#/")) stop("only local $ref is supported: ", ref)
    keys <- strsplit(substring(ref, 3L), "/", fixed = TRUE)[[1L]]
    node <- root
    for (k in keys) node <- node[[k]]
    if (is.null(node)) stop("unresolvable $ref: ", ref)
    node
}

# YAML scalars arrive as length-1 atomics, sequences as unnamed lists or
# atomic vectors, mappings as named lists
.schema_is_type <- function(type, x) {
    switch(type,
        string = is.character(x) && length(x) == 1L,
        number = is.numeric(x) && length(x) == 1L,
        integer = is.numeric(x) && length(x) == 1L && x == round(x),
        boolean = is.logical(x) && length(x) == 1L && !is.na(x),
        object = is.list(x) && (!length(x) ||
            (!is.null(names(x)) && all(nzchar(names(x))))),
        array = (is.list(x) && (!length(x) || is.null(names(x)))) ||
            (is.atomic(x) && !is.null(x)),
        null = is.null(x),
        FALSE
    )
}

.schema_describe <- function(x) {
    if (is.null(x)) return("null")
    if (is.list(x)) return(if (is.null(names(x))) "a sequence" else "a mapping")
    if (length(x) != 1L) return(sprintf("a sequence of %d %s", length(x), typeof(x)))
    sprintf("%s '%s'", typeof(x), as.character(x))
}

.schema_path <- function(path, name) {
    if (nzchar(path)) paste(path, name, sep = ".") else name
}

# index into a sequence, labelled with the element id when it has one
.schema_index <- function(path, i, el) {
    lab <- sprintf("%s[%d]", path, i)
    if (is.list(el) && is.character(el$id) && length(el$id) == 1L) {
        lab <- sprintf("%s(id=%s)", lab, el$id)
    }
    lab
}

# per-op argument rules the schema cannot express without oneOf
.mapping_ops_check <- function(m) {
    probs <- character()
    check_ops <- function(ops, path) {
        if (!is.list(ops)) return(invisible())
        for (i in seq_along(ops)) {
            o <- ops[[i]]
            if (!is.list(o) || !is.character(o$op)) next
            rule <- .OP_ARGS[[o$op]]
            if (is.null(rule)) next
            p <- sprintf("%s[%d](op=%s)", path, i, o$op)
            args <- setdiff(names(o), "op")
            miss <- setdiff(rule$required, args)
            extra <- setdiff(args, c(rule$required, rule$optional))
            if (length(miss)) {
                probs <<- c(probs, sprintf("%s: missing argument '%s'", p, miss))
            }
            if (length(extra)) {
                probs <<- c(probs, sprintf("%s: op '%s' does not take '%s'", p,
                    o$op, extra))
            }
        }
    }
    for (i in seq_along(m$parts)) {
        part <- m$parts[[i]]
        if (!is.list(part)) next
        pp <- .schema_index("parts", i, part)
        check_ops(part$ops, paste0(pp, ".ops"))
        for (j in seq_along(part$alternatives)) {
            alt <- part$alternatives[[j]]
            if (!is.list(alt)) next
            check_ops(alt$ops, paste0(.schema_index(paste0(pp, ".alternatives"), j, alt), ".ops"))
        }
    }
    probs
}
