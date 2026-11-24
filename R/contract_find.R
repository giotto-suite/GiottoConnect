#' @title Read a mapping file
#' @name readMapping
#' @description Read a YAML mapping file (or a mapping already parsed into a
#' list) and validate it against the mapping schema. All schema problems are
#' reported together, each with the YAML path it concerns.
#' @param x path to a YAML mapping file, the id of a mapping bundled with
#' GiottoConnect (e.g. `"seurat-assay5"`), or a parsed mapping list
#' @returns a `giottoMapping` (the parsed mapping as a list)
#' @examples
#' f <- system.file("mappings", package = "GiottoConnect")
#' list.files(f, pattern = "\\.yaml$")
#' @export
readMapping <- function(x) {
    file <- NA_character_
    if (is.character(x) && length(x) == 1L) {
        if (!file.exists(x)) {
            bundled <- system.file("mappings", paste0(x, ".yaml"),
                package = "GiottoConnect"
            )
            if (!nzchar(bundled)) {
                stop(wrap_txt(sprintf(
                    "No mapping file or bundled mapping named '%s'", x
                )), call. = FALSE)
            }
            x <- bundled
        }
        file <- normalizePath(x)
        m <- .read_yaml_strict(file)
    } else if (is.list(x)) {
        m <- unclass(x)
        file <- attr(x, "file") %null% NA_character_
    } else {
        stop(wrap_txt("`x` must be a file path, a mapping id or a parsed mapping"),
            call. = FALSE)
    }
    probs <- .mapping_schema_check(m)
    if (length(probs)) {
        src <- if (is.na(file)) "mapping" else basename(file)
        stop(sprintf("%s does not match the mapping schema:\n%s", src,
            paste0("  - ", probs, collapse = "\n")), call. = FALSE)
    }
    structure(m, class = "giottoMapping", file = file)
}

# YAML 1.1 reads `on`, `off`, `y`, `n`, `yes`, `no` as booleans (also as
# keys); mappings only treat `true`/`false` as booleans
.read_yaml_strict <- function(file) {
    keep <- function(x) if (tolower(x) %in% c("true", "false")) as.logical(x) else x
    yaml::read_yaml(file, handlers = list("bool#yes" = keep, "bool#no" = keep))
}

#' @export
print.giottoMapping <- function(x, ...) {
    h <- x$mapping
    cat(sprintf("<giottoMapping> %s v%s [%s]\n", h$id, h$version, h$status))
    cat(sprintf("  source: %s %s (%s), built with %s\n", h$source$ecosystem,
        h$source$class, h$source$package, h$source$built_with %null% "?"))
    cat(sprintf("  target: %s contract %s\n", h$target$contract, h$target$version))
    cat(sprintf("  params: %s\n", paste(names(x$params), collapse = ", ")))
    parts <- vapply(x$parts, function(p) sprintf("%s (%s)", p$id, p$part), character(1))
    cat(sprintf("  parts : %s\n", paste(parts, collapse = ", ")))
    n_loss <- length(x$lossy) + sum(vapply(x$parts, function(p) length(p$lossy), integer(1)))
    cat(sprintf("  declared losses: %d\n", n_loss))
    invisible(x)
}

#' @title Summarise a mapping
#' @name summary.giottoMapping
#' @description List a mapping's parameters, parts and every declared loss.
#' @param object a `giottoMapping`
#' @param ... unused
#' @returns invisibly, a list with `parts` and `losses` data.frames
#' @export
summary.giottoMapping <- function(object, ...) {
    parts <- do.call(rbind, lapply(object$parts, function(p) data.frame(
        id = p$id, part = p$part, from = p$from,
        via = p$hook %null% paste(vapply(p$ops, `[[`, character(1), "op"), collapse = ","),
        optional = isTRUE(p$optional),
        to = paste(names(p$to), unlist(p$to), sep = "=", collapse = " "),
        stringsAsFactors = FALSE
    )))
    loss_rows <- function(ls, part_id) lapply(ls, function(l) data.frame(
        part_id = part_id, kind = l$kind, what = l$path %null% l$field %null% "",
        note = l$note, stringsAsFactors = FALSE
    ))
    losses <- do.call(rbind, c(
        unlist(lapply(object$parts, function(p) loss_rows(p$lossy, p$id)), recursive = FALSE),
        loss_rows(object$lossy, "(mapping)")
    ))
    print(object)
    cat("\nParts:\n")
    print(parts, row.names = FALSE, right = FALSE)
    cat("\nDeclared losses:\n")
    if (is.null(losses)) cat("  none\n") else print(losses, row.names = FALSE, right = FALSE)
    invisible(list(parts = parts, losses = losses))
}

#' @title Find the mapping for an object
#' @name findMapping
#' @description Choose the mapping that applies to a profiled object. Mappings
#' are filtered by their `applies_to` declaration (class, package version and
#' structural predicates after parameter binding); among those that apply, a
#' fingerprint match is preferred, then the more specific mapping (more
#' `applies_to` predicates satisfied), then the one whose source paths
#' (optional parts included) are best covered; without an exact match the status is `"repair_needed"`.
#' Fingerprints cover the mapping's required parts: an optional part that is
#' absent from the object does not change it.
#' @param profile `giottoProfile` of the object
#' @param mappings list of mappings (`giottoMapping`, files or ids). Defaults
#' to all bundled mappings.
#' @param ecosystem optional ecosystem filter, e.g. `"seurat"`
#' @param target optional target filter (`mapping.target.contract`), e.g.
#' `"seurat"` for mappings from giotto objects to Seurat
#' @returns list with `status` (`"match"`, `"repair_needed"` or `"none"`),
#' `mapping`, `params`, `fingerprint`, `candidates` (one row per mapping
#' considered) and `reasons`
#' @export
findMapping <- function(profile, mappings = NULL, ecosystem = NULL, target = NULL) {
    checkmate::assert_class(profile, "giottoProfile")
    if (is.null(mappings)) {
        mappings <- list.files(system.file("mappings", package = "GiottoConnect"),
            pattern = "\\.yaml$", full.names = TRUE)
    }
    mappings <- lapply(mappings, function(m) if (inherits(m, "giottoMapping")) m else readMapping(m))
    if (!is.null(ecosystem)) {
        mappings <- Filter(function(m) m$mapping$source$ecosystem == ecosystem, mappings)
    }
    if (!is.null(target)) {
        mappings <- Filter(function(m) identical(m$mapping$target$contract, target), mappings)
    }
    rows <- list()
    res <- list()
    for (m in mappings) {
        why <- .applies_to(m, profile)
        params <- if (is.null(why$params)) list(values = list()) else why$params
        fp <- if (why$ok) .mapping_fingerprint(m, profile, params$values) else NA_character_
        # how much of what the mapping can read is in the object
        paths <- .pf_bind(.mapping_paths(m, optional = TRUE), params$values)
        cover <- mean(paths %in% profile$nodes$path)
        status <- if (!why$ok) "none" else if (identical(fp, m$mapping$fingerprint)) "match" else "repair_needed"
        rows[[length(rows) + 1L]] <- data.frame(id = m$mapping$id, version = m$mapping$version,
            status = status, coverage = cover,
            specificity = length(m$mapping$source$applies_to$predicates),
            reason = why$reason, stringsAsFactors = FALSE)
        res[[length(res) + 1L]] <- list(mapping = m, params = params, fingerprint = fp)
    }
    cand <- if (length(rows)) do.call(rbind, rows) else data.frame(id = character(),
        version = character(), status = character(), coverage = numeric(),
        specificity = integer(), reason = character())
    ok <- which(cand$status != "none")
    if (!length(ok)) {
        return(list(status = "none", mapping = NULL, params = NULL, fingerprint = NA_character_,
            candidates = cand, reasons = if (nrow(cand)) paste(cand$id, cand$reason, sep = ": ") else "no mappings"))
    }
    # an exact match first, then the more specific mapping (more satisfied
    # applies_to predicates), then the best covered one
    best <- ok[order(cand$status[ok] != "match", -cand$specificity[ok], -cand$coverage[ok])][1L]
    c(list(status = cand$status[best]), res[[best]],
        list(candidates = cand, reasons = cand$reason[best]))
}

# class, package version and structural predicates
.applies_to <- function(m, profile) {
    src <- m$mapping$source
    no <- function(r) list(ok = FALSE, reason = r, params = NULL)
    if (!identical(profile$source$class, src$class)) {
        return(no(sprintf("class %s, mapping needs %s", profile$source$class, src$class)))
    }
    at <- src$applies_to
    if (!is.null(at$package_version) && !is.na(profile$source$version)) {
        if (!.version_ok(profile$source$version, at$package_version)) {
            return(no(sprintf("%s %s does not satisfy %s", src$package, profile$source$version,
                at$package_version)))
        }
    }
    params <- .bind_params(m, profile)
    for (pr in at$predicates) {
        path <- .pf_bind(pr$path, params$values)
        if (grepl("{", path, fixed = TRUE)) {
            miss <- intersect(params$unbound, gsub("[{}]", "", regmatches(path,
                gregexpr("\\{[^}]+\\}", path))[[1L]]))
            return(no(sprintf("param '%s' could not be bound (%s)", miss[1L],
                .param_describe(m$params[[miss[1L]]]))))
        }
        r <- profile$nodes[profile$nodes$path == path, , drop = FALSE]
        if (!nrow(r)) return(no(sprintf("no node at %s", path)))
        if (!is.null(pr$class) && r$class != pr$class) {
            return(no(sprintf("%s is %s, mapping needs %s", path, r$class, pr$class)))
        }
        if (!is.null(pr$kind) && r$kind != pr$kind) {
            return(no(sprintf("%s is %s, mapping needs %s", path, r$kind, pr$kind)))
        }
    }
    list(ok = TRUE, reason = "applies", params = params)
}

.param_describe <- function(pr) {
    w <- pr$where
    cond <- c(
        if (!is.null(w$class)) sprintf("class %s", w$class),
        if (!is.null(w$role)) sprintf("role %s", paste(unlist(w$role), collapse = "+")),
        if (!is.null(w$has_child_role)) sprintf("data that is %s", paste(unlist(w$has_child_role), collapse = "+")),
        if (!is.null(w$name_matches)) sprintf("name ~ %s", w$name_matches),
        if (!is.null(pr$class)) sprintf("class %s", pr$class),
        if (!is.null(pr$value)) sprintf("named %s", pr$value)
    )
    sprintf("%s %s%s", pr$list, pr$rule, if (length(cond)) paste0(": ", paste(cond, collapse = ", ")) else "")
}

# constraints like ">= 5.0.0" or ">= 5.0.0, < 6"
.version_ok <- function(version, constraint) {
    v <- numeric_version(version)
    for (c in trimws(strsplit(constraint, ",", fixed = TRUE)[[1L]])) {
        op <- regmatches(c, regexpr("^(>=|<=|==|!=|>|<)", c))
        ref <- numeric_version(trimws(sub("^(>=|<=|==|!=|>|<)", "", c)))
        if (!length(op)) op <- "=="
        if (!do.call(op, list(v, ref))) return(FALSE)
    }
    TRUE
}

# profiles stored next to a mapping, for drift detection
.profile_write <- function(profile, file) {
    jsonlite::write_json(unclass(profile), file, auto_unbox = TRUE, pretty = TRUE,
        na = "null", digits = NA)
    invisible(file)
}

.profile_read <- function(file) {
    p <- jsonlite::read_json(file, simplifyVector = TRUE)
    p$nodes <- as.data.frame(p$nodes, stringsAsFactors = FALSE)
    structure(p, class = "giottoProfile")
}
