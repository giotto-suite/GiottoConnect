# Repair a mapping for an object whose structure drifted from the mapping's
# stored profile. Order: re-bind params (renamed assays, layers, graphs need
# nothing else), then rules (moved/renamed nodes, storage changes, declared
# alternatives), then ranked candidates for anything still broken (input for
# an LLM or a human, see repairMapping(llm =) and writeMappingRequest()).
# The engine never repairs on its own: a repaired mapping is a new version
# with status patched-unreviewed.

# listing verbs whose values are object names embedded in other names
.OBJECT_NAME_LISTS <- c("Assays", "Images", "Reductions", "Graphs", "Neighbors")

#' @title Repair a mapping for a drifted object
#' @name repairMapping
#' @description Adapt a mapping to an object whose structure differs from the
#' one the mapping was verified on. Parameters are re-bound first; then rules
#' rewrite moved or renamed source paths, add a `cast` when only the storage
#' of a node changed, and switch a part to one of its declared alternatives
#' when its source is gone. Whatever is still broken gets up to five ranked
#' candidate paths. A repaired mapping is validated (and, when `x` is the
#' object, converted and read back) before it is accepted.
#' @param mapping a `giottoMapping`, mapping file or bundled mapping id
#' @param x the drifted object, or its `giottoProfile`
#' @param llm optional LLM backend used for parts the rules cannot repair
#' @param old_profile profile the mapping was verified on. Defaults to the
#' `<id>.profile.json` stored next to the mapping file.
#' @param file optional path to write the repaired mapping to (YAML)
#' @returns a `giottoRepair`: list with `status` (`"not_needed"`,
#' `"repaired"` or `"failed"`), `mapping` (the possibly patched mapping),
#' `changes`, `unresolved`, `candidates`, `validation` and `readback`
#' @export
repairMapping <- function(mapping, x, llm = NULL, old_profile = NULL, file = NULL) {
    mapping <- if (inherits(mapping, "giottoMapping")) mapping else readMapping(mapping)
    is_obj <- !inherits(x, "giottoProfile")
    new <- if (is_obj) profileObject(x) else x
    old <- old_profile %null% .mapping_stored_profile(mapping)

    m <- mapping
    pn <- .bind_params(m, new)
    po <- if (is.null(old)) NULL else .bind_params(m, old)
    diff <- if (is.null(old)) NULL else diffProfiles(old, new,
        params_old = .object_params(m, po), params_new = .object_params(m, pn))

    changes <- list()
    note <- function(part_id, rule, detail) {
        changes[[length(changes) + 1L]] <<- data.frame(part_id = part_id, rule = rule,
            detail = detail, stringsAsFactors = FALSE)
    }
    unresolved <- character()

    for (i in seq_along(m$parts)) {
        p <- m$parts[[i]]
        node <- m$source_nodes[[p$from]]
        cur <- .pf_bind(node$path, pn$values)
        ok_new <- .part_bound(p, node, pn) && cur %in% new$nodes$path
        oldpath <- if (is.null(po)) NA_character_ else .pf_bind(node$path, po$values)
        in_old <- !is.null(old) && oldpath %in% old$nodes$path

        if (ok_new) {
            cast <- if (in_old) .rule_cast(old, new, oldpath, cur, p) else NULL
            if (!is.null(cast)) {
                m$parts[[i]]$ops <- c(list(cast$op), p$ops)
                note(p$id, "cast", cast$detail)
            }
            next
        }
        if (!in_old && isTRUE(p$optional) && !length(p$alternatives)) next
        if (!in_old && isTRUE(p$optional)) {
            # absent before as well: only switch if an alternative now binds
            alt <- .rule_alternative(m, p, new)
            if (!is.null(alt)) {
                m <- .apply_alternative(m, i, alt)
                note(p$id, "alternative", alt$detail)
            }
            next
        }

        mv <- if (is.null(diff)) NULL else .rule_moved(diff, node, oldpath, po, pn, m)
        if (!is.null(mv)) {
            m$source_nodes[[p$from]] <- mv$node
            note(p$id, mv$rule, mv$detail)
            next
        }
        alt <- .rule_alternative(m, p, new)
        if (!is.null(alt)) {
            m <- .apply_alternative(m, i, alt)
            note(p$id, "alternative", alt$detail)
            next
        }
        if (isTRUE(p$optional)) {
            # optional data that is gone from the object: skipped, not repaired
            note(p$id, "absent", sprintf("%s is gone and nothing replaces it; optional part skipped",
                oldpath))
            next
        }
        unresolved <- c(unresolved, p$id)
    }
    if (!is.null(diff)) {
        # new nodes the (possibly patched) mapping now reads are not "unmapped"
        used <- unlist(.bind_params(m, new)$values)
        added <- diff$new_path[diff$change == "added"]
        added <- added[!sub(".*/", "", added) %in% used]
        for (pth in added) {
            note(NA_character_, "added", sprintf("%s is new and unmapped (no change)", pth))
        }
    }
    candidates <- .repair_candidates(m, unresolved, old, new, po, pn)
    llm_log <- NULL
    if (length(unresolved) && !is.null(llm)) {
        lp <- .repair_with_llm(m, unresolved, candidates, diff, llm)
        m <- lp$mapping
        for (r in seq_len(nrow(lp$applied))) note(lp$applied$part_id[r], "llm", lp$applied$detail[r])
        unresolved <- setdiff(unresolved, lp$applied$part_id)
        llm_log <- lp$log
    }
    changes <- if (length(changes)) do.call(rbind, changes) else
        data.frame(part_id = character(), rule = character(), detail = character())

    patched <- any(changes$rule %in% c("moved", "renamed_key", "cast", "alternative", "llm"))
    pn2 <- .bind_params(m, new)
    validation <- validateMapping(m, new, params = pn2, checks = "structural")
    readback <- NULL
    ok <- !length(unresolved) && !any(validation$level == "FAIL")
    if (ok && is_obj) {
        g <- tryCatch(suppressWarnings(convertWithMapping(x, m, profile = new)),
            error = function(e) e)
        if (inherits(g, "error")) {
            ok <- FALSE
            validation <- rbind(validation, data.frame(check = "convert", part_id = NA,
                level = "FAIL", message = conditionMessage(g)))
        } else {
            readback <- readbackCheck(x, g, m)
            ok <- isTRUE(attr(readback, "pass"))
        }
    }

    by <- if (any(changes$rule == "llm")) "llm" else "rules"
    status <- if (!ok) "failed" else if (patched) "repaired" else "not_needed"
    if (status == "repaired") {
        m <- .mapping_bump(m, by = by, change = paste(
            sprintf("%s: %s", changes$part_id, changes$rule)[changes$rule != "added"], collapse = "; "))
        m$mapping$fingerprint <- .mapping_fingerprint(m, new, pn2$values)
    } else if (status == "failed") {
        m$mapping$status <- "failed"
    }
    if (!is.null(file) && status != "not_needed") .mapping_write(m, file)

    structure(list(status = status, mapping = m, changes = changes, unresolved = unresolved,
        candidates = candidates, validation = validation, readback = readback,
        diff = diff, llm_log = llm_log), class = "giottoRepair")
}

#' @export
print.giottoRepair <- function(x, ...) {
    h <- x$mapping$mapping
    cat(sprintf("<giottoRepair> %s -> %s v%s [%s]\n", x$status, h$id, h$version, h$status))
    if (nrow(x$changes)) {
        ch <- x$changes
        ch$part_id[is.na(ch$part_id)] <- ""
        print(ch, row.names = FALSE, right = FALSE)
    }
    if (length(x$unresolved)) {
        cat("unresolved:", paste(x$unresolved, collapse = ", "), "\n")
        if (nrow(x$candidates)) print(x$candidates, row.names = FALSE, right = FALSE)
    }
    if (!is.null(x$readback)) {
        cat(sprintf("readback: %s\n", if (isTRUE(attr(x$readback, "pass"))) "PASS" else "FAIL"))
    }
    invisible(x)
}


# rules ####

# storage changed but data did not: dense <-> sparse, integer <-> double
.rule_cast <- function(old, new, oldpath, newpath, part) {
    a <- old$nodes[old$nodes$path == oldpath, ]
    b <- new$nodes[new$nodes$path == newpath, ]
    if (!nrow(a) || !nrow(b)) return(NULL)
    roles <- strsplit(b$role_hints %na% "", ",", fixed = TRUE)[[1L]]
    int_ok <- "integer_valued" %in% roles
    if (a$kind == "sparse_matrix" && b$kind == "dense_matrix" && part$part == "expression") {
        return(list(op = list(op = "cast", to = "sparse"),
            detail = sprintf("%s is stored dense (was sparse); cast to sparse", newpath)))
    }
    if (!identical(a$dtype, b$dtype) && all(c(a$dtype, b$dtype) %in% c("integer", "double"))) {
        if (a$dtype == "integer" && !int_ok) return(NULL)
        return(list(op = list(op = "cast", to = a$dtype),
            detail = sprintf("%s is %s (was %s); cast to %s", newpath, b$dtype, a$dtype, a$dtype)))
    }
    NULL
}

# a removed subtree that reappears elsewhere with the same structure. Only
# nodes whose accessor targets a path (`get.on`) can follow a move; accessors
# that address data by name (LayerData, Embeddings) are handled by params.
.rule_moved <- function(diff, node, oldpath, po, pn, m) {
    if (is.null(node$get$on)) return(NULL)
    obj <- .object_params(m, po)
    told <- .pf_template(oldpath, obj)
    hit <- diff[diff$change %in% c("moved", "renamed_key") &
        (diff$path == told | startsWith(told, paste0(diff$path, "/"))), , drop = FALSE]
    if (nrow(hit) != 1L) return(NULL)
    new_t <- paste0(hit$new_path, substring(told, nchar(hit$path) + 1L))
    # put the mapping's own placeholders back for the parts that are params
    new_t <- .pf_template(.pf_bind(new_t, .object_params(m, pn)), .object_params(m, pn))
    node$path <- new_t
    if (!is.null(node$get$on) && startsWith(.pf_bind(node$get$on, po$values), hit$path)) {
        node$get$on <- paste0(hit$new_path, substring(.pf_bind(node$get$on, po$values),
            nchar(hit$path) + 1L))
    }
    list(node = node, rule = hit$change,
        detail = sprintf("%s -> %s", oldpath, .pf_bind(new_t, pn$values)))
}

.rule_alternative <- function(m, p, new) {
    for (alt in p$alternatives) {
        mm <- m
        mm$params <- c(m$params, alt$params)
        b <- .bind_params(mm, new)
        if (length(intersect(names(alt$params), b$unbound))) next
        node <- m$source_nodes[[alt$from]]
        path <- .pf_bind(node$path, b$values)
        if (grepl("{", path, fixed = TRUE) || !path %in% new$nodes$path) next
        bound <- unlist(b$values[names(alt$params)])
        return(list(alt = alt, detail = sprintf("switched to alternative '%s' (%s)", alt$id,
            paste(names(bound), bound, sep = " = ", collapse = ", "))))
    }
    NULL
}

# The alternative becomes the part's source. The outgoing source is kept as
# an alternative (named after its source node), together with the params only
# it used, so objects with the earlier structure still repair back.
.apply_alternative <- function(m, i, a) {
    p <- m$parts[[i]]
    alt <- a$alt
    node <- m$source_nodes[[p$from]]
    own <- .part_placeholders(p, node)
    others <- unlist(lapply(m$parts[-i], function(q) .part_placeholders(q, m$source_nodes[[q$from]])))
    own <- intersect(setdiff(own, others), names(m$params))
    stash <- list(id = sprintf("%s_from_%s", p$id, p$from), from = p$from, hook = p$hook,
        params = if (length(own)) m$params[own] else NULL, ops = p$ops, lossy = p$lossy)
    stash <- stash[!vapply(stash, is.null, logical(1))]

    p$from <- alt$from
    p$hook <- alt$hook
    p$ops <- alt$ops
    p$lossy <- alt$lossy
    if (!is.null(alt$to)) p$to <- utils::modifyList(p$to, alt$to)
    p$alternatives <- c(Filter(function(x) x$id != alt$id && x$id != stash$id, p$alternatives),
        list(stash))
    p <- p[!vapply(p, is.null, logical(1))]
    m$parts[[i]] <- p
    m$params[own] <- NULL
    for (nm in names(alt$params)) m$params[[nm]] <- alt$params[[nm]]
    m
}

.part_placeholders <- function(p, node) {
    txt <- c(node$path, node$get$on, unlist(node$get$args), unlist(p$ops),
        names(unlist(lapply(p$ops, `[[`, "map"))))
    txt <- txt[is.character(txt)]
    unique(gsub("[{}]", "", unlist(regmatches(txt, gregexpr("\\{[^}]+\\}", txt)))))
}


# candidates ####

# up to five new-profile paths per unresolved part, ranked by same kind
# (required; same class for S4), matching size and value roles, name
# similarity and a shared parent
.repair_candidates <- function(m, unresolved, old, new, po, pn, n = 5L) {
    empty <- data.frame(part_id = character(), path = character(), kind = character(),
        score = numeric(), stringsAsFactors = FALSE)
    if (!length(unresolved)) return(empty)
    out <- list()
    for (pid in unresolved) {
        p <- m$parts[[match(pid, vapply(m$parts, `[[`, character(1), "id"))]]
        node <- m$source_nodes[[p$from]]
        ref <- NULL
        if (!is.null(old) && !is.null(po)) {
            ref <- old$nodes[old$nodes$path == .pf_bind(node$path, po$values), , drop = FALSE]
        }
        target_name <- sub(".*/", "", .pf_bind(node$path, (po %null% pn)$values))
        target_parent <- .pf_parent(.pf_bind(node$path, (po %null% pn)$values))
        nd <- new$nodes
        if (!is.null(ref) && nrow(ref)) {
            nd <- nd[nd$kind == ref$kind, , drop = FALSE]
            # S4 nodes are only interchangeable with the same class
            if (ref$kind == "S4") nd <- nd[nd$class == ref$class, , drop = FALSE]
        }
        nd <- nd[nzchar(nd$path) & !endsWith(nd$path, "/*"), , drop = FALSE]
        if (!nrow(nd)) next
        roles <- function(h, pat) {
            r <- strsplit(h %na% "", ",", fixed = TRUE)[[1L]]
            r[grepl(pat, r)]
        }
        size_pat <- "_eq_"
        value_pat <- "^(integer_valued|nonneg|binary_valued)$"
        has_ref <- !is.null(ref) && nrow(ref)
        size_ref <- if (has_ref) roles(ref$role_hints, size_pat) else character()
        value_ref <- if (has_ref) roles(ref$role_hints, value_pat) else character()
        # shape: overlap of size roles; content: share of the source's value roles
        role <- vapply(nd$role_hints, function(h) {
            r <- roles(h, size_pat)
            if (!length(size_ref) && !length(r)) 1 else
                length(intersect(r, size_ref)) / max(1L, length(union(r, size_ref)))
        }, numeric(1), USE.NAMES = FALSE)
        value <- vapply(nd$role_hints, function(h) {
            if (!length(value_ref)) 1 else mean(value_ref %in% roles(h, value_pat))
        }, numeric(1), USE.NAMES = FALSE)
        nm <- sub(".*/", "", nd$path)
        sim <- 1 - drop(utils::adist(tolower(target_name), tolower(nm))) /
            pmax(nchar(target_name), nchar(nm), 1L)
        parent <- as.numeric(.pf_parent(nd$path) == target_parent)
        score <- round(2 * role + 2 * value + sim + 0.5 * parent, 3)
        o <- order(-score)[seq_len(min(n, nrow(nd)))]
        out[[length(out) + 1L]] <- data.frame(part_id = pid, path = nd$path[o], kind = nd$kind[o],
            score = score[o], stringsAsFactors = FALSE)
    }
    if (length(out)) do.call(rbind, out) else empty
}


# helpers ####

.object_params <- function(m, params) {
    if (is.null(params)) return(list())
    keep <- names(m$params)[vapply(m$params, function(pr) pr$list %in% .OBJECT_NAME_LISTS, logical(1))]
    params$values[intersect(keep, names(params$values))]
}

.mapping_stored_profile <- function(mapping) {
    f <- attr(mapping, "file")
    if (is.null(f) || is.na(f)) return(NULL)
    pf <- sub("\\.ya?ml$", ".profile.json", f)
    if (file.exists(pf)) .profile_read(pf) else NULL
}

.mapping_bump <- function(m, by, change) {
    old <- m$mapping$version
    v <- as.integer(strsplit(old, ".", fixed = TRUE)[[1L]])
    v[3L] <- v[3L] + 1L
    m$mapping$version <- paste(v, collapse = ".")
    m$mapping$status <- "patched-unreviewed"
    m$mapping$history <- c(m$mapping$history, list(list(date = format(Sys.Date()),
        from_version = old, change = change, by = by)))
    m
}

# booleans as true/false: YAML 1.1 yes/no would not read back
.mapping_write <- function(m, file) {
    tf <- function(x) structure(ifelse(x, "true", "false"), class = "verbatim")
    yaml::write_yaml(unclass(m), file, handlers = list(logical = tf))
    invisible(file)
}

# LLM patches for unresolved parts (contract_llm.R)
.repair_with_llm <- function(m, unresolved, candidates, diff, llm) {
    stop(wrap_txt("LLM repair is not available yet"), call. = FALSE)
}
