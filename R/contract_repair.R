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
#' @param llm optional [llmBackend()]. Parts the rules cannot repair, and
#' declared alternatives whose params do not bind by rule, are offered to it as
#' closed choices (see [llmPatch()]). Its choices pass the same validation,
#' conversion and readback as rule repairs.
#' @param old_profile profile the mapping was verified on. Defaults to the
#' `<id>.profile.json` stored next to the mapping file.
#' @param file optional path to write the repaired mapping to (YAML)
#' @param request_file optional path for a [writeMappingRequest()] handoff,
#' written when the repair fails
#' @returns a `giottoRepair`: list with `status` (`"not_needed"`,
#' `"repaired"` or `"failed"`), `mapping` (the possibly patched mapping),
#' `changes`, `unresolved`, `candidates`, `validation` and `readback`
#' @export
repairMapping <- function(mapping, x, llm = NULL, old_profile = NULL, file = NULL,
    request_file = NULL) {
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
    pending <- list()

    for (i in seq_along(m$parts)) {
        p <- m$parts[[i]]
        # template parts expand to whatever the object has; nothing to repair
        if (!is.null(p$for_each)) next
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
        pending <- c(pending, .alt_pending(m, i, new))
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
        # report new data (matrices, tables), not new columns or settings
        added <- diff$new_path[diff$change == "added" &
            diff$kind %in% c("sparse_matrix", "dense_matrix", "data.frame", "S4")]
        added <- added[!sub(".*/", "", added) %in% used & !startsWith(added, "commands/")]
        for (pth in added) {
            note(NA_character_, "added", sprintf("%s is new and unmapped (no change)", pth))
        }
    }
    candidates <- .repair_candidates(m, unresolved, old, new, po, pn)
    llm_log <- NULL
    if (!is.null(llm) && (length(unresolved) || length(pending))) {
        items <- .llm_items(m, unresolved, pending, candidates, new, pn)
        if (length(items)) {
            req <- .llm_request(m, items, new, diff)
            ans <- llmPatch(req, llm)
            llm_log <- ans$log
            ap <- .apply_llm_patches(m, ans$patches, items)
            m <- ap$mapping
            for (r in seq_len(nrow(ap$applied))) {
                note(ap$applied$part_id[r], "llm", ap$applied$detail[r])
            }
            resolved <- unique(ap$applied$part_id)
            unresolved <- setdiff(unresolved, resolved)
            changes <- Filter(function(ch) !(ch$rule == "absent" && ch$part_id %in% resolved), changes)
        }
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

    res <- structure(list(status = status, mapping = m, changes = changes, unresolved = unresolved,
        candidates = candidates, validation = validation, readback = readback,
        diff = diff, llm_log = llm_log, profile = new), class = "giottoRepair")
    if (status == "failed" && !is.null(request_file)) writeMappingRequest(res, request_file)
    res
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


# LLM items and patches ####

# declared alternatives of part i whose source exists but whose own params
# did not bind by rule: candidates for an LLM param choice
.alt_pending <- function(m, i, new) {
    p <- m$parts[[i]]
    out <- list()
    for (alt in p$alternatives) {
        if (!length(alt$params)) next
        mm <- m
        mm$params <- c(m$params, alt$params)
        b <- .bind_params(mm, new)
        miss <- intersect(names(alt$params), b$unbound)
        node <- m$source_nodes[[alt$from]]
        if (!length(miss) || grepl("{", .pf_bind(node$path, b$values), fixed = TRUE)) next
        if (!.pf_bind(node$path, b$values) %in% new$nodes$path) next
        out[[length(out) + 1L]] <- list(i = i, part_id = p$id, alt = alt, unbound = miss,
            values = b$values)
    }
    out
}

# closed-choice items: a source path per unresolved part, a value per
# unbound param of a pending alternative
.llm_items <- function(m, unresolved, pending, candidates, new, pn) {
    items <- list()
    ids <- vapply(m$parts, `[[`, character(1), "id")
    for (pid in unresolved) {
        cand <- candidates[candidates$part_id == pid, , drop = FALSE]
        if (!nrow(cand)) next
        p <- m$parts[[match(pid, ids)]]
        node <- m$source_nodes[[p$from]]
        items[[length(items) + 1L]] <- list(item = paste0(pid, ".from"), kind = "from", part_id = pid,
            part = p$part, options = cand$path, option_paths = cand$path,
            scores = stats::setNames(cand$score, cand$path),
            question = sprintf("Part `%s` (%s) used to read `%s`, which is gone. Which node now holds it?",
                pid, p$part, .pf_bind(node$path, pn$values)))
    }
    lists <- .LIST_PATHS[[m$mapping$source$ecosystem]]
    for (pd in pending) {
        p <- m$parts[[pd$i]]
        for (prm in pd$unbound) {
            pr <- pd$alt$params[[prm]]
            parent <- gsub("{on}", .pf_bind(pr$on %null% "", pd$values), lists[[pr$list]], fixed = TRUE)
            kids <- new$nodes[.pf_parent(new$nodes$path) == parent & new$nodes$path != parent &
                !endsWith(new$nodes$path, "/*"), , drop = FALSE]
            if (!is.null(pr$where$kind)) kids <- kids[kids$kind == pr$where$kind, , drop = FALSE]
            if (!nrow(kids)) next
            used_in <- paste(vapply(pd$alt$ops, function(o) {
                if (!is.null(o$map)) paste(names(o$map), unlist(o$map), sep = " -> ", collapse = ", ") else o$op
            }, character(1)), collapse = "; ")
            items[[length(items) + 1L]] <- list(item = sprintf("%s.%s.%s", p$id, pd$alt$id, prm),
                kind = "param", part_id = p$id, part = p$part, alt = pd$alt$id, param = prm,
                options = sub(".*/", "", kids$path), option_paths = kids$path, scores = numeric(),
                question = sprintf(paste0("Part `%s` (%s) has no source; its alternative `%s` reads `%s` ",
                    "and needs param `%s` (used as: %s). Which entry is `%s`?"),
                    p$id, p$part, pd$alt$id, parent, prm, used_in, prm))
        }
    }
    items
}

.apply_llm_patches <- function(m, patches, items) {
    applied <- data.frame(part_id = character(), detail = character(), stringsAsFactors = FALSE)
    if (!nrow(patches)) return(list(mapping = m, applied = applied))
    by_item <- stats::setNames(items, vapply(items, `[[`, character(1), "item"))
    ids <- function() vapply(m$parts, `[[`, character(1), "id")
    say <- function(pt) sprintf("%s = %s (%s): %s", pt$item, pt$value, pt$confidence, pt$rationale)

    # source paths, read structurally
    for (r in which(vapply(patches$item, function(x) by_item[[x]]$kind, character(1)) == "from")) {
        pt <- patches[r, ]
        it <- by_item[[pt$item]]
        i <- match(it$part_id, ids())
        node_id <- paste0(it$part_id, "_patched")
        m$source_nodes[[node_id]] <- list(path = pt$value, get = list(fn = "path", on = pt$value))
        m$parts[[i]]$from <- node_id
        applied <- rbind(applied, data.frame(part_id = it$part_id, detail = say(pt)))
    }
    # params of an alternative: all of them must be chosen
    param_items <- Filter(function(it) it$kind == "param", items)
    groups <- split(param_items, vapply(param_items, function(it) paste(it$part_id, it$alt), character(1)))
    for (g in groups) {
        chosen <- patches[patches$item %in% vapply(g, `[[`, character(1), "item"), , drop = FALSE]
        if (nrow(chosen) != length(g)) next
        i <- match(g[[1]]$part_id, ids())
        p <- m$parts[[i]]
        alt <- Filter(function(a) a$id == g[[1]]$alt, p$alternatives)[[1]]
        for (it in g) {
            v <- chosen$value[chosen$item == it$item]
            pr <- alt$params[[it$param]]
            alt$params[[it$param]] <- list(list = pr$list, on = pr$on, rule = "name_is", value = v)
        }
        m <- .apply_alternative(m, i, list(alt = alt, detail = ""))
        for (r in seq_len(nrow(chosen))) {
            applied <- rbind(applied, data.frame(part_id = g[[1]]$part_id, detail = say(chosen[r, ])))
        }
    }
    list(mapping = m, applied = applied)
}

#' @title Write a mapping repair request for a human or a stronger model
#' @name writeMappingRequest
#' @description Write a markdown handoff for a repair that failed: what the
#' mapping expected, what changed in the object, which parts are unresolved
#' with their ranked candidates, the validation and readback results, and the
#' contract lines for the affected parts.
#' @param repair a `giottoRepair` from [repairMapping()]
#' @param file output path (`.md`)
#' @returns `file`, invisibly
#' @export
writeMappingRequest <- function(repair, file) {
    checkmate::assert_class(repair, "giottoRepair")
    m <- repair$mapping
    ct <- giottoContract()
    tab <- function(df) {
        if (is.null(df) || !nrow(df)) return("(none)")
        df <- as.data.frame(unclass(df), stringsAsFactors = FALSE)
        df[] <- lapply(df, function(v) gsub("|", "\\|", as.character(v), fixed = TRUE))
        paste(c(paste0("| ", paste(names(df), collapse = " | "), " |"),
            paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
            apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))), collapse = "\n")
    }
    parts <- Filter(function(p) p$id %in% c(repair$unresolved, repair$changes$part_id), m$parts)
    contract <- vapply(parts, function(p) {
        cp <- ct$parts[[p$part]]
        sprintf("- `%s` -> %s via `%s` (keys: %s)%s", p$id, p$part, cp$constructor,
            paste(cp$nesting, collapse = "/"), if (isTRUE(p$optional)) ", optional" else "")
    }, character(1))
    val <- repair$validation[repair$validation$level != "OK", , drop = FALSE]
    rb <- if (is.null(repair$readback)) NULL else repair$readback[repair$readback$status != "MATCH", ]
    txt <- c(
        sprintf("# Mapping repair request: %s v%s", m$mapping$id, m$mapping$version), "",
        sprintf("Status: **%s**. Generated %s.", repair$status, format(Sys.time(), "%Y-%m-%d %H:%M")), "",
        "## Unresolved parts",
        if (length(repair$unresolved)) paste0("- `", repair$unresolved, "`") else "(none)", "",
        "## Ranked candidates", tab(repair$candidates), "",
        "## Changes applied so far", tab(repair$changes), "",
        "## Structural changes in the object", tab(repair$diff), "",
        "## Validation (non-OK)", tab(val), "",
        "## Readback (non-MATCH)", tab(rb), "",
        "## Contract for the affected parts", if (length(contract)) contract else "(none)", "",
        "## How to answer",
        paste("Edit the mapping YAML: point the part's `from` at a source node (add one under",
            "`source_nodes` with `path` and an accessor `get`), or declare an alternative with",
            "params. Then run `repairMapping()` or `convertWithMapping()` + `readbackCheck()`.")
    )
    writeLines(txt, file)
    invisible(file)
}
