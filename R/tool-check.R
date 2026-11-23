# Conformance check for methods added to the Suite (adr/0008 in GiottoClass):
# compares a giotto object before and after a run*() call. Each failed check
# says what to fix, so people and agents can iterate on it.

#' @title Check a method integration
#' @name checkGiottoMethod
#' @description
#' Compares a giotto object before and after running a method's `run*()`
#' wrapper and checks the rules of GiottoClass ADR 0008:
#' * only the declared slots changed (via [GiottoClass::manifestDiff()]);
#' * something was stored;
#' * stored dimension reductions and spatial enrichments cover every cell of
#'   their spatial unit, and new cell metadata columns are not all `NA`;
#' * the last history entry carries the method's param
#'   (`update_giotto_params(attachments = list(param = param))`);
#' * the help page of `fn` has references and asks users to cite them.
#'
#' Methods whose results fit no slot (e.g. an interaction table) return them
#' instead: pass `slots = "none"` and the returned object as `result`. Then
#' the object must be unchanged and `attr(result, "param")` must hold the
#' param.
#' @param before,after giotto objects before and after the run
#' @param slots slots the method may change, e.g.
#'   `c("cell_metadata", "spatial_enrichment")`, or `"none"` for methods that
#'   return their results
#' @param result the returned result, with `slots = "none"
#' @param param_class class of the param expected in the history; `NULL`
#'   accepts any param
#' @param fn,package name of the `run*()` function and its package, for the
#'   documentation check; skipped when `NULL`
#' @returns A `giottoMethodCheck`: a `data.table` with columns `check`,
#'   `pass` and `message`, printed as a report.
#' @examples
#' \dontrun{
#' g <- GiottoData::loadGiottoMini("visium")
#' g2 <- Giotto::runNeighborhoodPCA(g, feats_to_use = "hvf")
#' checkGiottoMethod(g, g2, slots = "dimension_reduction",
#'     param_class = "neighborhoodPCAParam", fn = "runNeighborhoodPCA",
#'     package = "Giotto")
#' }
#' @md
#' @export
checkGiottoMethod <- function(before, after, slots, param_class = NULL,
    fn = NULL, package = NULL, result = NULL) {
    checkmate::assert_class(before, "giotto")
    checkmate::assert_class(after, "giotto")
    checkmate::assert_character(slots, min.len = 1L)
    out <- list()
    add <- function(check, pass, message = "") {
        out[[length(out) + 1L]] <<- data.table::data.table(check = check,
            pass = isTRUE(pass), message = if (isTRUE(pass)) "" else message)
    }

    d <- manifestDiff(objManifest(before), objManifest(after))$detail
    changed <- c(names(d$added), names(d$modified), d$removed)
    families <- unique(sub("\\..*$", "", changed))

    returns <- identical(slots, "none")
    if (returns) {
        add("object unchanged", !length(changed),
            paste0("slots = \"none\" but changed: ", paste(changed, collapse = ", ")))
    } else {
        add("something stored", length(c(names(d$added), names(d$modified))) > 0L,
            "no slot changed: the run*() wrapper must write its results with the setters")
    }
    extra <- setdiff(families, c(slots, "parameters"))
    add("only declared slots changed", !length(extra),
        paste0("undeclared slot(s) changed: ", paste(extra, collapse = ", "),
            ". Write only the slots of the result-to-slot map, or declare them in `slots`"))
    add("nothing removed", !length(d$removed),
        paste0("removed: ", paste(d$removed, collapse = ", ")))

    # coverage of stored results
    for (k in names(d$added)) {
        entry <- d$added[[k]]
        fam <- sub("\\..*$", "", k)
        keys <- entry$keys
        n_ok <- tryCatch(.check_coverage(after, fam, keys), error = function(e) NA)
        if (is.na(n_ok)) next
        add(paste("covers all cells:", k), n_ok,
            "the stored object does not have one row per cell of its spatial unit")
    }
    for (k in names(d$modified)) {
        if (!startsWith(k, "cell_metadata")) next
        new_cols <- .new_meta_columns(before, after, k)
        allna <- vapply(new_cols, function(cn) all(is.na(cn$values)), logical(1))
        add(paste("new metadata filled:", k), !any(allna),
            paste0("all-NA column(s): ", paste(names(new_cols)[allna], collapse = ", "),
                ". Join results by cell_ID (addCellMetadata(by_column = TRUE))"))
    }

    # provenance
    h <- after@parameters
    p <- if (returns) attr(result, "param") else if (length(h)) {
        attr(h[[length(h)]], "attachments")$param
    }
    is_param <- !is.null(p) && (if (is.null(param_class)) {
        any(vapply(c("analyzeParam", "processParam", "reduceParam",
            "filterParam"), function(cl) methods::is(p, cl), logical(1)))
    } else methods::is(p, param_class))
    add(if (returns) "param on the result" else "param in history", is_param,
        paste0(if (returns) "set attr(result, \"param\") <- param before returning" else
            paste0("end the wrapper with update_giotto_params(gobject, description = ",
                "\"_<name>\", attachments = list(param = param))"),
        if (!is.null(param_class)) paste0(" with a ", param_class) else ""))

    # credit
    if (!is.null(fn) && !is.null(package)) {
        rd <- .check_rd_text(fn, package)
        add("help page has references", rd$has_refs,
            "add @references with the original paper (and \\doi{})")
        add("help page asks to cite", grepl("cite", rd$text, ignore.case = TRUE),
            "say in @description: 'If you use it, please cite the original paper'")
    }

    res <- data.table::rbindlist(out)
    class(res) <- c("giottoMethodCheck", class(res))
    res
}

#' @export
print.giottoMethodCheck <- function(x, ...) {
    x <- data.table::as.data.table(unclass(x))
    cat(sprintf("checkGiottoMethod: %d/%d checks pass\n", sum(x$pass), nrow(x)))
    for (i in seq_len(nrow(x))) {
        cat(if (x$pass[i]) "  ok   " else "  FAIL ", x$check[i], "\n", sep = "")
        if (!x$pass[i]) cat("        -> ", x$message[i], "\n", sep = "")
    }
    invisible(x)
}

#' @rdname checkGiottoMethod
#' @param ... arguments of `checkGiottoMethod()`
#' @details `expect_giotto_method()` is the \pkg{testthat} expectation: it
#'   fails with the report when any check fails.
#' @export
expect_giotto_method <- function(...) {
    package_check("testthat", repository = "CRAN")
    res <- checkGiottoMethod(...)
    msg <- paste(utils::capture.output(print(res)), collapse = "\n")
    testthat::expect(all(res$pass), msg)
    invisible(res)
}

.check_coverage <- function(g, fam, keys) {
    su <- keys$spat_unit
    n <- length(spatIDs(g, spat_unit = su))
    if (fam == "dimension_reduction") {
        m <- getDimReduction(g, spat_unit = su, feat_type = keys$feat_type,
            reduction = keys$reduction, reduction_method = keys$reduction_method,
            name = keys$name, output = "matrix")
        return(nrow(m) == n && !anyNA(m))
    }
    if (fam == "spatial_enrichment") {
        e <- getSpatialEnrichment(g, spat_unit = su, feat_type = keys$feat_type,
            name = keys$name, output = "data.table")
        return(nrow(e) == n)
    }
    NA
}

.new_meta_columns <- function(before, after, key) {
    parts <- strsplit(key, ".", fixed = TRUE)[[1]]
    su <- parts[2]; ft <- parts[3]
    b <- names(pDataDT(before, spat_unit = su, feat_type = ft))
    a <- pDataDT(after, spat_unit = su, feat_type = ft)
    new <- setdiff(names(a), b)
    stats::setNames(lapply(new, function(cn) list(values = a[[cn]])), new)
}

.check_rd_text <- function(fn, package) {
    db <- tools::Rd_db(package)
    hit <- Filter(function(rd) {
        fn %in% unlist(lapply(rd, function(x)
            if (identical(attr(x, "Rd_tag"), "\\alias")) as.character(x)))
    }, db)
    if (!length(hit)) return(list(has_refs = FALSE, text = ""))
    rd <- hit[[1]]
    tags <- vapply(rd, function(x) attr(x, "Rd_tag") %null% "", character(1))
    list(has_refs = "\\references" %in% tags,
        text = paste(utils::capture.output(tools::Rd2txt(rd)), collapse = " "))
}
