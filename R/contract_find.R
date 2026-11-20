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
        m <- yaml::read_yaml(file)
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
