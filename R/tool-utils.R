# Shared helpers for the external tool connections (R/tool-<name>.R).
# Tools are never imported: R tools are checked with package_check() when
# their function runs, Python tools get an environment at that point through
# reticulate::py_require() (adr/0008 in GiottoClass).

# Path of a Python tool's command line program, in an environment that
# reticulate builds with uv on first use and caches afterwards. Must be
# called before Python starts in the session, or with a session whose Python
# already provides the tool.
.tool_python_cli <- function(requirement, module, cli,
    python_version = NULL) {
    package_check("reticulate", repository = "CRAN")
    if (utils::packageVersion("reticulate") < "1.41") {
        stop(wrap_txt("[GiottoConnect] reticulate >= 1.41 is needed to",
            "install", requirement, "on first use"), call. = FALSE)
    }
    if (!reticulate::py_available(initialize = FALSE)) {
        reticulate::py_require(requirement, python_version = python_version)
    }
    if (!reticulate::py_module_available(module)) {
        stop(wrap_txt("[GiottoConnect] Python already started in this",
            "session without", requirement, "- restart R and run this",
            "function before other Python code"), call. = FALSE)
    }
    exe <- file.path(dirname(reticulate::py_exe()), cli)
    if (.Platform$OS.type == "windows") {
        exe <- file.path(dirname(reticulate::py_exe()), "Scripts",
            paste0(cli, ".exe"))
    }
    if (!file.exists(exe)) {
        stop("[GiottoConnect] '", cli, "' not found next to ",
            reticulate::py_exe(), call. = FALSE)
    }
    exe
}

# Run a command line program; stop with the end of its log when it fails.
.tool_run <- function(exe, args, log_file, verbose = NULL) {
    vmsg(.v = verbose, "running:", basename(exe), paste(args, collapse = " "))
    status <- system2(exe, args, stdout = log_file, stderr = log_file)
    if (!identical(as.integer(status), 0L)) {
        tail_log <- utils::tail(readLines(log_file, warn = FALSE), 15L)
        stop("[GiottoConnect] ", basename(exe), " failed (status ", status,
            "):\n", paste(tail_log, collapse = "\n"), call. = FALSE)
    }
    invisible(log_file)
}

# A journal article for methodCitation(). Keep it in sync with @references.
.cite_article <- function(author, title, journal, year, volume, pages, doi) {
    utils::bibentry("Article",
        author = author, title = title, journal = journal, year = year,
        volume = volume, pages = pages, doi = doi
    )
}
