contract_script <- function(name) {
    f <- system.file("scripts", name, package = "GiottoConnect")
    if (!nzchar(f)) f <- test_path("..", "..", "inst", "scripts", name)
    f
}

seurat_fixture <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            env <- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = env)
            cache <<- env$build_seurat_fixture()
        }
        cache
    }
})

# run R code in a fresh process (base R only, no callr)
run_fresh_r <- function(code) {
    f <- tempfile(fileext = ".R")
    writeLines(code, f)
    out <- system2(file.path(R.home("bin"), "Rscript"), f, stdout = TRUE, stderr = TRUE)
    attr(out, "status") <- attr(out, "status") %||% 0L
    out
}

withr_tempdir <- function(env = parent.frame()) {
    d <- tempfile("gc_test_")
    withr_defer <- function() unlink(d, recursive = TRUE)
    do.call(on.exit, list(substitute(withr_defer()), add = TRUE), envir = env)
    d
}
