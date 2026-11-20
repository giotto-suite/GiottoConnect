minimal_mapping <- function() {
    list(
        mapping = list(
            id = "test", version = "0.1.0", status = "draft",
            schema_version = "0.1.0",
            source = list(ecosystem = "seurat", class = "Seurat",
                package = "SeuratObject", applies_to = list()),
            target = list(contract = "giotto", version = "0.7.31")
        ),
        params = list(
            assay = list(list = "Assays", rule = "first")
        ),
        source_nodes = list(
            counts = list(path = "assays/{assay}/layers/counts",
                get = list(fn = "LayerData", args = list(layer = "counts")))
        ),
        parts = list(
            list(id = "expr_raw", part = "expression", from = "counts",
                to = list(spat_unit = "cell", feat_type = "rna", name = "raw"),
                ops = list(list(op = "cast", to = "integer")))
        )
    )
}

test_that("a valid mapping passes, from a list and from YAML", {
    m <- readMapping(minimal_mapping())
    expect_s3_class(m, "giottoMapping")
    expect_output(print(m), "expr_raw \\(expression\\)")

    f <- tempfile(fileext = ".yaml")
    yaml::write_yaml(minimal_mapping(), f)
    m2 <- readMapping(f)
    expect_equal(attr(m2, "file"), normalizePath(f))
    expect_equal(m2$parts[[1]]$id, "expr_raw")
})

test_that("schema problems are reported with their YAML path", {
    bad <- minimal_mapping()
    bad$mapping$status <- "approved"
    bad$parts[[1]]$part <- "expressions"
    bad$parts[[1]]$ops[[1]] <- list(op = "cast")
    bad$parts[[1]]$extra <- TRUE
    bad$source_nodes$counts$get$fn <- "eval"
    bad$mapping$version <- NULL
    err <- tryCatch(readMapping(bad), error = conditionMessage)

    expect_match(err, "mapping.status: 'approved' is not one of", fixed = TRUE)
    expect_match(err, "parts[1](id=expr_raw).part: 'expressions'", fixed = TRUE)
    expect_match(err, "parts[1](id=expr_raw).ops[1](op=cast): missing argument 'to'", fixed = TRUE)
    expect_match(err, "parts[1](id=expr_raw): unknown field 'extra'", fixed = TRUE)
    expect_match(err, "source_nodes.counts.get.fn: 'eval'", fixed = TRUE)
    expect_match(err, "mapping: missing required field 'version'", fixed = TRUE)
})

test_that("types and op arguments are checked", {
    bad <- minimal_mapping()
    bad$parts[[1]]$optional <- "yes"
    bad$parts[[1]]$ops <- list(list(op = "rename", map = list(a = "b"), cols = "a"))
    err <- tryCatch(readMapping(bad), error = conditionMessage)
    expect_match(err, "optional: expected boolean", fixed = TRUE)
    expect_match(err, "op 'rename' does not take 'cols'", fixed = TRUE)

    bad <- minimal_mapping()
    bad$parts <- "expr_raw"
    expect_error(readMapping(bad), "parts\\[1\\]: expected object")
})

test_that("readMapping resolves bundled ids and rejects unknown ones", {
    expect_error(readMapping("no-such-mapping"), "No mapping file or bundled mapping")
})
