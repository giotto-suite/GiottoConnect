h5ad_fixture <- local({
    cache <- list()
    function(legacy = FALSE) {
        key <- if (legacy) "legacy" else "modern"
        if (is.null(cache[[key]])) {
            env <- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = env)
            cache[[key]] <<- env$build_h5ad_fixture(tempfile(fileext = ".h5ad"), legacy = legacy)
        }
        cache[[key]]
    }
})

test_that("the h5ad profile reads structure, sizes and roles from the file", {
    skip_if_not_installed("hdf5r")
    fx <- h5ad_fixture()
    p <- adapterProfile(GiottoConnect:::.adapter_for("anndata"), fx$path)
    expect_s3_class(p, "giottoProfile")
    expect_equal(c(p$n_cells, p$n_feat), c(200L, 60L))
    nd <- p$nodes
    row <- function(path) nd[nd$path == path, ]
    expect_equal(row("X")$kind, "dense_matrix")
    expect_equal(row("X")$shape, "200x60")
    expect_equal(row("layers/counts")$class, "csr_matrix")
    expect_match(row("layers/counts")$role_hints, "integer_valued")
    expect_false(grepl("integer_valued", row("X")$role_hints))
    expect_equal(row("obs")$kind, "data.frame")
    expect_equal(row("obs/cell_type")$kind, "factor")
    expect_equal(row("uns/neighbors/params/n_neighbors")$value, "6")
})

test_that("a scanpy file converts through the bundled mapping", {
    skip_if_not_installed("hdf5r")
    fx <- h5ad_fixture()
    f <- findMapping(adapterProfile(GiottoConnect:::.adapter_for("anndata"), fx$path))
    expect_equal(f$mapping$mapping$id, "anndata-spatial")
    expect_equal(f$status, "match")
    expect_equal(f$params$values$counts, "layers/counts")
    expect_equal(f$params$values$normalized, "X")

    g <- anndataToGiotto(fx$path)
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    expect_equal(unname(as.matrix(raw)), unname(t(as.matrix(fx$counts))))
    # duplicated gene names are made unique
    expect_true(all(c("gene01--1", "gene01--2") %in% rownames(raw)))
    norm <- GiottoClass::getExpression(g, values = "normalized", output = "matrix")
    expect_equal(unname(as.matrix(norm)), unname(t(fx$norm)))

    md <- GiottoClass::pDataDT(g)
    expect_true(is.factor(md$cell_type))
    expect_true(is.na(md$cell_type[3L]))

    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    expect_equal(sl$sdimx, fx$xy[, 1])
    expect_equal(sl$sdimy, -fx$xy[, 2])

    dr <- GiottoClass::list_dim_reductions(g)
    expect_setequal(dr$name, c("pca", "umap"))
    pca <- GiottoClass::getDimReduction(g, reduction_method = "pca", name = "pca", output = "matrix")
    expect_equal(unname(pca[fx$cells, ]), fx$pca)

    nn <- GiottoClass::getNearestNetwork(g, nn_type = "kNN", name = "kNN.pca", output = "data.table")
    d <- as.matrix(stats::dist(fx$pca))
    expect_equal(nrow(nn), 200L * 5L)
    expect_equal(nn$distance, d[cbind(match(nn$from, fx$cells), match(nn$to, fx$cells))])
    expect_equal(nn$weight, 1 / (1 + nn$distance))

    l <- mappingProvenance(g)$losses
    expect_match(l$note[l$part_id %in% "nn" & l$field == "name"], "X_pca, 10 PCs")
    expect_true(attr(readbackCheck(fx$path, g), "pass"))
})

test_that("anndata 0.7 files (legacy categoricals, counts in X) convert", {
    skip_if_not_installed("hdf5r")
    fx <- h5ad_fixture(legacy = TRUE)
    g <- suppressWarnings(anndataToGiotto(fx$path))
    expect_equal(mappingProvenance(g)$params$counts, "X")
    expect_null(mappingProvenance(g)$params$normalized)
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    expect_equal(unname(as.matrix(raw)), unname(t(as.matrix(fx$counts))))
    md <- GiottoClass::pDataDT(g)
    expect_equal(as.character(md$cell_type), as.character(fx$type))
    expect_true(attr(readbackCheck(fx$path, g), "pass"))
})

test_that("an h5ad converts onto a GiottoDisk project", {
    skip_if_not_installed("hdf5r")
    skip_if_not_installed("GiottoDisk")
    fx <- h5ad_fixture()
    g <- anndataToGiotto(fx$path, source = withr_tempdir())
    expect_true(attr(readbackCheck(fx$path, g), "pass"))
    gm <- anndataToGiotto(fx$path)
    for (v in c("raw", "normalized")) {
        e <- GiottoClass::getExpression(g, values = v, output = "exprObj")[]
        expect_s4_class(e, "parquetExprStore")
        a <- as.matrix(GiottoClass::getExpression(gm, values = v, output = "matrix"))
        expect_equal(as.matrix(GiottoConnect:::.rb_as_matrix(e)[rownames(a), colnames(a)]), a,
            ignore_attr = TRUE)
    }
})

test_that("matrices stream from the file when the target is on disk", {
    skip_if_not_installed("hdf5r")
    skip_if_not_installed("GiottoDisk")
    skip_if_not(exists("h5adInput", envir = asNamespace("GiottoDisk")))
    fx <- h5ad_fixture()
    m <- readMapping("anndata-spatial")
    p <- adapterProfile(GiottoConnect:::.adapter_for(m), fx$path)
    params <- GiottoConnect:::.bind_params(m, p)
    m <- GiottoConnect:::.expand_for_each(m, p, params)
    prep <- function(stream) {
        GiottoConnect:::.prepare_parts(fx$path, m, params,
            GiottoConnect:::.new_ctx(m, params, p), stream = stream)
    }
    s <- prep(TRUE)
    expect_s4_class(s$expr_raw$data, "exprInput")
    expect_s4_class(s$expr_normalized$data, "exprInput")
    expect_s3_class(s$cell_meta$data, "data.frame")
    expect_true(methods::is(prep(FALSE)$expr_raw$data, "Matrix"))
})

test_that("mappings may only use the listings and accessors of their adapter", {
    m <- readMapping(system.file("mappings", "anndata-spatial.yaml", package = "GiottoConnect"))
    m$params$counts$list <- "Assays"
    m$source_nodes$obs$get$fn <- "LayerData"
    probs <- GiottoConnect:::.mapping_adapter_check(m)
    expect_true(any(grepl("list 'Assays' is not a listing of the anndata adapter", probs)))
    expect_true(any(grepl("'LayerData' is not an accessor of the anndata adapter", probs)))
    expect_true("anndata" %in% adapterRegistry()$name)
})

test_that("params bind paths, exclude roles and for_each strips name prefixes", {
    nodes <- data.frame(
        path = c("", "X", "layers", "layers/counts", "obsm", "obsm/X_pca", "obsm/spatial"),
        kind = c("list", "dense_matrix", "list", "sparse_matrix", "list", "dense_matrix", "dense_matrix"),
        class = NA_character_, dtype = NA_character_, shape = NA_character_,
        role_hints = c(NA, "nonneg", NA, "integer_valued,nonneg", NA, NA, NA),
        value = NA_character_, stringsAsFactors = FALSE)
    profile <- structure(list(nodes = nodes), class = "giottoProfile")
    m <- list(mapping = list(source = list(ecosystem = "anndata")), params = list(
        counts = list(list = "Matrices", rule = "first_where", as = "path",
            where = list(role = "integer_valued")),
        norm = list(list = "Matrices", rule = "first_where", as = "path",
            where = list(kind = list("sparse_matrix", "dense_matrix"), not_role = "integer_valued"))
    ))
    b <- GiottoConnect:::.bind_params(m, profile)
    expect_equal(b$values, list(counts = "layers/counts", norm = "X"))

    m$parts <- list(list(id = "dimred_{red_label}", part = "dim_reduction", from = "red",
        for_each = list(list = "Obsm", as = "red", strip = "^X_", exclude = "spatial"),
        to = list(name = "{red_label}")))
    m$source_nodes <- list(red = list(path = "obsm/{red}", get = list(fn = "h5adObsm")))
    e <- GiottoConnect:::.expand_for_each(m, profile, b)
    expect_equal(vapply(e$parts, `[[`, "", "id"), "dimred_pca")
    expect_equal(e$parts[[1]]$to$name, "pca")
    expect_equal(e$source_nodes[[e$parts[[1]]$from]]$path, "obsm/X_pca")
})

test_that("counts in raw/X bring their own feature table", {
    skip_if_not_installed("hdf5r")
    env <- new.env()
    sys.source(contract_script("build_fixtures.R"), envir = env)
    fx <- env$build_h5ad_fixture(tempfile(fileext = ".h5ad"), raw_extra = 5L)
    expect_no_warning(g <- anndataToGiotto(fx$path))
    expect_equal(mappingProvenance(g)$params$counts, "raw/X")
    fm <- GiottoClass::fDataDT(g)
    # raw/var rows: the var genes plus the raw-only ones, var columns joined
    expect_equal(nrow(fm), 65L)
    expect_true(all(sprintf("extra%02d", 1:5) %in% fm$feat_ID))
    expect_true(all(is.na(fm$highly_variable[startsWith(fm$feat_ID, "extra")])))
    expect_false(anyNA(fm$gene_ids))
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    expect_equal(nrow(raw), 65L)
    expect_true(attr(readbackCheck(fx$path, g), "pass"))
})

test_that("a file without counts converts its normalized values", {
    skip_if_not_installed("hdf5r")
    env <- new.env()
    sys.source(contract_script("build_fixtures.R"), envir = env)
    fx <- env$build_h5ad_fixture(tempfile(fileext = ".h5ad"))
    h5 <- hdf5r::H5File$new(fx$path, mode = "r+")
    h5$link_delete("layers/counts")
    h5$close_all()
    g <- anndataToGiotto(fx$path)
    expect_equal(GiottoClass::list_expression(g)$name, "normalized")
    expect_true("expr_raw" %in% mappingProvenance(g)$absent)
    expect_equal(nrow(GiottoClass::fDataDT(g)), 60L)
})
