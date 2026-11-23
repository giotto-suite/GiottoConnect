export_fixture <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) cache <<- seuratToGiotto(seurat_fixture())
        cache
    }
})

test_that("giottoToAnnData writes a scanpy-layout .h5ad that reads back", {
    skip_if_not_installed("anndataR")
    skip_if_not_installed("hdf5r")
    g <- export_fixture()
    f <- tempfile(fileext = ".h5ad")
    expect_no_warning(out <- giottoToAnnData(g, f))
    expect_true(file.exists(f))
    expect_equal(mappingProvenance(f)$mapping, "giotto-anndata")
    expect_true(attr(readbackCheck(g, f), "pass"))

    ad <- anndataR::read_h5ad(f)
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    norm <- GiottoClass::getExpression(g, values = "normalized", output = "matrix")
    expect_equal(dim(ad$X), rev(dim(raw)))
    # scanpy layout: normalized in X, counts in a layer
    expect_equal(as.matrix(Matrix::t(ad$layers$counts)), as.matrix(raw[colnames(ad$X), rownames(ad$X)]),
        ignore_attr = TRUE)
    expect_equal(as.matrix(Matrix::t(ad$X)), as.matrix(norm[colnames(ad$X), rownames(ad$X)]),
        ignore_attr = TRUE)
    expect_true(all(c("spatial", "X_pca", "X_umap") %in% names(ad$obsm)))
})

test_that("Giotto -> h5ad -> Giotto keeps parts and network names", {
    skip_if_not_installed("anndataR")
    skip_if_not_installed("hdf5r")
    g <- export_fixture()
    f <- tempfile(fileext = ".h5ad")
    giottoToAnnData(g, f)
    g2 <- anndataToGiotto(f)
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    raw2 <- GiottoClass::getExpression(g2, values = "raw", output = "matrix")
    expect_equal(as.matrix(raw2)[rownames(raw), colnames(raw)], as.matrix(raw))
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    sl2 <- GiottoClass::getSpatialLocations(g2, output = "data.table")
    expect_equal(sl2$sdimy[match(sl$cell_ID, sl2$cell_ID)], sl$sdimy)
    expect_equal(GiottoClass::list_nearest_networks(g2)$name, GiottoClass::list_nearest_networks(g)$name)
    expect_setequal(GiottoClass::list_dim_reductions(g2)$name, GiottoClass::list_dim_reductions(g)$name)
})

test_that("a chosen feature type is exported", {
    skip_if_not_installed("anndataR")
    SO <- asNamespace("SeuratObject")
    s <- seurat_fixture()
    neg <- matrix(stats::rpois(4 * ncol(s), 1), 4, ncol(s),
        dimnames = list(sprintf("NegProbe%02d", 1:4), colnames(s)))
    s[["ControlProbe"]] <- SO$CreateAssay5Object(counts = methods::as(neg, "CsparseMatrix"))
    g <- seuratToGiotto(s)
    f <- tempfile(fileext = ".h5ad")
    giottoToAnnData(g, f, feat_type = "NegControlProbe")
    ad <- anndataR::read_h5ad(f)
    expect_setequal(colnames(ad$X), sprintf("NegProbe%02d", 1:4))
    expect_true(attr(readbackCheck(g, f), "pass"))
})

test_that("a GiottoDisk-backed object exports", {
    skip_if_not_installed("anndataR")
    skip_if_not_installed("GiottoDisk")
    g <- seuratToGiotto(seurat_fixture(), source = withr_tempdir())
    f <- tempfile(fileext = ".h5ad")
    giottoToAnnData(g, f)
    expect_true(attr(readbackCheck(g, f), "pass"))
})

test_that("the Visium mini round-trips its networks and reductions", {
    skip_if_not_installed("anndataR")
    skip_if_not_installed("GiottoData")
    g <- suppressWarnings(GiottoData::loadGiottoMini("visium"))
    f <- tempfile(fileext = ".h5ad")
    suppressWarnings(giottoToAnnData(g, f))
    expect_true(attr(readbackCheck(g, f), "pass"))
    g2 <- suppressWarnings(anndataToGiotto(f))
    expect_setequal(GiottoClass::list_nearest_networks(g2)$name, c("sNN.pca", "custom_NN"))
    expect_setequal(GiottoClass::list_spatial_networks(g2)$name, c("Delaunay_network", "spatial_network"))
    expect_setequal(GiottoClass::list_expression(g2)$name, c("raw", "normalized", "scaled"))
    w <- GiottoClass::getNearestNetwork(g, nn_type = "sNN", name = "sNN.pca", output = "data.table")
    w2 <- GiottoClass::getNearestNetwork(g2, nn_type = "sNN", name = "sNN.pca", output = "data.table")
    expect_equal(w2$weight[order(w2$from, w2$to)], w$weight[order(w$from, w$to)])
})
