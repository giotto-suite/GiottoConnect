spe_export_env <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            cache <<- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = cache)
        }
        cache
    }
})

test_that("giottoToSpatialExperiment writes a SpatialExperiment that converts back", {
    skip_if_not_installed("SpatialExperiment")
    g <- spe_export_env()$build_giotto_fixture()
    x <- giottoToSpatialExperiment(g, sfe = "no")
    expect_s4_class(x, "SpatialExperiment")
    expect_false(methods::is(x, "SpatialFeatureExperiment"))
    expect_equal(mappingProvenance(x)$mapping, "giotto-spe")
    expect_setequal(SummarizedExperiment::assayNames(x), c("counts", "logcounts"))
    expect_setequal(SingleCellExperiment::reducedDimNames(x), c("PCA", "UMAP"))

    # image space: y = -sdimy
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    xy <- SpatialExperiment::spatialCoords(x)
    expect_equal(unname(xy[sl$cell_ID, 2]), -sl$sdimy)

    d <- GiottoClass::getDimReduction(g, reduction_method = "pca", name = "pca", output = "dimObj")
    expect_equal(unname(attr(SingleCellExperiment::reducedDim(x, "PCA"), "rotation")), unname(d@misc$loadings))
    h <- SingleCellExperiment::colPair(x, "kNN.pca")
    expect_equal(S4Vectors::metadata(h)$giotto$type, "kNN")
    expect_true(attr(readbackCheck(g, x), "pass"))

    # and back, with the same names
    g2 <- spatialExperimentToGiotto(x)
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    raw2 <- GiottoClass::getExpression(g2, values = "raw", output = "matrix")
    expect_equal(as.matrix(raw2[rownames(raw), colnames(raw)]), as.matrix(raw))
    expect_setequal(GiottoClass::list_dim_reductions(g2)$name, c("pca", "umap"))
    a <- as.data.frame(GiottoClass::getNearestNetwork(g, nn_type = "kNN", name = "kNN.pca", output = "data.table"))
    b <- as.data.frame(GiottoClass::getNearestNetwork(g2, nn_type = "kNN", name = "kNN.pca", output = "data.table"))
    m <- match(paste(a$from, a$to), paste(b$from, b$to))
    expect_false(anyNA(m))
    expect_equal(b$weight[m], a$weight)
})

test_that("polygons, transcripts and spatial networks make a SpatialFeatureExperiment", {
    skip_if_not_installed("SpatialFeatureExperiment")
    g <- spe_export_env()$build_giotto_sfe_fixture()
    x <- giottoToSpatialExperiment(g)
    expect_s4_class(x, "SpatialFeatureExperiment")
    expect_true("cellSeg" %in% SpatialFeatureExperiment::colGeometryNames(x))
    expect_equal(SpatialFeatureExperiment::rowGeometryNames(x), "txSpots")
    expect_equal(SpatialFeatureExperiment::colGraphNames(x), "knn4")
    # spatialCoords are the Giotto locations, not the polygon centroids
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    expect_equal(unname(SpatialExperiment::spatialCoords(x)[sl$cell_ID, 1]), sl$sdimx)
    expect_true(attr(readbackCheck(g, x), "pass"))

    g2 <- spatialExperimentToGiotto(x)
    a <- GiottoClass::getSpatialNetwork(g, name = "knn4", output = "networkDT")
    b <- GiottoClass::getSpatialNetwork(g2, name = "knn4", output = "networkDT")
    expect_equal(nrow(b), nrow(a))
    expect_equal(sort(b$weight), sort(a$weight))
    p1 <- GiottoClass::getFeatureInfo(g, feat_type = "rna", return_giottoPoints = TRUE)
    p2 <- GiottoClass::getFeatureInfo(g2, feat_type = "rna", return_giottoPoints = TRUE)
    expect_equal(as.vector(table(p2[]$feat_ID)), as.vector(table(p1[]$feat_ID)))
    expect_equal(nrow(GiottoClass::getPolygonInfo(g2, polygon_name = "cell", return_giottoPolygon = TRUE)[]),
        nrow(GiottoClass::getPolygonInfo(g, polygon_name = "cell", return_giottoPolygon = TRUE)[]))
})

test_that("one feature type is converted at a time", {
    skip_if_not_installed("SpatialExperiment")
    g <- spe_export_env()$build_giotto_sfe_fixture()
    x <- giottoToSpatialExperiment(g, feat_type = "negprobes", sfe = "no")
    expect_setequal(rownames(x), c("Negative1", "Negative2"))
    expect_true(attr(readbackCheck(g, x), "pass"))
})

test_that("a GiottoDisk-backed object converts to a SpatialExperiment", {
    skip_if_not_installed("SpatialExperiment")
    skip_if_not_installed("GiottoDisk")
    g <- spe_export_env()$build_giotto_fixture(source = withr_tempdir())
    x <- giottoToSpatialExperiment(g, sfe = "no")
    expect_true(attr(readbackCheck(g, x), "pass"))
})
