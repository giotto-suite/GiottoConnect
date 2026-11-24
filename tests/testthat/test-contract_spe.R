spe_env <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            cache <<- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = cache)
        }
        cache
    }
})

spe_mapping_of <- function(x) {
    findMapping(adapterProfile(GiottoConnect:::.adapter_for("spatialexperiment"), x),
        ecosystem = "spatialexperiment")
}

test_that("a SpatialExperiment converts through spe-basic", {
    skip_if_not_installed("SpatialExperiment")
    x <- spe_env()$build_spe_fixture()
    f <- spe_mapping_of(x)
    expect_equal(f$mapping$mapping$id, "spe-basic")
    expect_equal(f$status, "match")

    expect_no_warning(g <- spatialExperimentToGiotto(x))
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    counts <- SummarizedExperiment::assay(x, "counts")
    expect_equal(as.matrix(raw)[rownames(counts), colnames(counts)], as.matrix(counts))
    expect_setequal(GiottoClass::list_expression(g)$name[GiottoClass::list_expression(g)$feat_type == "rna"],
        c("raw", "normalized"))

    # the control altExp is its own feature type
    neg <- GiottoClass::getExpression(g, feat_type = "NegControlProbe", values = "raw", output = "matrix")
    expect_equal(nrow(neg), 4L)

    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    xy <- SpatialExperiment::spatialCoords(x)
    expect_equal(sl$sdimx[match(colnames(x), sl$cell_ID)], unname(xy[, 1]))
    expect_equal(sl$sdimy[match(colnames(x), sl$cell_ID)], -unname(xy[, 2]))

    d <- GiottoClass::getDimReduction(g, reduction_method = "pca", name = "pca", output = "dimObj")
    pca <- SingleCellExperiment::reducedDim(x, "PCA")
    expect_equal(unname(d@misc$loadings), unname(attr(pca, "rotation")))
    expect_equal(d@misc$eigenvalues, attr(pca, "varExplained"))
    expect_setequal(GiottoClass::list_dim_reductions(g)$name, c("pca", "umap"))
    expect_true(attr(readbackCheck(x, g), "pass"))
})

test_that("VisiumIO layouts use spe-visium and convert pixels to microns", {
    skip_if_not_installed("SpatialExperiment")
    x <- spe_env()$build_spe_visium_fixture(n_side = 6L, um_per_px = 0.25)
    f <- spe_mapping_of(x)
    expect_equal(f$mapping$mapping$id, "spe-visium")
    g <- spatialExperimentToGiotto(x)
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    # neighbouring bins are one 16 um bin apart
    expect_equal(sort(unique(round(diff(sort(unique(sl$sdimx))), 6))), 16)
    l <- mappingProvenance(g)$losses
    expect_match(l$note[l$part_id %in% "locs" & l$field == "units"], "0.25000 um per px")
})

test_that("a SpatialFeatureExperiment brings geometries, transcripts and graphs", {
    skip_if_not_installed("SpatialFeatureExperiment")
    x <- spe_env()$build_sfe_fixture()
    f <- spe_mapping_of(x)
    expect_equal(f$mapping$mapping$id, "sfe")
    expect_no_warning(g <- spatialExperimentToGiotto(x))

    # control rows of the main matrix are split into their feature types
    ex <- GiottoClass::list_expression(g)
    expect_setequal(unique(ex$feat_type), c("rna", "negprobes", "SystemControl"))
    expect_equal(nrow(GiottoClass::getExpression(g, values = "raw", output = "matrix")), 40L)
    expect_setequal(rownames(GiottoClass::getExpression(g, feat_type = "negprobes", values = "raw",
        output = "matrix")), c("Negative1", "Negative2"))

    expect_equal(GiottoClass::list_spatial_info(g)$spat_info, "cell")
    expect_setequal(GiottoClass::list_feature_info(g)$feat_info, c("rna", "negprobes", "SystemControl"))
    np <- GiottoClass::getFeatureInfo(g, feat_type = "negprobes", return_giottoPoints = TRUE)
    expect_equal(nrow(np[]), 12L)

    sn <- GiottoClass::getSpatialNetwork(g, name = "knn4", output = "networkDT")
    expect_equal(sn$weight, 1 / sn$distance)
    expect_false(any(duplicated(paste(pmin(sn$from, sn$to), pmax(sn$from, sn$to)))))
    expect_true(attr(readbackCheck(x, g), "pass"))
})

test_that("SpatialExperiments convert onto a GiottoDisk project", {
    skip_if_not_installed("SpatialExperiment")
    skip_if_not_installed("GiottoDisk")
    x <- spe_env()$build_spe_fixture()
    g <- spatialExperimentToGiotto(x, source = withr_tempdir())
    expect_true(attr(readbackCheck(x, g), "pass"))
})

test_that("select_rows keeps or drops rows by name", {
    m <- matrix(1:6, 3, dimnames = list(c("Negative1", "geneA", "SystemControl2"), c("a", "b")))
    keep <- GiottoConnect:::.apply_op(m, list(op = "select_rows", pattern = "^(Negative|SystemControl)"))
    expect_equal(rownames(keep), c("Negative1", "SystemControl2"))
    drop <- GiottoConnect:::.apply_op(m, list(op = "select_rows", pattern = "^Negative", invert = TRUE))
    expect_equal(rownames(drop), c("geneA", "SystemControl2"))
})

test_that("spot spacing is measured against every spot, whatever their order", {
    xy <- as.matrix(expand.grid(x = seq(0, 990, by = 10), y = seq(0, 990, by = 10)))
    shuffled <- xy[sample(nrow(xy)), ]
    expect_equal(GiottoConnect:::.median_nn_spacing(shuffled), 10)
})
