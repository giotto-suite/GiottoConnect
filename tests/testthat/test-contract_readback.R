test_that("readbackCheck passes for in-memory and disk conversions", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    rb <- readbackCheck(s, seuratToGiotto(s))
    expect_s3_class(rb, "giottoReadback")
    expect_true(attr(rb, "pass"))
    expect_true(all(rb$status == "MATCH"))
    expect_setequal(unique(rb$part_id), c("expr_raw", "expr_normalized", "cell_meta",
        "feat_meta", "locs", "polygons", "dimred_pca", "dimred_umap", "nn"))
    expect_output(print(rb), "PASS")

    skip_if_not_installed("GiottoDisk")
    d <- withr_tempdir()
    rb_disk <- readbackCheck(s, seuratToGiotto(s, source = d))
    expect_true(attr(rb_disk, "pass"))
})

test_that("readbackCheck flags changed content as UNDECLARED", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    g <- seuratToGiotto(s)
    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    raw[1, ] <- raw[1, ] + 1
    g2 <- GiottoClass::setGiotto(g, GiottoClass::createExprObj(raw, name = "raw"), verbose = FALSE)
    sl <- GiottoClass::getSpatialLocations(g2)
    sl[]$sdimy <- -sl[]$sdimy
    g2 <- GiottoClass::setGiotto(g2, sl, verbose = FALSE)
    rb <- readbackCheck(s, g2)
    expect_false(attr(rb, "pass"))
    bad <- rb[rb$status == "UNDECLARED", ]
    expect_setequal(paste(bad$part_id, bad$aspect), c("expr_raw values", "locs values"))
})

test_that("differences covered by a recorded loss are DECLARED", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    g <- seuratToGiotto(s)
    p <- GiottoClass::getPolygonInfo(g, name = "cell", return_giottoPolygon = TRUE)
    p[] <- terra::shift(p[], dx = 1)
    g2 <- GiottoClass::setGiotto(g, p, verbose = FALSE)
    rb <- readbackCheck(s, g2)
    geo <- rb[rb$part_id == "polygons" & rb$aspect == "geometry", ]
    expect_equal(geo$status, "DECLARED")
    expect_true(attr(rb, "pass"))
    expect_error(readbackCheck(s, GiottoClass::giotto()), "no conversion record")
})
