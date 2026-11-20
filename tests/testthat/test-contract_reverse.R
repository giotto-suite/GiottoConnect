giotto_fixture <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) cache <<- seuratToGiotto(seurat_fixture())
        cache
    }
})

test_that("giotto objects select the Giotto-to-Seurat mapping", {
    skip_if_not_installed("SeuratObject")
    f <- findMapping(profileObject(giotto_fixture()), ecosystem = "giotto")
    expect_equal(f$mapping$mapping$id, "giotto-seurat5")
    expect_equal(f$status, "match")
    expect_equal(f$params$values$spat_unit, "cell")
    expect_equal(f$params$values$feat_type, "rna")
})

test_that("giottoToSeurat converts and reads back, in memory and from GiottoDisk", {
    skip_if_not_installed("SeuratObject")
    g <- giotto_fixture()
    s2 <- giottoToSeurat(g)
    expect_s4_class(s2, "Seurat")
    rb <- readbackCheck(g, s2)
    expect_true(attr(rb, "pass"))
    expect_true(all(rb$status == "MATCH"))
    prov <- mappingProvenance(s2)
    expect_equal(prov$mapping, "giotto-seurat5")
    expect_true("expr_raw" %in% prov$parts)

    skip_if_not_installed("GiottoDisk")
    gd <- seuratToGiotto(seurat_fixture(), source = withr_tempdir())
    sd <- giottoToSeurat(gd)
    expect_true(attr(readbackCheck(gd, sd), "pass"))
    expect_equal(as.matrix(SeuratObject::LayerData(sd, layer = "counts")),
        as.matrix(SeuratObject::LayerData(s2, layer = "counts")))
})

test_that("Seurat -> Giotto -> Seurat returns the original data", {
    skip_if_not_installed("SeuratObject")
    SO <- asNamespace("SeuratObject")
    s <- seurat_fixture()
    s2 <- giottoToSeurat(seuratToGiotto(s))
    cells <- colnames(s)
    expect_setequal(colnames(s2), cells)
    for (layer in c("counts", "data")) {
        a <- SO$LayerData(s, layer = layer)
        b <- SO$LayerData(s2, layer = layer)
        expect_equal(as.matrix(b[rownames(a), cells]), as.matrix(a), ignore_attr = TRUE)
    }
    md <- s[[]]
    expect_equal(s2[[]][cells, names(md)], md[cells, ], ignore_attr = TRUE)

    c1 <- SO$GetTissueCoordinates(s[["fov"]], which = "centroids")
    c2 <- SO$GetTissueCoordinates(s2[["fov"]][["centroids"]])
    c2 <- c2[match(c1$cell, c2$cell), ]
    expect_equal(c2$x, c1$x)
    expect_equal(c2$y, c1$y)
    g1 <- SO$GetTissueCoordinates(s[["fov"]][["segmentation"]])
    g2 <- SO$GetTissueCoordinates(s2[["fov"]][["segmentation"]])
    expect_equal(tapply(g2$x, g2$cell, max)[cells], tapply(g1$x, g1$cell, max)[cells])
    expect_equal(tapply(g2$y, g2$cell, min)[cells], tapply(g1$y, g1$cell, min)[cells])

    for (r in c("pca", "umap")) {
        expect_equal(unname(SO$Embeddings(s2, r)[cells, ]), unname(SO$Embeddings(s, r)[cells, ]))
    }
    # the kNN graph comes back as Seurat's binary adjacency with self edges
    expect_equal(as.matrix(s2[["rna_nn"]][cells, cells]), as.matrix(s[["Xenium_nn"]][cells, cells]),
        ignore_attr = TRUE)
})

test_that("a Giotto object that did not come from Seurat converts", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("GiottoData")
    skip_on_cran()
    g <- suppressWarnings(GiottoData::loadGiottoMini("visium"))
    s <- suppressWarnings(giottoToSeurat(g))
    expect_s4_class(s, "Seurat")
    expect_equal(ncol(s), length(GiottoClass::spatIDs(g)))
    expect_true(attr(suppressWarnings(readbackCheck(g, s)), "pass"))
})
