test_that("adapters and targets are registered and dispatch", {
    reg <- adapterRegistry()
    expect_setequal(reg$name[reg$role == "adapter"], c("seurat", "giotto", "anndata", "spatialexperiment"))
    expect_setequal(reg$name[reg$role == "target"], c("giotto", "seurat", "anndata", "spatialexperiment"))

    sa <- GiottoConnect:::.adapter_for("seurat")
    ga <- GiottoConnect:::.adapter_for("giotto")
    expect_s4_class(sa, "giottoConnectAdapter")
    expect_equal(sa@package, "SeuratObject")
    expect_equal(GiottoConnect:::adapterListPaths(sa)[["Layers"]], "assays/{on}/layers")
    expect_equal(GiottoConnect:::adapterListPaths(ga)[["Expressions"]], "expression/{on}")
    expect_true("LayerData" %in% GiottoConnect:::adapterAccessors(sa))
    expect_false("LayerData" %in% GiottoConnect:::adapterAccessors(ga))
    expect_s4_class(GiottoConnect:::.target_for("seurat"), "giottoConnectTarget")
    expect_equal(GiottoConnect:::.adapter_for("anndata")@package, "hdf5r")
    expect_error(GiottoConnect:::.adapter_for("spatialdata"), "No adapter for ecosystem 'spatialdata'")
})

test_that("mappings may only call accessors of their adapter", {
    m <- unclass(readMapping("seurat-assay5"))
    m$source_nodes$counts$get$fn <- "getExpression"
    expect_error(readMapping(m),
        "source_nodes.counts.get.fn: 'getExpression' is not an accessor of the seurat adapter")
    m2 <- unclass(readMapping("seurat-assay5"))
    m2$mapping$source$ecosystem <- "spatialdata"
    expect_error(readMapping(m2), "no adapter for 'spatialdata'")
    m3 <- unclass(readMapping("seurat-assay5"))
    m3$mapping$target$contract <- "sce"
    expect_error(readMapping(m3), "no target 'sce'")
})

test_that("the giotto adapter states cell and feature counts", {
    skip_if_not_installed("SeuratObject")
    g <- seuratToGiotto(seurat_fixture())
    p <- GiottoConnect:::adapterProfile(GiottoConnect:::.adapter_for("giotto"), g)
    expect_equal(p$size_from, "user")
    expect_equal(c(p$n_cells, p$n_feat), c(200L, 50L))
})
