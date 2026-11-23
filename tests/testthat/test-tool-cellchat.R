# CellChat connection. The offline tests cover the param and the micron scale
# of the spatial mode; the CellChat run is opt-in (GIOTTOCONNECT_TOOL_TESTS).

test_that("cellchatParam stores CellChat's options", {
    p <- cellchatParam(species = "human", type = "truncatedMean",
        center_distance = 100, contact.knn.k = 6)
    expect_s4_class(p, "analyzeParam")
    expect_equal(p$species, "human")
    expect_equal(p$type, "truncatedMean")
    expect_equal(p$contact.knn.k, 6)
    expect_error(cellchatParam(species = "rat"))
    expect_error(cellchatParam(thresh = 2))
})

test_that("the spatial mode derives ratio and tol from the coordinates", {
    # square grid with 137 units between neighbours, 100 microns apart
    locs <- as.matrix(expand.grid(x = 0:4 * 137, y = 0:4 * 137))
    rownames(locs) <- paste0("c", seq_len(nrow(locs)))
    sf <- GiottoConnect:::.cellchat_spatial_factors(locs,
        cellchatParam(center_distance = 100))
    expect_equal(sf$ratio, 100 / 137)
    expect_equal(sf$tol, 50)
    sf <- GiottoConnect:::.cellchat_spatial_factors(locs,
        cellchatParam(ratio = 1, tol = 5))
    expect_equal(c(sf$ratio, sf$tol), c(1, 5))
    expect_error(GiottoConnect:::.cellchat_spatial_factors(locs,
        cellchatParam()), "micron")
})

test_that("runCellChat returns CellChat's interactions on the visium mini", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("CellChat")
    skip_if_not_installed("GiottoData")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("visium")
    res <- runCellChat(g, cell_type = "leiden_clus",
        param = cellchatParam(species = "mouse", type = "truncatedMean",
            center_distance = 100, tol = 32.5, nboot = 20))
    expect_s3_class(res, "data.table")
    expect_true(all(c("source", "target", "ligand", "receptor", "prob",
        "pval", "pathway_name") %in% names(res)))
    expect_gt(nrow(res), 0L)
    expect_true(all(res$pval <= 0.05))
    expect_s4_class(attr(res, "cellchat"), "CellChat")
    expect_equal(attr(res, "cellchat")@options$datatype, "spatial")
    expect_s4_class(attr(res, "param"), "cellchatParam")
})
