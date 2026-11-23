# runBanksy() (moved from Giotto). Banksy runs are skipped when it is absent.

.banksy_fixture <- function() {
    set.seed(3)
    n <- 400
    cells <- sprintf("c_%03d", seq_len(n))
    x <- stats::runif(n, 0, 100)
    y <- stats::runif(n, 0, 100)
    region <- as.integer(x > 50) + 1L
    rates <- cbind(c(8, 2, 5, 5, 1, 6), c(2, 8, 5, 5, 6, 1))
    m <- matrix(stats::rpois(6 * n, rates[, region]), 6, n,
        dimnames = list(paste0("g", 1:6), cells))
    emb <- matrix(stats::rnorm(n * 6), n, 6,
        dimnames = list(cells, paste0("Dim.", 1:6)))
    emb[, 1] <- emb[, 1] + 4 * region
    g <- GiottoClass::createGiottoObject(expression = m,
        spatial_locs = data.frame(sdimx = x, sdimy = y, cell_ID = cells),
        verbose = FALSE)
    g <- GiottoClass::setExpression(g, GiottoClass::createExprObj(log1p(m),
        name = "normalized"), name = "normalized", verbose = FALSE)
    g <- GiottoClass::setDimReduction(g, GiottoClass::createDimObj(
        coordinates = emb, name = "pca", spat_unit = "cell",
        feat_type = "rna", method = "pca"), spat_unit = "cell",
        feat_type = "rna", reduction = "cells", reduction_method = "pca",
        name = "pca", verbose = FALSE)
    list(g = g, m = log1p(m), emb = emb, x = x, y = y, region = region)
}

test_that("banksyParam checks its arguments", {
    expect_equal(banksyParam(lambda = 0.5)$lambda, 0.5)
    expect_error(banksyParam(lambda = 2))
    expect_error(banksyParam(npcs = 0))
})

test_that("runBanksy on expression equals Banksy called directly", {
    skip_if_not_installed("Banksy")
    skip_if_not_installed("SpatialExperiment")
    fx <- .banksy_fixture()
    g <- runBanksy(fx$g, feats_to_use = NULL, npcs = 5, resolution = 0.3,
        cluster_method = "banksy", verbose = FALSE)
    pcs <- GiottoClass::getDimReduction(g, reduction_method = "banksy",
        name = "banksy", output = "matrix")
    se <- SpatialExperiment::SpatialExperiment(assays = list(e = fx$m),
        spatialCoords = cbind(sdimx = fx$x, sdimy = fx$y))
    se <- Banksy::computeBanksy(se, assay_name = "e", k_geom = 15,
        verbose = FALSE)
    se <- Banksy::runBanksyPCA(se, lambda = 0.8, npcs = 5, seed = 1234)
    se <- Banksy::clusterBanksy(se, lambda = 0.8, npcs = 5,
        k_neighbors = 30, resolution = 0.3, seed = 1234)
    ref_pcs <- SingleCellExperiment::reducedDim(se, "PCA_M0_lam0.8")
    expect_equal(unname(pcs[colnames(fx$m), ]), unname(ref_pcs),
        tolerance = 1e-10, ignore_attr = TRUE)
    cd <- SummarizedExperiment::colData(se)
    ref <- as.character(cd[[grep("^clust_", colnames(cd))]])
    md <- GiottoClass::pDataDT(g)
    expect_identical(md$banksy[match(colnames(fx$m), md$cell_ID)], ref)
})

test_that("signed input warns and the default clustering is Giotto's", {
    skip_if_not_installed("Banksy")
    skip_if_not_installed("Giotto")
    fx <- .banksy_fixture()
    expect_warning(
        g <- runBanksy(fx$g, input = "dim_reduction", npcs = 5,
            resolution = 0.3, verbose = FALSE),
        "absolute value")
    md <- GiottoClass::pDataDT(g)
    expect_false(anyNA(md$banksy))
    expect_true("banksy" %in% GiottoClass::list_nearest_networks(g)$name)
    expect_error(runBanksy(fx$g, input = "dim_reduction",
        dimensions_to_use = 50:60), "no stored dimension")
})

test_that("runBanksy conforms to the extension contract", {
    skip_if_not_installed("Banksy")
    skip_if_not_installed("Giotto")
    fx <- .banksy_fixture()
    g <- runBanksy(fx$g, feats_to_use = NULL, npcs = 5, verbose = FALSE)
    expect_giotto_method(fx$g, g,
        slots = c("dimension_reduction", "nn_network", "cell_metadata"),
        param_class = "banksyParam")
})
