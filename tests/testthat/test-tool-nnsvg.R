# nnSVG connection. The offline tests cover the param and the feature
# selection; the nnSVG runs are opt-in (GIOTTOCONNECT_TOOL_TESTS).

test_that("nnsvgParam stores nnSVG's options", {
    p <- nnsvgParam(n_neighbors = 5, order = "Sum_coords")
    expect_s4_class(p, "analyzeParam")
    expect_identical(p$n_neighbors, 5L)
    expect_equal(p$order, "Sum_coords")
    expect_identical(nnsvgParam()$order, "AMMD")
    expect_identical(nnsvgParam()$seed, 1L)
    expect_null(nnsvgParam(seed = NULL)$seed)
    expect_error(nnsvgParam(n_neighbors = 0))
    expect_error(nnsvgParam(order = "random"))
})

test_that("feats_to_use selects by metadata column or IDs", {
    fm <- data.table::data.table(feat_ID = c("a", "b", "c"),
        hvf = c("yes", "no", "yes"), keep = c(FALSE, TRUE, TRUE))
    f <- GiottoConnect:::.nnsvg_feats
    expect_equal(f(NULL, fm, c("a", "b", "c")), c("a", "b", "c"))
    expect_equal(f("hvf", fm, c("a", "b", "c")), c("a", "c"))
    expect_equal(f("keep", fm, c("a", "b", "c")), c("b", "c"))
    expect_equal(f(c("c", "b", "z"), fm, c("a", "b", "c")), c("b", "c"))
    expect_error(f("z", fm, c("a", "b", "c")), "no features")
})

test_that("mean, var and spcov follow nnSVG's logcounts definitions", {
    x <- Matrix::Matrix(c(0, 1, 2, 3, 1, 1, 4, 0), nrow = 2, sparse = TRUE)
    res <- cbind(sigma.sq = c(1, 4), mean = NA, var = NA, spcov = NA)
    out <- GiottoConnect:::.nnsvg_add_moments(res, x)
    expect_equal(out[, "mean"], Matrix::rowMeans(x))
    expect_equal(out[, "var"], apply(as.matrix(x), 1, stats::var))
    expect_equal(out[, "spcov"], sqrt(c(1, 4)) / Matrix::rowMeans(x))
})

test_that("analyzeData is reproducible and keeps the caller's random state", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("nnSVG")
    set.seed(42)
    locs <- matrix(stats::runif(200), ncol = 2)
    x <- rbind(g1 = locs[, 1] * 3 + stats::rnorm(100, sd = 0.3),
        g2 = stats::rnorm(100))
    set.seed(7)
    before <- .Random.seed
    a <- analyzeData(x, nnsvgParam(n_neighbors = 5), spatial_locs = locs)
    expect_identical(.Random.seed, before)
    b <- analyzeData(x, nnsvgParam(n_neighbors = 5), spatial_locs = locs)
    expect_equal(a$LR_stat, b$LR_stat)
    expect_equal(a$feat_ID, c("g1", "g2"))
    expect_lt(a$padj[1], 0.05)
})

test_that("runNnSVG stores nnSVG's statistics on the visium mini", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("nnSVG")
    skip_if_not_installed("GiottoData")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("visium")
    g2 <- runNnSVG(g, feats_to_use = "hvf", param = nnsvgParam(n_neighbors = 5))
    fm <- fDataDT(g2)
    expect_true(all(c("nnSVG_LR_stat", "nnSVG_rank", "nnSVG_pval",
        "nnSVG_padj", "nnSVG_prop_sv", "nnSVG_sigma_sq") %in% names(fm)))
    tested <- fm$hvf == "yes"
    expect_false(anyNA(fm$nnSVG_LR_stat[tested]))
    expect_true(all(is.na(fm$nnSVG_LR_stat[!tested])))
    expect_setequal(fm$nnSVG_rank[tested], seq_len(sum(tested)))
    chk <- checkGiottoMethod(g, g2, slots = "feat_metadata",
        param_class = "nnsvgParam", fn = "runNnSVG", package = "GiottoConnect")
    expect_true(all(chk$pass))
})
