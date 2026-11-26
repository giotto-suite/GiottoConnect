# SPARK-X connection. The offline tests cover the param; the SPARK-X runs are
# opt-in (GIOTTOCONNECT_TOOL_TESTS).

test_that("sparkxParam stores SPARK-X's options", {
    p <- sparkxParam(option = "single", n_cores = 2)
    expect_s4_class(p, "analyzeParam")
    expect_equal(p$option, "single")
    expect_identical(p$n_cores, 2L)
    expect_identical(sparkxParam()$option, "mixture")
    expect_error(sparkxParam(option = "gaussian"))
    expect_error(sparkxParam(n_cores = 0))
})

test_that("analyzeData finds a spatial gene", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("SPARK")
    set.seed(42)
    locs <- matrix(stats::runif(400), ncol = 2)
    x <- rbind(g1 = stats::rpois(200, 1 + 10 * locs[, 1]),
        g2 = stats::rpois(200, 3))
    res <- analyzeData(x, sparkxParam(), spatial_locs = locs)
    expect_equal(res$feat_ID, c("g1", "g2"))
    expect_lt(res$padj[1], 1e-6)
    expect_gt(res$pval[2], 0.01)
})

test_that("runSparkX stores SPARK-X's p-values on the visium mini", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("SPARK")
    skip_if_not_installed("GiottoData")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("visium")
    g2 <- runSparkX(g, feats_to_use = "hvf")
    fm <- fDataDT(g2)
    expect_true(all(c("sparkx_pval", "sparkx_padj") %in% names(fm)))
    tested <- fm$hvf == "yes"
    expect_false(anyNA(fm$sparkx_pval[tested]))
    expect_true(all(is.na(fm$sparkx_pval[!tested])))
    chk <- checkGiottoMethod(g, g2, slots = "feat_metadata",
        param_class = "sparkxParam", fn = "runSparkX", package = "GiottoConnect")
    expect_true(all(chk$pass))
})
