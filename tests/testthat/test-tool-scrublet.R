# Scrublet connection. The offline tests cover the param; the Scrublet runs
# are opt-in (GIOTTOCONNECT_TOOL_TESTS) and need Python (installed by
# reticulate::py_require() on first use).

test_that("scrubletParam stores Scrublet's options", {
    p <- scrubletParam(expected_doublet_rate = 0.1, n_prin_comps = 20)
    expect_s4_class(p, "analyzeParam")
    expect_equal(p$expected_doublet_rate, 0.1)
    expect_identical(p$n_prin_comps, 20L)
    expect_identical(scrubletParam()$seed, 1234L)
    expect_error(scrubletParam(expected_doublet_rate = 2))
    expect_error(scrubletParam(n_prin_comps = 0))
})

test_that("runScrublet stores doublet scores on the visium mini", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("reticulate", minimum_version = "1.41")
    skip_if_not_installed("GiottoData")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("visium")
    # spots, not cells: Scrublet sets no threshold, so only scores are stored
    expect_warning(g1 <- runScrublet(g), "could not set a doublet score")
    cm <- pDataDT(g1)
    expect_true("scrublet_score" %in% names(cm))
    expect_false("scrublet_doublet" %in% names(cm))
    expect_false(anyNA(cm$scrublet_score))
    expect_true(all(cm$scrublet_score >= 0 & cm$scrublet_score <= 1))

    # a threshold gives the calls; the same seed gives the same scores
    g2 <- runScrublet(g, param = scrubletParam(threshold = 0.2))
    cm2 <- pDataDT(g2)
    expect_type(cm2$scrublet_doublet, "logical")
    expect_equal(cm2$scrublet_doublet, cm2$scrublet_score > 0.2)
    expect_equal(cm2$scrublet_score, cm$scrublet_score)

    chk <- checkGiottoMethod(g, g2, slots = "cell_metadata",
        param_class = "scrubletParam", fn = "runScrublet",
        package = "GiottoConnect")
    expect_true(all(chk$pass))
})
