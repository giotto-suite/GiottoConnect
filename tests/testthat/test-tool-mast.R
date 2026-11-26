# MAST connection. The offline tests cover the param and the argument checks;
# the MAST runs are opt-in (GIOTTOCONNECT_TOOL_TESTS).

test_that("mastParam stores MAST's options", {
    p <- mastParam(method = "glm", ebayes = FALSE)
    expect_s4_class(p, "analyzeParam")
    expect_equal(p$method, "glm")
    expect_false(p$ebayes)
    expect_equal(mastParam()$method, "bayesglm")
    expect_true(mastParam()$ebayes)
    expect_error(mastParam(method = "lm"))
})

test_that("runMAST needs both sides of a comparison", {
    g <- methods::new("giotto")
    expect_error(runMAST(g, cluster_column = "clus", group_1 = 1),
        "both `group_1` and `group_2`")
})

test_that("analyzeData tests one group against another", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("MAST")
    set.seed(1)
    n <- 60
    groups <- rep(c("a", "b"), each = n / 2)
    x <- rbind(
        up = c(stats::rnorm(n / 2, 3, 0.3), stats::rnorm(n / 2, 1, 0.3)),
        flat = stats::rnorm(n, 2, 0.3)
    )
    colnames(x) <- paste0("c", seq_len(n))
    names(groups) <- colnames(x)
    res <- analyzeData(x, mastParam(), groups = groups, test = "a", ref = "b")
    expect_equal(res$feats[1], "up")
    expect_gt(res$coef[res$feats == "up"], 1)
    expect_lt(res$fdr[res$feats == "up"], 1e-6)
    expect_equal(unique(res$cluster), "a_vs_b")
    expect_error(analyzeData(x, mastParam(), groups = groups, test = "a",
        ref = "z"), "no cells in group 'z'")
})

test_that("runMAST labels comparisons and keeps the best features", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("MAST")
    skip_if_not_installed("GiottoData")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("visium")
    pw <- runMAST(g, cluster_column = "leiden_clus", group_1 = 1,
        group_2 = c(2, 3))
    expect_equal(unique(pw$cluster), "1_vs_2_3")
    expect_equal(nrow(pw), nrow(fDataDT(g)))
    expect_false(is.unsorted(pw$fdr))
    chk <- checkGiottoMethod(g, g, slots = "none", result = pw,
        param_class = "mastParam", fn = "runMAST", package = "GiottoConnect")
    expect_true(all(chk$pass))

    ova <- runMAST(g, cluster_column = "leiden_clus", subset_clusters = 1:3,
        min_feats = 5)
    expect_setequal(unique(ova$cluster),
        c("1_vs_others", "2_vs_others", "3_vs_others"))
    expect_true(all(tapply(ova$ranking, ova$cluster, max) >= 5))
    expect_false(anyNA(ova$feats))
    expect_true(all(ova$ranking <= 5 | (ova$fdr < 0.001 & ova$coef > 1)))
})
