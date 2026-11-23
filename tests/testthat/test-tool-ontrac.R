# ONTraC connection. The offline tests cover the command line and the file
# reader; the full run needs Python and is opt-in (GIOTTOCONNECT_TOOL_TESTS).

test_that("ontracParam builds ONTraC's command line options", {
    p <- ontracParam(k = 5, n_neighbors = 30, "trajectory-construct" = "BF")
    a <- GiottoConnect:::.ontrac_args(p)
    expect_equal(a[match("--k-clusters", a) + 1L], "5")
    expect_equal(a[match("--n-neighbors", a) + 1L], "30")
    expect_equal(a[match("--trajectory-construct", a) + 1L], "BF")
    expect_false(any(c("--version", "--out_dir") %in% a))
    expect_error(ontracParam(k = 9))
})

test_that("ONTraC's output files become a per-cell table", {
    out <- tempfile("ontrac_test_"); dir.create(out)
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    dir.create(file.path(out, "NT")); dir.create(file.path(out, "GNN"))
    ids <- c("a", "b", "c")
    data.table::fwrite(data.table::data.table(Cell_ID = ids, x = 1:3, y = 1:3,
        Niche_NTScore = c(.1, .5, .9), Cell_NTScore = c(.2, .4, .8)),
        file.path(out, "NT", "NTScore.csv.gz"))
    data.table::fwrite(data.table::data.table(Cell_ID = ids,
        NicheCluster_0 = c(.9, .2, .1), NicheCluster_1 = c(.1, .8, .9)),
        file.path(out, "GNN", "cell_level_niche_cluster.csv.gz"))
    data.table::fwrite(data.table::data.table(c(0.3, 0.7)),
        file.path(out, "NT", "niche_cluster_score.csv.gz"), col.names = FALSE)
    data.table::fwrite(data.table::data.table(c(0, .5), c(.5, 0)),
        file.path(out, "GNN", "consolidate_out_adj.csv.gz"), col.names = FALSE)

    res <- GiottoConnect:::.ontrac_read(out)
    expect_equal(res$cell_ID, ids)
    expect_equal(res$NTScore, c(.2, .4, .8))
    expect_equal(res$niche_cluster,
        c("NicheCluster_0", "NicheCluster_1", "NicheCluster_1"))
    expect_equal(attr(res, "niche_clusters")$NTScore, c(.3, .7))
    expect_equal(dim(attr(res, "connectivity")), c(2L, 2L))
})

test_that("runONTraC writes cell metadata and a spatial enrichment", {
    skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "",
        "set GIOTTOCONNECT_TOOL_TESTS to run external tools")
    skip_if_not_installed("GiottoData")
    options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
        giotto.no_python_warn = TRUE)
    g <- GiottoData::loadGiottoMini("vizgen")
    g <- runONTraC(g, cell_type = "leiden_clus", spat_unit = "aggregate",
        param = ontracParam(n_neighbors = 30, n_local = 10, epochs = 100))
    md <- GiottoClass::pDataDT(g, spat_unit = "aggregate")
    expect_true(all(c("ontrac_NTScore", "ontrac_niche_cluster") %in% names(md)))
    expect_false(anyNA(md$ontrac_NTScore))
    enr <- GiottoClass::getSpatialEnrichment(g, spat_unit = "aggregate", name = "ontrac",
        output = "spatEnrObj")
    expect_setequal(enr[]$cell_ID, md$cell_ID)
    expect_true(is.matrix(enr@misc$connectivity))
    h <- g@parameters[[length(g@parameters)]]
    expect_s4_class(attr(h, "attachments")$param, "ontracParam")
})
