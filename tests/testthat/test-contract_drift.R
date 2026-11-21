drift <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            env <- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = env)
            cache <<- list(variants = env$build_drift_variants(), expected = env$DRIFT_EXPECTED)
        }
        cache
    }
})

outcome <- function(x) {
    f <- findMapping(profileObject(x), ecosystem = "seurat")
    if (f$status == "none") return(list(find = f, repair = NULL))
    list(find = f, repair = repairMapping(f$mapping, x))
}

test_that("params absorb renamed assays, layers and graphs (D1, D2, D7, D8)", {
    skip_if_not_installed("SeuratObject")
    v <- drift()$variants
    for (d in c("D1", "D2", "D7", "D8")) {
        o <- outcome(v[[d]])
        expect_equal(o$find$status, "match", label = d)
        expect_equal(o$repair$status, "not_needed", label = d)
        expect_true(attr(o$repair$readback, "pass"), label = d)
    }
    expect_equal(findMapping(profileObject(v$D2))$params$values$layer_raw, "raw_counts")
    expect_equal(findMapping(profileObject(v$D7))$params$values$nn_graph, "knn_graph")
})

test_that("a moved coordinate source is read through the declared alternative (D3)", {
    skip_if_not_installed("SeuratObject")
    s <- drift()$variants$D3
    o <- outcome(s)
    # the required structure is unchanged; the conversion switches by itself
    expect_equal(o$find$status, "match")
    g0 <- seuratToGiotto(s)
    prov <- mappingProvenance(g0)
    expect_match(prov$alternatives, "^locs: switched to alternative 'locs_from_meta'")
    expect_true("polygons" %in% prov$absent)
    sl0 <- GiottoClass::getSpatialLocations(g0, output = "data.table")
    expect_equal(sl0$sdimx[match(colnames(s), sl0$cell_ID)], s$x_um, ignore_attr = TRUE)
    expect_true(attr(readbackCheck(s, g0), "pass"))

    # repair stamps the same switch into a new mapping version
    r <- o$repair
    expect_equal(r$status, "repaired")
    expect_equal(r$mapping$mapping$version, "0.1.1")
    expect_equal(r$mapping$mapping$status, "patched-unreviewed")
    expect_equal(r$mapping$mapping$history[[1]]$by, "rules")
    expect_setequal(r$changes$rule, c("alternative", "absent"))

    g <- seuratToGiotto(s, mapping = r$mapping)
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    sl <- sl[match(colnames(s), sl$cell_ID), ]
    expect_equal(sl$sdimx, s$x_um, ignore_attr = TRUE)
    expect_equal(sl$sdimy, -s$y_um, ignore_attr = TRUE)
    expect_equal(findMapping(profileObject(s), mappings = list(r$mapping))$status, "match")

    # the original source is kept as an alternative, so the patch repairs back
    back <- repairMapping(r$mapping, seurat_fixture())
    expect_equal(back$status, "repaired")
    expect_true(attr(back$readback, "pass"))
})

test_that("a storage change gets a cast (D4)", {
    skip_if_not_installed("SeuratObject")
    r <- outcome(drift()$variants$D4)$repair
    expect_equal(r$status, "repaired")
    expect_equal(r$changes$rule, "cast")
    expect_equal(r$mapping$parts[[1]]$ops[[1]], list(op = "cast", to = "sparse"))
    expect_true(attr(r$readback, "pass"))
})

test_that("an extra reduction is converted by the template (D5) and mm coordinates warn (D6)", {
    skip_if_not_installed("SeuratObject")
    v <- drift()$variants
    r5 <- outcome(v$D5)$repair
    expect_equal(r5$status, "not_needed")
    expect_false(any(r5$validation$check == "unmapped" & grepl("tsne", r5$validation$message)))
    expect_true("tsne" %in% GiottoClass::list_dim_reductions(seuratToGiotto(v$D5))$name)

    expect_warning(g6 <- seuratToGiotto(v$D6), "looks like mm")
    v6 <- validateMapping("seurat-assay5", profileObject(v$D6), result = g6, checks = "semantic")
    expect_equal(v6$level[v6$check == "coord_units"], "WARN")
})

test_that("v3/v4 assays use the seurat-assay mapping (D9)", {
    skip_if_not_installed("SeuratObject")
    s <- drift()$variants$D9
    f <- findMapping(profileObject(s), ecosystem = "seurat")
    expect_equal(f$mapping$mapping$id, "seurat-assay")
    # optional parts (FOV, reductions, graphs) are absent, not changed
    expect_equal(f$status, "match")
    r <- repairMapping(f$mapping, s)
    expect_equal(r$status, "not_needed")
    expect_true(all(r$changes$rule[r$changes$rule != "added"] == "absent"))
    expect_true(attr(r$readback, "pass"))
    g <- seuratToGiotto(s)
    expect_equal(mappingProvenance(g)$mapping, "seurat-assay")
})

test_that("objects without required data are refused, not mis-converted", {
    skip_if_not_installed("SeuratObject")
    SO <- asNamespace("SeuratObject")
    norm <- log1p(SO$LayerData(seurat_fixture(), layer = "counts"))
    s <- SO$CreateSeuratObject(SO$CreateAssay5Object(data = norm), assay = "RNA")
    f <- findMapping(profileObject(s), ecosystem = "seurat")
    expect_equal(f$mapping$mapping$id, "seurat-assay5")
    expect_true("layer_raw" %in% f$params$unbound)
    expect_error(suppressWarnings(seuratToGiotto(s)), "structural validation failed")
    expect_equal(repairMapping(f$mapping, s)$status, "failed")

    # an object of a class no mapping covers
    none <- findMapping(profileObject(list(a = 1)), ecosystem = "seurat")
    expect_equal(none$status, "none")
    expect_error(seuratToGiotto(list(a = 1)), "No bundled mapping applies")
})

test_that("a required part that cannot be repaired fails with candidates", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    a <- s[["Xenium"]]
    a$counts <- NULL
    s[["Xenium"]] <- suppressWarnings(a)
    s@misc$stash <- SeuratObject::LayerData(seurat_fixture(), layer = "counts")
    r <- repairMapping("seurat-assay5", s)
    expect_equal(r$status, "failed")
    expect_equal(r$mapping$mapping$status, "failed")
    expect_true("expr_raw" %in% r$unresolved)
    expect_equal(r$candidates$path[r$candidates$part_id == "expr_raw"][1], "misc/stash")
})

test_that("repaired mappings are written as readable YAML", {
    skip_if_not_installed("SeuratObject")
    f <- tempfile(fileext = ".yaml")
    r <- repairMapping("seurat-assay5", drift()$variants$D4, file = f)
    m <- readMapping(f)
    expect_equal(m$mapping$version, "0.1.1")
    expect_equal(m$params$layer_raw$on, "{assay}")
})
