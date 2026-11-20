fixture_giotto <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) cache <<- seuratToGiotto(seurat_fixture())
        cache
    }
})

test_that("the bundled mapping reads, summarises and declares its losses", {
    m <- readMapping("seurat-assay5")
    expect_equal(m$mapping$status, "verified")
    expect_match(m$mapping$fingerprint, "^[0-9a-f]{16}$")
    # `on:` keys survive YAML 1.1 boolean parsing
    expect_equal(m$params$layer_raw$on, "{assay}")
    s <- expect_output(summary(m), "Declared losses")
    expect_true(all(c("dropped", "estimated") %in% s$losses$kind))
})

test_that("hookRegistry lists the hooks and their packages", {
    h <- hookRegistry()
    expect_setequal(h$hook, c("seurat_graph_to_nn", "seurat_visium_to_um",
        "anndata_distances_to_nn", "giotto_polygons_to_vertices", "giotto_nn_to_graph"))
    expect_match(h$packages[h$hook == "seurat_visium_to_um"], "Seurat")
})

test_that("ops transform tables and matrices", {
    op <- GiottoConnect:::.apply_op
    df <- data.frame(x = 1:2, y = c(3, 4), row.names = c("a", "b"))
    expect_named(op(df, list(op = "select_cols", cols = "y")), "y")
    expect_named(op(df, list(op = "exclude_cols", cols = "y")), "x")
    expect_named(op(df, list(op = "rename", map = list(x = "sdimx"))), c("sdimx", "y"))
    expect_equal(op(df, list(op = "scale", factor = -1, cols = "y"))$y, c(-3, -4))
    expect_equal(op(df, list(op = "rownames_to_col", col = "cell_ID"))$cell_ID, c("a", "b"))
    same <- cbind(df, cell_ID = c("a", "b"))
    expect_equal(names(op(same, list(op = "rownames_to_col", col = "cell_ID"))), c("cell_ID", "x", "y"))
    other <- cbind(df, cell_ID = c("p", "q"))
    out <- op(other, list(op = "rownames_to_col", col = "cell_ID"))
    expect_equal(out$cell_ID, c("a", "b"))
    expect_equal(out$cell_ID.source, c("p", "q"))
    expect_equal(op(df, list(op = "set", field = "z", value = "k"))$z, c("k", "k"))
    expect_error(op(df, list(op = "rename", map = list(q = "r"))), "missing column")
    m <- matrix(c(1, 2, 3, 4), 2)
    expect_equal(op(m, list(op = "transpose")), t(m))
    expect_s4_class(op(m, list(op = "cast", to = "sparse")), "CsparseMatrix")
    expect_equal(storage.mode(op(m, list(op = "cast", to = "integer"))), "integer")
    expect_error(op(m + 0.5, list(op = "cast", to = "integer")), "not integer-valued")
})

test_that("params bind on the profile without the source package", {
    skip_if_not_installed("SeuratObject")
    m <- readMapping("seurat-assay5")
    b <- GiottoConnect:::.bind_params(m, profileObject(seurat_fixture()))
    expect_equal(b$values$assay, "Xenium")
    expect_equal(b$values$layer_raw, "counts")
    expect_equal(b$values$segmentation, "segmentation")
    expect_equal(b$values$nn_graph, "Xenium_nn")
    expect_length(b$unbound, 0L)
})

test_that("findMapping matches the fixture and absorbs renamed assays", {
    skip_if_not_installed("SeuratObject")
    f <- findMapping(profileObject(seurat_fixture()), ecosystem = "seurat")
    expect_equal(f$status, "match")
    expect_equal(f$mapping$mapping$id, "seurat-assay5")

    env <- new.env()
    sys.source(contract_script("build_fixtures.R"), envir = env)
    f2 <- findMapping(profileObject(env$build_seurat_fixture(assay = "Spatial")))
    expect_equal(f2$status, "match")
    expect_equal(f2$params$values$assay, "Spatial")

    none <- findMapping(profileObject(list(a = 1)))
    expect_equal(none$status, "none")
})

test_that("seuratToGiotto converts every mapped part in memory", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    g <- fixture_giotto()
    expect_s4_class(g, "giotto")

    raw <- GiottoClass::getExpression(g, values = "raw", output = "matrix")
    src <- SeuratObject::LayerData(s, "counts")
    expect_equal(as.matrix(raw[rownames(src), colnames(src)]), as.matrix(src),
        ignore_attr = TRUE)

    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    ct <- SeuratObject::GetTissueCoordinates(s[["fov"]], which = "centroids")
    sl <- sl[match(ct$cell, sl$cell_ID), ]
    expect_equal(sl$sdimx, ct$x)
    expect_equal(sl$sdimy, -ct$y)

    expect_length(GiottoClass::spatIDs(GiottoClass::getPolygonInfo(g, name = "cell",
        return_giottoPolygon = TRUE)), 200L)
    nn <- GiottoClass::list_nearest_networks(g)
    expect_equal(nn$nn_type, "kNN")
    expect_equal(nn$name, "kNN.pca")
    e <- GiottoClass::getNearestNetwork(g, nn_type = "kNN", name = "kNN.pca", output = "data.table")
    expect_false(any(e$from == e$to))
    emb <- SeuratObject::Embeddings(s, "pca")
    expect_equal(e$distance, sqrt(rowSums((emb[e$from, ] - emb[e$to, ])^2)), ignore_attr = TRUE)
    expect_equal(e$weight, 1 / (1 + e$distance))

    fm <- GiottoClass::fDataDT(g)
    expect_false(any(c("nCount_Xenium", "seurat_clusters") %in% names(fm)))
})

test_that("conversions are recorded with their losses", {
    skip_if_not_installed("SeuratObject")
    prov <- mappingProvenance(fixture_giotto())
    expect_equal(prov$mapping, "seurat-assay5")
    expect_equal(prov$fingerprint, readMapping("seurat-assay5")$mapping$fingerprint)
    k <- prov$losses[prov$losses$kind == "recovered_from_parameters", ]
    expect_equal(unique(k$part_id), "nn")
    expect_match(k$note[k$field == "name"], "FindNeighbors.Xenium.pca on 'pca'; k = 9 neighbours")
    expect_match(k$note[k$field == "distance"], "recomputed in pca dims 1:10 \\(euclidean\\)")
    expect_true("assays/Xenium/layers/scale.data" %in% prov$losses$field)
    expect_null(mappingProvenance(GiottoClass::giotto()))
})

test_that("seuratToGiotto writes disk-backable parts to a GiottoDisk project", {
    skip_if_not_installed("SeuratObject")
    skip_if_not_installed("GiottoDisk")
    s <- seurat_fixture()
    d <- withr_tempdir()
    g <- seuratToGiotto(s, source = d)
    store <- GiottoClass::getExpression(g, values = "raw")[]
    expect_s4_class(store, "parquetExprStore")
    src <- SeuratObject::LayerData(s, "counts")
    m <- methods::as(store, "dgCMatrix")
    expect_equal(as.matrix(m[rownames(src), colnames(src)]), as.matrix(src), ignore_attr = TRUE)
    expect_true(inherits(GiottoClass::getPolygonInfo(g, name = "cell",
        return_giottoPolygon = TRUE)[], "dataStore"))
    expect_true(inherits(GiottoClass::getSpatialLocations(g)@coordinates, "data.table"))
})

test_that("validateMapping reports unmapped data and missing sources", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    m <- readMapping("seurat-assay5")
    v <- validateMapping(m, profileObject(s))
    expect_s3_class(v, "giottoValidation")
    expect_false(any(v$level == "FAIL"))
    expect_false(any(v$level == "WARN"))

    s2 <- s
    s2@misc$extra <- matrix(1:4, 2)
    s2@reductions$pca <- NULL
    v2 <- validateMapping(m, profileObject(s2))
    expect_true(any(v2$check == "unmapped" & grepl("misc/extra", v2$message)))
    # reductions come from a for_each template: a missing one is simply absent
    expect_false("dimred_pca" %in% v2$part_id)
    expect_true("dimred_umap" %in% v2$part_id)

    m3 <- m
    m3$source_nodes$counts$path <- "assays/{assay}/layers/nope"
    v3 <- validateMapping(m3, profileObject(s))
    expect_equal(v3$level[v3$check == "source" & v3$part_id == "expr_raw"], "FAIL")
    expect_error(convertWithMapping(s, m3), "structural validation failed")
})

test_that("semantic validation checks the converted object", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    v <- validateMapping("seurat-assay5", profileObject(s), result = fixture_giotto(),
        checks = "semantic")
    expect_setequal(unique(v$check), c("raw_counts", "feat_ids", "dims", "cell_ids", "coord_units",
        "dim_reduction", "nn_network"))
    expect_true(all(v$level == "OK"))
})

test_that("networks are named from Seurat's FindNeighbors record when it exists", {
    skip_if_not_installed("SeuratObject")
    expect_equal(GiottoClass::list_nearest_networks(fixture_giotto())$name, "kNN.pca")

    s <- seurat_fixture()
    s@commands <- list()
    g <- seuratToGiotto(s)
    nn <- GiottoClass::list_nearest_networks(g)
    expect_equal(nn$name, "Xenium_nn")
    expect_equal(nn$nn_type, "kNN")
    l <- mappingProvenance(g)$losses
    expect_true(any(l$kind == "ambiguous" & l$field == "name"))
    expect_match(l$note[l$field == "k"], "k = 9 recovered from out-degree")
    expect_true(all(c("distance", "weight") %in% l$field[l$kind == "dropped"]))
    e <- GiottoClass::getNearestNetwork(g, nn_type = "kNN", name = "Xenium_nn", output = "data.table")
    expect_true(all(is.na(e$distance)))
})

test_that("Visium objects select the Visium mapping and convert pixels to microns", {
    skip_if_not_installed("Seurat")
    env <- new.env()
    sys.source(contract_script("build_fixtures.R"), envir = env)
    s <- env$build_seurat_visium_fixture(n_side = 6L, um_per_px = 0.25)
    f <- findMapping(profileObject(s), ecosystem = "seurat")
    expect_equal(f$mapping$mapping$id, "seurat-visium")
    expect_equal(f$status, "match")
    cand <- f$candidates
    expect_lt(cand$coverage[cand$id == "seurat-assay5"], cand$coverage[cand$id == "seurat-visium"])

    g <- seuratToGiotto(s)
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    expect_equal(sort(unique(round(diff(sort(unique(sl$sdimx))), 6))), 16)
    expect_true(all(sl$sdimy < 0))
    l <- mappingProvenance(g)$losses
    expect_match(l$note[l$part_id == "locs" & l$field == "units"], "0.25000 um per px")
    expect_equal(GiottoClass::list_nearest_networks(g)$name, "kNN.pca")
})

test_that("v3/v4 objects convert through seurat-assay and read back", {
    skip_if_not_installed("SeuratObject")
    env <- new.env()
    sys.source(contract_script("build_fixtures.R"), envir = env)
    s <- env$build_seurat_v3_fixture()
    f <- findMapping(profileObject(s), ecosystem = "seurat")
    expect_equal(f$mapping$mapping$id, "seurat-assay")
    expect_equal(f$status, "match")
    g <- seuratToGiotto(s)
    rb <- readbackCheck(s, g)
    expect_true(attr(rb, "pass"))
    expect_true(all(rb$status == "MATCH"))
    expect_equal(GiottoClass::list_nearest_networks(g)$name, "kNN.pca")
})

test_that("the assay param binds the default assay only", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    s[["Extra"]] <- SeuratObject::CreateAssayObject(counts = SeuratObject::LayerData(s, layer = "counts"))
    p <- profileObject(s)
    expect_equal(p$nodes$value[p$nodes$path == "active.assay"], "Xenium")
    f <- findMapping(p, ecosystem = "seurat")
    expect_equal(f$mapping$mapping$id, "seurat-assay5")
    expect_equal(f$candidates$status[f$candidates$id == "seurat-assay"], "none")

    SeuratObject::DefaultAssay(s) <- "Extra"
    f2 <- findMapping(profileObject(s), ecosystem = "seurat")
    expect_equal(f2$mapping$mapping$id, "seurat-assay")
    expect_equal(f2$params$values$assay, "Extra")
})

test_that("for_each expands one part per matching child", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    s[["harmony"]] <- SeuratObject::CreateDimReducObject(
        embeddings = SeuratObject::Embeddings(s, "pca")[, 1:3], key = "harmony_", assay = "Xenium")
    s[["spatial"]] <- SeuratObject::CreateDimReducObject(
        embeddings = `colnames<-`(SeuratObject::Embeddings(s, "umap"), c("SPATIAL_1", "SPATIAL_2")),
        key = "SPATIAL_", assay = "Xenium")
    m <- readMapping("seurat-assay5")
    p <- profileObject(s)
    x <- GiottoConnect:::.expand_for_each(m, p, GiottoConnect:::.bind_params(m, p))
    ids <- vapply(x$parts, `[[`, character(1), "id")
    expect_setequal(grep("^dimred_", ids, value = TRUE), c("dimred_pca", "dimred_umap", "dimred_harmony"))
    expect_equal(x$source_nodes$reduction__harmony$get$args$reduction, "harmony")
    expect_length(GiottoConnect:::.mapping_schema_check(unclass(x)), 0L)
    g <- seuratToGiotto(s)
    expect_true("harmony" %in% GiottoClass::list_dim_reductions(g)$name)
    expect_false("spatial" %in% GiottoClass::list_dim_reductions(g)$name)
    expect_true(attr(readbackCheck(s, g), "pass"))
})

test_that("declared drops are listed only when the object holds them", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    g <- seuratToGiotto(s)
    l <- mappingProvenance(g)$losses
    has_scale <- "scale.data" %in% SeuratObject::Layers(s)
    expect_identical(any(grepl("scale.data", l$field)), has_scale)
})

test_that("a spatial mapping without any coordinates warns", {
    skip_if_not_installed("SeuratObject")
    s <- seurat_fixture()
    s@images <- list()
    expect_warning(g <- seuratToGiotto(s, mapping = "seurat-assay5"), "no spatial locations found")
    expect_true(all(c("locs", "polygons") %in% mappingProvenance(g)$absent))
})

