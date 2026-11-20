fixtures_env <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            cache <<- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = cache)
        }
        cache
    }
})

convert_and_check <- function(x, expect_id) {
    f <- findMapping(profileObject(x), ecosystem = "seurat")
    expect_equal(f$mapping$mapping$id, expect_id)
    expect_equal(f$status, "match")
    g <- seuratToGiotto(x)
    expect_true(attr(readbackCheck(x, g), "pass"))
    g
}

test_that("VisiumV1 objects convert, with v3 and v5 assays", {
    skip_if_not_installed("Seurat")
    env <- fixtures_env()
    for (cls in c("Assay", "Assay5")) {
        x <- env$build_seurat_visium_v1_fixture(assay_class = cls, um_per_px = 0.885)
        g <- convert_and_check(x, "seurat-visium-v1")
        sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
        # neighbouring spots within a row are one 100 um pitch apart
        row1 <- sl[abs(sl$sdimy - max(sl$sdimy)) < 1e-6, ]
        expect_equal(unique(round(diff(sort(row1$sdimx)), 6)), 100)
        l <- mappingProvenance(g)$losses
        expect_match(l$note[l$part_id == "locs" & l$field == "units"], "100 um spot pitch")
    }
})

test_that("Slide-seq objects convert", {
    skip_if_not_installed("Seurat")
    x <- fixtures_env()$build_seurat_slideseq_fixture()
    g <- convert_and_check(x, "seurat-slideseq")
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    ct <- SeuratObject::GetTissueCoordinates(x[["puck"]])
    sl <- sl[match(ct$cells, sl$cell_ID), ]
    expect_equal(sl$sdimx, ct$x)
    expect_equal(sl$sdimy, -ct$y)
    expect_setequal(names(sl), c("sdimx", "sdimy", "cell_ID"))
})

test_that("the family mapping wins over the generic one", {
    skip_if_not_installed("Seurat")
    x <- fixtures_env()$build_seurat_visium_v1_fixture()
    x@reductions <- list()
    x@graphs <- list()
    f <- findMapping(profileObject(x), ecosystem = "seurat")
    expect_equal(f$mapping$mapping$id, "seurat-visium-v1")
    # the generic mapping applies too, but reads less of the object
    cand <- f$candidates
    expect_true("seurat-assay" %in% cand$id[cand$status != "none"])
    expect_lt(cand$coverage[cand$id == "seurat-assay"], cand$coverage[cand$id == "seurat-visium-v1"])
    g <- seuratToGiotto(x)
    expect_true("nn" %in% mappingProvenance(g)$absent)
})

test_that("Atera objects (Seurat 6 LoadAtera) use the FOV mapping", {
    f <- path.expand("~/Documents/Datasets/giottoconnect_cache/seurat6/seurat_atera_tiny.rds")
    skip_if_not(file.exists(f), "tiny Atera Seurat object not cached")
    skip_if_not_installed("SeuratObject")
    x <- readRDS(f)
    m <- findMapping(profileObject(x), ecosystem = "seurat")$mapping
    expect_equal(m$mapping$id, "seurat-assay5")
    g <- suppressWarnings(seuratToGiotto(x))
    expect_true(attr(readbackCheck(x, g), "pass"))
})
