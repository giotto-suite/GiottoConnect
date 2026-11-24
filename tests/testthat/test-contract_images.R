img_env <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            cache <<- new.env()
            sys.source(contract_script("build_fixtures.R"), envir = cache)
        }
        cache
    }
})

# a giotto object with an image whose red / green channels encode x / y, so a
# misplaced image shows as a colour difference at the spots
giotto_with_image <- function() {
    g <- img_env()$build_giotto_fixture()
    sl <- GiottoClass::getSpatialLocations(g, output = "data.table")
    e <- c(min(sl$sdimx) - 50, max(sl$sdimx) + 50, min(sl$sdimy) - 50, max(sl$sdimy) + 50)
    nr <- 120L
    nc <- 150L
    a <- array(0, c(nr, nc, 3L))
    a[, , 1] <- matrix(seq(0, 255, length.out = nc), nr, nc, byrow = TRUE)
    a[, , 2] <- matrix(seq(0, 255, length.out = nr), nr, nc)
    a[, , 3] <- 128
    r <- terra::rast(a, extent = terra::ext(e))
    im <- GiottoClass::createGiottoLargeImage(r, name = "image", use_rast_ext = TRUE, verbose = FALSE)
    GiottoClass::setGiotto(g, im, verbose = FALSE)
}

spot_colours <- function(g, name) {
    sl <- as.data.frame(GiottoClass::getSpatialLocations(g, output = "data.table"))
    sl <- sl[order(sl$cell_ID), ]
    r <- GiottoConnect:::.giotto_image_rgb(GiottoClass::getGiottoImage(g, name = name))
    as.matrix(terra::extract(r, as.matrix(sl[, c("sdimx", "sdimy")])))
}

test_that("images are resampled onto a grid from the origin", {
    r <- terra::rast(array(1, c(10L, 20L, 3L)), extent = terra::ext(100, 300, -250, -150))
    out <- GiottoConnect:::.image_to_origin(r)
    expect_equal(out$scale, 0.1)
    expect_equal(dim(out$image), c(25L, 30L, 3L))
    # padding between the origin and the image is white
    expect_equal(out$image[1, 1, 1], 1)
    expect_error(GiottoConnect:::.image_to_origin(terra::rast(array(1, c(2L, 2L, 3L)),
        extent = terra::ext(-10, -1, -5, -1))), "outside")
})

test_that("the rule 'first' can prefer a name", {
    kids <- data.frame(path = c("images/alignment", "images/image"), class = "giottoLargeImage")
    pick <- GiottoConnect:::.param_pick
    expect_equal(pick(list(rule = "first", prefer = "^image$"), kids, NULL), "image")
    expect_equal(pick(list(rule = "first", prefer = "^hires$"), kids, NULL), "alignment")
})

test_that("an image goes to Seurat as a VisiumV2 and comes back in place", {
    skip_if_not_installed("Seurat")
    g <- giotto_with_image()
    s <- giottoToSeurat(g)
    expect_s4_class(s[["fov"]], "VisiumV2")
    # spots are drawn at the spot size, as for Seurat's own Visium reader
    expect_equal(s[["fov"]][["centroids"]]@radius, s[["fov"]]@scale.factors$spot)
    expect_true(attr(readbackCheck(g, s), "pass"))
    g2 <- seuratToGiotto(s)
    expect_equal(GiottoClass::list_images(g2)$name, "fov")
    expect_lt(mean(abs(spot_colours(g2, "fov") - spot_colours(g, "image"))), 0.02)
})

test_that("an image goes to SpatialExperiment imgData and comes back in place", {
    skip_if_not_installed("SpatialExperiment")
    g <- giotto_with_image()
    x <- giottoToSpatialExperiment(g, sfe = "no")
    expect_equal(SpatialExperiment::imgData(x)$image_id, "image")
    expect_true(attr(readbackCheck(g, x), "pass"))
    g2 <- spatialExperimentToGiotto(x)
    expect_lt(mean(abs(spot_colours(g2, "image") - spot_colours(g, "image"))), 0.02)
})

test_that("an image goes to AnnData uns/spatial and comes back in place", {
    skip_if_not_installed("anndataR")
    g <- giotto_with_image()
    f <- tempfile(fileext = ".h5ad")
    giottoToAnnData(g, f)
    expect_true(attr(readbackCheck(g, f), "pass"))
    ad <- anndataR::read_h5ad(f)
    expect_length(dim(ad$uns$spatial$image$images$hires), 3L)
    g2 <- anndataToGiotto(f)
    expect_lt(mean(abs(spot_colours(g2, "image") - spot_colours(g, "image"))), 0.02)
})

test_that("Seurat Visium images are placed in the units of the spots", {
    skip_if_not_installed("Seurat")
    s <- img_env()$build_seurat_visium_v1_fixture()
    g <- seuratToGiotto(s)
    im <- GiottoClass::getGiottoImage(g, name = "slice1")
    sf <- s[["slice1"]]@scale.factors
    # 8 lowres pixels at the lowres scale, in microns (the fixture has 0.885
    # um per full-resolution pixel)
    expect_equal(unname(terra::ext(im@raster_object)[2]), 8 / sf$lowres * 0.885, tolerance = 1e-3)
})

test_that("a SpatialFeatureExperiment holds the image as a SpatRasterImage", {
    skip_if_not_installed("SpatialFeatureExperiment")
    g <- giotto_with_image()
    x <- giottoToSpatialExperiment(g, sfe = "yes")
    im <- SpatialExperiment::imgData(x)$data[[1L]]
    expect_s4_class(im, "SpatRasterImage")
    expect_true(attr(readbackCheck(g, x), "pass"))
    # SFE subsetting crops the image
    expect_no_error(x[, 1:20])
    g2 <- spatialExperimentToGiotto(x)
    expect_lt(mean(abs(spot_colours(g2, "image") - spot_colours(g, "image"))), 0.02)
})
