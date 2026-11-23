test_that("the bundled contract describes every mappable part", {
    ct <- giottoContract("bundled")
    expect_s3_class(ct, "giottoContract")
    expect_setequal(names(ct$parts), c("expression", "cell_metadata", "feat_metadata",
        "spatial_locs", "spatial_info", "dim_reduction", "nn_network", "spatial_network",
        "feature_info", "images"))
    expect_equal(ct$orientation$expression, "features_x_cells")
    expect_equal(ct$parts$expression$nesting, c("spat_unit", "feat_type", "name"))
    expect_equal(ct$parts$expression$class, "exprObj")
    disk <- vapply(ct$parts, function(p) isTRUE(p$disk_backable), logical(1))
    expect_equal(names(disk)[disk], c("expression", "spatial_info", "nn_network", "spatial_network",
        "feature_info"))
    expect_output(print(ct), "\\[disk\\]")
})

test_that("the live contract matches the bundled one for the same GiottoClass", {
    skip_if_not_installed("GiottoClass")
    ver <- as.character(utils::packageVersion("GiottoClass"))
    bundled <- giottoContract("bundled")
    skip_if(bundled$giottoclass_version != ver, "no bundled contract for this GiottoClass")
    live <- giottoContract("live")
    strip <- function(p) jsonlite::fromJSON(jsonlite::toJSON(p, auto_unbox = TRUE))
    expect_equal(strip(live$parts), strip(bundled$parts))
    expect_silent(giottoContract("auto"))
})
