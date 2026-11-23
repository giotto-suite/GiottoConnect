# Images between Giotto and the other ecosystems. Giotto places an image by
# its extent in spatial coordinates; Seurat (Visium), SpatialExperiment and
# AnnData (scanpy) place an RGB array by one scale factor (image pixels per
# coordinate unit) with the origin at the image's top-left corner (y down).

# RGB array (H x W x 3, values 0..1 or 0..255, or an R raster) -> SpatRaster
# of 0..255 values covering x in [0, W / scale * unit], y in
# [-H / scale * unit, 0] (Giotto: y up)
.rgb_to_giotto_raster <- function(a, scale, unit = 1) {
    if (inherits(a, "raster") || is.character(a)) {
        m <- as.matrix(a)
        rgb <- grDevices::col2rgb(m)
        a <- array(t(rgb), c(nrow(m), ncol(m), 3L))
    } else {
        a <- as.array(a)
        if (length(dim(a)) == 2L) a <- array(a, c(dim(a), 3L))
        a <- a[, , seq_len(min(3L, dim(a)[3L])), drop = FALSE]
        if (max(a, na.rm = TRUE) <= 1) a <- a * 255
    }
    r <- terra::rast(round(a))
    terra::ext(r) <- c(0, ncol(r) / scale * unit, -nrow(r) / scale * unit, 0)
    r
}

# Giotto large image -> 3-layer SpatRaster of 0..1 values (grey images are
# repeated over the three channels)
.giotto_image_rgb <- function(img) {
    if (!methods::is(img, "giottoLargeImage")) {
        stop(wrap_txt("only giottoLargeImage images can be exported"), call. = FALSE)
    }
    r <- img@raster_object
    r <- if (terra::nlyr(r) >= 3L) r[[1:3]] else r[[c(1L, 1L, 1L)]]
    mx <- suppressWarnings(as.numeric(img@max_intensity))
    if (!length(mx) || is.na(mx) || mx <= 0) mx <- max(terra::global(r, "max", na.rm = TRUE)[[1]])
    # integer images are 8- or 16-bit
    if (isTRUE(img@is_int)) mx <- if (mx <= 255) 255 else 65535
    r / mx
}

# Resample an RGB SpatRaster onto a grid from the origin (x >= 0, y <= 0), so
# that pixel (row i, column j) covers x / scale in [j - 1, j] and -y / scale in
# [i - 1, i]; areas without image are white. The long side has at most
# `max_px` pixels. Returns list(image = H x W x 3 array of 0..1, scale).
.image_to_origin <- function(r, max_px = 2000L) {
    e <- as.vector(terra::ext(r))
    if (e[["xmax"]] <= 0 || e[["ymin"]] >= 0) {
        stop(wrap_txt("the image lies outside x >= 0, y <= 0"), call. = FALSE)
    }
    native <- 1 / terra::res(r)[1L]
    scale <- min(native, max_px / max(e[["xmax"]], -e[["ymin"]]))
    nc <- ceiling(e[["xmax"]] * scale)
    nr <- ceiling(-e[["ymin"]] * scale)
    grid <- terra::rast(ncols = nc, nrows = nr, nlyrs = 3L,
        extent = terra::ext(0, nc / scale, -nr / scale, 0))
    a <- terra::as.array(terra::resample(r, grid, method = "bilinear"))
    a[is.na(a)] <- 1
    list(image = pmin(pmax(a, 0), 1), scale = scale)
}

# spot diameter in coordinate units, from the centre spacing (Visium: 55 um
# spots, 100 um apart)
.spot_diameter <- function(xy) 0.55 * .median_nn_spacing(as.matrix(xy))


# hooks ####

# Giotto image -> list(image, scale) on a grid from the origin
.hook_giotto_image_to_array <- function(node, params, ctx) {
    r <- .giotto_image_rgb(node)
    e <- as.vector(terra::ext(r))
    out <- .image_to_origin(r)
    shift <- e[["xmin"]] > 0 || e[["ymax"]] < 0
    ctx$add_loss("estimated", "image", sprintf(
        "%d x %d px at %.4g px per unit%s", dim(out$image)[2L], dim(out$image)[1L], out$scale,
        if (shift) "; the area between the origin and the image is padded white" else ""))
    out
}

# Seurat VisiumV1 / VisiumV2 lowres image, in the units of the locs part
.hook_seurat_visium_image <- function(node, params, ctx) {
    sf <- methods::slot(node, "scale.factors")
    um <- .visium_um_per_px(node, params, ctx)$um
    .rgb_to_giotto_raster(methods::slot(node, "image"), sf[["lowres"]], um)
}

# first SpatialExperiment image (hires preferred), in the units of the locs
# part: microns for VisiumIO layouts, as stored otherwise
.hook_spe_image <- function(node, params, ctx, to_um = FALSE) {
    SPE <- asNamespace("SpatialExperiment")
    d <- SPE$imgData(ctx$x)
    i <- match("hires", d$image_id)
    if (is.na(i)) i <- 1L
    unit <- if (to_um) .spe_um_per_px(ctx$x, SPE$spatialCoords(ctx$x))$um else 1
    im <- d$data[[i]]
    r <- if (methods::is(im, "SpatRasterImage")) {
        # SpatialFeatureExperiment: placed by its extent in spatialCoords
        # units, y down; Giotto: y up
        g <- methods::as(im, "SpatRaster")
        e <- as.vector(terra::ext(g))
        g <- terra::flip(g, direction = "vertical")
        g <- g[[seq_len(min(3L, terra::nlyr(g)))]]
        terra::ext(g) <- c(e[["xmin"]], e[["xmax"]], -e[["ymax"]], -e[["ymin"]]) * unit
        g
    } else {
        .rgb_to_giotto_raster(SPE$imgRaster(ctx$x, sample_id = d$sample_id[i], image_id = d$image_id[i]),
            d$scaleFactor[i], unit)
    }
    .hook_result(r, to = list(name = d$image_id[i]))
}

# scanpy uns/spatial/<library>: hires (else lowres) image and its scale factor
.hook_anndata_spatial_image <- function(node, params, ctx) {
    res <- intersect(c("hires", "lowres"), names(node$images))[1L]
    if (is.na(res)) stop(wrap_txt("uns/spatial has no hires or lowres image"), call. = FALSE)
    sf <- node$scalefactors[[paste0("tissue_", res, "_scalef")]]
    .rgb_to_giotto_raster(node$images[[res]], as.numeric(sf))
}
