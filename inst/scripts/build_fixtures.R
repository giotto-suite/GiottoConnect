# Small deterministic Seurat v5 test objects, built in code.
# Sourced by the tests (tests/testthat/helper-contract-fixtures.R), by
# inst/scripts/verify_mapping.R and by examples/seurat_comparison/.
#   build_seurat_fixture()         Assay5 + FOV (centroids, segmentation): SeuratObject
#   build_seurat_visium_fixture()  Assay5 + VisiumV2 image (Visium HD bins): needs Seurat

build_seurat_fixture <- function(n_cells = 200L, n_feat = 50L, assay = "Xenium",
                                 seed = 1L) {
    stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
    SO <- asNamespace("SeuratObject")
    set.seed(seed)
    s <- .fx_expression(n_cells, n_feat, assay)
    cells <- colnames(s)

    # FOV: centroids + one square segmentation per cell, microns in 0-2000
    xy <- data.frame(x = stats::runif(n_cells, 0, 2000),
        y = stats::runif(n_cells, 0, 2000), cell = cells)
    sq <- do.call(rbind, lapply(seq_len(n_cells), function(i) {
        data.frame(x = xy$x[i] + c(-5, 5, 5, -5), y = xy$y[i] + c(-5, -5, 5, 5),
            cell = cells[i])
    }))
    s[["fov"]] <- SO$CreateFOV(
        list(centroids = SO$CreateCentroids(xy), segmentation = SO$CreateSegmentation(sq)),
        type = c("centroids", "segmentation"), assay = assay
    )
    .fx_reductions_graphs(s, assay)
}

# Visium HD 16 um bins on a square grid, coordinates in full-resolution
# pixels as Load10X_Spatial(bin.size = 16) stores them
build_seurat_visium_fixture <- function(n_side = 12L, n_feat = 50L, bin_um = 16L,
                                        um_per_px = 0.2737, seed = 1L) {
    stopifnot(requireNamespace("Seurat", quietly = TRUE))
    SO <- asNamespace("SeuratObject")
    set.seed(seed)
    assay <- sprintf("Spatial.%03dum", bin_um)
    n <- n_side^2
    s <- .fx_expression(n, n_feat, assay,
        cells = sprintf("s_%03dum_%05d_%05d-1", bin_um, rep(seq_len(n_side), n_side),
            rep(seq_len(n_side), each = n_side)))
    bin_px <- bin_um / um_per_px
    coords <- data.frame(
        imagecol = 500 + (rep(seq_len(n_side), n_side) - 1) * bin_px,
        imagerow = 300 + (rep(seq_len(n_side), each = n_side) - 1) * bin_px,
        row.names = colnames(s)
    )
    sf <- Seurat::scalefactors(spot = bin_px, fiducial = 20 * bin_px, hires = 0.08, lowres = 0.008)
    fov <- SO$CreateFOV(coords[, c("imagecol", "imagerow")], type = "centroids",
        radius = sf[["spot"]], assay = assay, key = "slice1_")
    img <- methods::new("VisiumV2", boundaries = fov@boundaries, molecules = fov@molecules,
        assay = fov@assay, key = fov@key, image = array(0.5, c(8L, 8L, 3L)),
        scale.factors = sf, coords_x_orientation = "horizontal")
    s[[sprintf("slice1.%03dum", bin_um)]] <- img
    .fx_reductions_graphs(s, assay)
}

# counts / data / scale.data layers and cell metadata
.fx_expression <- function(n_cells, n_feat, assay,
                           cells = sprintf("cell_%03d", seq_len(n_cells))) {
    SO <- asNamespace("SeuratObject")
    feats <- sprintf("gene%02d", seq_len(n_feat))
    counts <- matrix(stats::rpois(n_cells * n_feat, 2), n_feat, n_cells,
        dimnames = list(feats, cells))
    counts <- methods::as(counts, "CsparseMatrix")
    norm <- log1p(counts)
    scaled <- t(scale(t(as.matrix(norm))))
    scaled[is.na(scaled)] <- 0

    s <- SO$CreateSeuratObject(SO$CreateAssay5Object(counts = counts), assay = assay)
    s <- SO$`LayerData<-`(s, layer = "data", value = norm)
    s <- SO$`LayerData<-`(s, layer = "scale.data", value = scaled)
    s$seurat_clusters <- factor(sample(0:3, n_cells, replace = TRUE))
    s
}

# pca (10 dims), umap (2 dims), binary kNN (k = 10 incl. self), weighted SNN,
# and the FindNeighbors command record Seurat writes for them
.fx_reductions_graphs <- function(s, assay, k = 10L) {
    SO <- asNamespace("SeuratObject")
    cells <- colnames(s)
    n_cells <- length(cells)
    norm <- SO$LayerData(s, layer = "data", assay = assay)
    pcs <- stats::prcomp(t(as.matrix(norm)), rank. = 10L)
    emb <- pcs$x
    dimnames(emb) <- list(cells, paste0("PC_", 1:10))
    load <- pcs$rotation
    colnames(load) <- paste0("PC_", 1:10)
    s[["pca"]] <- SO$CreateDimReducObject(embeddings = emb, loadings = load,
        stdev = pcs$sdev[1:10], key = "PC_", assay = assay)
    um <- matrix(stats::rnorm(n_cells * 2L), n_cells,
        dimnames = list(cells, c("umap_1", "umap_2")))
    s[["umap"]] <- SO$CreateDimReducObject(embeddings = um, key = "umap_", assay = assay)

    d <- as.matrix(stats::dist(emb[, 1:5]))
    nn <- Matrix::sparseMatrix(
        i = rep(seq_len(n_cells), each = k),
        j = as.vector(apply(d, 1, function(r) order(r)[seq_len(k)])),
        x = 1, dims = c(n_cells, n_cells), dimnames = list(cells, cells)
    )
    snn <- nn %*% Matrix::t(nn)
    snn@x <- snn@x / k
    graphs <- paste0(assay, c("_nn", "_snn"))
    s[[graphs[1]]] <- SO$as.Graph(nn)
    s[[graphs[2]]] <- SO$as.Graph(snn)
    cmd <- sprintf("FindNeighbors.%s.pca", assay)
    s@commands[[cmd]] <- methods::new("SeuratCommand", name = cmd,
        time.stamp = as.POSIXct("2026-01-01", tz = "UTC"), assay.used = assay,
        call.string = "FindNeighbors(object, reduction = \"pca\")",
        params = list(reduction = "pca", dims = 1:10, assay = assay, k.param = k,
            compute.SNN = TRUE, graph.name = graphs))
    s
}
