# Small deterministic Seurat v5 test objects, built in code with SeuratObject.
# Sourced by the tests (tests/testthat/helper-contract-fixtures.R) and by the
# comparison scripts in examples/seurat_comparison/.

build_seurat_fixture <- function(n_cells = 200L, n_feat = 50L, assay = "Xenium",
                                 seed = 1L) {
    stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
    SO <- asNamespace("SeuratObject")
    set.seed(seed)
    cells <- sprintf("cell_%03d", seq_len(n_cells))
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

    # FOV: centroids + one square segmentation per cell, microns in 0-2000
    xy <- data.frame(x = stats::runif(n_cells, 0, 2000),
        y = stats::runif(n_cells, 0, 2000), cell = cells)
    sq <- do.call(rbind, lapply(seq_len(n_cells), function(i) {
        data.frame(x = xy$x[i] + c(-5, 5, 5, -5), y = xy$y[i] + c(-5, -5, 5, 5),
            cell = cells[i])
    }))
    fov <- SO$CreateFOV(
        list(centroids = SO$CreateCentroids(xy), segmentation = SO$CreateSegmentation(sq)),
        type = c("centroids", "segmentation"), assay = assay
    )
    s[["fov"]] <- fov

    # reductions: pca (10 dims) and umap (2 dims)
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

    # graphs: binary kNN (k = 10, self included) and weighted SNN
    k <- 10L
    d <- as.matrix(stats::dist(emb[, 1:5]))
    nn <- Matrix::sparseMatrix(
        i = rep(seq_len(n_cells), each = k),
        j = as.vector(apply(d, 1, function(r) order(r)[seq_len(k)])),
        x = 1, dims = c(n_cells, n_cells), dimnames = list(cells, cells)
    )
    snn <- nn %*% Matrix::t(nn)
    snn@x <- snn@x / k
    s[[paste0(assay, "_nn")]] <- SO$as.Graph(nn)
    s[[paste0(assay, "_snn")]] <- SO$as.Graph(snn)
    s
}
