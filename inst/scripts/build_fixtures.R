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

# giotto object made from the Seurat fixture (in memory, or on GiottoDisk
# when `source` is a directory)
build_giotto_fixture <- function(source = NULL) {
    GiottoConnect::seuratToGiotto(build_seurat_fixture(), source = source)
}

# Seurat v3/v4 layout: a v3 Assay (counts/data/scale.data slots) with the
# same FOV, reductions, graphs and FindNeighbors record as the v5 fixture
build_seurat_v3_fixture <- function(n_cells = 200L, n_feat = 50L, assay = "Xenium",
                                    seed = 1L) {
    stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
    SO <- asNamespace("SeuratObject")
    set.seed(seed)
    s <- .fx_expression(n_cells, n_feat, assay, assay_class = "Assay")
    cells <- colnames(s)
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

# VisiumV1 spots on a hexagonal grid (100 um pitch), Seurat v4-style v3 Assay
build_seurat_visium_v1_fixture <- function(n_side = 12L, n_feat = 50L, um_per_px = 0.885,
                                           assay_class = "Assay", seed = 1L) {
    stopifnot(requireNamespace("Seurat", quietly = TRUE))
    set.seed(seed)
    n <- n_side^2
    s <- .fx_expression(n, n_feat, "Spatial", cells = sprintf("SPOT%04d-1", seq_len(n)),
        assay_class = assay_class)
    pitch_px <- 100 / um_per_px
    row <- rep(seq_len(n_side) - 1L, n_side)
    col <- rep(seq_len(n_side) - 1L, each = n_side)
    coords <- data.frame(tissue = 1L, row = row, col = col,
        imagerow = 500 + row * pitch_px * sqrt(3) / 2,
        imagecol = 500 + col * pitch_px + (row %% 2) * pitch_px / 2,
        row.names = colnames(s))
    sf <- Seurat::scalefactors(spot = 0.55 * pitch_px, fiducial = 0.9 * pitch_px, hires = 0.2, lowres = 0.06)
    img <- methods::new("VisiumV1", image = array(0.5, c(8L, 8L, 3L)), scale.factors = sf,
        coordinates = coords, spot.radius = 0.01, assay = "Spatial", key = "slice1_")
    s[["slice1"]] <- img
    .fx_reductions_graphs(s, "Spatial")
}

# Slide-seq beads read with Seurat's own ReadSlideSeq()
build_seurat_slideseq_fixture <- function(n_cells = 200L, n_feat = 50L, seed = 1L) {
    stopifnot(requireNamespace("Seurat", quietly = TRUE))
    set.seed(seed)
    s <- .fx_expression(n_cells, n_feat, "Spatial", cells = sprintf("BEAD%04d", seq_len(n_cells)))
    f <- tempfile(fileext = ".csv")
    writeLines(c("barcodes,xcoord,ycoord", sprintf("%s,%.3f,%.3f", colnames(s),
        stats::runif(n_cells, 0, 3000), stats::runif(n_cells, 0, 3000))), f)
    img <- Seurat::ReadSlideSeq(f, assay = "Spatial")
    s[["puck"]] <- img[colnames(s)]
    .fx_reductions_graphs(s, "Spatial")
}

# counts / data / scale.data layers and cell metadata
.fx_expression <- function(n_cells, n_feat, assay,
                           cells = sprintf("cell_%03d", seq_len(n_cells)),
                           assay_class = c("Assay5", "Assay")) {
    assay_class <- match.arg(assay_class)
    SO <- asNamespace("SeuratObject")
    feats <- sprintf("gene%02d", seq_len(n_feat))
    counts <- matrix(stats::rpois(n_cells * n_feat, 2), n_feat, n_cells,
        dimnames = list(feats, cells))
    counts <- methods::as(counts, "CsparseMatrix")
    norm <- log1p(counts)
    scaled <- t(scale(t(as.matrix(norm))))
    scaled[is.na(scaled)] <- 0

    if (assay_class == "Assay5") {
        s <- SO$CreateSeuratObject(SO$CreateAssay5Object(counts = counts), assay = assay)
        s <- SO$`LayerData<-`(s, layer = "data", value = norm)
        s <- SO$`LayerData<-`(s, layer = "scale.data", value = scaled)
    } else {
        s <- SO$CreateSeuratObject(SO$CreateAssayObject(counts = counts), assay = assay)
        s <- SO$SetAssayData(s, layer = "data", new.data = norm, assay = assay)
        s <- SO$SetAssayData(s, layer = "scale.data", new.data = scaled, assay = assay)
    }
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

# Structural drift variants of the FOV fixture, each with the outcome the
# mapping machinery should reach (see docs/contract-prototype).
build_drift_variants <- function() {
    stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
    SO <- asNamespace("SeuratObject")
    base <- build_seurat_fixture()
    v <- list()

    # D1 assay renamed
    v$D1 <- build_seurat_fixture(assay = "Spatial")

    # D2 layer counts -> raw_counts
    s <- base
    a <- s[["Xenium"]]
    a$raw_counts <- SO$LayerData(a, layer = "counts")
    a$counts <- NULL
    s[["Xenium"]] <- a
    v$D2 <- s

    # D3 coordinates moved into meta.data, FOV removed
    s <- base
    ct <- SO$GetTissueCoordinates(s[["fov"]], which = "centroids")
    s$x_um <- ct$x[match(colnames(s), ct$cell)]
    s$y_um <- ct$y[match(colnames(s), ct$cell)]
    s[["fov"]] <- NULL
    v$D3 <- s

    # D4 counts stored as a dense matrix
    s <- base
    a <- s[["Xenium"]]
    a$counts <- as.matrix(SO$LayerData(a, layer = "counts"))
    s[["Xenium"]] <- a
    v$D4 <- s

    # D5 extra reduction
    s <- base
    s[["tsne"]] <- SO$CreateDimReducObject(
        embeddings = `colnames<-`(SO$Embeddings(s, "umap") * 2, c("tSNE_1", "tSNE_2")),
        key = "tSNE_", assay = "Xenium")
    v$D5 <- s

    # D6 coordinates in mm
    s <- base
    xy <- SO$GetTissueCoordinates(s[["fov"]], which = "centroids")
    seg <- SO$GetTissueCoordinates(s[["fov"]][["segmentation"]])
    xy[, 1:2] <- xy[, 1:2] / 1000
    seg[, 1:2] <- seg[, 1:2] / 1000
    s[["fov"]] <- NULL
    s[["fov"]] <- SO$CreateFOV(
        list(centroids = SO$CreateCentroids(xy), segmentation = SO$CreateSegmentation(seg)),
        type = c("centroids", "segmentation"), assay = "Xenium")
    v$D6 <- s

    # D7 graph renamed
    s <- base
    s[["knn_graph"]] <- s[["Xenium_nn"]]
    s[["Xenium_nn"]] <- NULL
    v$D7 <- s

    # D8 slot the mapping does not use removed
    s <- base
    s[["Xenium_snn"]] <- NULL
    v$D8 <- s

    # D9 v3 Assay instead of Assay5
    s <- SO$CreateSeuratObject(SO$CreateAssayObject(counts = SO$LayerData(base, layer = "counts")),
        assay = "Xenium")
    s[["fov"]] <- base[["fov"]]
    v$D9 <- s

    v
}

DRIFT_EXPECTED <- data.frame(
    variant = paste0("D", 1:9),
    change = c("assay renamed Xenium -> Spatial", "layer counts -> raw_counts",
        "coordinates moved to meta.data x_um/y_um, FOV removed", "counts stored as a dense matrix",
        "extra reduction tsne", "coordinates in mm", "graph Xenium_nn -> knn_graph",
        "unused graph Xenium_snn removed", "v3 Assay (counts only) instead of Assay5"),
    expected = c("params", "params", "rules (alternative)", "rules (cast)", "params, converted by the for_each template",
        "params, WARN looks like mm", "params", "params (fingerprint unchanged)",
        "seurat-assay mapping, optional parts absent"),
    stringsAsFactors = FALSE
)

# Drift that rules alone cannot resolve: input for the LLM repair benchmark.
# truth: the option a correct repair picks (NA: no correct option exists)
build_llm_cases <- function() {
    stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
    SO <- asNamespace("SeuratObject")
    base <- build_seurat_fixture()
    counts <- SO$LayerData(base, layer = "counts")
    drop_counts <- function(s) {
        a <- s[["Xenium"]]
        a$counts <- NULL
        s[["Xenium"]] <- suppressWarnings(a)
        s
    }
    ct <- SO$GetTissueCoordinates(base[["fov"]], which = "centroids")
    ct <- ct[match(colnames(base), ct$cell), ]
    cases <- list()

    # L1 counts stashed under misc
    s <- drop_counts(base)
    s@misc$stash <- counts
    cases$L1 <- list(object = s, truth = c(expr_raw.from = "misc/stash"))

    # L2 two plausible integer matrices, the raw one has the less obvious name
    s <- drop_counts(base)
    corrected <- counts
    corrected@x <- pmax(0, corrected@x - 1)
    s@misc$counts_corrected <- Matrix::drop0(corrected)
    s@misc$umi_raw <- counts
    cases$L2 <- list(object = s, truth = c(expr_raw.from = "misc/umi_raw"))

    # L3 CosMx-style coordinate columns, FOV removed
    s <- base
    s$CenterX_local_px <- ct$x %% 100
    s$CenterY_local_px <- ct$y %% 100
    s$CenterX_global_px <- ct$x
    s$CenterY_global_px <- ct$y
    s$Area <- stats::runif(ncol(s), 50, 150)
    s[["fov"]] <- NULL
    cases$L3 <- list(object = s, truth = c(
        locs.locs_from_meta.x_col = "CenterX_global_px",
        locs.locs_from_meta.y_col = "CenterY_global_px"))

    # L4 coordinates only as a "spatial" reduction, FOV removed (rules)
    s <- base
    emb <- as.matrix(ct[, c("x", "y")])
    dimnames(emb) <- list(ct$cell, c("SPATIAL_1", "SPATIAL_2"))
    s[["spatial"]] <- SO$CreateDimReducObject(embeddings = emb, key = "SPATIAL_", assay = "Xenium")
    s[["fov"]] <- NULL
    cases$L4 <- list(object = s, truth = c(locs = "rules: locs_from_reduction"))

    # L5 counts gone and nothing replaces them
    cases$L5 <- list(object = drop_counts(base), truth = c(expr_raw.from = NA_character_))
    cases
}


# AnnData .h5ad ####
# Written with hdf5r following the anndata on-disk spec (no python).
#   build_h5ad_fixture(path)                  modern scanpy: dense normalized X,
#       layers/counts (csr), categorical obs, obsm X_pca/X_umap/spatial,
#       obsp distances + connectivities with a uns/neighbors record
#   build_h5ad_fixture(path, legacy = TRUE)   anndata 0.7: counts in X (csc),
#       raw/X copy, pre-0.8 categorical (`__categories`), obsm spatial only

build_h5ad_fixture <- function(path, n_cells = 200L, n_feat = 60L, k = 5L,
                               legacy = FALSE, seed = 1L) {
    stopifnot(requireNamespace("hdf5r", quietly = TRUE))
    set.seed(seed)
    cells <- sprintf("cell_%03d", seq_len(n_cells))
    feats <- sprintf("gene%02d", seq_len(n_feat))
    feats[n_feat] <- feats[1L]                     # one duplicated gene name
    counts <- Matrix::rsparsematrix(n_cells, n_feat, density = 0.2,
        rand.x = function(n) as.double(stats::rpois(n, 4) + 1))
    xy <- cbind(stats::runif(n_cells, 0, 2000), stats::runif(n_cells, 0, 2000))
    pca <- matrix(stats::rnorm(n_cells * 10L), n_cells, 10L)
    type <- factor(sample(c("T cell", "B cell", "tumor"), n_cells, TRUE))
    type[3L] <- NA

    if (file.exists(path)) unlink(path)
    h5 <- hdf5r::H5File$new(path, mode = "w")
    on.exit(h5$close_all())
    enc <- function(node, type, version = "0.1.0") {
        hdf5r::h5attr(node, "encoding-type") <- type
        hdf5r::h5attr(node, "encoding-version") <- version
    }
    # obs x var matrices; hdf5r writes dims reversed, so write t()
    dense <- function(where, name, m) {
        where[[name]] <- t(m)
        if (!legacy) enc(where[[name]], "array", "0.2.0")
    }
    sparse <- function(where, name, m, csr = TRUE) {
        g <- where$create_group(name)
        cm <- methods::as(if (csr) Matrix::t(m) else m, "CsparseMatrix")
        g[["data"]] <- if (all(cm@x == round(cm@x))) as.integer(cm@x) else cm@x
        g[["indices"]] <- cm@i
        g[["indptr"]] <- cm@p
        enc(g, if (csr) "csr_matrix" else "csc_matrix")
        hdf5r::h5attr(g, "shape") <- dim(m)
    }
    frame <- function(name, ids, cols) {
        g <- h5$create_group(name)
        g[["_index"]] <- ids
        for (cn in names(cols)) {
            v <- cols[[cn]]
            if (is.factor(v) && legacy) {
                codes <- as.integer(v) - 1L
                codes[is.na(codes)] <- -1L
                g[[cn]] <- codes
                if (!g$exists("__categories")) g$create_group("__categories")
                g[["__categories"]][[cn]] <- levels(v)
                hdf5r::h5attr(g[[cn]], "categories") <- "__categories"
            } else if (is.factor(v)) {
                cg <- g$create_group(cn)
                codes <- as.integer(v) - 1L
                codes[is.na(codes)] <- -1L
                cg[["codes"]] <- codes
                cg[["categories"]] <- levels(v)
                enc(cg, "categorical", "0.2.0")
                hdf5r::h5attr(cg, "ordered") <- FALSE
            } else {
                g[[cn]] <- v
                if (!legacy) enc(g[[cn]], if (is.character(v)) "string-array" else "array", "0.2.0")
            }
        }
        enc(g, "dataframe", "0.2.0")
        hdf5r::h5attr(g, "_index") <- "_index"
        hdf5r::h5attr(g, "column-order") <- if (length(cols)) names(cols) else character()
        g
    }
    norm <- log1p(as.matrix(counts) / Matrix::rowSums(counts) * 1e4)
    if (legacy) {
        sparse(h5, "X", counts, csr = FALSE)
        raw <- h5$create_group("raw")
        sparse(raw, "X", counts, csr = FALSE)
        frame("obs", cells, list(n_counts = Matrix::rowSums(counts), cell_type = type))
        frame("var", feats, list())
        h5$link_move_from(h5, "var", "raw/var")
        frame("var", feats, list())
        obsm <- h5$create_group("obsm")
        dense(obsm, "spatial", xy)
        h5$create_group("uns")
    } else {
        enc(h5, "anndata")
        dense(h5, "X", norm)
        layers <- h5$create_group("layers")
        enc(layers, "dict")
        sparse(layers, "counts", counts)
        frame("obs", cells, list(n_counts = Matrix::rowSums(counts), cell_type = type,
            center_x = xy[, 1], center_y = xy[, 2]))
        frame("var", feats, list(highly_variable = rep(c(1L, 0L), length.out = n_feat)))
        obsm <- h5$create_group("obsm")
        enc(obsm, "dict")
        dense(obsm, "X_pca", pca)
        dense(obsm, "X_umap", pca[, 1:2])
        dense(obsm, "spatial", xy)
        # kNN in PCA space (30 dims would be scanpy's; 10 here), self excluded
        d <- as.matrix(stats::dist(pca))
        nn <- t(apply(d, 1L, function(r) order(r)[2:(k + 1L)]))
        dist <- Matrix::sparseMatrix(i = rep(seq_len(n_cells), each = k), j = as.vector(t(nn)),
            x = d[cbind(rep(seq_len(n_cells), each = k), as.vector(t(nn)))], dims = c(n_cells, n_cells))
        obsp <- h5$create_group("obsp")
        enc(obsp, "dict")
        sparse(obsp, "distances", dist)
        sparse(obsp, "connectivities", methods::as(dist > 0, "dMatrix") * 0.5)
        uns <- h5$create_group("uns")
        enc(uns, "dict")
        nb <- uns$create_group("neighbors")
        enc(nb, "dict")
        nb[["connectivities_key"]] <- "connectivities"
        nb[["distances_key"]] <- "distances"
        enc(nb[["connectivities_key"]], "string", "0.2.0")
        enc(nb[["distances_key"]], "string", "0.2.0")
        pr <- nb$create_group("params")
        enc(pr, "dict")
        for (nm in c("n_neighbors", "n_pcs")) {
            pr[[nm]] <- if (nm == "n_neighbors") k + 1L else 10L
            enc(pr[[nm]], "numeric-scalar", "0.2.0")
        }
        pr[["method"]] <- "umap"
        pr[["metric"]] <- "euclidean"
        enc(pr[["method"]], "string", "0.2.0")
        enc(pr[["metric"]], "string", "0.2.0")
    }
    invisible(list(path = path, counts = counts, norm = norm, xy = xy, pca = pca,
        type = type, cells = cells, feats = feats))
}
