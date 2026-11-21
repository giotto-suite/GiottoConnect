# Build Seurat v5 objects from the GiottoData vendor-layout crops with
# Seurat's own loaders, then run the standard Seurat steps so every part
# (normalized data, pca, umap, nn/snn graphs) is present.
# Output: <cache>/seurat_<name>.rds
suppressMessages({
    library(Seurat)
    library(SeuratObject)
})
options(arrow.unsafe_metadata = TRUE)

crops <- path.expand("~/Library/Caches/org.R-project.R/R/GiottoData")
cache <- Sys.getenv("GIOTTOCONNECT_CACHE",
    path.expand("~/Documents/Datasets/giottoconnect_cache"))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

standard_steps <- function(s, n_pcs = 20L) {
    s <- NormalizeData(s, verbose = FALSE)
    s <- FindVariableFeatures(s, nfeatures = min(2000L, nrow(s)), verbose = FALSE)
    s <- ScaleData(s, verbose = FALSE)
    s <- RunPCA(s, npcs = n_pcs, verbose = FALSE)
    s <- FindNeighbors(s, dims = seq_len(n_pcs), verbose = FALSE)
    s <- RunUMAP(s, dims = seq_len(n_pcs), verbose = FALSE)
    s
}

build <- list(
    xenium = function() LoadXenium(file.path(crops, "xenium_mini_lung_ffpe"),
        fov = "fov", segmentations = "cell"),
    cosmx = function() LoadNanostring(file.path(crops, "cosmx_mini_colon_wtx"),
        assay = "Nanostring", fov = "fov"),
    visiumhd = function() Load10X_Spatial(file.path(crops, "visiumhd_mini_crc"),
        bin.size = 16)
)

for (nm in names(build)) {
    out <- file.path(cache, sprintf("seurat_%s.rds", nm))
    t0 <- Sys.time()
    s <- build[[nm]]()
    s <- standard_steps(s)
    saveRDS(s, out)
    message(sprintf("%-9s %5d cells x %5d features, %.1f s -> %s", nm, ncol(s), nrow(s),
        as.numeric(difftime(Sys.time(), t0, units = "secs")), out))
}
