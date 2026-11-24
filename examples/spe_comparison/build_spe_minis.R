# Build SpatialExperiment / SpatialFeatureExperiment objects from the GiottoData
# vendor-layout crops with the Bioconductor readers, then run the usual steps
# (logNormCounts, PCA, UMAP; a kNN spatial graph for the SFE).
# Output: <cache>/spe_<name>.rds
# SFE's readXenium() needs RBioFormats (Java); the Xenium SFE is not built.
suppressMessages({
    library(SpatialExperimentIO)
    library(SpatialFeatureExperiment)
    library(VisiumIO)
})
cache <- Sys.getenv("GIOTTOCONNECT_CACHE",
    path.expand("~/Documents/Datasets/giottoconnect_cache"))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
crop <- function(nm) suppressMessages(GiottoData::getSpatialDataset(nm))

standard_steps <- function(x, n_pcs = 20L) {
    set.seed(1)
    # cells without counts have no size factor
    x <- x[, Matrix::colSums(SummarizedExperiment::assay(x, "counts")) > 0]
    x <- scuttle::logNormCounts(x)
    x <- scater::runPCA(x, ncomponents = n_pcs)
    scater::runUMAP(x, dimred = "PCA")
}

build <- list(
    spe_xenium = function() readXeniumSXE(crop("xenium_mini_lung_ffpe")),
    spe_cosmx = function() readCosmxSXE(crop("cosmx_mini_colon_wtx")),
    spe_visiumhd = function() import(TENxVisiumHD(spacerangerOut = crop("visiumhd_mini_crc"),
        bin_size = "016", images = "lowres")),
    sfe_cosmx = function() {
        x <- suppressMessages(readCosMX(crop("cosmx_mini_colon_wtx"), add_molecules = TRUE))
        colGraph(x, "knn5") <- findSpatialNeighbors(x, method = "knearneigh", k = 5)
        x
    }
)

for (nm in names(build)) {
    out <- file.path(cache, sprintf("%s.rds", nm))
    t0 <- Sys.time()
    x <- standard_steps(build[[nm]]())
    saveRDS(x, out)
    message(sprintf("%-13s %5d cells x %5d features, %.1f s -> %s", nm, ncol(x), nrow(x),
        as.numeric(difftime(Sys.time(), t0, units = "secs")), out))
}
