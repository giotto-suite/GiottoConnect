# Readback check of seuratToGiotto() (memory and GiottoDisk) on the fixture
# and the GiottoData crops. Run from the GiottoConnect package root.
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE, width = 150)
suppressMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(SeuratObject)
})
cache <- Sys.getenv("GIOTTOCONNECT_CACHE",
    path.expand("~/Documents/Datasets/giottoconnect_cache"))
env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
datasets <- list(
    fixture = function() env$build_seurat_fixture(),
    fixture_visium = function() env$build_seurat_visium_fixture(),
    xenium = function() readRDS(file.path(cache, "seurat_xenium.rds")),
    cosmx = function() readRDS(file.path(cache, "seurat_cosmx.rds")),
    visiumhd = function() readRDS(file.path(cache, "seurat_visiumhd.rds"))
)
all <- list()
for (ds in names(datasets)) {
    s <- datasets[[ds]]()
    for (backend in c("memory", "disk")) {
        src <- if (backend == "disk") tempfile(paste0(ds, "_")) else NULL
        g <- suppressWarnings(seuratToGiotto(s, source = src))
        rb <- readbackCheck(s, g)
        cat(sprintf("\n## %s / %s (%s)\n", ds, backend, mappingProvenance(g)$mapping))
        print(rb)
        all[[length(all) + 1L]] <- cbind(dataset = ds, backend = backend,
            mapping = mappingProvenance(g)$mapping, as.data.frame(unclass(rb)))
    }
}
all <- do.call(rbind, all)
dir.create(file.path("examples", "seurat_comparison", "results"), showWarnings = FALSE)
write.csv(all, file.path("examples", "seurat_comparison", "results", "phase3_readback.csv"),
    row.names = FALSE)
cat("\nSummary:\n")
print(table(paste(all$dataset, all$backend), all$status))
