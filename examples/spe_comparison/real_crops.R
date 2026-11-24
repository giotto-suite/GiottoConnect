# spatialExperimentToGiotto() on SpatialExperiment / SFE objects read from the
# GiottoData crops with the Bioconductor readers (build_spe_minis.R): mapping
# choice, bound params, parts, readback and time, in memory and onto GiottoDisk.
#   Rscript examples/spe_comparison/real_crops.R
# Run from the GiottoConnect package root. Writes results/real_crops.csv.
options(giotto.has_conda = FALSE, giotto.use_conda = FALSE, giotto.no_python_warn = TRUE,
    lifecycle_verbosity = "quiet")
suppressMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(SpatialFeatureExperiment)
})
cache <- Sys.getenv("GIOTTOCONNECT_CACHE",
    path.expand("~/Documents/Datasets/giottoconnect_cache"))
objs <- c(spe_xenium = "readXeniumSXE", spe_cosmx = "readCosmxSXE",
    spe_visiumhd = "VisiumIO TENxVisiumHD (16 um)", sfe_cosmx = "SFE readCosMX + kNN graph")
secs <- function(t0) round(as.numeric(Sys.time() - t0, units = "secs"), 2)
rows <- list()
for (n in names(objs)) {
    x <- readRDS(file.path(cache, paste0(n, ".rds")))
    p <- adapterProfile(.adapter_for("spatialexperiment"), x)
    fm <- findMapping(p, ecosystem = "spatialexperiment")
    t0 <- Sys.time()
    g <- suppressWarnings(spatialExperimentToGiotto(x))
    t_mem <- secs(t0)
    rb <- readbackCheck(x, g)
    prov <- mappingProvenance(g)
    t0 <- Sys.time()
    gd <- suppressWarnings(spatialExperimentToGiotto(x, source = tempfile("spe_disk_")))
    t_disk <- secs(t0)
    v <- validateMapping(fm$mapping, p, result = g)
    rows[[n]] <- data.frame(
        object = n, reader = objs[[n]], class = class(x)[1], cells = ncol(x), features = nrow(x),
        mapping = fm$mapping$mapping$id, status = fm$status,
        feat_types = paste(unique(GiottoClass::list_expression(g)$feat_type), collapse = ", "),
        parts = paste(prov$parts, collapse = ", "),
        readback_match = sum(rb$status == "MATCH"), readback_declared = sum(rb$status == "DECLARED"),
        readback_undeclared = sum(rb$status == "UNDECLARED"),
        readback_disk_pass = isTRUE(attr(readbackCheck(x, gd), "pass")),
        warnings = paste(unique(v$check[v$level == "WARN"]), collapse = ", "),
        memory_s = t_mem, disk_s = t_disk, stringsAsFactors = FALSE)
    print(rows[[n]][, c("object", "mapping", "status", "feat_types", "readback_match",
        "readback_undeclared", "readback_disk_pass", "warnings")])
}
utils::write.csv(do.call(rbind, rows), "examples/spe_comparison/results/real_crops.csv", row.names = FALSE)
