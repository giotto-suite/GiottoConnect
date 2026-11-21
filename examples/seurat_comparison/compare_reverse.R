# The legacy giottoToSeuratV5 has been removed; run this at commit 83c26c8
# (results: results/phase7_reverse.csv).
# Giotto -> Seurat: mapping (giottoToSeurat) vs legacy giottoToSeuratV5, on
# giotto objects made from the fixture and the crops (in memory and on disk)
# and on the GiottoData Visium mini. Each result is checked against the
# original Seurat object where there is one (round trip).
# Run from the GiottoConnect package root.
options(giotto.has_conda = FALSE, giotto.use_conda = FALSE, giotto.no_python_warn = TRUE, width = 170,
    lifecycle_verbosity = "quiet", giottodisk.dgc_max_rows = Inf, giottodisk.dgc_max_cols = Inf)
suppressMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(SeuratObject)
})
cache <- Sys.getenv("GIOTTOCONNECT_CACHE", path.expand("~/Documents/Datasets/giottoconnect_cache"))
env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
sources <- list(
    fixture = function() env$build_seurat_fixture(),
    xenium = function() readRDS(file.path(cache, "seurat_xenium.rds")),
    cosmx = function() readRDS(file.path(cache, "seurat_cosmx.rds"))
)
legacy_v5 <- get("giottoToSeuratV5", asNamespace("GiottoConnect"))

check <- function(s0, s) {
    if (inherits(s, "error")) return(list(status = "error", detail = substr(conditionMessage(s), 1, 70)))
    a0 <- DefaultAssay(s0)
    out <- character()
    m0 <- LayerData(s0, layer = "counts", assay = a0)
    m <- LayerData(s, layer = "counts")
    cells <- intersect(colnames(m0), colnames(m))
    feats <- intersect(rownames(m0), rownames(m))
    out["counts"] <- if (length(cells) == ncol(m0) && length(feats) == nrow(m0) &&
        isTRUE(all.equal(as.numeric(m[feats, cells]), as.numeric(m0[feats, cells])))) "ok" else "WRONG"
    fov <- Filter(function(i) inherits(s[[i]], "FOV"), Images(s))
    ct0 <- GetTissueCoordinates(s0[[Images(s0)[1]]], which = "centroids")
    if (!length(fov)) out["locs"] <- "missing" else {
        b <- if ("centroids" %in% Boundaries(s[[fov[1]]])) "centroids" else Boundaries(s[[fov[1]]])[1]
        ct <- GetTissueCoordinates(s[[fov[1]]][[b]])
        mm <- merge(ct0, ct, by = "cell")
        out["locs"] <- if (nrow(mm) != nrow(ct0)) "WRONG ids" else
            if (isTRUE(all.equal(mm$x.x, mm$x.y)) && isTRUE(all.equal(mm$y.x, mm$y.y))) "ok" else "WRONG (axes)"
        out["segmentation"] <- if (any(vapply(Boundaries(s[[fov[1]]]), function(bb)
            inherits(s[[fov[1]]][[bb]], "Segmentation"), logical(1)))) "present" else "missing"
    }
    out["reductions"] <- paste(Reductions(s), collapse = "+")
    out["graphs"] <- paste(Graphs(s), collapse = "+")
    list(status = "ran", detail = paste(names(out), out, sep = ": ", collapse = "; "))
}

rows <- list()
for (ds in names(sources)) {
    s0 <- sources[[ds]]()
    for (backend in c("memory", "disk")) {
        g <- suppressWarnings(seuratToGiotto(s0, source = if (backend == "disk") tempfile() else NULL))
        for (impl in c("mapping", "legacy")) {
            t0 <- Sys.time()
            s <- tryCatch(suppressWarnings(suppressMessages(
                if (impl == "mapping") giottoToSeurat(g) else legacy_v5(g))), error = function(e) e)
            secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
            ck <- check(s0, s)
            rows[[length(rows) + 1L]] <- data.frame(dataset = ds, backend = backend, impl = impl,
                seconds = secs, status = ck$status, detail = ck$detail)
            cat(sprintf("%-8s %-6s %-8s %6.2fs %s\n", ds, backend, impl, secs, ck$detail))
        }
    }
}
# a Giotto object that did not come from Seurat
g <- suppressWarnings(GiottoData::loadGiottoMini("visium"))
for (impl in c("mapping", "legacy")) {
    s <- tryCatch(suppressWarnings(suppressMessages(
        if (impl == "mapping") giottoToSeurat(g) else legacy_v5(g))), error = function(e) e)
    d <- if (inherits(s, "error")) paste("error:", substr(conditionMessage(s), 1, 70)) else
        sprintf("%d cells x %d feats; images: %s; reductions: %s; graphs: %s", ncol(s), nrow(s),
            paste(Images(s), collapse = "+"), paste(Reductions(s), collapse = "+"), paste(Graphs(s), collapse = "+"))
    rows[[length(rows) + 1L]] <- data.frame(dataset = "visium_mini (GiottoData)", backend = "memory",
        impl = impl, seconds = NA, status = if (inherits(s, "error")) "error" else "ran", detail = d)
    cat(sprintf("visium_mini %-8s %s\n", impl, d))
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path("examples", "seurat_comparison", "results", "phase7_reverse.csv"), row.names = FALSE)
