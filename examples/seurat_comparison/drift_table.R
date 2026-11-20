# Drift variants D1-D9: what findMapping / repairMapping (rules only) do with
# each, whether the result converts and reads back, and what the legacy
# GiottoClass::seuratToGiottoV5 does with the same object.
# Run from the GiottoConnect package root.
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE, width = 160)
suppressMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(SeuratObject)
})
env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
v <- suppressWarnings(env$build_drift_variants())
truth <- env$build_seurat_fixture()

# does a giotto object hold the fixture's raw counts and correctly oriented locs?
legacy_check <- function(g, s) {
    raw <- tryCatch(GiottoClass::getExpression(g, values = "raw", output = "matrix"),
        error = function(e) NULL)
    src <- LayerData(truth, layer = "counts")
    raw_ok <- !is.null(raw) && setequal(colnames(raw), colnames(src)) &&
        isTRUE(all.equal(as.numeric(raw[rownames(src), colnames(src)]), as.numeric(src)))
    sl <- tryCatch(GiottoClass::getSpatialLocations(g, output = "data.table"), error = function(e) NULL)
    ct <- GetTissueCoordinates(truth[["fov"]], which = "centroids")
    locs <- if (is.null(sl)) "none" else {
        m <- merge(sl, ct, by.x = "cell_ID", by.y = "cell")
        if (!nrow(m)) "wrong ids" else if (isTRUE(all.equal(m$sdimx, m$x)) && isTRUE(all.equal(m$sdimy, -m$y))) "x,-y"
        else if (isTRUE(all.equal(m$sdimx, m$y)) && isTRUE(all.equal(m$sdimy, -m$x))) "y,-x (swapped)"
        else "other"
    }
    sprintf("raw %s, locs %s", if (raw_ok) "ok" else "WRONG/missing", locs)
}

rows <- list()
for (d in names(v)) {
    s <- v[[d]]
    p <- profileObject(s)
    f <- findMapping(p, ecosystem = "seurat")
    rep_status <- "-"
    rules <- "-"
    warns <- "-"
    readback <- "-"
    if (f$status != "none") {
        r <- repairMapping(f$mapping, s)
        rep_status <- r$status
        rr <- r$changes$rule[r$changes$rule != "added"]
        rules <- if (length(rr)) paste(unique(rr), collapse = ", ") else "-"
        readback <- if (is.null(r$readback)) "-" else if (isTRUE(attr(r$readback, "pass"))) "PASS" else "FAIL"
        if (r$status != "failed") {
            w <- character()
            g <- withCallingHandlers(seuratToGiotto(s, mapping = r$mapping),
                warning = function(cnd) {
                    w <<- c(w, conditionMessage(cnd))
                    invokeRestart("muffleWarning")
                })
            msgs <- unlist(strsplit(w, "\n"))
            msgs <- trimws(sub("^\\s*- ", "", msgs[grepl("^\\s*- \\[", msgs)]))
            warns <- if (length(msgs)) paste(unique(sub("\\].*", "]", msgs)), collapse = " ") else "-"
        }
    }
    a <- if (Assays(s)[1] %in% Assays(s)) DefaultAssay(s) else Assays(s)[1]
    leg <- tryCatch(suppressWarnings(suppressMessages(GiottoClass::seuratToGiottoV5(s,
        spatial_assay = a, dim_reduction = intersect(c("pca", "umap"), Reductions(s)),
        nn_network = intersect(paste0(a, "_nn"), Graphs(s)), verbose = FALSE))),
        error = function(e) e)
    legacy <- if (inherits(leg, "error")) paste("error:", substr(conditionMessage(leg), 1, 60)) else
        if (!inherits(leg, "giotto")) sprintf("returned a %s", class(leg)[1]) else legacy_check(leg, s)
    rows[[d]] <- data.frame(variant = d,
        change = env$DRIFT_EXPECTED$change[env$DRIFT_EXPECTED$variant == d],
        expected = env$DRIFT_EXPECTED$expected[env$DRIFT_EXPECTED$variant == d],
        findMapping = f$status, repair = rep_status, rules = rules, readback = readback,
        conversion_warnings = warns, legacy = legacy, stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path("examples", "seurat_comparison", "results", "phase4_drift.csv"),
    row.names = FALSE)
print(tab, row.names = FALSE, right = FALSE)
