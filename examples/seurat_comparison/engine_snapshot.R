# Snapshot of engine behaviour, to prove a refactor changes nothing.
#   Rscript engine_snapshot.R <out.rds>      record
#   Rscript engine_snapshot.R <out.rds> <baseline.rds>   record and compare
# Per object: mapping choice, status, fingerprint, bound params, a content
# hash of every converted part (memory and GiottoDisk), readback outcome.
# Run from the GiottoConnect package root.
args <- commandArgs(trailingOnly = TRUE)
options(giotto.has_conda = FALSE, giotto.use_conda = FALSE, giotto.no_python_warn = TRUE,
    arrow.unsafe_metadata = TRUE, giottodisk.dgc_max_rows = Inf, giottodisk.dgc_max_cols = Inf,
    lifecycle_verbosity = "quiet")
suppressMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(SeuratObject)
})
cache <- path.expand("~/Documents/Datasets/giottoconnect_cache")
env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
dv <- suppressWarnings(env$build_drift_variants())
seurat_objs <- list(
    fixture = env$build_seurat_fixture(), v3 = env$build_seurat_v3_fixture(),
    visium_fx = env$build_seurat_visium_fixture(), visium_v1_fx = env$build_seurat_visium_v1_fixture(),
    slideseq_fx = env$build_seurat_slideseq_fixture(),
    D3 = dv$D3, D4 = dv$D4, D5 = dv$D5, D9 = dv$D9,
    xenium = readRDS(file.path(cache, "seurat_xenium.rds")),
    cosmx = readRDS(file.path(cache, "seurat_cosmx.rds")),
    visiumhd = readRDS(file.path(cache, "seurat_visiumhd.rds")),
    atera = readRDS(file.path(cache, "seurat6/seurat_atera_tiny.rds"))
)
h <- function(x) digest::digest(x, algo = "xxhash64")
mat <- function(m) {
    m <- if (methods::is(m, "Matrix") || is.matrix(m)) m else methods::as(m, "dgCMatrix")
    m <- methods::as(m, "CsparseMatrix")
    m[order(rownames(m)), order(colnames(m)), drop = FALSE]
}

giotto_content <- function(g) {
    out <- list()
    ex <- GiottoClass::list_expression(g)
    for (i in seq_len(nrow(ex))) {
        e <- GiottoClass::getExpression(g, values = ex$name[i], spat_unit = ex$spat_unit[i],
            feat_type = ex$feat_type[i], output = "exprObj")[]
        # rna keeps the key it had before other feature types were mapped
        key <- if (ex$feat_type[i] == "rna") paste("expr", ex$name[i]) else
            paste("expr", ex$feat_type[i], ex$name[i])
        out[[key]] <- h(as.matrix(mat(e)))
    }
    md <- as.data.frame(GiottoClass::pDataDT(g))
    out$cell_meta <- h(md[order(md$cell_ID), sort(names(md))])
    fm <- as.data.frame(GiottoClass::fDataDT(g))
    out$feat_meta <- h(fm[order(fm$feat_ID), sort(names(fm))])
    sl <- tryCatch(as.data.frame(GiottoClass::getSpatialLocations(g, output = "data.table")), error = function(e) NULL)
    if (!is.null(sl)) out$locs <- h(sl[order(sl$cell_ID), c("cell_ID", "sdimx", "sdimy")])
    pn <- GiottoClass::list_spatial_info(g)
    for (nm in pn$spat_info) {
        p <- GiottoClass::getPolygonInfo(g, name = nm, return_giottoPolygon = TRUE)
        sv <- if (inherits(p[], "SpatVector")) p[] else GiottoClass::as.terra(p[])
        geo <- as.data.frame(terra::geom(sv))
        geo$id <- sv$poly_ID[geo$geom]
        out[[paste("poly", nm)]] <- h(geo[order(geo$id, geo$part, seq_len(nrow(geo))), c("id", "x", "y")])
    }
    dr <- GiottoClass::list_dim_reductions(g)
    for (i in seq_len(NROW(dr))) {
        m <- GiottoClass::getDimReduction(g, reduction_method = dr$dim_type[i], name = dr$name[i], output = "matrix")
        out[[paste("dimred", dr$name[i])]] <- h(unname(m[order(rownames(m)), ]))
    }
    nn <- GiottoClass::list_nearest_networks(g)
    for (i in seq_len(NROW(nn))) {
        e <- as.data.frame(GiottoClass::getNearestNetwork(g, nn_type = nn$nn_type[i], name = nn$name[i],
            output = "data.table"))[, c("from", "to", "weight", "distance")]
        out[[paste("nn", nn$nn_type[i], nn$name[i])]] <- h(e[order(e$from, e$to), ])
    }
    out
}

seurat_content <- function(s) {
    out <- list()
    for (a in Assays(s)) for (l in Layers(s[[a]])) out[[paste("layer", a, l)]] <- h(as.matrix(mat(LayerData(s, layer = l, assay = a))))
    md <- s[[]]
    out$meta <- h(md[order(rownames(md)), sort(names(md))])
    for (im in Images(s)) for (b in Boundaries(s[[im]])) {
        ct <- GetTissueCoordinates(s[[im]][[b]])
        out[[paste("geo", im, b)]] <- h(ct[order(ct$cell, seq_len(nrow(ct))), c("cell", "x", "y")])
    }
    for (r in Reductions(s)) { e <- Embeddings(s, r); out[[paste("red", r)]] <- h(unname(e[order(rownames(e)), ])) }
    for (gname in Graphs(s)) out[[paste("graph", gname)]] <- h(as.matrix(mat(s[[gname]])))
    out
}

record <- function(x, to) {
    p <- profileObject(x)
    f <- findMapping(p)
    row <- list(mapping = f$mapping$mapping$id %||% "-", status = f$status,
        fingerprint = f$fingerprint, params = f$params$values)
    if (f$status == "none") return(row)
    res <- suppressWarnings(if (to == "giotto") seuratToGiotto(x, mapping = f$mapping) else giottoToSeurat(x))
    row$content <- if (to == "giotto") giotto_content(res) else seurat_content(res)
    rb <- suppressWarnings(readbackCheck(x, res))
    row$readback <- table(factor(rb$status, levels = c("MATCH", "DECLARED", "UNDECLARED")))
    row$losses <- mappingProvenance(res)$losses[, c("part_id", "kind", "field", "note", "source")]
    if (to == "giotto") {
        gd <- suppressWarnings(seuratToGiotto(x, mapping = f$mapping, source = tempfile("snap_")))
        row$content_disk <- giotto_content(gd)
    }
    row
}

snap <- list()
for (n in names(seurat_objs)) {
    snap[[paste0("s2g:", n)]] <- record(seurat_objs[[n]], "giotto")
    cat(n, snap[[paste0("s2g:", n)]]$mapping, snap[[paste0("s2g:", n)]]$status, "\n")
}
giotto_objs <- list(
    fixture = seuratToGiotto(seurat_objs$fixture),
    fixture_disk = seuratToGiotto(seurat_objs$fixture, source = tempfile("snapg_")),
    xenium = suppressWarnings(seuratToGiotto(seurat_objs$xenium)),
    visium_mini = suppressWarnings(GiottoData::loadGiottoMini("visium"))
)
for (n in names(giotto_objs)) {
    snap[[paste0("g2s:", n)]] <- record(giotto_objs[[n]], "seurat")
    cat(n, snap[[paste0("g2s:", n)]]$mapping, snap[[paste0("g2s:", n)]]$status, "\n")
}
saveRDS(snap, args[[1]])
cat("recorded", length(snap), "objects to", args[[1]], "\n")

if (length(args) > 1) {
    base <- readRDS(args[[2]])
    diffs <- character()
    for (k in union(names(base), names(snap))) {
        a <- base[[k]]
        b <- snap[[k]]
        for (field in union(names(a), names(b))) {
            if (!identical(a[[field]], b[[field]])) diffs <- c(diffs, sprintf("%s $%s", k, field))
        }
    }
    if (length(diffs)) {
        cat("DIFFERENCES:\n", paste0("  ", diffs, collapse = "\n"), "\n")
        quit(status = 1)
    }
    cat("IDENTICAL to baseline:", length(snap), "objects\n")
}
