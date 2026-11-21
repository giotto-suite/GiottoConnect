# Part-by-part comparison: mapping engine (memory + GiottoDisk) vs the legacy
# GiottoClass::seuratToGiottoV5, each checked against the Seurat source read
# through SeuratObject accessors. Xenium is also compared with GiottoDisk's
# native importer on the same crop (independent oracle).
# Run from the GiottoConnect package root.
options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE,
    arrow.unsafe_metadata = TRUE, lifecycle_verbosity = "quiet",
    giottodisk.dgc_max_rows = Inf, giottodisk.dgc_max_cols = Inf)
suppressMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(SeuratObject)
})
cache <- Sys.getenv("GIOTTOCONNECT_CACHE",
    path.expand("~/Documents/Datasets/giottoconnect_cache"))
out_dir <- file.path("examples", "seurat_comparison", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
datasets <- list(
    fixture = function() env$build_seurat_fixture(),
    xenium = function() readRDS(file.path(cache, "seurat_xenium.rds")),
    cosmx = function() readRDS(file.path(cache, "seurat_cosmx.rds")),
    visiumhd = function() readRDS(file.path(cache, "seurat_visiumhd.rds"))
)

as_mat <- function(x) if (methods::is(x, "Matrix") || is.matrix(x)) x else methods::as(x, "dgCMatrix")
row <- function(status, detail = "") list(status = status, detail = detail)

checks <- list(
    expr_raw = function(g, s, a) {
        src <- LayerData(s, "counts", assay = a)
        m <- as_mat(getExpression(g, values = "raw", spat_unit = "cell", feat_type = "rna")[])
        if (!setequal(colnames(m), colnames(src))) return(row("wrong", "cell ids differ"))
        if (!setequal(rownames(m), rownames(src))) return(row("wrong", "feature ids differ"))
        same <- isTRUE(all.equal(as.numeric(m[rownames(src), colnames(src)]), as.numeric(src)))
        row(if (same) "correct" else "wrong", sprintf("%dx%d", nrow(m), ncol(m)))
    },
    expr_normalized = function(g, s, a) {
        if (!"data" %in% Layers(s[[a]])) return(row("n/a"))
        src <- LayerData(s, "data", assay = a)
        m <- as_mat(getExpression(g, values = "normalized", spat_unit = "cell", feat_type = "rna")[])
        same <- isTRUE(all.equal(as.numeric(m[rownames(src), colnames(src)]), as.numeric(src),
            tolerance = 1e-6))
        row(if (same) "correct" else "wrong")
    },
    cell_meta = function(g, s, a) {
        src <- s[[]]
        dt <- pDataDT(g, spat_unit = "cell", feat_type = "rna")
        miss <- setdiff(names(src), names(dt))
        if (length(miss)) return(row("wrong", paste("missing:", paste(miss, collapse = ","))))
        dt <- dt[match(rownames(src), dt$cell_ID)]
        bad <- names(src)[!vapply(names(src), function(cl) {
            identical(as.character(src[[cl]]), as.character(dt[[cl]]))
        }, logical(1))]
        row(if (length(bad)) "wrong" else "correct", paste(bad, collapse = ","))
    },
    feat_meta = function(g, s, a) {
        dt <- fDataDT(g, spat_unit = "cell", feat_type = "rna")
        cell_cols <- intersect(setdiff(names(s[[]]), "orig.ident"), names(dt))
        if (length(cell_cols)) return(row("wrong", paste("holds cell metadata:", paste(cell_cols, collapse = ","))))
        row(if (setequal(dt$feat_ID, Features(s[[a]]))) "correct" else "wrong", sprintf("%d features", nrow(dt)))
    },
    locs = function(g, s, a) {
        fov <- Filter(function(i) inherits(s[[i]], "FOV"), Images(s))
        if (!length(fov)) return(row("n/a", "no FOV"))
        ct <- GetTissueCoordinates(s[[fov[1]]], which = "centroids")
        if (inherits(s[[fov[1]]], "VisiumV2")) {
            # Visium HD: full-resolution pixels -> microns (bin size / spot scale factor)
            bin <- as.numeric(sub(".*\\.0*(\\d+)um$", "\\1", fov[1]))
            ct$x <- ct$x * bin / s[[fov[1]]]@scale.factors$spot
            ct$y <- ct$y * bin / s[[fov[1]]]@scale.factors$spot
        }
        sl <- getSpatialLocations(g, spat_unit = "cell", output = "data.table")
        m <- merge(sl, ct, by.x = "cell_ID", by.y = "cell")
        if (nrow(m) != nrow(ct)) return(row("wrong", sprintf("%d of %d cells", nrow(m), nrow(ct))))
        eq <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-8))
        orient <- if (eq(m$sdimx, m$x) && eq(m$sdimy, -m$y)) "x,-y" else
            if (eq(m$sdimx, m$y) && eq(m$sdimy, -m$x)) "y,-x (axes swapped)" else
            if (eq(m$sdimx, m$x) && eq(m$sdimy, m$y)) "x,y (not flipped)" else "other"
        row(if (orient == "x,-y") "correct" else "wrong", orient)
    },
    polygons = function(g, s, a) {
        fov <- Filter(function(i) inherits(s[[i]], "FOV"), Images(s))
        if (!length(fov)) return(row("n/a", "no FOV"))
        seg <- Filter(function(b) inherits(s[[fov[1]]][[b]], "Segmentation"), Boundaries(s[[fov[1]]]))
        if (!length(seg)) return(row("n/a", "no segmentation"))
        src <- GetTissueCoordinates(s[[fov[1]]][[seg[1]]])
        sb <- data.table::as.data.table(src)[, .(xmin = min(x), xmax = max(x), ymin = min(-y), ymax = max(-y)), by = cell]
        p <- getPolygonInfo(g, name = "cell", return_giottoPolygon = TRUE)
        sv <- if (inherits(p[], "SpatVector")) p[] else GiottoClass::as.terra(p[])
        geo <- data.table::as.data.table(terra::geom(sv))
        ids <- GiottoClass::spatIDs(p)
        gb <- geo[, .(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)), by = geom]
        gb$cell <- ids[gb$geom]
        mm <- merge(sb, gb, by = "cell")
        if (nrow(mm) != nrow(sb)) return(row("wrong", sprintf("%d of %d polygons", nrow(mm), nrow(sb))))
        d <- max(abs(c(mm$xmin.x - mm$xmin.y, mm$xmax.x - mm$xmax.y, mm$ymin.x - mm$ymin.y, mm$ymax.x - mm$ymax.y)))
        n_bad <- sum(apply(abs(cbind(mm$xmin.x - mm$xmin.y, mm$xmax.x - mm$xmax.y,
            mm$ymin.x - mm$ymin.y, mm$ymax.x - mm$ymax.y)), 1, max) > 1e-6)
        valid <- terra::is.valid(terra::vect(lapply(split(src, src$cell), function(v) cbind(v$x, -v$y)), type = "polygons"))
        st <- if (d < 1e-6) "correct" else if (n_bad <= sum(!valid)) "correct (repaired)" else "wrong"
        row(st, sprintf("%d polygons, %d differ (max %.3g), %d invalid in source", nrow(sb), n_bad, d, sum(!valid)))
    },
    pca = function(g, s, a) .check_dr(g, s, "pca"),
    umap = function(g, s, a) .check_dr(g, s, "umap"),
    nn = function(g, s, a) {
        gn <- paste0(a, "_nn")
        if (!gn %in% Graphs(s)) return(row("n/a"))
        nnl <- list_nearest_networks(g)
        if (is.null(nnl) || !nrow(nnl)) return(row("missing"))
        e <- getNearestNetwork(g, spat_unit = "cell", feat_type = "rna", nn_type = nnl$nn_type[1],
            name = nnl$name[1], output = "data.table")
        tr <- Matrix::summary(methods::as(s[[gn]], "CsparseMatrix"))
        ids <- rownames(s[[gn]])
        src <- unique(paste(ids[tr$i], ids[tr$j])[tr$i != tr$j])
        und <- function(f, t) ifelse(f < t, paste(f, t), paste(t, f))
        got <- unique(und(e$from, e$to))
        exp_und <- unique(und(ids[tr$i][tr$i != tr$j], ids[tr$j][tr$i != tr$j]))
        ok <- setequal(got, exp_und)
        row(if (ok) "correct" else "wrong", sprintf("type %s name %s; %d edges (source %d directed)",
            nnl$nn_type[1], nnl$name[1], nrow(e), length(src)))
    }
)
.check_dr <- function(g, s, r) {
    if (!r %in% Reductions(s)) return(row("n/a"))
    dr <- list_dim_reductions(g)
    hit <- dr[dr$name == r, ]
    if (!nrow(hit)) return(row("missing"))
    m <- getDimReduction(g, spat_unit = "cell", feat_type = "rna", reduction_method = hit$dim_type[1],
        name = r, output = "matrix")
    src <- Embeddings(s, r)
    ok <- isTRUE(all.equal(as.numeric(m[rownames(src), ]), as.numeric(src)))
    row(if (ok) "correct" else "wrong", sprintf("%dx%d", nrow(m), ncol(m)))
}

convert <- list(
    new_memory = function(s, a, d) GiottoConnect::seuratToGiotto(s),
    new_disk = function(s, a, d) GiottoConnect::seuratToGiotto(s, source = d),
    legacy = function(s, a, d) GiottoClass::seuratToGiottoV5(s, spatial_assay = a,
        dim_reduction = intersect(c("pca", "umap"), Reductions(s)),
        nn_network = intersect(paste0(a, "_nn"), Graphs(s)), verbose = FALSE)
)

res <- list()
timing <- list()
objs <- list()
for (ds in names(datasets)) {
    s <- datasets[[ds]]()
    a <- DefaultAssay(s)
    for (impl in names(convert)) {
        d <- tempfile(paste0(ds, "_"))
        t0 <- Sys.time()
        g <- tryCatch(suppressWarnings(suppressMessages(convert[[impl]](s, a, d))),
            error = function(e) e)
        secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        timing[[length(timing) + 1L]] <- data.frame(dataset = ds, impl = impl, seconds = round(secs, 2),
            error = if (inherits(g, "error")) conditionMessage(g) else "")
        if (inherits(g, "error")) {
            for (pt in names(checks)) res[[length(res) + 1L]] <- data.frame(dataset = ds, part = pt,
                impl = impl, status = "error", detail = substr(conditionMessage(g), 1, 80))
            next
        }
        objs[[paste(ds, impl)]] <- g
        for (pt in names(checks)) {
            r <- tryCatch(checks[[pt]](g, s, a), error = function(e) row("missing", substr(conditionMessage(e), 1, 80)))
            res[[length(res) + 1L]] <- data.frame(dataset = ds, part = pt, impl = impl,
                status = r$status, detail = r$detail)
        }
    }
}
parts <- do.call(rbind, res)
timing <- do.call(rbind, timing)
wide <- reshape(parts[, c("dataset", "part", "impl", "status")], idvar = c("dataset", "part"),
    timevar = "impl", direction = "wide")
names(wide) <- sub("^status\\.", "", names(wide))
write.csv(parts, file.path(out_dir, paste0(Sys.getenv("RESULT_TAG", "phase2"), "_parts.csv")), row.names = FALSE)
write.csv(timing, file.path(out_dir, paste0(Sys.getenv("RESULT_TAG", "phase2"), "_timing.csv")), row.names = FALSE)
print(wide, row.names = FALSE)
cat("\n")
print(timing, row.names = FALSE)
cat("\nDetails where new and legacy differ:\n")
for (i in seq_len(nrow(wide))) {
    w <- wide[i, ]
    if (identical(w$new_memory, w$legacy)) next
    det <- parts[parts$dataset == w$dataset & parts$part == w$part, c("impl", "status", "detail")]
    cat(sprintf("- %s / %s: %s\n", w$dataset, w$part,
        paste(sprintf("%s=%s (%s)", det$impl, det$status, det$detail), collapse = "; ")))
}

# independent oracle: GiottoDisk's own Xenium importer on the same crop
crops <- path.expand("~/Library/Caches/org.R-project.R/R/GiottoData")
gx <- suppressWarnings(suppressMessages(GiottoDisk::importXeniumDisk(
    file.path(crops, "xenium_mini_lung_ffpe"), backend = tempfile("oracle_"))$create_gobject()))
gn <- objs[["xenium new_disk"]]
o_ids <- GiottoClass::spatIDs(gx, spat_unit = "cell")
n_ids <- GiottoClass::spatIDs(gn, spat_unit = "cell")
oe <- as_mat(GiottoClass::getExpression(gx, values = "raw", spat_unit = "cell", feat_type = "rna")[])
ne <- as_mat(GiottoClass::getExpression(gn, values = "raw", spat_unit = "cell", feat_type = "rna")[])
fg <- intersect(rownames(oe), rownames(ne))
cc <- intersect(o_ids, n_ids)
ol <- GiottoClass::getSpatialLocations(gx, spat_unit = "cell", output = "data.table")
nl <- GiottoClass::getSpatialLocations(gn, spat_unit = "cell", output = "data.table")
ml <- merge(ol, nl, by = "cell_ID")
oracle <- data.frame(
    check = c("cells", "features", "raw counts on shared cells x genes", "locs max |dx|", "locs max |dy|",
        "polygons"),
    giottodisk_import = c(length(o_ids), nrow(oe), "", "", "",
        length(GiottoClass::getPolygonInfo(gx, name = "cell", return_giottoPolygon = TRUE)$poly_ID)),
    seuratToGiotto = c(length(n_ids), nrow(ne), "", "", "",
        length(GiottoClass::getPolygonInfo(gn, name = "cell", return_giottoPolygon = TRUE)$poly_ID)),
    agreement = c(sprintf("%d shared", length(cc)),
        sprintf("%d shared; only import: %s; only Seurat: %s", length(fg),
            paste(setdiff(rownames(oe), rownames(ne)), collapse = ","),
            paste(setdiff(rownames(ne), rownames(oe)), collapse = ",")),
        if (isTRUE(all.equal(as.numeric(oe[fg, cc]), as.numeric(ne[fg, cc])))) "identical" else "differ",
        sprintf("%.3f", max(abs(ml$sdimx.x - ml$sdimx.y))),
        sprintf("%.3f", max(abs(ml$sdimy.x - ml$sdimy.y))),
        "")
)
write.csv(oracle, file.path(out_dir, paste0(Sys.getenv("RESULT_TAG", "phase2"), "_xenium_oracle.csv")), row.names = FALSE)
cat("\nXenium: seuratToGiotto(LoadXenium()) vs GiottoDisk::importXeniumDisk():\n")
print(oracle, row.names = FALSE)
