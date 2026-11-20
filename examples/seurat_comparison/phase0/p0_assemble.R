options(giotto.has_conda = FALSE, giotto.no_python_warn = TRUE)
suppressMessages({ library(GiottoClass); library(SeuratObject); library(GiottoDisk) })
SP <- Sys.getenv("SP", tempdir())
s <- readRDS(file.path(SP, "fixture_seurat5.rds"))
cnt <- LayerData(s, "counts", assay = "Xenium"); nrm <- LayerData(s, "data", assay = "Xenium")
ct <- GetTissueCoordinates(s[["fov"]], which = "centroids")
sg <- GetTissueCoordinates(s[["fov"]][["segmentation"]])
G <- s[["Xenium_nn"]]; sm <- Matrix::summary(as(G, "CsparseMatrix"))
edges <- data.table::data.table(from = rownames(G)[sm$i], to = colnames(G)[sm$j], weight = sm$x, distance = 1 - sm$x)

parts <- function() list(
  createExprObj(cnt, name = "raw", spat_unit = "cell", feat_type = "rna"),
  createExprObj(nrm, name = "normalized", spat_unit = "cell", feat_type = "rna"),
  createCellMetaObj(data.table::data.table(cell_ID = rownames(s[[]]), s[[]]), spat_unit = "cell", feat_type = "rna"),
  createSpatLocsObj(data.table::data.table(sdimx = ct$x, sdimy = ct$y, cell_ID = ct$cell), name = "raw", spat_unit = "cell"),
  createGiottoPolygonsFromDfr(data.table::data.table(poly_ID = sg$cell, x = sg$x, y = sg$y), name = "cell", calc_centroids = FALSE),
  createDimObj(Embeddings(s, "pca"), name = "pca", method = "pca", spat_unit = "cell", feat_type = "rna"),
  createDimObj(Embeddings(s, "umap"), name = "umap", method = "umap", spat_unit = "cell", feat_type = "rna"),
  createNearestNetObj(name = "kNN.pca", network = edges, nn_type = "kNN", spat_unit = "cell", feat_type = "rna")
)
slot_backing <- function(g) {
  b <- function(x) if (inherits(x, "dataStore")) class(x)[1] else paste0("mem:", class(x)[1])
  c(expr_raw = b(getExpression(g, values = "raw")[]),
    polys = b(getPolygonInfo(g, return_giottoPolygon = TRUE)@spatVector),
    nn = b(getNearestNetwork(g, nn_type = "kNN", name = "kNN.pca")@network),
    locs = b(getSpatialLocations(g)@coordinates),
    meta = b(getCellMetadata(g)@metaDT))
}
# in memory
t0 <- Sys.time()
gm <- setGiotto(giotto(), parts(), verbose = FALSE)
cat("memory assemble s:", round(as.numeric(Sys.time() - t0), 2), "\n"); print(gm)
print(slot_backing(gm))
# as.list round trip
gr <- setGiotto(giotto(), as.list(gm), verbose = FALSE)
d <- manifestDiff(objManifest(gm, fingerprint = "sample"), objManifest(gr, fingerprint = "sample"))
cat("as.list round trip changed:", d$changed, "\n"); if (isTRUE(d$changed)) print(d$summary)
# disk
proj <- file.path(SP, "p0_disk_proj"); unlink(proj, recursive = TRUE)
src <- sourceCreate(proj, type = "gDirSource")
t0 <- Sys.time()
gd <- setGiotto(giotto(source = src), parts(), initialize = FALSE, verbose = FALSE)
gd <- methods::initialize(gd)
cat("disk assemble s:", round(as.numeric(Sys.time() - t0), 2), "\n")
print(slot_backing(gd))
md <- manifestDiff(objManifest(gm, fingerprint = "sample"), objManifest(gd, fingerprint = "sample"))
cat("memory vs disk manifest changed:", md$changed, "\n"); print(md$summary)
lm <- function(g, f) { m <- objManifest(g, fingerprint = "none"); f(m$slots) }
cat("dtype mem/disk:", lm(gm, function(x) x$expression$cell$rna$raw$dtype), "/", lm(gd, function(x) x$expression$cell$rna$raw$dtype), "\n")
cat("directed mem/disk:", lm(gm, function(x) x$nn_network$cell$rna$kNN$kNN.pca$directed), "/", lm(gd, function(x) x$nn_network$cell$rna$kNN$kNN.pca$directed), "\n")
st <- getExpression(gd, values = "raw")[]
cat("store class:", class(st), " dim:", dim(st), "\n")
cat("methods for as.matrix/[ :", paste(head(showMethods("[", classes = class(st)[1], printTo = FALSE), 4), collapse=" | "), "\n")
sub <- tryCatch(st[rownames(cnt)[1:5], colnames(cnt)[1:7]], error = function(e) paste("ERR", conditionMessage(e)))
cat("subset class:", class(sub), "\n"); print(sub)
mat <- tryCatch(as(st, "dgCMatrix"), error = function(e) tryCatch(as.matrix(st), error = function(e2) paste("ERR", conditionMessage(e2))))
cat("materialised class:", class(mat), "\n")
if (inherits(mat, c("Matrix","matrix"))) cat("disk raw readback equal:", isTRUE(all.equal(as.matrix(mat[rownames(cnt), colnames(cnt)]), as.matrix(cnt), check.attributes = FALSE)), "\n")
cat("project files:", length(list.files(proj, recursive = TRUE)), "\n")
fl <- function(g) { m <- objManifest(g, fingerprint = "none"); rapply(m$slots, function(z) z, how = "unlist") }
a <- fl(gm); b <- fl(gd); k <- union(names(a), names(b))
dd <- data.frame(key = k, mem = unname(a[k]), disk = unname(b[k]))
print(dd[grepl("dtype|directed|sparse|nnz", dd$key), ], row.names = FALSE)
