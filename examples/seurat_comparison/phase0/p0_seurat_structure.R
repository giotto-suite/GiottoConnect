suppressMessages(library(SeuratObject))
set.seed(1)
nc <- 200; nf <- 50
cells <- sprintf("cell_%03d", seq_len(nc)); feats <- sprintf("gene%02d", seq_len(nf))
m <- matrix(rpois(nc * nf, 2), nf, nc, dimnames = list(feats, cells))
m <- as(m, "CsparseMatrix")
s <- CreateSeuratObject(CreateAssay5Object(counts = m), assay = "Xenium")
LayerData(s, "data") <- log1p(m)
LayerData(s, "scale.data") <- as.matrix(scale(t(as.matrix(log1p(m)))))[, feats] |> t()
s$seurat_clusters <- factor(sample(0:3, nc, TRUE))
xy <- data.frame(x = runif(nc, 0, 2000), y = runif(nc, 0, 2000), cell = cells)
cen <- CreateCentroids(xy)
sq <- do.call(rbind, lapply(seq_len(nc), function(i)
  data.frame(x = xy$x[i] + c(-5, 5, 5, -5), y = xy$y[i] + c(-5, -5, 5, 5), cell = cells[i])))
seg <- CreateSegmentation(sq)
fov <- CreateFOV(list(centroids = cen, segmentation = seg),
                 type = c("centroids", "segmentation"), assay = "Xenium")
s[["fov"]] <- fov
pcs <- prcomp(t(as.matrix(log1p(m))), rank. = 10)
s[["pca"]] <- CreateDimReducObject(embeddings = `rownames<-`(pcs$x, cells) |> `colnames<-`(paste0("PC_", 1:10)),
  loadings = `colnames<-`(pcs$rotation, paste0("PC_", 1:10)), stdev = pcs$sdev[1:10], key = "PC_", assay = "Xenium")
s[["umap"]] <- CreateDimReducObject(embeddings = `dimnames<-`(matrix(rnorm(nc * 2), nc), list(cells, c("umap_1", "umap_2"))),
  key = "umap_", assay = "Xenium")
# kNN graph k=10 (binary) and an snn (weighted)
d <- as.matrix(dist(pcs$x[, 1:5]))
nn <- Matrix::sparseMatrix(i = rep(seq_len(nc), each = 10),
  j = as.vector(apply(d, 1, function(r) order(r)[1:10])), x = 1, dims = c(nc, nc), dimnames = list(cells, cells))
s[["Xenium_nn"]] <- as.Graph(nn)
snn <- nn %*% Matrix::t(nn); snn@x <- snn@x / 10
s[["Xenium_snn"]] <- as.Graph(snn)
saveRDS(s, file.path(Sys.getenv("SP", tempdir()), "fixture_seurat5.rds"))

cat("== Version:", as.character(Version(s)), " class attr pkg:", attr(class(s), "package"), "\n")
cat("== Seurat slots:", paste(slotNames(s), collapse = ", "), "\n")
a <- s[["Xenium"]]
cat("== Assay5 class/slots:", class(a), "|", paste(slotNames(a), collapse = ", "), "\n")
cat("   layers:", paste(names(a@layers), collapse = ","), " counts class:", class(a@layers$counts),
    " dimnames NULL:", is.null(dimnames(a@layers$counts)[[1]]), "\n")
cat("   cells slot:", class(a@cells), "dim", dim(a@cells), " colnames:", paste(colnames(a@cells), collapse = ","), "\n")
cat("   LayerData named:", !is.null(rownames(LayerData(s, "counts", assay = "Xenium"))), "\n")
cat("== FOV slots:", paste(slotNames(fov), collapse = ", "), "| boundaries:", paste(names(s[["fov"]]@boundaries), collapse = ","),
    "| Boundaries():", paste(Boundaries(s[["fov"]]), collapse = ","), "\n")
cat("   Centroids class/slots:", class(cen), "|", paste(slotNames(cen), collapse = ", "), " coords:", class(cen@coords), dim(cen@coords), "\n")
cat("   Segmentation class/is sp:", class(seg), is(seg, "SpatialPolygons"), "|", paste(slotNames(seg), collapse = ", "), "\n")
cat("   GetTissueCoordinates(fov) cols:", paste(names(GetTissueCoordinates(s[["fov"]])), collapse = ","), "\n")
cat("   GetTissueCoordinates(fov, which=centroids) cols:", paste(names(GetTissueCoordinates(s[["fov"]], which = "centroids")), collapse = ","), "\n")
sg <- GetTissueCoordinates(s[["fov"]][["segmentation"]])
cat("   segmentation coords:", nrow(sg), "rows, cols", paste(names(sg), collapse = ","), "\n")
g <- s[["Xenium_nn"]]
cat("== Graph class/super:", class(g), is(g, "dgCMatrix"), "| slots extra:", paste(setdiff(slotNames(g), slotNames("dgCMatrix")), collapse = ","),
    "| assay.used:", g@assay.used, "| nnz/n:", length(g@x) / nrow(g), " binary:", all(g@x == 1), "\n")
cat("   snn binary:", all(s[["Xenium_snn"]]@x == 1), "\n")
r <- s[["pca"]]
cat("== DimReduc slots:", paste(slotNames(r), collapse = ", "), "| emb:", dim(Embeddings(r)), "| loadings:", dim(Loadings(r)), "| stdev:", length(Stdev(r)), "\n")
cat("== listing verbs: Assays", Assays(s), "| Layers", Layers(s), "| Reductions", Reductions(s),
    "| Graphs", Graphs(s), "| Neighbors", length(Neighbors(s)), "| Images", Images(s), "\n")
cat("== meta.data cols:", paste(names(s[[]]), collapse = ","), "\n")
# raw attribute walk without methods: what does attributes() expose?
cat("== attributes(s) names:", paste(names(attributes(s)), collapse = ","), "\n")
cat("== attributes(Assay5) names:", paste(names(attributes(a)), collapse = ","), "\n")
