suppressMessages({ library(Seurat); library(SeuratObject) }); options(arrow.unsafe_metadata = TRUE)
C <- path.expand("~/Library/Caches/org.R-project.R/R/GiottoData")
print(names(formals(LoadXenium))); print(names(formals(ReadXenium)))
x <- tryCatch(LoadXenium(file.path(C, "xenium_mini_lung_ffpe"), fov = "fov", segmentations = "cell"), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
if (!is.null(x)) { cat("boundaries:", Boundaries(x[["fov"]]), "\n"); saveRDS(x, file.path(Sys.getenv("SP", tempdir()), "seurat_xenium_mini.rds")) }
md <- data.table::fread(file.path(C, "cosmx_mini_colon_wtx/S0_metadata_file.csv.gz"), select = c("fov","cell_ID"))
cat("cosmx metadata cells:", nrow(md), " fovs:", paste(sort(unique(md$fov)), collapse=","), "\n")
