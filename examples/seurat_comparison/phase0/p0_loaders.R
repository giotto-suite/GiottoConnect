suppressMessages({ library(Seurat); library(SeuratObject) })
options(arrow.unsafe_metadata = TRUE)
C <- path.expand("~/Library/Caches/org.R-project.R/R/GiottoData")
SP <- Sys.getenv("SP", tempdir())
desc <- function(s, label) {
  cat("\n####", label, "\n")
  if (is.null(s)) return(invisible())
  cat("cells x feats:", ncol(s), "x", nrow(s), "| assays:", Assays(s), "| layers:", Layers(s),
      "| images:", Images(s), "\n")
  for (im in Images(s)) cat("  image", im, "class:", class(s[[im]]),
      "| boundaries:", tryCatch(paste(Boundaries(s[[im]]), collapse = ","), error = function(e) "-"), "\n")
  saveRDS(s, file.path(SP, paste0("seurat_", label, ".rds")))
}
try1 <- function(label, f) {
  t0 <- Sys.time()
  s <- tryCatch(f(), error = function(e) { cat("\n####", label, "\nERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(s)) { desc(s, label); cat("time:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n") }
  s
}
x <- try1("xenium_mini", function() LoadXenium(file.path(C, "xenium_mini_lung_ffpe"), fov = "fov"))
cx <- try1("cosmx_mini", function() LoadNanostring(file.path(C, "cosmx_mini_colon_wtx"), assay = "Nanostring", fov = "fov"))
vh <- try1("visiumhd_mini", function() Load10X_Spatial(file.path(C, "visiumhd_mini_crc"), bin.size = 16))
