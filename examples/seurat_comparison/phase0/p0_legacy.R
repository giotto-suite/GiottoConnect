options(giotto.has_conda = FALSE, giotto.use_conda = FALSE)
suppressMessages({ library(GiottoClass); library(SeuratObject) })
SP <- Sys.getenv("SP", tempdir())
s <- readRDS(file.path(SP, "fixture_seurat5.rds"))
gc_dir <- "/Users/rubendries/Documents/Software/GiottoSuite/GiottoConnect/R"
e <- new.env(parent = asNamespace("GiottoClass"))
for (f in c("seurat.R", "spatialexperiment.R")) sys.source(file.path(gc_dir, f), envir = e)

run <- function(label, f) {
  cat("\n#### ", label, "\n")
  t0 <- Sys.time()
  res <- tryCatch(f(), error = function(err) { cat("ERROR:", conditionMessage(err), "\n"); NULL })
  cat("time:", format(round(difftime(Sys.time(), t0, units = "secs"), 2)), "\n")
  if (inherits(res, "giotto")) {
    man <- objManifest(res, fingerprint = "sample")
    saveRDS(man, file.path(SP, paste0("legacy_manifest_", gsub("[^A-Za-z0-9]", "_", label), ".rds")))
    print(res)
    cat("feat_metadata cols:", paste(names(pDataDT(res, spat_unit = "cell", feat_type = "rna") |> invisible()), collapse = ","), "\n")
    cat("fDataDT cols:", tryCatch(paste(names(fDataDT(res)), collapse = ","), error = function(e) conditionMessage(e)), "\n")
    cat("nn names:", tryCatch(paste(list_nearest_networks(res)$name, collapse = ","), error = function(e) "none"), "\n")
    cat("polys:", tryCatch(paste(list_spatial_info(res)$spat_info, collapse = ","), error = function(e) "none"), "\n")
  } else cat("returned class:", class(res), "\n")
  invisible(res)
}
r1 <- run("GiottoConnect seuratToGiottoV5", function()
  e$seuratToGiottoV5(s, spatial_assay = "Xenium", dim_reduction = c("pca", "umap"),
                     nn_network = c("Xenium_nn"), verbose = FALSE))
r2 <- run("GiottoClass seuratToGiottoV5", function()
  GiottoClass::seuratToGiottoV5(s, spatial_assay = "Xenium", dim_reduction = c("pca", "umap"),
                                nn_network = c("Xenium_nn"), verbose = FALSE))
cat("\n#### copies identical (deparse)?",
    identical(deparse(e$seuratToGiottoV5), deparse(GiottoClass::seuratToGiottoV5)), "\n")
