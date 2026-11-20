# anndataToGiotto on real .h5ad files: mapping choice, bound params, parts,
# readback, losses and time, in memory and onto GiottoDisk.
#   Rscript examples/anndata_comparison/real_files.R
# Run from the GiottoConnect package root. Writes results/real_files.csv.
options(giotto.has_conda = FALSE, giotto.use_conda = FALSE, giotto.no_python_warn = TRUE,
    lifecycle_verbosity = "quiet")
suppressMessages(pkgload::load_all(".", quiet = TRUE))
ds <- path.expand("~/Documents/Datasets")
files <- c(
    seqfish_squidpy = "squidpy_examples/seqfish.h5ad",
    merfish_squidpy = "squidpy_examples/merfish.h5ad",
    slideseqv2_squidpy = "squidpy_examples/slideseqv2.h5ad",
    visium_kuppe_P12 = "Visium/Kuppe2022_MI_zenodo/Visium_RZ_BZ_P12.h5ad",
    seqfish_aki = "SeqFISH/AKI_seqFISH_allCS_subsampled_dr_25000.h5ad"
)
secs <- function(t0) round(as.numeric(Sys.time() - t0, units = "secs"), 2)
rows <- list()
for (n in names(files)) {
    f <- file.path(ds, files[[n]])
    t0 <- Sys.time()
    p <- adapterProfile(.adapter_for("anndata"), f)
    t_profile <- secs(t0)
    fm <- findMapping(p, ecosystem = "anndata")
    t0 <- Sys.time()
    g <- suppressWarnings(anndataToGiotto(f))
    t_mem <- secs(t0)
    rb <- readbackCheck(f, g)
    prov <- mappingProvenance(g)
    t0 <- Sys.time()
    gd <- suppressWarnings(anndataToGiotto(f, source = tempfile("anndata_disk_")))
    t_disk <- secs(t0)
    rbd <- readbackCheck(f, gd)
    v <- validateMapping(fm$mapping, p, result = g)
    rows[[n]] <- data.frame(
        file = n, size_mb = round(file.size(f) / 2^20), cells = p$n_cells, features = p$n_feat,
        mapping = fm$mapping$mapping$id, status = fm$status,
        params = paste(names(prov$params), unlist(prov$params), sep = "=", collapse = "; "),
        parts = paste(prov$parts, collapse = ", "),
        readback_match = sum(rb$status == "MATCH"), readback_declared = sum(rb$status == "DECLARED"),
        readback_undeclared = sum(rb$status == "UNDECLARED"),
        readback_disk_pass = isTRUE(attr(rbd, "pass")),
        warnings = paste(unique(v$check[v$level == "WARN"]), collapse = ", "),
        profile_s = t_profile, memory_s = t_mem, disk_s = t_disk,
        stringsAsFactors = FALSE)
    print(rows[[n]])
}
out <- do.call(rbind, rows)
utils::write.csv(out, "examples/anndata_comparison/results/real_files.csv", row.names = FALSE)
