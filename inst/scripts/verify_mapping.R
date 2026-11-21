# Verify a bundled mapping on its fixture and stamp it: converts in memory and
# on disk, requires validation without FAIL and a passing readback, then
# writes the fingerprint,
# status and build date into the YAML and stores the fixture profile next to
# it (<id>.profile.json) for drift detection.
# Run from the GiottoConnect package root:
#   Rscript inst/scripts/verify_mapping.R seurat-assay5
args <- commandArgs(trailingOnly = TRUE)
id <- if (length(args)) args[[1]] else "seurat-assay5"
pkgload::load_all(".", quiet = TRUE)
env <- new.env()
sys.source(file.path("inst", "scripts", "build_fixtures.R"), envir = env)
builders <- list(
    seurat5_fixture_small = env$build_seurat_fixture,
    seurat5_visium_fixture = env$build_seurat_visium_fixture,
    seurat_v3_fixture = env$build_seurat_v3_fixture,
    giotto_fixture_small = env$build_giotto_fixture,
    seurat_visium_v1_fixture = env$build_seurat_visium_v1_fixture,
    seurat_slideseq_fixture = env$build_seurat_slideseq_fixture,
    anndata_fixture = function() env$build_h5ad_fixture(tempfile(fileext = ".h5ad"))$path
)

yaml_file <- file.path("inst", "mappings", paste0(id, ".yaml"))
m <- readMapping(yaml_file)
x <- builders[[m$mapping$fixtures[[1]]]]()
p <- adapterProfile(GiottoConnect:::.adapter_for(m), x)
g <- convertWithMapping(x, m, profile = p)
if (identical(m$mapping$target$contract, "giotto")) {
    # the same mapping onto a GiottoDisk project
    res <- convertWithMapping(x, m, profile = p, source = tempfile("verify_"))
    v <- validateMapping(m, p, result = res)
} else {
    # a GiottoDisk-backed source through the same mapping
    xd <- builders[[m$mapping$fixtures[[1]]]](source = tempfile("verify_"))
    res <- convertWithMapping(xd, m)
    v <- validateMapping(m, profileObject(xd), result = res)
}
if (any(v$level == "FAIL")) stop("validation failed; mapping not stamped")
rb <- readbackCheck(x, g, m)
if (!isTRUE(attr(rb, "pass"))) {
    print(rb[rb$status == "UNDECLARED", ])
    stop("readback failed; mapping not stamped")
}

fp <- mappingProvenance(g)$fingerprint
txt <- readLines(yaml_file)
set_field <- function(txt, field, value, after) {
    i <- grep(sprintf("^  %s:", field), txt)
    line <- sprintf("  %s: %s", field, value)
    if (length(i)) txt[i[1]] <- line else {
        j <- grep(sprintf("^  %s:", after), txt)[1]
        txt <- append(txt, line, after = j)
    }
    txt
}
txt <- set_field(txt, "status", "verified", "version")
txt <- set_field(txt, "built", sprintf('"%s"', Sys.Date()), "status")
txt <- set_field(txt, "fingerprint", sprintf('"%s"', fp), "target")
writeLines(txt, yaml_file)
GiottoConnect:::.profile_write(p, file.path("inst", "mappings", paste0(id, ".profile.json")))
message(sprintf("%s verified: fingerprint %s (%d OK, %d WARN)", id, fp,
    sum(v$level == "OK"), sum(v$level == "WARN")))
