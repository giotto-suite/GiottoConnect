# Integration test: full GiottoConnect pipeline with mini MERFISH dataset
#
# Run manually from R: source("tests/test_mini_merfish.R")
#
# Requirements:
#   - conda (miniconda or anaconda) installed and on PATH
#   - GiottoData package installed
#   - Internet access (to download mini MERFISH dataset on first run)
#   - GiottoConnect installed or loaded via devtools::load_all()

library(GiottoConnect)
library(GiottoData)

# ---------------------------------------------------------------------------
# 0. Set up Python environment (skip if already done in this session)
# ---------------------------------------------------------------------------
setup_giotto_connect_env()

# ---------------------------------------------------------------------------
# 1. Load mini MERFISH Giotto object
# ---------------------------------------------------------------------------
message("Loading mini MERFISH Giotto object ...")
gobject <- GiottoData::loadGiottoMini("merfish")

# Inspect available slots to verify obs_set_paths below
message("\nAvailable spatial units / feature types:")
print(GiottoClass::list_spatial_locations(gobject))

message("\nAvailable dimension reductions:")
print(GiottoClass::list_dim_reductions(gobject))

message("\nAvailable cell metadata columns:")
meta <- GiottoClass::getCellMetadata(gobject, output = "data.table")
print(names(meta))

# ---------------------------------------------------------------------------
# 2. Run full pipeline
# ---------------------------------------------------------------------------
output_dir <- tempfile("gc_test_")
dir.create(output_dir)
message("\nOutput directory: ", output_dir)

result_dir <- giotto_to_vitessce_site(
    gobject         = gobject,
    output_dir      = output_dir,
    repo_name       = "merfish-vitessce-demo",
    github_username = "testuser",
    dataset_name    = "mini_merfish",
    data_mode       = "github",
    verbose         = TRUE
)

# ---------------------------------------------------------------------------
# 3. Verify output files
# ---------------------------------------------------------------------------
message("\n--- Verifying output files ---")

expected_files <- c(
    "src/App.jsx",
    "src/main.jsx",
    "src/index.css",
    "src/vitessce_config.json",
    "index.html",
    "vite.config.js",
    "package.json",
    ".github/workflows/deploy.yml",
    ".gitignore"
)

all_ok <- TRUE
for (f in expected_files) {
    path <- file.path(result_dir, f)
    if (file.exists(path)) {
        message("  OK : ", f)
    } else {
        message("  MISSING: ", f)
        all_ok <- FALSE
    }
}

# ---------------------------------------------------------------------------
# 4. Verify zarr structure
# ---------------------------------------------------------------------------
message("\n--- Verifying zarr store ---")
zarr_dir <- file.path(result_dir, "public", "data", "mini_merfish.zarr")

if (dir.exists(zarr_dir)) {
    message("  OK : public/data/mini_merfish.zarr/ exists")
} else {
    message("  MISSING: public/data/mini_merfish.zarr/")
    all_ok <- FALSE
}

zmeta_path <- file.path(zarr_dir, ".zmetadata")
if (file.exists(zmeta_path)) {
    message("  OK : .zmetadata found (zarr metadata consolidated)")
} else {
    message("  MISSING: .zmetadata (run zarr.consolidate_metadata())")
    all_ok <- FALSE
}

# ---------------------------------------------------------------------------
# 5. Validate Vitessce config JSON
# ---------------------------------------------------------------------------
message("\n--- Validating Vitessce config ---")
config_path <- file.path(result_dir, "src", "vitessce_config.json")

if (file.exists(config_path)) {
    config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)

    if (!is.null(config$version)) {
        message("  OK : config$version = ", config$version)
    } else {
        message("  WARNING: config$version not found")
    }

    if (!is.null(config$datasets) && length(config$datasets) >= 1L) {
        message("  OK : ", length(config$datasets), " dataset(s) in config")
    } else {
        message("  WARNING: no datasets found in config")
    }
} else {
    message("  MISSING: vitessce_config.json")
    all_ok <- FALSE
}

# ---------------------------------------------------------------------------
# 6. Check vite.config.js has correct base
# ---------------------------------------------------------------------------
vite_cfg <- readLines(file.path(result_dir, "vite.config.js"), warn = FALSE)
if (any(grepl("base: '/merfish-vitessce-demo/'", vite_cfg))) {
    message("  OK : vite.config.js base path is correct")
} else {
    message("  WARNING: vite.config.js base path may be incorrect")
}

# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------
message("\n")
if (all_ok) {
    message("=== All checks PASSED ===")
} else {
    message("=== Some checks FAILED — see above ===")
}

message("\nTo test the site locally:")
message("  cd ", result_dir)
message("  npm install && npm run dev")
message("\nThen open http://localhost:5173/merfish-vitessce-demo/ in your browser.")
