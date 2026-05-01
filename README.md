# GiottoConnect

**GiottoConnect** is an R package that serves as a hub for Giotto interoperability functions and connections to external visualization and analysis tools.

The first major feature is a **Vitessce viewer pipeline** that converts a Giotto spatial single-cell object into a static, interactive website deployable via GitHub Pages — no server required.

## Pipeline

```
Giotto object
    ↓  export_giotto_to_zarr()       Stage 1
SpatialData-Zarr store
    ↓  create_vitessce_config()      Stage 2
vitessce_config.json
    ↓  scaffold_vitessce_site()      Stage 3
GitHub Pages repo (Vite + React + Vitessce)
    ↓  git push + enable Pages
Live interactive viewer
```

Or in one call:

```r
giotto_to_vitessce_site(gobject, output_dir, repo_name, github_username)
```

## Installation

```r
# Install from GitHub (once published)
remotes::install_github("drieslab/GiottoConnect")
```

### Python dependencies

GiottoConnect calls Python via `reticulate`. Set up the conda environment once:

```r
library(GiottoConnect)
setup_giotto_connect_env()   # creates "giotto-connect" conda env
```

This installs: `anndata`, `zarr`, `vitessce`, `spatialdata`, `numpy`.

## Quick Start

```r
library(GiottoConnect)
library(GiottoData)

# 1. Set up Python env (once per machine)
setup_giotto_connect_env()

# 2. Load a Giotto object
gobject <- GiottoData::loadGiottoMini("merfish")

# 3. Run the full pipeline
result_dir <- giotto_to_vitessce_site(
    gobject         = gobject,
    output_dir      = "~/vitessce-output",
    repo_name       = "merfish-viewer",
    github_username = "myuser",
    dataset_name    = "mini_merfish"
)
```

Follow the printed instructions to push to GitHub and enable Pages.

## Function Reference

| Function | Description |
|---|---|
| `setup_giotto_connect_env()` | Set up the Python conda environment |
| `export_giotto_to_zarr()` | Stage 1: Giotto → SpatialData-Zarr |
| `create_vitessce_config()` | Stage 2: Generate Vitessce JSON config |
| `scaffold_vitessce_site()` | Stage 3: Scaffold Vite/React GitHub Pages repo |
| `giotto_to_vitessce_site()` | Convenience wrapper for all three stages |

## Deployment Modes

**`data_mode = "github"` (default):** The zarr data is copied into `public/data/` in the repo and served via GitHub Pages. Suitable for datasets < ~500 MB.

**`data_mode = "remote"`:** The zarr data is hosted externally (e.g., AWS S3, GCS, Figshare). The repo contains only the React app; the config points to the remote URL.

```r
giotto_to_vitessce_site(
    ...,
    data_mode     = "remote",
    data_base_url = "https://my-bucket.s3.amazonaws.com/data/"
)
```

## Viewer Features (v1)

- Spatial view: cell centroids / polygon shapes
- UMAP scatterplot
- Observation-set panel (cluster labels, cell types, etc.)
- Linked selection: click cells in one view to highlight them in all views

## Architecture Notes

- Stage 1 uses `GiottoClass::giottoToSpatialData()` as the canonical conversion function. Additions needed for future slots (expression matrix, images) should be implemented in GiottoClass, not here.
- Template files in `inst/templates/` use `{{TOKEN}}` placeholders replaced via `gsub()`. glue/sprintf are intentionally avoided because JSX source uses `{expression}` syntax throughout.
- The `SpatialDataWrapper` from the vitessce Python package is used in Stage 2 (not `AnnDataWrapper`) since Stage 1 produces a full SpatialData store.

## Future Extensions

- Expression matrix with heatmap and gene list views
- Morphology images (OME-Zarr tiling)
- Multi-sample support
- Additional interoperability targets beyond Vitessce
