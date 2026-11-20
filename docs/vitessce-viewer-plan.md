# giotto-vitessce-viewer

## Project Goal

Build an R package (or set of modular R functions) that converts a Giotto object into a static Vitessce-powered website deployable via GitHub Pages. The pipeline is: **Giotto → AnnData Zarr → Vitessce view config + React app → GitHub Pages repo**.

## Architecture Overview

The workflow has three discrete stages, each implemented as one or more exported R functions:

### Stage 1: Export Giotto to AnnData-Zarr

Convert the relevant slots of a Giotto object into an AnnData-Zarr store that Vitessce can consume.

**Function:** `export_giotto_to_zarr(gobject, output_dir, ...)`

- Use `GiottoClass::giottoToSpatialData()` or the interoperability path documented at <https://giottosuite.com/articles/interoperability_spatialdata.html> to convert to SpatialData/AnnData via reticulate.
- Alternatively, use `vitessceAnalysisR::seurat_to_anndata_zarr()` as a reference pattern — but the preferred path is Giotto → SpatialData (Python spatialdata package) → write zarr, since SpatialData is a first-class Vitessce input.
- **Minimum data to export (v1):**
  - Cell centroids → `obsm/X_spatial` (n × 2 array)
  - Cell polygon geometries → SpatialData shapes element or `obsm/X_segmentations`
  - UMAP embedding → `obsm/X_umap` (n × 2 array)
  - Cell metadata (cluster labels, any obs columns) → `obs/` columns
- **Not yet needed:** expression matrix (`X`), images, multi-resolution OME-NGFF. Design the function signature so these can be toggled on later.
- Run `zarr.consolidate_metadata()` on the output store (Vitessce requires `.zmetadata`).
- Use `vitessce.data_utils.optimize_adata()` if going through AnnData to ensure compatible dtypes.

### Stage 2: Generate Vitessce View Config

Build the JSON view config that tells Vitessce where the data lives and which views to render.

**Function:** `create_vitessce_config(zarr_path, dataset_name, base_url, ...)`

- Use the `vitessce` Python package (via reticulate) or `vitessceR` to construct the config programmatically.
- The `base_url` should be the final GitHub Pages URL where the zarr data will be served, e.g. `https://<username>.github.io/<repo>/data/`.
- **Views to include (v1):**
  - `SCATTERPLOT` with `mapping = "UMAP"`
  - `SPATIAL` for the centroid/polygon spatial view
  - `OBS_SETS` for cluster label selection
- Coordinate the SCATTERPLOT and SPATIAL views so cell selections sync between them.
- Export the config as a JSON file (`vitessce_config.json`).
- Design for future extensibility: expression heatmap, feature list, layer controller, etc.

### Stage 3: Scaffold the GitHub Pages Site

Generate a ready-to-push GitHub repository containing a Vite + React app that embeds Vitessce.

**Function:** `scaffold_vitessce_site(config_json, zarr_dir, output_dir, repo_name, github_username, ...)`

- Create the following directory structure inside `output_dir`:
  ```
  <repo_name>/
  ├── public/
  │   └── data/
  │       └── <dataset>.zarr/    # copied or symlinked from Stage 1 output
  ├── src/
  │   ├── App.jsx                # imports Vitessce, loads config
  │   ├── my-view-config.js      # the exported JSON config
  │   ├── main.jsx
  │   └── index.css              # Vitessce container styling
  ├── .github/
  │   └── workflows/
  │       └── deploy.yml         # GitHub Actions workflow for gh-pages deploy
  ├── vite.config.js             # base set to /<repo_name>/
  ├── package.json               # vitessce, react@18, vite dependencies
  ├── index.html
  └── README.md
  ```
- The GitHub Actions workflow should build with `npm run build` and deploy to GitHub Pages using `actions/upload-pages-artifact` + `actions/deploy-pages`.
- Follow the official Vitessce tutorial pattern: <https://vitessce.io/docs/tutorial-gh-pages/>
- The zarr data goes into `public/data/` so Vite copies it to the build output and it is served alongside the app on GitHub Pages (no external hosting like S3 needed).
- **Important:** Data files in the GitHub repo are served via GitHub Pages, so URLs in the view config must point to `https://<username>.github.io/<repo>/data/<dataset>.zarr` (not raw GitHub URLs).

### Convenience Wrapper

**Function:** `giotto_to_vitessce_site(gobject, output_dir, repo_name, github_username, dataset_name, ...)`

- Orchestrates Stage 1 → 2 → 3 in sequence.
- Prints instructions for the user on how to push to GitHub and enable Pages.

## Development Rules

- **R-first:** Use R packages wherever possible. Only fall back to Python (via reticulate) when the Python tool is clearly better or more advanced — e.g., `spatialdata` for writing SpatialData-Zarr, `vitessce` Python for config generation if `vitessceR` lacks features. Document every Python dependency and why R wasn't sufficient.
- **data.table:** When you need a data.frame-like structure (e.g., extracting cell metadata, building obs tables), use `data.table`, not base data.frames or tibbles.
- **Giotto integration:** Consider which of the new functions belong inside Giotto Suite (GiottoClass or GiottoUtils) vs. staying in this standalone package. Functions that do generic Giotto-object-to-Zarr conversion are strong candidates for upstreaming into GiottoClass. Document this reasoning in code comments with `# TODO(giotto-upstream): ...` where applicable.
- **Two data-hosting modes:** Design all functions to support at least two deployment modes:
  1. **Self-contained (GitHub-hosted):** Zarr data lives in `public/data/` inside the repo and is served via GitHub Pages. Suitable for small datasets (< ~500 MB). This is the default for the prototype.
  2. **External hosting:** Zarr data is hosted on an established external service (e.g., AWS S3, Google Cloud Storage, Figshare, Zenodo, or an institutional server). The view config points to the remote URL. The GitHub Pages site contains only the React app, not the data. Implement this as a `data_mode` parameter (e.g., `"github"` vs `"remote"`) with a `data_base_url` argument for the remote case. The scaffold function should omit copying zarr data into the repo when `data_mode = "remote"`.

## Technical Constraints

- **Language bridge:** R is the primary interface. Python dependencies (spatialdata, vitessce, anndata, zarr) are called via `reticulate`. Ensure a consistent Python environment — recommend using `reticulate::virtualenv_create()` or conda.
- **Node.js:** The scaffold step generates static files but does NOT run `npm install` or `npm run build` from R. The user (or GitHub Actions) handles that. The R function just writes the file tree.
- **File size:** GitHub repos have a soft limit of ~1 GB and individual files at 100 MB. For the prototype with a mini Giotto object this is fine. For larger datasets, use `data_mode = "remote"` (see Development Rules above).
- **Vitessce version:** Use `vitessce@3.x` (latest stable). Schema version `1.0.16`.
- **React version:** React 18 (required for Vitessce compatibility).

## Prototype Test Case

Use a **mini MERFISH Giotto object** for development and testing. This can be created via:
```r
# Example: use the mini seqFISH or MERFISH dataset from GiottoData
# or subset an existing MERFISH object to ~500 cells
library(GiottoData)
gobject <- GiottoData::loadGiottoMini("merfish")
```

The prototype should demonstrate:
1. Spatial view showing cell polygons/centroids colored by cluster
2. UMAP scatterplot colored by the same cluster labels
3. Linked selection between spatial and UMAP views
4. Deployed and functional on GitHub Pages

## Future Extensions (do not implement now, but design for)

- Multi-sample support (sample selector or multi-dataset configs)
- Expression matrix with heatmap and gene list views
- Morphology images (OME-NGFF / OME-Zarr tiling)
- SpatialData as a first-class intermediate (beyond just AnnData)
- Automatic `git init` + `git push` from R
- Custom theming / branding of the viewer site
- Support for additional Giotto data types (e.g., spatial networks, enrichment results)

## Key References

- Giotto → SpatialData interop: <https://giottosuite.com/articles/interoperability_spatialdata.html>
- Vitessce GitHub Pages tutorial: <https://vitessce.io/docs/tutorial-gh-pages/>
- Vitessce data file types: <https://vitessce.io/docs/data-file-types/>
- vitessceR export to static files: <http://r-docs.vitessce.io/articles/export_files.html>
- vitessce-python data export: <http://python-docs.vitessce.io/notebooks/data_export_files.html>
- Vitessce view config API (Python): <https://vitessce.github.io/vitessce-python/api_config.html>
- AnnDataWrapper parameters: <https://vitessce.github.io/vitessce-python/api_data.html>
- GitHub discussion on serving data from the repo itself: <https://github.com/vitessce/vitessce/discussions/2162>

## Repository Structure (this project)

```
giotto-vitessce-viewer/
├── CLAUDE.md              # this file
├── R/
│   ├── export_zarr.R      # Stage 1
│   ├── create_config.R    # Stage 2
│   ├── scaffold_site.R    # Stage 3
│   └── giotto_to_site.R   # convenience wrapper
├── inst/
│   └── templates/         # template files for the Vite/React app
│       ├── App.jsx
│       ├── main.jsx
│       ├── index.css
│       ├── index.html
│       ├── vite.config.js
│       ├── package.json
│       └── deploy.yml
├── tests/
│   └── test_mini_merfish.R
└── README.md
```
