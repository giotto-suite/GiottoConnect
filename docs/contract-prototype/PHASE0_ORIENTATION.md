# Phase 0: orientation

2026-11-19, branch `contract-prototype` (off `dev_connect`).

**Stack.** Scratch library `rlibs/suite_dev`:
- GiottoClass 0.7.30 (`gsource` 9122ab7)
- GiottoDisk 0.0.0.8 (`dev` 99d2c8b)
- Giotto 4.3.6
- sedonadb 0.5

From the user library: SeuratObject 5.4.0, Seurat 5.5.0, duckdb 1.5.1.

**Probes.** Scripts are in `examples/seurat_comparison/phase0/`.

## 1. Legacy Seurat → Giotto path

| Copy | Function | Lines | Status on GiottoClass 0.7.30 |
|---|---|---|---|
| GiottoConnect `R/seurat.R` | `seuratToGiottoV5` | 367 | **Fails**: `could not find function "set_expression_values"` (removed from GiottoClass) |
| GiottoConnect `R/seurat.R` | `giottoToSeuratV5` | 386 | not run |
| GiottoClass `R/interoperability.R` | `seuratToGiottoV4/V5`, `giottoToSeuratV4/V5`, `.createGiottoPolygonFromSeurat` | (4,359-line file) | V5 runs; it is a different copy from GiottoConnect's |

Both packages export the same names, so whichever loads last masks the other. That the
GiottoConnect copy broke on a GiottoClass internal rename is the maintenance problem the
prototype addresses.

Output of the GiottoClass copy on the synthetic fixture (`p0_legacy_check.R`):

| Part | Result |
|---|---|
| raw, normalized expression | correct |
| PCA / UMAP | correct |
| polygons | 200 polygons, 1,000 vertices: correct |
| spatial locs | **wrong**: `sdimx = y`, `sdimy = -x`. This applies the VisiumV1 row/col convention to FOV data, so the result is rotated 90° |
| feature metadata | **wrong**: the cell metadata is cbind'd positionally onto features (GiottoClass warns "falling back to positional cbind") |
| NN network | stored, but `nn_type` defaults to sNN for a binary kNN graph |
| `scale.data` | read and dropped silently |

## 2. SeuratObject 5.4.0 structure (`p0_seurat_structure.R`)

- **`Seurat`** slots: `assays, meta.data, active.assay, active.ident, graphs, neighbors,
  reductions, images, project.name, misc, version, commands, tools`.
  - `x@version` gives the package version.
  - `attributes()` exposes every slot without SeuratObject loaded.
- **`Assay5`** slots: `layers, cells, features, default, assay.orig, meta.data, misc, key`.
  - **Layer matrices have NULL dimnames.** Names live in `cells`/`features`, which are
    `LogMap` objects (a logical matrix of name × layer).
  - `LayerData()` returns named matrices, so **the source side must read through accessors,
    not slot paths**.
- **`FOV`** slots: `molecules, boundaries, coords_x_orientation, assay, key, misc`.
  `Boundaries()` lists the boundary names.
- **`Centroids`**: a `coords` matrix plus `cells`, `nsides`, `radius`, `theta`, `bbox`.
- **`Segmentation`** extends sp `SpatialPolygons` and also carries an `sf.data` slot.
- **`GetTissueCoordinates()`** returns `data.frame(x, y, cell)` for centroids (`which = "centroids"`)
  and for a segmentation (closed rings, 5 rows per square). **No sp/sf needed.**
- **`Graph`** extends `dgCMatrix` and adds `assay.used`, which is empty unless set.
  - The `_nn` graph is binary with nnz/n = k.
  - `_snn` is weighted, so k cannot be recovered from it.
- **`DimReduc`** slots: `cell.embeddings, feature.loadings, feature.loadings.projected,
  assay.used, global, stdev, jackstraw, misc, key`.
- **Listing verbs** (`Assays`, `Layers`, `Reductions`, `Graphs`, `Neighbors`, `Images`,
  `Boundaries`) all work. They are what param binding needs.

## 3. Seurat loaders on GiottoData crops (`p0_loaders.R`, `p0_xen2.R`)

| Crop | Loader | Cells × feats | Structure notes |
|---|---|---|---|
| `xenium_mini_lung_ffpe` | `LoadXenium(segmentations = "cell")` | 1,531 × 377 | **4 assays** (Xenium + BlankCodeword, ControlCodeword, ControlProbe). Boundaries are `centroids` + **`segmentations`** (plural). Without `segmentations=` there are centroids only. Reads the h5-only crop fine. |
| `cosmx_mini_colon_wtx` | `LoadNanostring` | 1,752 × 20,378 | Boundaries are `centroids` + **`segmentation`** (singular) |
| `visiumhd_mini_crc` | `Load10X_Spatial(bin.size = 16)` | 386 × 18,132 | Assay `Spatial.016um`. Image is **`VisiumV2`** (not FOV) with `centroids` only |

So boundary names, assay names and image classes vary by loader. Boundary and assay names
are params chosen by class/role; image class selects the mapping family.

**Crop discrepancy.** The GiottoData description for the CosMx crop says 9 FOVs and
16,766 cells, but the cached files contain **one FOV (71) and 1,752 cells**: the metadata
and `S0_fov_positions_file` (one row) agree with each other. This is a GiottoData follow-up.

## 4. Orientation and the independent oracle (`p0_orient.R`)

Example cell `afihnpcm-1` in the Xenium crop:

| Reader | x | y |
|---|---|---|
| `cells.parquet` | 19.260 | 364.696 |
| `LoadXenium` (`flip.xy = FALSE`, the default) | 19.260 | 364.696 |
| GiottoDisk `importXeniumDisk` | 19.169 | **−364.595** |
| legacy `seuratToGiottoV5` | y | −x (swapped) |

- **Giotto convention is `sdimx = x`, `sdimy = −y`.** The mapping expresses this as a
  `scale` op (factor −1 on y), and the legacy axis swap is a bug.
- The ~0.1 µm offset in GiottoDisk's centroids (probably computed from polygons) will be a
  tolerance in the oracle comparison, not a loss.
- Compared with `LoadXenium`, GiottoDisk's importer also produces:
  - nucleus polygons and locs
  - control feature types (`NegControlProbe`, `NegControlCodeword`, `UnassignedCodeword`)
  - 376 rna features against Seurat's 377 (to explain in Phase 2)

## 5. Assembly with existing tools, in memory and on disk (`p0_assemble.R`)

The fixture parts were built with the constructors (`createExprObj` ×2, `createCellMetaObj`,
`createSpatLocsObj`, `createGiottoPolygonsFromDfr`, `createDimObj` ×2,
`createNearestNetObj`) and assembled two ways:
- **Memory**: `setGiotto(giotto(), parts)`, 0.38 s.
- **Disk**: `setGiotto(giotto(source = sourceCreate(dir)), parts, initialize = FALSE)` then
  `initialize()`, 0.31 s. This is the same pattern as GiottoDisk's own
  `sourceWrite(gDirSource, giotto)`.

Results:
- **Disk-backed**: expression → `parquetExprStore`, polygons → `parquetGeomStore`,
  NN network → `parquetEdgeStore`.
- **In memory, as expected**: spatial locs, cell metadata, dim reductions.
- **Disk readback**: `as(store, "dgCMatrix")` equals the source counts exactly. Subsetting a
  store is lazy.
- **`as.list` → `setGiotto(giotto(), …)` round trip**: `manifestDiff` reports only "object
  identity changed" (the uid), so content is identical.
- **Manifest memory vs disk**: `dtype` reports the storage class (`dgCMatrix` vs
  `parquetExprStore`), `directed` is NA for edge stores, and disk fingerprints hash
  handles. So **`manifestDiff` flags backend changes as modifications**. The readback check
  must compare content (via store subset → matrix), not manifest fields.

## 6. Corrections to the original spec (all applied in v2)

1. The converters exist in two copies (GiottoConnect + GiottoClass); one of them no longer runs.
2. GiottoClass is already in Imports; checkmate and jsonlite are already imported; rlang
   is not; digest is a GiottoClass import.
3. Seurat v5 layer names are not reachable by slot path, so accessors are needed.
4. Centroids and segmentation need no hooks (`GetTissueCoordinates`), so one hook remains.
5. `as.list`/`setGiotto` round-trip.
6. Giotto has no slot → class → nesting table. The contract must generate one.
7. The GiottoClass manifest + `manifestDiff` exist and provide the target vocabulary. They
   are not backend-neutral for dtype/fingerprint.

## 7. Environment notes

- GiottoDisk `dev` head (f0f431f) does not install in the scratch library: its
  `VisiumDiskReader` needs a `VisiumReader` class from a newer Giotto. The installed 99d2c8b
  covers the Xenium, CosMx and Visium HD importers used here.
- A stale `.git/index.lock` (empty, dated Aug 22, no git process running) was removed from
  GiottoConnect to create the branch.

## Next (Phase 1)

- Profiler, fingerprint, diff, schema and validator, contract generation.
- The fixture builder in `inst/scripts/build_fixtures.R`, based on `p0_seurat_structure.R`.
- The legacy oracle snapshot is regenerated from that builder.
