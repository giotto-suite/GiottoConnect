# GiottoConnect 0.1.0

## Converters
- `seuratToGiotto()` and `anndataToGiotto()` convert with a declarative mapping
  chosen for the object's structure, not per object version. Bundled mappings
  cover Seurat v3/v4 and v5 assays with FOV, Visium (V1, V2, HD) or Slide-seq
  images, and AnnData `.h5ad` files from scanpy/squidpy (read with hdf5r, no
  python).
- Converted: expression (every assay or `feature_types` category as its own
  feature type), cell and feature metadata, spatial locations, polygons,
  transcripts (Seurat FOV molecules), dimension reductions with PCA loadings and
  eigenvalues, nearest-neighbour networks and squidpy spatial networks.
- `source =` builds a GiottoDisk-backed object; AnnData matrices stream from the
  file.
- `mappingProvenance()` records the mapping, bound parameters, losses, absent
  parts and fallback sources used; `readbackCheck()` compares the result with
  the source.
- Lower-level tools: `findMapping()`, `convertWithMapping()`,
  `validateMapping()`, `repairMapping()`, `profileObject()`, `giottoContract()`.

## Removed
- `seuratToGiottoV4()`, `seuratToGiottoV5()`, `giottoToSeuratV4()`,
  `giottoToSeuratV5()` (use `seuratToGiotto()` / `giottoToSeurat()`).
- `giottoToAnnData()` and `giottoToSpatialData()`; `anndataToGiotto()` no longer
  needs python.
