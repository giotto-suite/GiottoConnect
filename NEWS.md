# GiottoConnect 0.1.0

## Converters
- `giottoToAnnData()` writes one spatial unit and feature type of a giotto
  object (in memory or GiottoDisk-backed) as an `.h5ad` in scanpy's layout,
  without python (anndataR): normalized values in `X`, counts in
  `layers/counts`, other values as layers, reductions in `obsm` (PCA loadings in
  `varm/PCs`), nn and spatial networks in `obsp` with a `uns` record.
  `anndataToGiotto()` reads them back under their Giotto names, and also reads
  every scanpy/squidpy network and extra layer of other files.
- `giottoToSeurat()` takes `spat_unit` / `feat_type` and exports every
  reduction (with PCA loadings and stdev), every nearest-neighbour network as a
  Graph, and transcripts as FOV molecules.
- `giottoToSpatialExperiment()` converts one spatial unit and feature type to
  a SpatialExperiment, or a SpatialFeatureExperiment when there are polygons,
  transcripts or spatial networks: `counts` / `logcounts` and other assays,
  scater-named reductions with PCA loadings, networks as colPairs (and SFE
  spatialGraphs), polygons as `cellSeg`, transcripts as `txSpots`.
  `spatialExperimentToGiotto()` reads them back under their Giotto names.
- Tissue images go both ways. Seurat Visium images, SpatialExperiment
  `imgData` and scanpy `uns/spatial` images become Giotto images placed in the
  spot coordinates; `giottoToSeurat()` writes a VisiumV2 (so
  `SpatialFeaturePlot()` shows the tissue), `giottoToSpatialExperiment()` an
  `imgData` image and `giottoToAnnData()` a scanpy `uns/spatial` library.
- `seuratToGiotto()` reads every nearest-neighbour Graph, not only the first.
- `spatialExperimentToGiotto()` converts SpatialExperiment and
  SpatialFeatureExperiment objects through the same mapping engine: the layouts
  of SpatialExperimentIO (Xenium, CosMx, MERSCOPE), VisiumIO (Visium, Visium HD;
  pixels converted to microns) and SpatialFeatureExperiment (cell and nucleus
  segmentations, transcripts, spatial graphs). Control probes, as altExps or
  inside the main matrix, become their own feature types; scater PCA loadings
  and variances are kept.
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

## External tools
- `runSparkX()` runs SPARK-X (SPARK package, GitHub) for spatially variable
  genes on counts and stores its combined and Benjamini-Yekutieli adjusted
  p-values as feature metadata (`sparkxParam()` for `analyzeData()`). Use it
  in place of Giotto's `spark()`, which runs the original SPARK model.
- Banksy, BayesSpace, CellChat, nnSVG and ONTraC declare their papers
  (`GiottoClass::methodCitation()`): the citation is shown once per session
  when the tool runs and `GiottoClass::objCitations()` lists it.
  `checkGiottoMethod()` checks that a param declares one (`cite = FALSE` for
  Suite-native methods that wrap no published method).
- `runMAST()` runs MAST for marker features between clusters, either one
  comparison (`group_1` against `group_2`) or each cluster against the
  others, and returns MAST's hurdle test per feature (`mastParam()` for
  `analyzeData()`). It replaces Giotto's `findMastMarkers()` and
  `findMastMarkers_one_vs_all()` with the same results. MAST is in Suggests.

## Dependencies
- reticulate moves to Suggests: only the Bento and SpatialData converters and
  the python tools (ONTraC) use it, and they check for it when they run.
- Requires GiottoClass 0.8.0 (`methodCitation()`).

## Removed
- `gefToGiotto()`: Giotto's `importStereoSeq()` (and GiottoDisk's
  `importStereoSeqDisk()`) read Stereo-seq GEF files. rhdf5 is no longer
  suggested.
- `seuratToGiottoV4()`, `seuratToGiottoV5()`, `giottoToSeuratV4()`,
  `giottoToSeuratV5()` (use `seuratToGiotto()` / `giottoToSeurat()`).
- `giottoToSpatialData()`, and the reticulate-based `giottoToAnnData()`
  (replaced by the one above); `anndataToGiotto()` no longer needs python.
- The Vitessce site pipeline (`export_giotto_to_zarr()`, `create_vitessce_config()`,
  `scaffold_vitessce_site()`, `giotto_to_vitessce_site()`,
  `setup_giotto_connect_env()`).
