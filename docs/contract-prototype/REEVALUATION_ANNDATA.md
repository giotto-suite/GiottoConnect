# Re-evaluation after the second ecosystem (AnnData)

## What AnnData needed

| Piece | Size | Notes |
|---|---|---|
| `R/contract_adapter_anndata.R` | ~400 lines | profile walker over HDF5 + six element accessors (ported from GiottoDisk `.h5ad_*`, plus pre-0.8 categoricals) |
| `inst/mappings/anndata-spatial.yaml` | ~120 lines | one mapping for scanpy/squidpy layouts, anndata 0.7 to current |
| hook `anndata_distances_to_nn` | ~50 lines | obsp distances → kNN edges; name from the scanpy neighbors record |
| `anndataToGiotto()` | ~20 lines | same shape as `seuratToGiotto()`; replaces the reticulate converter |
| engine | 16 lines | `as: path`, `where.not_role`, `where.kind` lists, `for_each.strip`; all generic |
| schema / checks | small | listing verbs checked against the adapter instead of a Seurat/Giotto enum |

Targets, validation, readback and repair did not change (readback now checks the
adapter's package). `engine_snapshot.R`: Seurat/Giotto results are **identical to
the baseline on 17/17 objects**. Tests: 416 run. All contract tests pass. The
7 failures (SpatialExperiment and sf/terra/stars converters) already fail on HEAD.

## Real files (`examples/anndata_comparison/results/real_files.csv`)

| File | Cells × genes | Parts mapped | Readback (memory) | Disk readback | Memory / disk |
|---|---|---|---|---|---|
| squidpy seqfish (anndata 0.7) | 19,416 × 351 | counts, meta, locs, umap | 16 MATCH | pass | 3.1 s / 2.5 s |
| squidpy merfish | 73,655 × 161 | normalized, meta, locs | 12 MATCH | pass | 1.8 s / 2.1 s |
| squidpy slideseqv2 | 41,786 × 4,000 | normalized, meta, locs, pca, umap, kNN | 24 MATCH | pass | 3.5 s / 4.1 s |
| Kuppe Visium P12 | 3,392 × 14,514 | counts (raw/X), normalized, meta, locs, pca, umap | 23 MATCH | pass | 2.9 s / 3.5 s |
| AKI seqFISH (537 MB) | 225,000 × 988 | counts (layers), normalized, meta, locs, pca, umap, kNN.pca | 27 MATCH | pass | 10.6 s / 23.1 s |

- There are no DECLARED or UNDECLARED readback differences.
- Profiling takes 0.4–0.7 s per file and reads no matrix data.
- Repair also runs through the adapter. On Kuppe it gave 2 rule fixes and 23 MATCH.

## What held
- **The adapter seam:** the input does not have to be an R object. A file path
  plus an adapter that profiles and reads it was enough.
- **The node table and listing verbs** describe HDF5 just as well as S4 slots.
- **One mapping covers** anndata 0.7 to current, CSR/CSC/dense matrices, counts in
  `X`, `layers/*` or `raw/X`, and files without counts.
- **The package boundary is settled.** GiottoConnect decodes AnnData elements and
  assembles the object. GiottoDisk keeps only the streaming inputs `h5adInput()` and
  `h5adEdgeInput()`; its temporary loaders and `importH5adDisk()` are removed in
  GiottoDisk PR #172. The Stereo-seq readers read GEF, not
  h5ad, so they don't overlap. The removed loaders are a reference for parts not
  mapped yet: PCA loadings and eigenvalues, 3-D coordinates, and squidpy
  spatial-network names (GiottoDisk 7be5447, `R/h5ad.R`).

## Weak points and proposals

1. **Disk conversion did not stream: fixed.** The engine read every part into R, and
   GiottoDisk's writer for in-memory matrices works in one piece, outside the memory
   budget. Now, when the target is a GiottoDisk project, an adapter can return a
   streaming input for expression parts (`adapterStream()`; `h5adInput()` for
   AnnData). The Giotto target writes it straight into the project. Other ecosystems
   are unaffected: the snapshot is identical on 17/17. AKI to disk, both routes
   with `giottodisk.chunk_ram_frac = 0.02`:

   | Route | Peak RSS | Memory footprint | R heap during build | Time |
   |---|---|---|---|---|
   | read into R (before) | 16.9 GB | 13.3 GB | 13.1 GB | 29 s |
   | stream (now) | 7.8 GB | 3.1 GB | 3.5 GB | 32 s |

   At the default budget (25% of free RAM) on this 128 GB machine, both routes
   peak around 16 GB, because GiottoDisk sizes its batches from free RAM.
   Still open: data already in R (Seurat, SCE, BPCells) needs the deferred chunked
   writer in GiottoDisk; obsp networks are still read into R (6.75M edges, ~0.25 GB).
2. **Fingerprints flagged absent optional parts as drift: fixed.**
   - Fingerprints now cover only a mapping's required parts and predicates, and
     integer vs float storage counts as the same type.
   - Conversion switches by itself to a declared alternative when a part's source
     is absent. Drift case D3 (coordinates moved to `meta.data`) now converts with
     its locations; before, it lost them.
   - Absent optional parts and the alternatives used are recorded in the
     provenance (`absent`, `alternatives`). A spatial mapping without any
     coordinates warns.
   - Mapping choice ranks by coverage of all source paths, optional ones included,
     so a family mapping beats the generic one.
   - Result: all 5 real h5ad files match, and the Seurat snapshot differs only as
     intended (D3 gains locations; 4 objects go from `repair_needed` to `match`;
     no mapping choice changes).
3. **Declared losses were listed when their path was absent: fixed.** Only
   drops of data the object holds are listed.
4. **Feature metadata always came from `var`: fixed.** It now follows the counts.
   For `raw/X` counts it is `raw/var` (all its genes), with the `var` columns joined
   by feature ID (`NA` for raw-only genes); otherwise it is `var`.
   - A new `adapterReads()` names the extra paths an accessor reads, so `raw/var`
     no longer shows as unmapped.
   - Accessor args that refer only to an unbound optional param are dropped, so a
     file without counts still gets its feature table.
5. **Not mapped yet:**
   - squidpy spatial networks (`obsp/spatial_*`);
   - obsm data frames (deconvolution);
   - 3-D coordinates;
   - `varm` loadings;
   - Visium scale factors / images in `uns/spatial`.

   Spatial networks need a contract part, which Seurat would share.
6. **Units:** merfish and Kuppe store coordinates rescaled to [0, 1].
   Validation flags them, but the units cannot be recovered from these files.
7. **Reverse direction:** `giottoToAnnData()` is still the reticulate version. A
   Giotto → h5ad target (hdf5r writer) would follow the Seurat target pattern.

## Follow-up: modalities
- **Seurat:** every assay other than the main one becomes its own feature type
  (counts, data and feature metadata; Assay5 and v3). The Xenium/Atera controls
  take the names Giotto's readers use (`ControlProbe` → `NegControlProbe`,
  `ControlCodeword` → `NegControlCodeword`, `BlankCodeword` → `UnassignedCodeword`),
  and `ADT` becomes `protein`. Real Xenium and Atera objects gain their control
  feature types; all existing content is byte-identical.
- **AnnData:** the categories of `var/feature_types` (scanpy `read_10x_h5`) split the
  same way: "Gene Expression" is `rna`, and the other categories get the same Giotto
  names. None of the real files has `feature_types`, so this is tested on a fixture.
  Feature-type subsets are read into R, not streamed.
- Found on the way: GiottoDisk could not read back a store written from an
  all-zero matrix (the Atera crop's `BlankCodeword`/`ControlCodeword`). The fix is
  on GiottoDisk branch `fix/empty-store-write` (local, not pushed).

