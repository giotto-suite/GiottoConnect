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
2. **Fingerprints flag absent optional parts as drift.** 4 of the 5 real files are
   `repair_needed` (and warn) although all convert correctly. The causes are an
   absent optional part (normalized, nn) and counts stored as float. **Proposal:**
   - treat unbound optional parts as neutral in the fingerprint;
   - compare integer-valued storage by role, not dtype.

   Seurat fingerprints would change, so this needs a separate re-stamp with the
   snapshot proof.
3. **Mapping-level declared losses are listed even when their path is absent**
   (both ecosystems). **Proposal:** keep only those whose path is present.
4. **Feature metadata always comes from `var`**, even when the counts come from
   `raw/X` with their own `raw/var`. In slideseqv2, `raw` has 17,733 genes vs 4,000
   in `var`. That file has no integer `raw`, so it was not hit. **Proposal:** bind
   the feature table from the matrix that was chosen.
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
