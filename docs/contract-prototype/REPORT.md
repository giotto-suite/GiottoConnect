# Mapping-based interoperability: prototype report

GiottoConnect, branch `contract-prototype` (off `dev_connect`, 11 commits, not pushed),
2026-11-20.

**Scope:** Seurat → Giotto, in memory and on GiottoDisk.

**Stack:**
- GiottoClass 0.7.31 (`gsource`)
- GiottoDisk 0.0.0.8 (`dev`)
- Giotto 4.3.6 (`gsource`)
- SeuratObject 5.4.0
- Seurat 5.5.0 and 6.0.0

Phase write-ups: `PHASE0_ORIENTATION.md`, `PHASE2_CHECKPOINT.md`, `PHASE3_READBACK.md`,
`PHASE4_CHECKPOINT.md` and `PHASE5_LLM.md`. Scripts and result tables are in
`examples/seurat_comparison/`.

## Summary

- **The new route converts every mapped part correctly** on the fixtures and the three
  GiottoData crops (Xenium, CosMx, Visium HD), in memory and on disk, under Seurat 5.5 and
  6.0: 35 of 35 parts. The legacy `seuratToGiottoV5` gets 15 of 36 wrong, and one of its
  two copies no longer runs.
- **Drift is mostly absorbed by parameters.** All nine drift variants reach the expected
  outcome without an LLM, and Seurat 6.0.0 objects needed no mapping change.
- **Supporting an object family is a YAML file, not a function.** Seurat v3/v4 support is
  a second mapping that differs from the v5 one in about ten lines.
- **An LLM helps only where choices depend on naming, and can be wrong in ways the checks
  cannot see.** It stays opt-in, selector-only, and behind human review.

## Code size

| | lines (code, excl. comments/blank) |
|---|---|
| Engine, ecosystem-neutral (`R/contract_*.R`) | 2,513 |
| of which: profile 485, repair 423, engine 361, validate 230, LLM 197, readback 181, contract 173, find 170, schema 154, hooks 139 | |
| Seurat mappings (YAML): `seurat-assay5` 130, `seurat-assay` 130, `seurat-visium` 88 | 348 |
| Seurat hooks (graph → NN network, Visium pixels → µm) | ~120 |
| `seuratToGiotto()` wrapper | 20 |
| **Legacy Seurat → Giotto** in GiottoConnect: `seuratToGiottoV5` | 367 (+ `.createGiottoPolygonFromSeurat`) |
| Legacy pair incl. `giottoToSeuratV5` | 753 (V4 pair on top) |

The engine is larger than one converter, but it is written once:
- every later ecosystem adds YAML, hooks and an adapter;
- every later Seurat structure adds YAML only.

The legacy code also covers images, molecules and spatial networks, which the mappings do
not yet.

## Dependencies added

- **Imports:** `yaml`, `digest` (already a GiottoClass import).
- **Suggests:** `GiottoDisk`, `httr2`.
- No Python, no `jsonvalidate`, no `ellmer`.
- New code calls SeuratObject accessors. Seurat is needed only by the Visium hook, which
  reads `VisiumV2` (defined in Seurat).

## Public API (all documented, 269 tests, 1 opt-in live test skipped)

| Area | Functions |
|---|---|
| Converting | `seuratToGiotto()`, `convertWithMapping()` |
| Mappings | `findMapping()`, `readMapping()`, `validateMapping()` |
| Checking | `readbackCheck()`, `mappingProvenance()` |
| Profiling | `profileObject()`, `fingerprint()`, `diffProfiles()` |
| Contract and hooks | `giottoContract()`, `hookRegistry()` |
| Repair | `repairMapping()`, `writeMappingRequest()`, `llmBackend()`, `llmPatch()` |

## Mappings, parts and declared losses

Mappings are chosen by object structure, not by technology. FOV-based loaders (Xenium,
CosMx, MERSCOPE) share one mapping.

| mapping | applies to | parts | alternatives | declared losses |
|---|---|---|---|---|
| `seurat-assay5` | default assay is `Assay5`, FOV image | expression raw/normalized, cell + feature metadata, centroids, polygons, pca, umap, kNN network | locs from `meta.data` columns, locs from a `spatial` reduction | units assumed µm; invalid polygons repaired; `scale.data`, PCA loadings, sNN graph dropped |
| `seurat-visium` | default assay `Assay5`, `VisiumV2` image | same, locs via pixels → µm | | tissue image, `scale.data`, loadings, sNN dropped |
| `seurat-assay` | default assay is a v3/v4 `Assay` | as `seurat-assay5` | as `seurat-assay5` | as `seurat-assay5`, plus `var.features` |

**Recorded per conversion (in the object history):**
- **Network naming:** the network is named `kNN.<reduction>`, with k recovered from the
  FindNeighbors record, and distances are recomputed in the recorded space.
  - Without a record, the name is `ambiguous` and distances are dropped.
- **Visium scale:** the µm/px factor (equal to Space Ranger's `microns_per_pixel` on the
  Visium HD crop).

**Hooks: 2.** `seurat_graph_to_nn` and `seurat_visium_to_um`. Everything else is SeuratObject
accessors plus 10 ops: `select_cols`, `exclude_cols`, `rename`, `transpose`, `scale`,
`cast`, `rownames_to_col`, `set`, `set_colnames`, and `table` as a cast target.

## Correctness

### Against the source (`results/phase2_parts.csv`, `results/seurat6_parts.csv`)

| part | new (memory and disk) | legacy |
|---|---|---|
| expression raw/normalized, cell metadata, pca, umap | correct | correct |
| feature metadata | correct | **wrong**: holds the cell metadata |
| spatial locations | correct (FOV: x, −y; Visium HD: µm) | **wrong**: axes swapped (y, −x) |
| polygons | correct; 2 of 1531 Xenium polygons repaired, declared | **wrong**: no `poly_ID`, not linked to cells |
| nn network | correct: kNN, no self edges, real distances | **wrong**: typed sNN, self edges, no usable weights |

**Readback** (`readbackCheck()`, `results/phase3_readback.csv`): PASS on every dataset and
backend. The only non-MATCH is the 2 declared polygon repairs.

**Independent oracles: GiottoDisk's native importers on the same crops.**
- **Xenium:** identical raw counts on shared genes; centroids within 0.55 µm (GiottoDisk
  computes them from polygons).
- **Visium HD:** identical locations on all 386 bins.

**Real object:** `xe_beforeRCTD_sample_cores.rds` (400,216 cells, Seurat 5.1, no image).
It repairs by rule (locations from its `spatial` reduction), reads back PASS, and
converts in seconds.

## Drift (`results/phase4_drift.csv`)

| | change | outcome | legacy |
|---|---|---|---|
| D1 | assay renamed | params | axes swapped |
| D2 | layer `counts` → `raw_counts` | params | axes swapped (finds the layer by regex) |
| D3 | coordinates moved to `meta.data`, FOV removed | rules: alternative | **error** |
| D4 | counts stored dense | rules: cast | axes swapped |
| D5 | extra reduction | params, WARN unmapped | axes swapped |
| D6 | coordinates in mm | params, WARN looks like mm | no warning |
| D7 | graph renamed | params | axes swapped |
| D8 | unused graph removed | params (fingerprint unchanged) | axes swapped |
| D9 | v3 `Assay` | `seurat-assay` mapping | axes swapped |

Every new-route result passes the readback. Seurat 6.0.0 objects change only in metadata
the mappings do not read: `findMapping` → `match`, readback PASS, identical part results.

## LLM repair (`results/phase5_llm.csv`)

| | situation | baseline (no model) | qwen3.8:27b (local Ollama) | Haiku 4.5 (subagent replay) |
|---|---|---|---|---|
| L1 | counts moved to `misc` | declines | correct | correct |
| L2 | `counts_corrected` vs `umi_raw` | declines | correct | correct |
| L3 | CosMx local vs global coordinate columns | declines | **wrong, accepted** (local) | correct (global) |
| L4 | coordinates in a `spatial` reduction | rules | rules | rules |
| L5 | no correct option | refuses | refuses | refuses (its wrong pick was rejected by the gate) |

- **Per call:** 400–800 input and 120–280 output tokens, 7–18 s on the local 27B model.
- **Haiku 4.5** was measured through Claude Code subagents: exact prompts, answers replayed
  through the gate (`llm_replay.R`). An API run needs `ANTHROPIC_API_KEY`.
- **The L3 safety point:** internally consistent but wrong choices pass validation and
  readback, so `patched-unreviewed` must stay a hard stop before `verified`.

## Speed (single runs, seconds)

| | new memory | new disk | legacy |
|---|---|---|---|
| Xenium crop (1,531 cells) | 0.20 | 0.35 | 7.3 |
| CosMx crop (1,752 × 20,378) | 0.54 | 0.56 | **154** |
| Visium HD crop (386 bins) | 0.11 | 0.28 | 0.37 |

The serialized RSS benchmark and the full Xenium lung run (162k cells) are still to do
(Phase 8).

## Bugs found along the way

**Legacy converters**
- GiottoConnect's `seuratToGiottoV5` no longer runs on GiottoClass ≥ 0.7.30 (calls the
  removed `set_expression_values`). It differs from GiottoClass's copy, and both packages
  export the same names.
- GiottoClass's `seuratToGiottoV5`:
  - axes are swapped for FOV data;
  - feature metadata is the cell metadata;
  - polygons carry no cell IDs;
  - kNN graphs are stored as sNN with self edges;
  - `scale.data` is dropped silently;
  - CosMx takes 150 s because molecules are fetched one gene at a time.
- `giottoToSeuratV5`: writes v3 `Assay` for extra feature types, and passes centroids as
  the segmentation.

**Other packages**
- **GiottoDisk** `importXeniumDisk` drops all-zero genes (`TCL1A` in the crop). Seurat
  keeps them.
- **GiottoDisk** `importVisiumHDDisk` errors (`invalid filename argument`) when given the
  crop root instead of `binned_outputs/`.
- **GiottoDisk** `dev` needs Giotto `gsource` ≥ 611c419 (`VisiumReader`). The install
  order matters.
- **GiottoData:** the cached CosMx crop has 1 FOV / 1,752 cells; its description says
  9 FOVs / 16,766 cells.

**Found and fixed in the prototype (selection)**
- YAML 1.1 booleans (`on:` keys, `yes`/`no`).
- Locale-dependent fingerprint sort.
- S4 dispatch in the profiler loading SeuratObject.
- NA network weights (found by the readback).
- Assay binding that depended on counts.
- A duplicate `cell_ID` column on a real object.

## What felt wrong or too limited in the schema

- **No iteration over sources.** Reduction names (`pca`, `umap`) are fixed, so an extra
  `tsne` is only reported. Xenium's control assays cannot become feature types without
  one part per assay. A part template ("for each reduction / each extra assay") is the
  most needed addition.
- **Two addresses per source node.** `path` (for profile, fingerprint, repair) and `get`
  (for extraction) can disagree. Rules can only follow moves for path-read nodes, and
  patches fall back to structural `path` reads. This is correct, but it is the main
  conceptual cost of reading through accessors.
- **Name rules are regexes.** They miss vendor naming (CosMx `CenterX_global_px`), which
  is exactly where the LLM step is needed and where it can be wrong.
- **The op list grew.** `set_colnames` and the `table` cast were added for reductions.
  Each op addition is a schema change, so keep growth deliberate.
- **Units and orientation are losses, not data.** The contract has no unit or orientation
  field, so "assumed µm, y down" is declared rather than checked. A `units` key on
  spatial parts would let validation compare instead of guess.
- **Some Seurat semantics sit in generic rules.** The `value_of active.assay` binding and
  the listing-path table work, but belong in the ecosystem adapter (see below).

## What AnnData / SpatialExperiment will stress

- **AnnData** (next, per spec §8): the source is a file, not an R object. The adapter has
  to:
  - profile and read `.h5ad` through hdf5r;
  - stream large matrices into GiottoDisk via `GiottoDisk::h5adInput()`;
  - handle cells × features orientation (`transpose`);
  - decode categoricals (including the pre-0.8 `__categories` layout);
  - read `obsm` coordinates and `obsp` graphs.

  `obsp/distances` gives real distances, better than Seurat.
- **SpatialExperiment:** the getters are SummarizedExperiment generics (`assay()`,
  `colData`, `spatialCoords`, `reducedDims`). DelayedArray/HDF5Array assays are currently
  "stop" classes in the profiler, so the adapter must profile them without materialising.
- **Both need the adapter seam** (spec §8): `get`, the listing paths, `open` and
  streaming inputs move out of the engine into a per-ecosystem class. Seurat moves first;
  its readback checks prove no behaviour change.

## Open decisions

1. **Retire the legacy converters**, in GiottoConnect and the GiottoClass duplicates, in a
   separate round with your review. The new route now covers everything they got right.
   They still cover images, molecules and spatial networks, which the mappings do not.
2. **Adapter seam before AnnData,** or AnnData directly on the current engine?
3. **Part templates** (for each reduction / extra assay): add to the schema now or after
   the second ecosystem?
4. **GiottoDisk chunked matrix input** (`matrixColumnInput`): a separate GiottoDisk PR,
   needed before the full-size (162k-cell) cost run.
5. **LLM step:** keep opt-in, measure Haiku (needs a key), and decide whether to add
   domain hints to the prompt (e.g. whole-section vs per-FOV coordinates).
6. **Missing mapping families:** VisiumV1, Slide-seq, Atera (Seurat 6 `LoadAtera`, no crop
   yet).
7. **Housekeeping:**
   - DESCRIPTION `Remotes` still point to `drieslab/`;
   - `GiottoConnect/CLAUDE.md` still describes the old vitessce project;
   - the legacy test files fail on `dev_connect` as well;
   - `profileObject()` could later move into GiottoClass next to the manifest walker.
8. **Remaining phases:** Phase 7 (Giotto → Seurat as a mapping) and Phase 8 (the full
   benchmark, including adaptability fix sizes and RSS).
