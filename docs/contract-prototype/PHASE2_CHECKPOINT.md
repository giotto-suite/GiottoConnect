# Checkpoint 2: first mapping and engine

2026-11-19, branch `contract-prototype`. Stack as in PHASE0_ORIENTATION.md.

## The mapping

- **File**: [`inst/mappings/seurat-assay5.yaml`](../../inst/mappings/seurat-assay5.yaml),
  106 non-comment lines, `status: verified`, fingerprint `b529eda4ec3f9d20`.
- **Stored profile**: `seurat-assay5.profile.json`, kept next to it for drift detection.
- **Stamping**: `inst/scripts/verify_mapping.R` converts the fixture in memory and on disk,
  requires validation with no FAIL, then writes the fingerprint, status and date.

`summary(readMapping("seurat-assay5"))`, as recorded in
`examples/seurat_comparison/results/phase2_validation.txt`:

| part | from | via | to |
|---|---|---|---|
| expr_raw | `LayerData(assay, layer_raw)` | | expression cell/rna/raw |
| expr_normalized | `LayerData(assay, layer_norm)` | | expression cell/rna/normalized |
| cell_meta | `x[[]]` | rownames_to_col | cell_metadata cell/rna |
| feat_meta | `x[[assay]][[]]` | rownames_to_col | feat_metadata cell/rna |
| locs | `GetTissueCoordinates(centroids)` | rename, scale y × −1 | spatial_locs cell/raw |
| polygons | `GetTissueCoordinates(segmentation)` | rename, scale y × −1 | spatial_info cell |
| pca, umap | `Embeddings()` | | dim_reduction |
| nn | `x[[nn_graph]]` | hook `seurat_graph_to_nn` | nn_network kNN / kNN.pca |

**Parameters** are bound on the profile with no SeuratObject calls: `assay` (first
Assay5 with integer counts), `layer_raw` (first integer non-negative layer),
`layer_norm`, `fov`, `centroids`/`segmentation` (by boundary class, so Xenium's
`segmentations` binds too), and `nn_graph` (`_nn$` or binary-valued).

**Declared losses (7)**:

| part | kind | what | note |
|---|---|---|---|
| locs | estimated | units | assumed micron |
| polygons | estimated | geometry | invalid source rings are repaired by `make_valid` |
| nn | dropped | distance | Seurat graphs store no distances |
| nn | ambiguous | name | the graph's reduction is not recorded |
| mapping | dropped | `scale.data` | |
| mapping | dropped | PCA loadings and stdev | |
| mapping | dropped | the `_snn` graph | |

Per conversion, the hook also records `k` (recovered from out-degree).

**Hooks: 1** (`seurat_graph_to_nn`).

## Validation output (Xenium crop)

Structural checks: 7 params bound, every source path found, every target valid. They
raise **6 WARN for unmapped data**: the BlankCodeword, ControlCodeword and ControlProbe
assays (counts + data).

Semantic checks on the result: **10 OK, 0 WARN**.
- raw counts are integer and non-negative
- 1531 cells and 377 features, as in the source
- 100% cell-ID agreement across metadata, locs and polygons
- coordinate range 424 µm
- pca/umap have 1531 rows
- kNN: 29,089 edges, all over cells

`findMapping()` returns `match` for the fixture, for the fixture with the assay renamed
to `Spatial`, and for the real Xenium crop.

## Comparison with `seuratToGiottoV5`

Each part is checked against the Seurat source read through SeuratObject accessors
(`examples/seurat_comparison/compare_parts.R`; results in `results/phase2_parts.csv`).
The legacy converter is GiottoClass's copy, because the GiottoConnect copy no longer runs.

| part | new (memory) | new (GiottoDisk) | legacy |
|---|---|---|---|
| expression raw / normalized | correct ×4 | correct ×4 | correct ×4 |
| cell metadata | correct ×4 | correct ×4 | correct ×4 |
| feature metadata | correct ×4 | correct ×4 | **wrong ×4**: holds the cell metadata |
| spatial locs | correct ×3, Visium HD not mapped | same | **wrong ×4**: axes swapped (`y, −x`) |
| polygons | correct ×3 | correct ×3 | **wrong ×3**: no `poly_ID`, so no link to cells |
| pca, umap | correct ×4 | correct ×4 | correct ×4 |
| nn network | correct ×4 (kNN, self edges dropped) | correct ×4 | **wrong ×4**: typed sNN, includes self edges, named `<assay>_nn` |

The 4 datasets are the fixture, the Xenium crop, the CosMx crop and the Visium HD crop.
On Xenium, the new mapping's polygons are correct but 2 of the 1531 cells differ: they
are the cells whose source rings are invalid, repaired as declared.

**Visium HD.** Its image is `VisiumV2`, not FOV, so `locs` is unbound and skipped (the
part is optional). This is the expected point for a second mapping family. The legacy
converter produces swapped coordinates for it.

**Time to convert** (s, single run):

| dataset | new memory | new disk | legacy |
|---|---|---|---|
| fixture | 0.51 | 1.69 | 0.83 |
| Xenium | 0.19 | 0.74 | 8.17 |
| CosMx | 0.52 | 0.53 | **149.4** |
| Visium HD | 0.08 | 0.30 | 0.74 |

Legacy CosMx is slow because it fetches molecules one gene at a time. The new mapping
does not map molecules yet.

## Independent oracle: GiottoDisk `importXeniumDisk()` on the same crop

| | importXeniumDisk | seuratToGiotto(LoadXenium()) |
|---|---|---|
| cells | 1531 | 1531 (all shared) |
| features | 376 | 377 |
| raw counts, shared cells × genes | | identical |
| centroids | | max diff 0.28 µm (x), 0.55 µm (y) |
| polygons | 1531 | 1531 |

- The extra gene in the Seurat route is `TCL1A`, which has zero counts in the crop.
  GiottoDisk's importer drops all-zero genes; Seurat keeps them.
- The centroid offsets are GiottoDisk computing centroids from polygons, against
  Xenium's own cell centroids.

## Code size

- **Engine (generic, reused for every ecosystem)**: 1,542 code lines in `R/contract_*.R`.
- **Seurat-specific code**:
  - the 106-line YAML
  - the 20-line `seuratToGiotto()` wrapper
  - the hook (~35 lines)
- **Legacy Seurat pair in GiottoConnect**: 753 lines.

The legacy code also covers images, molecules and spatial networks, which the mapping
does not yet.

## Bugs found and fixed while building this

| Area | Bug | Fix |
|---|---|---|
| YAML | `on:` keys and `y` values parsed as booleans (YAML 1.1) | only `true`/`false` are booleans |
| Fingerprint | the hash depended on `LC_COLLATE` | radix sort |
| Contract | constructor data argument was assumed to be first (wrong for `createNearestNetObj`) | now declared per part |
| Validation | GiottoDisk refuses to materialise >100 × >100 | sample < 100 cells |
| Validation | used the deprecated `getPolygonInfo(polygon_name =)` | switched to `name =` |

## Open points

- Control assays (Xenium) are WARN-only. GiottoDisk maps them as extra feat_types; the
  mapping could do the same with one part per extra assay.
- **VisiumV2 needs a second mapping** (`seurat-visium.yaml`): coordinates via
  `GetTissueCoordinates` on the image, plus the scale factors.
- **Not mapped yet**: molecules, images, spatial networks, PCA loadings.

## Decisions after the checkpoint (Ruben, 2026-11-19)

1. **Visium gets its own mapping.** `inst/mappings/seurat-visium.yaml` (verified, fingerprint
   `a18c0d2cb8a248fa`) covers `VisiumV2` images. Mappings follow the Seurat object
   structure, not the technology: FOV-based loaders (Xenium, CosMx, MERSCOPE) share
   `seurat-assay5`.
   - **New hook `seurat_visium_to_um`**: pixels → microns. Visium HD uses bin size / spot
     scale factor; standard Visium uses a 100 µm spot pitch over the measured spacing.
   - **Visium HD crop**: `findMapping` selects it (`match`, while `seurat-assay5` would
     need repair). The scale is 0.27370 µm/px, exactly Space Ranger's
     `microns_per_pixel`. Locations are identical to `GiottoDisk::importVisiumHDDisk(bin =
     16, micron = TRUE)` on all 386 bins (dx = dy = 0).
2. **Xenium control assays stay WARN-only.**
3. **Networks are named `<kNN|sNN>.<reduction>` only when Seurat's `FindNeighbors` record
   lists the graph.** The record also gives k (`k.param` − 1, since Seurat counts the cell
   itself). Without a record, the Seurat graph name is kept and an `ambiguous` loss is
   recorded. On the real crops this gives `kNN.pca` with k = 19.

Hooks: **2**.
