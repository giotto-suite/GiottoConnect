# GiottoConnect: mapping-based interoperability (prototype spec v2)

Seurat ↔ Giotto, built as a declarative mapping engine on top of existing GiottoClass,
GiottoDisk and SeuratObject tools. Work phase by phase and stop at every **CHECKPOINT**.

v2 amends the original CONTRACT_PROTOTYPE.md (2026-11-19). The main changes:
- reuse existing tools instead of new parallel machinery (§1)
- read sources through package accessors instead of raw slot paths
- on-disk Giotto (GiottoDisk) as the primary target
- unversioned user commands
- a benchmark against the legacy converters

---

## 0. Goal

The legacy converters (`seuratToGiottoV4/V5`, `giottoToSeuratV4/V5`) are hand-written
field remapping, one function per ecosystem and version. They break when upstream changes
and they lose data silently. They are replaced by:

1. **Mapping files** (YAML, validated against `inst/schema/mapping.schema.json`) that say
   how a family of source objects maps onto Giotto. They are deterministic, with **no LLM
   at runtime**.
2. **A Giotto contract**: the GiottoClass manifest vocabulary plus a generated
   part → constructor → nesting-keys table.
3. **Profiler + fingerprint**: a structural description of any object without its package
   loaded, used to choose a mapping and to detect drift.
4. **Repair**: parameter binding first, then rules, then declared alternatives, then
   optionally an LLM restricted to enumerated candidates.

User commands are per ecosystem and **never carry version names**: `seuratToGiotto()`,
later `giottoToSeurat()`. Object versions are handled by choosing a mapping file.

Backends: `seuratToGiotto(x, source = NULL)` gives an in-memory object;
`seuratToGiotto(x, source = "<dir>")` gives a GiottoDisk project (the future default).

### Success criteria
- A Seurat v5 spatial object (expression raw + normalized, metadata, centroids,
  segmentation, PCA/UMAP, NN graph) converts to Giotto, both in memory and on disk, using
  only mapping + engine.
- The readback check (§4.10) passes with no UNDECLARED differences on the synthetic fixture
  and the GiottoData minis.
- Drift variants D1–D9 behave as in §6.
- The benchmark (Phase 8) shows the new approach correct, more flexible under drift, and
  cheaper to adapt than the legacy converters.

### Non-goals (prototype)
- Do not modify GiottoClass or GiottoDisk. Needed changes there are proposed as separate
  PRs (e.g. a chunked matrix input for GiottoDisk).
- Do not remove the legacy converters yet. They are the benchmark baseline. Removal
  (GiottoConnect and the GiottoClass duplicates) is a separate round after Phase 8.
- No other ecosystems (SCE, SpatialExperiment, AnnData, SpatialData) yet.
- No LLM call anywhere in the conversion path.

---

## 1. Reuse map

| Need | Existing tool | Where |
|---|---|---|
| Read Seurat data (named matrices, FOV coords) | `LayerData`, `Embeddings`, `Loadings`, `Stdev`, `GetTissueCoordinates`, `[[` | SeuratObject (exported) |
| List Seurat names for param binding | `Assays`, `Layers`, `Reductions`, `Graphs`, `Neighbors`, `Images`, `Boundaries` | SeuratObject |
| Build Giotto parts | `createExprObj`, `createCellMetaObj`, `createFeatMetaObj`, `createSpatLocsObj`, `createDimObj`, `createNearestNetObj`, `createGiottoPolygonsFromDfr` | GiottoClass `R/create.R` |
| Assemble | `setGiotto(giotto(), <list of subobjects>)` (initialises once) | GiottoClass `R/methods-setGiotto.R` |
| Disk backend | `gobject@source` + backend-aware setters → `GiottoDisk::sourceWrite` | GiottoClass `gsource`, GiottoDisk |
| Describe / diff Giotto objects | `objManifest`, `objManifest_json`, `manifestDiff`, `inst/schema/giotto-manifest-0.1.0.json` | GiottoClass |
| Hash | `digest::digest(algo = "xxhash64")` (same as manifest) | digest |
| Package checks, messages | `package_check`, `wrap_txt`, `vmsg`, `%null%` | GiottoUtils |
| Independent oracle on vendor data | `importXeniumDisk`, `importCosMxDisk`, `importVisiumHDDisk` | GiottoDisk |
| Test data | `getSpatialDataset()` vendor-layout crops, `loadGiottoMini()` | GiottoData |

New code is only what has no equivalent:
- the Seurat-side profiler
- the mapping schema, validator and engine
- the part→constructor table
- repair
- the LLM backends

---

## 2. Ground rules

- Branch `contract-prototype` off `dev_connect`; PRs target `dev_connect`. Commit per phase.
- New code: `R/contract_*.R`. Tests: `tests/testthat/test-contract_*.R`. Data: `inst/contract/`,
  `inst/mappings/`, `inst/schema/`. Comparisons live outside the package in
  `examples/seurat_comparison/`, which is build-ignored.
- Exported names: `seuratToGiotto()`, `profileObject()`, `fingerprint()`, `giottoContract()`,
  `readMapping()`, `findMapping()`, `validateMapping()`, `convertWithMapping()`,
  `readbackCheck()`, `diffProfiles()`, `repairMapping()`, `writeMappingRequest()`.
  Never `profile()`, `convert()` or `validate()`.
- **Dependencies**:
  - Imports: add `yaml` and `digest`.
  - Suggests: add `GiottoDisk` and `httr2`.
  - SeuratObject and Seurat are already in Suggests. New code uses SeuratObject only; Seurat
    only in example build scripts.
  - Forbidden: `jsonvalidate`, `ellmer`, Python.
- Tests needing optional packages use `skip_if_not_installed()`. LLM tests use
  `skip_if(Sys.getenv("GIOTTO_LLM_PROVIDER") == "")`. Never commit API keys.
- Fail loudly: no `try(silent = TRUE)`.
- Development stack: scratch library `~/Documents/Software/rlibs/suite_dev`, with GiottoClass
  `gsource` ≥ 0.7.28, GiottoDisk `dev` and sedonadb, plus duckdb, SeuratObject and Seurat
  from the user library.
- Verify every stated fact against code or a live object; record deviations in phase reports.

---

## 3. Layout

```
R/
  contract_profile.R     profileObject(), print, fingerprint(), diffProfiles()
  contract_schema.R      hand-written validator driven by mapping.schema.json
  contract_giotto.R      giottoContract(): manifest schema + part/constructor table
  contract_find.R        readMapping(), findMapping(), summary/print for mappings
  contract_engine.R      param binding, source_nodes, ops, convertWithMapping()
  contract_hooks.R       hookRegistry() + seurat_graph_to_nn
  contract_validate.R    validateMapping()
  contract_readback.R    readbackCheck()
  contract_repair.R      rules, alternatives, candidates, repairMapping()
  contract_llm.R         anthropic / ollama backends
  seurat.R               seuratToGiotto() thin wrapper (replaces the deprecated stub)
inst/
  schema/mapping.schema.json
  contract/giotto-contract-<GiottoClass version>.json
  mappings/seurat-assay5.yaml (+ .profile.json)
  scripts/build_fixtures.R, scripts/regenerate_contract.R
examples/seurat_comparison/   build_seurat_minis.R, benchmark scripts, results
docs/contract-prototype/      this file, PHASE0_ORIENTATION.md, REPORT.md
```

---

## 4. Core specifications

### 4.1 Paths
Structural paths are `/`-separated (Seurat names contain dots): `assays/{assay}/layers/counts`,
`meta.data`, `images/{fov}/boundaries/centroids`. A segment is an S4 slot, list element or
environment entry. `{name}` is a parameter bound at match time, not drift.

### 4.2 Profile (`profileObject()`)
- An S3 `giottoProfile` with `source` (class, package, version; Seurat uses `x@version`),
  `profiled_at`, `profiler_version`, `n_cells`/`n_feat` (and which node they came from), plus
  a node table.
- Node table columns: `path, kind, class, dtype, shape, n_children, role_hints, sample`.
  - Field names follow the GiottoClass manifest leaf vocabulary (`shape`, `dtype`,
    `sparse`, `nnz`, `columns`).
  - `kind` is one of S4, list, environment, data.frame, dense_matrix, sparse_matrix,
    vector, factor, null, external_pointer, function, other.
  - `role_hints` are cheap and comma-separated: `rows_eq_ncells`, `cols_eq_ncells`,
    `rows_eq_nfeat`, `cols_eq_nfeat`, `integer_valued`, `nonneg`, `binary_valued`,
    `two_numeric_cols`, `has_dimnames`.
- Walks via `attributes()`/`slotNames` and works without the source package loaded.
- Depth limit 8; at most 50 children walked per list, with truncation recorded.
- Stop classes are recorded but not descended:
  - `sp` polygon internals
  - external pointers
  - DelayedArray / BPCells / HDF5Array / terra (each with a note naming the package
    needed to profile it)
- Matrix role hints are computed on at most 10,000 non-zero values.

### 4.3 Fingerprint
- `digest::digest(algo = "xxhash64")` over sorted `(templated_path, kind, dtype, ndim)`.
- Mapping scope (default) covers nodes the mapping references plus their parents;
  `scope = "full"` covers all nodes.
- Parameter values in paths are replaced by `{param}` first.
- Sizes and values are excluded.
- `findMapping()` computes one fingerprint per candidate mapping.

### 4.4 Mapping file
```yaml
mapping:
  id: seurat-assay5                # no package version in the id
  version: 0.1.0
  status: verified                 # verified | patched-unreviewed | failed | draft
  built: 2026-11-19
  schema_version: 0.1.0
  source:
    ecosystem: seurat
    class: Seurat
    package: SeuratObject
    built_with: "5.4.0"
    applies_to:                    # filter used by findMapping()
      package_version: ">= 5.0.0"
      predicates:
        - { path: "assays/{assay}", class: Assay5 }
  target: { contract: giotto, version: "0.7.30" }
  fingerprint: "<hash>"
  fixtures: [seurat5_fixture_small]
  history: []

params:
  assay:      { list: Assays, rule: first_where, where: { has_layer_role: integer_valued } }
  layer_raw:  { list: Layers, on: "{assay}", rule: first_where, where: { role: [integer_valued, nonneg] } }
  layer_norm: { list: Layers, on: "{assay}", rule: name_is, value: data, optional: true }
  fov:        { list: Images, rule: first_of_class, class: FOV, optional: true }
  nn_graph:   { list: Graphs, rule: first_where, where: { name_matches: "_nn$", or_role: binary_valued }, optional: true }

source_nodes:
  counts:       { path: "assays/{assay}/layers/{layer_raw}",
                  get: { fn: LayerData, args: { assay: "{assay}", layer: "{layer_raw}" } } }
  normalized:   { path: "assays/{assay}/layers/{layer_norm}",
                  get: { fn: LayerData, args: { assay: "{assay}", layer: "{layer_norm}" } } }
  cell_meta:    { path: "meta.data", get: { fn: "[[" } }
  centroids:    { path: "images/{fov}/boundaries/centroids",
                  get: { fn: GetTissueCoordinates, on: "images/{fov}", args: { which: centroids } } }
  segmentation: { path: "images/{fov}/boundaries/segmentation",
                  get: { fn: GetTissueCoordinates, on: "images/{fov}/boundaries/segmentation" } }
  pca:          { path: "reductions/pca", get: { fn: Embeddings, args: { reduction: pca } } }
  umap:         { path: "reductions/umap", get: { fn: Embeddings, args: { reduction: umap } } }
  nn:           { path: "graphs/{nn_graph}", get: { fn: "[[", args: { i: "{nn_graph}" } } }

parts:
  - id: expr_raw
    part: expression
    from: counts
    to: { spat_unit: cell, feat_type: rna, name: raw }
  - id: locs
    part: spatial_locs
    from: centroids
    ops:
      - { op: rename, map: { x: sdimx, y: sdimy, cell: cell_ID } }
    to: { spat_unit: cell, name: raw }
    alternatives:                  # tried in order by repair rules (§4.11)
      - id: locs_from_meta
        from: cell_meta
        ops:
          - { op: select_cols, cols: ["{x_col}", "{y_col}"] }
          - { op: rownames_to_col, col: cell_ID }
          - { op: rename, map: { "{x_col}": sdimx, "{y_col}": sdimy } }
    lossy:
      - { kind: estimated, field: units, note: "FOV carries no unit; assuming micron" }
  # ... cell_metadata, spatial_info, dim_reduction (pca, umap), nn_network

lossy:
  - { kind: dropped, path: "assays/{assay}/layers/scale.data", note: "recomputable; no Giotto slot" }
```

- `get.fn` is a schema enum of SeuratObject exported getters. The engine calls
  `getExportedValue("SeuratObject", fn)` and never evaluates strings.
- If `on:` is given, the getter is applied to that node; otherwise to the object.
- Allowed ops (adding one is a schema change):

  | op | args |
  |---|---|
  | `select_cols` | `cols` |
  | `exclude_cols` | `cols` |
  | `rename` | `map` |
  | `transpose` | — |
  | `scale` | `factor`, `cols` |
  | `cast` | `to`: integer, double, character, sparse |
  | `rownames_to_col` | `col` |
  | `set` | `field`, `value` |

- Loss kinds: `dropped`, `estimated`, `ambiguous`, `recovered_from_parameters`.
- Parts: `expression`, `cell_metadata`, `feat_metadata`, `spatial_locs`, `spatial_info`,
  `dim_reduction`, `nn_network`.

### 4.5 Schema validation
Hand-written R checks driven by `mapping.schema.json`, supporting only `type`, `required`,
`enum`, `properties`, `items` and `additionalProperties: false`. Errors name the YAML path.

### 4.6 Hooks
- One initial hook, `seurat_graph_to_nn`: a Graph → edge table, with `k` recovered as
  nnz / n for binary kNN graphs (not `_snn`) and recorded as a `recovered_from_parameters`
  loss.
- Centroids and segmentation need no hook (accessor + `rename`).
- `hookRegistry()` lists each hook, its packages (checked via `package_check`) and a
  description. Report the count.

### 4.7 Giotto contract (`giottoContract()`)
- `inst/scripts/regenerate_contract.R` writes `inst/contract/giotto-contract-<ver>.json`,
  containing:
  - the GiottoClass version
  - a reference to the manifest schema
  - per part: the subobject class, the constructor and its `formals()`, nesting keys
    (`spat_unit`, `feat_type`, `name`, plus `nn_type`/`reduction`/`method` where relevant),
    required args, orientation (expression = features × cells), and `disk_backable` (true
    if the part's setter writes through `sourceWrite`; today expression, polygons, NN/spatial
    networks and points)
- `giottoContract()` uses the live GiottoClass when its version matches, otherwise the
  bundled JSON, with a warning.

### 4.8 Engine (`convertWithMapping(x, mapping, profile = NULL, source = NULL)`)
1. Profile `x` if needed, bind params, resolve `source_nodes` via `get`.
2. Per part: run ops or the hook, then the GiottoClass constructor for that part.
3. Assemble `g <- giotto()`.
   - If `source` is given: `GiottoDisk::sourceCreate()` and set `g@source` (GiottoDisk checked
     via `package_check`).
   - Then `setGiotto(g, parts)`, ordered as expression and metadata → locs → polygons →
     reductions → networks.
   - Disk-backable parts are written by the setters; the others stay in memory.
4. Provenance (mapping id, version, fingerprint, status, declared losses) goes in subobject
   `misc` and an object attribute.

Large inputs: prototype minis go through `sourceWrite(memoryMatrix)`. A chunked
`matrixColumnInput` (dgCMatrix / BPCells column blocks into GiottoDisk's `exprInput`
iterator) is proposed as a separate GiottoDisk PR before the full-size run.

### 4.9 Validation (`validateMapping(mapping, profile)`)
Returns a tidy table `(check, part_id, level OK/WARN/FAIL, message)`.

**Structural checks**
- The mapping is schema-valid.
- Params bind.
- Each node's `path` exists.
- Each `to` is valid under the contract.
- Hooks are registered and their packages are available.
- WARN for any substantial source node that is neither mapped nor declared lossy.

**Semantic checks**, on `objManifest()` of the result
- Raw expression is integer and non-negative.
- Cell IDs overlap ≥ 99% across expression, metadata, locs and polygons.
- Dimensions agree.
- Coordinates in µm with range < 50 → WARN "looks like mm"; range > 1e6 → WARN.
- dim_reduction rows equal n_cells.
- NN network nodes ⊆ cells.

### 4.10 Readback check (`readbackCheck(x, result, mapping)`)
- For every mapped part, read the Giotto content back through Giotto getters. For disk
  stores, read sampled cells through the store subset.
- Compare it with the source accessor output on shape, IDs and a sampled-value fingerprint.
- Classify each difference as MATCH / DECLARED / UNDECLARED; any UNDECLARED → FAIL.
- Phase 7 adds the reverse mapping, making it a true round trip.

### 4.11 Drift and repair
- `diffProfiles(old, new)` follows `manifestDiff`'s flatten → added/removed/modified and adds
  `moved` and `renamed_key`. It also reports `kind_changed` and `dtype_changed`.
- `repairMapping(mapping, new_profile, llm = NULL)` tries, in order:
  1. **Params**: re-bind. Renamed assay, layer or graph are absorbed here.
  2. **Rules**:
     - a single `moved`/`renamed_key` with the same kind/dtype/shape → rewrite the path
     - only `added` → WARN "unmapped new data"
     - integer↔double with integer values → add `cast`
     - primary node missing → try the part's declared `alternatives`, binding their params
       by role hints (e.g. two numeric meta.data columns)
  3. **Candidates**: up to 5 per broken part, ranked by kind (required), shape/role match,
     `adist` name similarity and shared parent.
  4. **LLM (optional)**: sees broken parts + diff rows + candidates.
     - Its patch can set `from`, `dimnames.*` or `alternative`, with `new_value` restricted
       to the enumerated candidates or alternative ids.
     - Low confidence is rejected.
  5. **Validate**, then run the readback check.
- On success: bump the patch version, set `status: patched-unreviewed`, append to `history`
  (`by: params|rules|llm`), refresh the fingerprint.
- On failure: `status: failed`, and `writeMappingRequest()` writes a markdown handoff
  (contract excerpt, profile, diff, schema).
- `findMapping(profile)`: filters on `applies_to`, then matches the fingerprint, else
  returns the closest mapping with "repair needed", or "no mapping applies".

---

## 5. LLM backends (`contract_llm.R`)

- `GIOTTO_LLM_PROVIDER` is anthropic or ollama; `GIOTTO_LLM_MODEL` sets the model.
  `ANTHROPIC_API_KEY` and `OLLAMA_HOST` (default `http://localhost:11434`) are read from the
  environment only.
- **anthropic**: `httr2` POST to `/v1/messages` with headers `anthropic-version: 2023-06-01`.
  - Structured output via a single tool whose `input_schema` is the per-request patch schema
    (enums = that request's part ids and candidates), with `tool_choice` forcing that tool.
  - Default model: the current Haiku (check the docs for the id).
  - `temperature = 0`.
- **ollama**: `/api/chat` with `format` = JSON schema, `stream = false`, temperature 0.
- Every call is logged (provider, model, prompt size, latency, tokens, result) as an
  attribute. The key is never logged.

---

## 6. Fixtures and drift variants

The base fixture is built in code by `inst/scripts/build_fixtures.R` with SeuratObject:
- 200 cells × 50 features, `rpois` with a fixed seed
- `counts`, `data` and `scale.data` layers
- meta.data including `seurat_clusters`
- an FOV with centroids and square segmentations, in µm, range 0–2000
- pca (10 dims) and umap (2 dims)
- `{assay}_nn` and `{assay}_snn` graphs
- the assay is named `Xenium`

| variant | change | expected |
|---|---|---|
| D1 | assay `Xenium` → `Spatial` | params |
| D2 | layer `counts` → `raw_counts` | params (role `integer_valued, nonneg`) |
| D3 | coords moved to meta.data `x_um`,`y_um`, FOV removed | rules via declared alternative |
| D4 | counts stored as double with integer values | rules add `cast` |
| D5 | extra reduction `tsne` | no change, WARN unmapped |
| D6 | coords in mm (÷1000) | semantic WARN "looks like mm" |
| D7 | graph `{assay}_nn` → `knn_graph` | params (`binary_valued`) |
| D8 | unrelated slot removed | fingerprint unchanged |
| D9 | v3 `Assay` instead of `Assay5` | `findMapping` picks another mapping or reports none applies; never mis-converts |

---

## 7. Phases

| Phase | Content | Stop |
|---|---|---|
| 0 | Orientation, read-only. Map the legacy code path. Check the dev stack; check `setGiotto` with `source` writes to disk; check `as.list`/`setGiotto` round trip. Probe SeuratObject structures and accessors. Check Seurat loaders on the GiottoData crops. Record the legacy converter output on the fixture. Write `PHASE0_ORIENTATION.md`. | **CHECKPOINT 0** |
| 1 | Profiler, fingerprint, diff, schema + validator, contract generation, base fixture | |
| 2 | `seurat-assay5.yaml`, the hook, engine (memory + disk), validation, `seuratToGiotto()`, summary/print | **CHECKPOINT 2** |
| 3 | Readback check; fix or declare every difference | |
| 4 | D1–D9, `findMapping`, `repairMapping` without LLM | **CHECKPOINT 4** |
| 5 | LLM backends; run the remaining variants; `writeMappingRequest()` | |
| 6 | Report (`REPORT.md`) | |
| 7 | Giotto → Seurat mapping on the same engine (manifest as profile; `CreateSeuratObject`, `CreateFOV`, `CreateDimReducObject` as constructors); `giottoToSeurat()` | |
| 8 | Benchmark vs legacy converters on the minis | **CHECKPOINT 8**, then full size |

### Benchmark (Phase 8, `examples/seurat_comparison/`)

**Axes**
- correctness (manifestDiff vs truth + readback)
- coverage (converted / dropped, silent vs declared)
- flexibility (D1–D9 outcomes)
- adaptability (lines/files of the minimal fix actually written, R vs YAML)
- cost (wall time + per-process RSS, serialized, `caffeinate -dimsu`)

Every comparison runs in memory and on disk.

**Independent oracle**: for vendor-layout minis, compare against GiottoDisk's native
importer on the same folder.

**Data, in order (agree each step before the next)**
1. The synthetic fixture.
2. GiottoData `xenium_mini_lung_ffpe`.
3. GiottoData `cosmx_mini_colon_wtx`.
4. GiottoData `visiumhd_mini_crc` (VisiumV2 images → a different mapping family).
5. `loadGiottoMini()` and disk objects for the reverse direction.
6. The full Xenium lung run (162k cells), for cost only.

Not yet covered: classic Visium (VisiumV1) and vendor-layout MERSCOPE.

### Report contents
- Lines of code (R vs YAML) new vs legacy.
- Dependencies added.
- Hooks.
- Parts and declared losses.
- Legacy bugs found.
- The drift table.
- LLM log summary.
- Schema limitations, and what AnnData/SpatialExperiment would stress (note GiottoDisk's
  `h5adInput` family as reusable building blocks).
- Open decisions.
