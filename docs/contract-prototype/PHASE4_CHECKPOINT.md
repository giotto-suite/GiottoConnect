# Checkpoint 4: drift detection and rule-based repair

2026-11-19, branch `contract-prototype`.

## How drift is handled

1. **`findMapping()`** filters mappings by `applies_to`, then compares the mapping-scoped
   fingerprint: `match`, `repair_needed` or `none` (with a reason such as "param 'assay'
   could not be bound (Assays first_where: class Assay5, data that is integer_valued+nonneg)").
2. **Parameters are re-bound on every object.** Renamed assays, layers and graphs need no
   repair.
3. **`repairMapping(mapping, x)`** compares the object with the profile stored next to the
   mapping (`<id>.profile.json`) and applies rules. It needs no LLM and no SeuratObject
   calls; only the listing-path table is Seurat-specific.
   - **cast**: a source node's storage changed but its values did not (dense ↔ sparse,
     integer ↔ double). A `cast` op is added.
   - **moved / renamed_key**: a removed subtree reappears with the same structure. The
     path is rewritten, but only for nodes whose accessor targets a path (`get.on`).
     Accessors that address data by name (`LayerData`, `Embeddings`) follow renames
     through params instead.
   - **alternative**: the source of a part is gone and one of its declared alternatives
     binds. The part switches source. The outgoing source becomes an alternative,
     carrying the params only it used, so the patch can repair back.
   - **absent**: an optional part whose data is gone. Noted; the part is skipped.
   - **added**: new nodes are reported. No change is made.
   - Anything else is **unresolved**. It gets up to 5 candidate paths, ranked by kind
     (and class for S4), size and value roles, name similarity and a shared parent. These
     candidates are the input for the LLM step (Phase 5).
4. **The patch is validated before it is accepted.** Structural checks run, and when the
   object is given it is converted and read back. On success the patch version is bumped,
   the status is set to `patched-unreviewed`, `history` gets `by: rules`, and the
   fingerprint is re-stamped. Otherwise the status is `failed`.
   `repairMapping(file =)` writes the patched YAML.

The engine never repairs on its own. A verified mapping stays deterministic, and repair
produces a new version.

## Drift table (`examples/seurat_comparison/drift_table.R`, `results/phase4_drift.csv`)

| | change | expected | findMapping | repair | rules | readback | conversion warnings | legacy `seuratToGiottoV5` |
|---|---|---|---|---|---|---|---|---|
| D1 | assay Xenium → Spatial | params | match | not needed | | PASS | | raw ok, locs swapped |
| D2 | layer counts → raw_counts | params | match | not needed | | PASS | | raw ok, locs swapped |
| D3 | coordinates moved to meta.data, FOV removed | rules (alternative) | repair_needed | **repaired** | alternative, absent | PASS | params/sources for the gone FOV | **error** (subscript out of bounds) |
| D4 | counts stored dense | rules (cast) | repair_needed | **repaired** | cast | PASS | | raw ok, locs swapped |
| D5 | extra reduction tsne | params, WARN unmapped | match | not needed | | PASS | `[unmapped]` tsne | raw ok, locs swapped |
| D6 | coordinates in mm | params, WARN mm | match | not needed | | PASS | `[coord_units]` looks like mm | no warning |
| D7 | graph Xenium_nn → knn_graph | params | match | not needed | | PASS | | raw ok, locs swapped |
| D8 | unused graph removed | params | match | not needed | | PASS | | raw ok, locs swapped |
| D9 | v3 `Assay` instead of `Assay5` | no mapping applies | none | | | | | raw ok, locs swapped |

All nine reach the expected outcome. No variant needs the LLM step.

- **D3** reads coordinates from `x_um`/`y_um` (bound by name pattern), which equal the
  original centroids with y flipped. The patch (`v0.1.1`) repairs back to the FOV source on
  the original fixture.
- **D9** is refused, not mis-converted. A `seurat-assay` (v3/v4) mapping would cover it.
  The legacy converter accepts it, as it accepts everything, but with swapped axes.
- **Legacy on D2.** It finds `raw_counts` only because SeuratObject matches layer names by
  regex (`counts` matches `raw_counts`).

## Also tested

- A required part with no source fails repair, with ranked candidates. With the counts
  layer moved to `misc/stash`, `misc/stash` ranks first and the `data` layer second.
- Patched YAML is written with `true`/`false`, because `yaml::write_yaml` defaults to the
  YAML 1.1 `yes`/`no`, which the strict reader rejects.

## Fixes made while building this

| Bug | Fix |
|---|---|
| A part counted as bound when its path was bound, even with an unbound placeholder in its ops | boundness now checks path, accessor and op arguments (engine and repair) |
| `rename` map keys were not bound (`"{x_col}": sdimx`) | keys are bound like values |
| The moved rule could rewrite paths that an accessor never reads | restricted to `get.on` nodes |
| Candidate ranking ignored value roles, so the `data` layer outranked integer counts | value roles scored separately |
