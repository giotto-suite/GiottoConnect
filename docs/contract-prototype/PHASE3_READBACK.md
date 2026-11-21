# Phase 3: readback check

2026-11-19, branch `contract-prototype`.

## What `readbackCheck()` checks

`readbackCheck(x, result)` re-runs the mapping's preparation step on the source.
"Prepared" is the accessor output after hooks and ops, plus the final, bound target keys.
It then reads every part back from the giotto object through Giotto getters. Disk-backed
expression is read for a sample of up to 99 cells through the store.

It compares:
- placement under the target keys
- IDs
- shape
- columns
- values (expression, tables, reductions)
- extent per polygon
- edges, weights and distances for networks

Each difference is MATCH, DECLARED (covered by a declared or recorded loss) or
UNDECLARED, and the check fails on any UNDECLARED.

This replaces the round trip through `giottoToSeuratV5` from the original spec.
That converter is legacy; a true round trip comes with the Giotto → Seurat mapping (Phase 7).

## Results (`examples/seurat_comparison/readback_all.R`, `results/phase3_readback.csv`)

| dataset | mapping | memory | GiottoDisk |
|---|---|---|---|
| fixture (FOV) | seurat-assay5 | PASS, 30 MATCH | PASS, 30 MATCH |
| fixture (Visium HD) | seurat-visium | PASS, 27 MATCH | PASS, 27 MATCH |
| Xenium crop | seurat-assay5 | PASS, 29 MATCH, 1 DECLARED | PASS, 29 MATCH, 1 DECLARED |
| CosMx crop | seurat-assay5 | PASS, 30 MATCH | PASS, 30 MATCH |
| Visium HD crop | seurat-visium | PASS, 27 MATCH | PASS, 27 MATCH |

On Xenium, the DECLARED item is 2 of 1531 polygons whose invalid source rings were
repaired (`polygons: estimated geometry`).

## What the readback found and how it was fixed

The first run had two UNDECLARED differences in every dataset.

1. **Network weights were all NA.** GiottoClass derives kNN weights as 1 / (1 + distance)
   and the hook passed NA distances, so the stored network had no usable weights.
   - **Fix:** when Seurat's FindNeighbors record lists the graph, the hook recomputes each
     edge's distance in the recorded reduction, dims and metric (`annoy.metric`,
     `l2.norm`), then applies Giotto's weight formula. This is recorded as
     `recovered_from_parameters`.
   - **Check on Xenium:** 97.8% of Seurat's edges fall inside the exact 19-NN radius
     recomputed in pca 1:20. The rest is Annoy's approximate search, so the space is right.
   - **Without a record,** distance and weight are recorded as `dropped`.
2. **The polygon comparator used data.table syntax in a package that is not data.table
   aware.** Rewritten in base R.

## Other changes

- **The engine is split into preparation and building**, so the readback compares against
  exactly what the constructors received.
- **Placeholders in `to` (e.g. `name: "{nn_graph}"`) are now bound.** They previously only
  worked because the hook overrides the name.
