# Phase 5: LLM-assisted repair

2026-11-20, branch `contract-prototype`.

## Design

The LLM only *selects*. Repair runs in tiers (params, rules, then LLM). Only items the
rules leave open are offered to the model, as closed choices:
- **`from`**: a candidate source path for a part whose source is gone. Applied as a
  structural `path` read; validation flags it for review.
- **`param`**: a value for an unbound param of a declared alternative, e.g. which metadata
  column is `x_col`. Applied by switching the part to that alternative.

How the reply is constrained and checked:
- Each request carries a JSON schema built for it: `anyOf` one object per item, with
  `item` as a `const`, `value` an `enum` of that item's options, and `confidence`.
  Structured outputs enforce it on the Anthropic API (`output_config.format`, Haiku 4.5 by
  default) and Ollama (`format`).
- Replies are checked again: off-list values, unknown items and `low` confidence are
  dropped, and the rationale is cut at 300 characters.
- The prompt holds paths, kinds, classes, shapes, role hints and names, but no data
  values. It is about 400–800 tokens.

The gate is the same as for rule repairs: structural validation, conversion,
`readbackCheck()`. A new semantic check, `feat_ids`, rejects expression parts whose
features are not the object's features. On success: `patched-unreviewed` with
`history.by: llm`. On failure: `failed`, and `writeMappingRequest()` writes a markdown
handoff.

`convertWithMapping()` never reaches a backend (tested). Keys are read from the
environment at call time; `httr2` redacts them.

## Benchmark (`examples/seurat_comparison/llm_repair.R`, `results/phase5_llm.csv`)

Declared before running: keep the LLM step only if it beats the baseline. The baseline
takes the top-ranked candidate path when it leads by ≥ 0.5, and makes no param choices.

| case | drift | truth | baseline | qwen3.8:27b (Ollama, local) | Haiku 4.5 (subagent) |
|---|---|---|---|---|---|
| L1 | counts moved to `misc/stash` | `misc/stash` | not repaired (lead 0.33) | **correct**, 479 in / 117 out tokens, 7 s | **correct** |
| L2 | `misc/counts_corrected` + `misc/umi_raw` | `umi_raw` | not repaired (the wrong `counts_corrected` ranks first) | **correct**, 7 s | **correct** |
| L3 | CosMx columns: local and global `CenterX/Y_*_px`, FOV removed | global | not repaired | **WRONG ACCEPTED**: picked the local columns | **correct** (global) |
| L4 | coordinates only as a `spatial` reduction | rules | correct (rules) | correct (rules) | (rules) |
| L5 | counts gone, nothing replaces them | refuse | refused | refused (its pick was dropped by the check) | refused: it picked the `data` layer (medium confidence); the gate rejected it (not integer counts) |

**How Haiku was run.** No API key is set on this machine, so each case's exact prompt and
schema went to a Claude Code Haiku subagent in a fresh context. Its JSON answers were
replayed through the same check and gate (`llm_replay.R`, rows `backend = replay`).
- The output format was enforced by the checker, not by API structured outputs.
- Token counts are not comparable.
- `llm_repair.R` runs Haiku through the API once `ANTHROPIC_API_KEY` is set.

### Reading the results

- **The LLM beats the baseline where choosing needs meaning** (L1, L2): it read "stash"
  and "umi_raw" correctly.
- **L3 is the safety finding.** Local per-FOV coordinates are internally consistent, so
  validation and readback cannot tell them from global ones; the wrong choice reached
  `patched-unreviewed`.
  - Mitigations: human review is required before `verified`; a domain hint in the prompt
    ("whole-section, not per-FOV, coordinates"); or a rule for CosMx naming. The prompt
    was deliberately not changed after the run.
- **L4 needed no LLM.** A declared `locs_from_reduction` alternative handles it. It also
  repairs the real `xe_beforeRCTD_sample_cores.rds` (400,216 cells, Seurat 5.1, no
  image): readback PASS, converted in seconds.

## Fixes found while building this

| Issue | Fix |
|---|---|
| Requiring integer data in the `assay` param made every assay-dependent part unresolved when only the counts were gone | dropped; mappings bind the default assay by class, and missing counts surface as an unbound `layer_raw` |
| A param no part reads (after a patch) failed validation | now WARN |
| `rownames_to_col` duplicated an existing `cell_ID` column (real object) | reused when it equals the row names, else kept as `cell_ID.source` |
| The repair report listed every new column and command setting | now only new matrices and tables |
| New op `set_colnames` and `cast: table` | used for 2-d reductions |
| New accessor `path` | structural read, for patches only |

## Recommendation

Keep the LLM step, as opt-in and selector-only.
- **It wins where naming carries the meaning.** Haiku got 3/3 where the baseline got 0/3,
  and the local 27B model got 2/3.
- **Model choice matters.** The local model accepted per-FOV coordinates; Haiku used the
  FOV removal as the cue for global ones.
- **The gate catches wrong picks it can measure.** Haiku's L5 pick (normalized values as
  counts) was rejected.
- **The gate cannot see semantically wrong but consistent choices** (the local model's
  L3), so `patched-unreviewed` must stay a hard stop before a patched mapping is
  verified.
