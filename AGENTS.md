# GiottoConnect

Interoperability hub of the Giotto Suite: conversions between Giotto and other
ecosystems (Seurat, AnnData, SpatialData, Bento, Stereo-seq GEF).

Seurat ↔ Giotto and AnnData → Giotto run on the **mapping engine**, a declarative and
ecosystem-neutral converter driven by YAML mapping files. The other ecosystems still use
hand-written converters and move to the engine one adapter at a time (SpatialExperiment
next; see `docs/contract-prototype/REEVALUATION_ANNDATA.md`).

## Documentation map

| Doc | Role |
|---|---|
| `AGENTS.md` (this file) | Code navigation, constraints, invariants. Read first when modifying code. |
| `docs/contract-prototype/CONTRACT_PROTOTYPE_v2.md` | Mapping engine spec: schema, engine, validation, repair, LLM step, ecosystem adapters (§8). |
| `docs/contract-prototype/REPORT.md` | Prototype results, open decisions. Phase write-ups sit next to it. |
| `docs/contract-prototype/FOLLOWUPS_OTHER_PACKAGES.md` | Findings in GiottoDisk / GiottoData still to assess. |
| `examples/seurat_comparison/` | Comparison and benchmark scripts with their result tables (Rbuildignored). |
| `DESCRIPTION` `Remotes:` | Authoritative upstream branch pins (see below). |

## Upstream branch pins

| Package | Branch | Required because | Drop the pin when |
|---|---|---|---|
| `GiottoClass` | `gsource` | `giotto(source =)` and the backend-aware setters (disk conversion), `objManifest()`, `update_giotto_params(attachments =)`. Floored at `>= 0.7.31`; `inst/contract/` bundles a contract per supported GiottoClass version (0.7.31, 0.8.0; regenerate with `inst/scripts/regenerate_contract.R`). | those are on `dev`. |
| `GiottoUtils` | `dev` | Suite convention. | `main` catches up. |
| `GiottoDisk` (Suggests) | `dev` | `sourceCreate()` / `sourceWrite()` for `seuratToGiotto(source =)`, `h5adInput()` for `anndataToGiotto(source =)`. Needs Giotto `gsource` installed first (`VisiumReader`). | it is released. |

## Package structure

```
R/
  seurat.R               seuratToGiotto() / giottoToSeurat(): thin wrappers over the mapping engine
  contract_adapter.R     adapters (seurat, giotto) and targets (giotto, seurat); adapterRegistry()
  contract_profile.R     profileObject(), fingerprint(), diffProfiles()
  contract_schema.R      mapping schema checker (inst/schema/mapping.schema.json)
  contract_giotto.R      giottoContract(): part -> constructor/setter/nesting, generated from GiottoClass
  contract_find.R        readMapping(), findMapping(), summary/print
  contract_engine.R      param binding, accessors, ops, convertWithMapping()
  contract_hooks.R       hookRegistry() and the hooks
  contract_validate.R    validateMapping() (structural + semantic)
  contract_readback.R    readbackCheck()
  contract_repair.R      repairMapping(), writeMappingRequest()
  contract_llm.R         llmBackend(), llmPatch() (opt-in, repair only)
  contract_target_seurat.R  Seurat target: assembly, key/semantic checks, readback comparators
  anndata.R, spatialdata.R, bento.R, stereoseq_gef.R
  tool-utils.R           shared helpers for external tool connections (py_require CLI path, run + log)
  tool-<name>.R          one external tool per file (adr/0008 in GiottoClass): <name>Param, analyzeData
                         method, run<Tool>(); the tool is checked or installed only when it runs
inst/
  mappings/<id>.yaml (+ <id>.profile.json)   one file per source structure, not per technology
  contract/giotto-contract-<GiottoClass version>.json
  schema/mapping.schema.json
  scripts/build_fixtures.R, verify_mapping.R, regenerate_contract.R
```

## Mapping engine: invariants

- **No LLM in the conversion path.** `convertWithMapping()` / `.prepare_parts()` never
  reach `contract_llm.R` (tested). LLMs only choose among options inside
  `repairMapping(llm =)`.
- **Verified mappings read through the source package's exported accessors** (`get.fn`
  is a schema enum). The structural `get.fn: path` read is for repair patches only;
  validation warns on it.
- **Mappings follow object structure.** A new loader or package release that keeps the
  structure needs no change. A new structure gets a new mapping file, never a new
  function.
- **Every patch is gated:** structural validation, then conversion, then `readbackCheck()`.
  Repair output is `patched-unreviewed` until a human marks it `verified`.
- **Losses are declared or recorded,** never silent. A readback difference without a
  loss is UNDECLARED and fails.
- **Disk backing is slot polymorphism:** `convertWithMapping(source =)` only attaches a
  GiottoDisk source; the GiottoClass setters write.

## Traps (each one cost a bug)

- **The profiler must not dispatch on S4 objects.** `inherits()`, `length()`, `dim()` and
  `is.matrix()` load the class's package (e.g. SeuratObject). Use
  `isS4`/`typeof`/`class`/`attr`/`attributes`/`unclass`.
- **YAML 1.1** reads `on`, `y`, `yes`, `no` as booleans. Read mappings with
  `.read_yaml_strict()`, write them with `.mapping_write()` (`true`/`false`).
- **Fingerprints sort with `method = "radix"`.** Locale-dependent sorting made the test
  suite and scripts disagree.
- **Seurat v5 layer matrices have no dimnames** (names live in `LogMap`s). Never read
  layers by path in a verified mapping.
- **GiottoDisk refuses to materialise > 100 × > 100.** Sample fewer than 100 cells when
  reading stores back.
- **Reverse templating (`.pf_template`) only for object names** (assays, reductions,
  images). A layer value like `data` matches `meta.data`.

## Adding an ecosystem

Everything ecosystem-specific sits behind two S4 classes in `R/contract_adapter.R`;
the engine only dispatches. `R/contract_adapter_anndata.R` is the example of an
adapter whose input is a file (.h5ad, read with hdf5r) rather than an R object.

- **Source side:** a `giottoConnectAdapter` subclass supplying `adapterProfile()` (how
  to walk the input; the default is `profileObject()`), `adapterListPaths()` (listing
  verb → profile path), `adapterAccessors()` (the functions a mapping may call),
  `adapterObjectLists()` (listings whose values are object names) and, when needed,
  `adapterRead()`. Register it in `.ADAPTERS` under the mapping's
  `source.ecosystem`.
- **Output side (only for X ← Giotto):** a `giottoConnectTarget` subclass supplying
  `targetBuild()`, `targetCheckKeys()`, `targetSemantic()` and `targetComparators()`.
  Register it in `.TARGETS` under the mapping's `target.contract`.
- Then mapping files, hooks if the op list cannot express a step, and a thin
  `<eco>ToGiotto()` / `giottoTo<Eco>()` wrapper. `readMapping()` rejects listing verbs
  and accessors the adapter does not supply.
- Add a fixture builder to `inst/scripts/build_fixtures.R` and stamp the mapping with
  `Rscript inst/scripts/verify_mapping.R <id>`.
- Prove refactors with `examples/seurat_comparison/engine_snapshot.R`: it records
  mapping choice, fingerprints, params, content hashes and readback for every fixture
  and crop, and compares against a saved baseline.

## Changing mappings or the contract

- After editing a mapping: `Rscript inst/scripts/verify_mapping.R <id>`. It converts the
  fixture in memory and on disk, requires no FAIL, and re-stamps the fingerprint and the
  stored profile.
- After a GiottoClass upgrade: `Rscript inst/scripts/regenerate_contract.R`, which writes
  `inst/contract/giotto-contract-<version>.json`.
- Adding an op, accessor or list verb is a schema change. Update
  `inst/schema/mapping.schema.json` and `.OP_ARGS`.

## Development workflow

- **Library:** develop against the scratch library
  `R_LIBS=~/Documents/Software/rlibs/suite_dev` (suite packages at their pinned branches).
- **Tests:** `testthat::test_local()`. The live LLM test is skipped unless
  `GIOTTO_LLM_PROVIDER` is set. `R CMD check` WARNINGs about loading come from
  GiottoClass's `.onLoad` (`getClass()` without `methods::`), not from this package.
- **Docs:** roxygen2 8.0.0 (suite pin).
- **Branches:** feature branches off `dev_connect`; PRs target `dev_connect`.
- **Results:** comparison results go in `examples/seurat_comparison/results/`, with time
  as a column. Run on the crops first; serialize timing runs under `caffeinate -dimsu`.
