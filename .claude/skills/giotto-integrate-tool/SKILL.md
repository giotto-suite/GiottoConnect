---
name: giotto-integrate-tool
description: Integrate an external spatial-omics method into the Giotto Suite — connect an outside tool (R or Python package) in GiottoConnect, or reimplement a published method natively in Giotto. Use when asked to add, wrap, connect, port or reimplement a tool or method (e.g. "add CellChat to Giotto", "connect ONTraC", "native version of BANKSY").
---

# Integrate a method into the Giotto Suite

Follow GiottoClass `adr/0008-extensions-attach-to-the-verbs.md`. Everything you
need is below; do not explore the suite to rediscover it.

## 1. Decide the kind (one minute)
- **Connection** (default): the tool keeps running its own code. One file
  `GiottoConnect/R/tool-<name>.R`. The tool is never in `Imports`: R tools go
  in `Suggests` (or only an install hint if GitHub-only) and are checked with
  `package_check()` when the function runs; Python tools use
  `.tool_python_cli()` (reticulate `py_require()`, no conda).
- **Native**: the method is a short composition of steps Giotto already has
  and scale matters. One file `Giotto/R/<method>.R`, a descriptive name (not
  the tool's), parity against the original.

## 2. Read the tool, briefly
Find: inputs (expression? coordinates? labels?), outputs, the main call or
CLI, its defaults, the paper (title, journal, DOI), licence, version.
For R packages print the bodies of the main functions; for Python, `--help`.

## 3. Write three pieces (copy `templates/connection.R` or `templates/native.R`)
1. `<name>Param()` + `setClass("<name>Param", contains = "<verb>Param")`.
   Verb by output: `analyzeData` (tables), `reduceData` (embedding),
   `clusterData` (labels), `processData` (same-shape values).
2. `setMethod("<verb>", signature(x = <data class>, param = "<name>Param"))`
   that **returns** the result (never a giotto object). Multi-output tools:
   one per-cell table, extra outputs as attributes. Several inputs (matrix +
   labels + coordinates): dispatch on the main one (`x = "allMatrix"`,
   features x cells) and take the others as named arguments of the method.
3. `run<Tool>(gobject, ...)`: defaults via `set_default_spat_unit()` /
   `set_default_feat_type()`, inputs via accessors, call the verb, write via
   setters, end with
   `update_giotto_params(gobject, description = "_<name>", attachments = list(param = param))`.

### Result -> slot
| Output | Write with |
|---|---|
| label / value per cell | `addCellMetadata(g, spat_unit, feat_type, new_metadata = dt, by_column = TRUE, column_cell_ID = "cell_ID")` |
| value per feature | `addFeatMetadata()` |
| cells x dims | `createDimObj(coordinates, name, spat_unit, feat_type, method, reduction = "cells")` + `setDimReduction(g, x, spat_unit, feat_type, reduction = "cells", reduction_method, name)` |
| cells x categories | `createSpatEnrObj(enrichment_data = dt_with_cell_ID, name, spat_unit, feat_type, method, misc)` + `setSpatialEnrichment(g, x, spat_unit, feat_type, name)` |
| features x cells | `createExprObj()` + `setExpression()` |
| describes a stored object (e.g. per-cluster scores of an enrichment) | that object's `misc` |
| anything else (interaction tables, tool objects) | **return it** from `run<Tool>()`, leave the giotto object unchanged, set `attr(res, "param") <- param`; never invent a spat_unit or feat_type to park it |

### Inputs
`getExpression(g, spat_unit, feat_type, values = "normalized", output = "matrix")`
(disk-backed: `output = "exprObj"` then `Giotto:::.read_expr_matrix(ex[], feats)`),
`getSpatialLocations(g, spat_unit, name = "raw", output = "data.table")` (columns
`cell_ID, sdimx, sdimy`), `pDataDT()` / `fDataDT()`, `getDimReduction(..., output = "matrix")`,
`createNetwork(locs, kNNNetworkParam(k, include_weight = FALSE, output = "data.table"))`
(directed kNN, `from, to, distance`), `reduceData(x, pcaParam("irlba", ncp, center = TRUE, scale = FALSE))`.

### R tools
`package_check("<pkg>", repository = "CRAN" | "Bioc" | "github:owner/repo")`
at the top of the method. CRAN/Bioc tools go in `Suggests`; GitHub-only
tools are not listed anywhere (install hint only) and are called with
`getExportedValue("<pkg>", "<fun>")`. Tests: `skip_if_not_installed("<pkg>")`.
Load a tool's lazy data with `utils::data(list = "<name>", package = "<pkg>")`.

### Python tools
`exe <- .tool_python_cli("<pkg>[extras]==<ver>", module = "<module>", cli = "<cli>", python_version = ">=3.10,<3.13")`,
then `.tool_run(exe, args, log_file)`. Write inputs to `tempfile()`; read
outputs with `data.table::fread()`. Pin the version in the param.

## 4. Credit (required)
Help page of the param **and** the wrapper (the check reads each page):
name the original method and paper, `@references` with `\doi{}` (on the
wrapper: `@inherit <name>Param references`), and the sentence "If you use
it, please cite the original paper". Native: state differences from the
original and the parity result.

## 5. Test and check
- `tests/testthat/test-tool-<name>.R`: offline tests (argument building, output
  parsing with tiny fake files) always; the real run behind
  `skip_if(Sys.getenv("GIOTTOCONNECT_TOOL_TESTS") == "")` + `skip_if_not_installed()`.
- Run on a mini and check conformance:
  `checkGiottoMethod(before, after, slots = c(...), param_class = "<name>Param", fn = "run<Tool>", package = "<pkg>")`;
  for methods that return: `checkGiottoMethod(g, g, slots = "none", result = res, ...)`.
  Fix every FAIL using its message.
- Parity: same input through the original tool directly; compare values
  (max abs diff) or labels (`aricode::ARI`). Scripts go in
  `~/Documents/Software/examples/<name>_comparison/`, not in the package.

## Environment facts (this machine)
- Library: `R_LIBS=/Users/rubendries/Documents/Software/rlibs/suite_dev`.
- Before loading minis: `options(giotto.has_conda = FALSE, giotto.use_conda = FALSE, giotto.no_python_warn = TRUE)`.
- Minis (`GiottoData::loadGiottoMini()`): `visium` (624 spots, `leiden_clus`, `hvf`),
  `vizgen` (use `spat_unit = "aggregate"`, 462 cells, `leiden_clus`), `cosmx`, `seqfish`
  (no normalized values: `normalizeGiotto()` first). `starmap`, `spatialgenomics` fail to load.
- Document: `roxygen2::roxygenise()` (8.0.0); then `git checkout --` any
  unrelated `man/*.Rd` or Collate reordering it produced. Install:
  `R CMD INSTALL --no-multiarch .` with the library above.
- GiottoConnect does not import data.table: always subset with a comma
  (`dt[i, ]`, `dt[, cols, with = FALSE]`) or use `data.table::` functions;
  `dt[match(...)]` falls through to `[.data.frame` and fails.
- A wrapper that calls another `run*()` wrapper: run the inner one inside
  `gwith_options(list(giotto.update_param = FALSE), ...)` and record one
  history entry yourself (nested `update_giotto_params()` fails on `...`).
- Work on a branch in a worktree; do not push.

## Worked examples
`GiottoConnect/R/tool-ontrac.R` (Python CLI, results to metadata + enrichment),
`GiottoConnect/R/tool-cellchat.R` (R package, returns its table),
`Giotto/R/neighborhood-pca.R` (native `reduceData`), `Giotto/R/hmrf.R`
(native `clusterData`, plus the `findSpatialDomains()` composition).

## Done means
Conformance all ok, parity reported, docs credit the authors, tests pass,
nothing changed outside the new file(s), NAMESPACE and man pages.
