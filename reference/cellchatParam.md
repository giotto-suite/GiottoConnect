# CellChat parameters

Parameters for
[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
that run CellChat (Jin et al. 2021) and its spatial mode from CellChat
v2 (Jin et al. 2025), which infers cell-cell communication between cell
groups from the expression of ligands, receptors and their cofactors
(CellChatDB), with spatial distance between the groups as a constraint
when coordinates are given. CellChat is an R package from GitHub,
install it with `remotes::install_github("jinworks/CellChat")`. **If you
use it, please cite the original paper** (see References).

## Usage

``` r
cellchatParam(
  species = c("mouse", "human", "zebrafish"),
  db_search = NULL,
  non_protein = FALSE,
  type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
  trim = 0.1,
  population_size = FALSE,
  variable_both = TRUE,
  thresh_p = 0.05,
  do_fast = NULL,
  nboot = 100L,
  seed = 1L,
  min_cells = 10L,
  thresh = 0.05,
  distance_use = TRUE,
  interaction_range = 250,
  scale_distance = 0.01,
  k_min = 10L,
  contact_dependent = TRUE,
  contact_range = 100,
  ratio = NULL,
  center_distance = NULL,
  tol = NULL,
  ...
)

# S4 method for class 'allMatrix,cellchatParam'
analyzeData(x, param, cell_type, spatial_locs = NULL, ...)
```

## Arguments

- species:

  CellChatDB to use: `"mouse"`, `"human"` or `"zebrafish"`

- db_search:

  CellChatDB annotations to keep (`subsetDB()`); `NULL` keeps all
  protein signaling ("Secreted Signaling", "ECM-Receptor", "Cell-Cell
  Contact"), CellChat's recommendation

- non_protein:

  logical, also keep "Non-protein Signaling"

- type, trim:

  how average expression per group is computed (CellChat's
  `computeCommunProb()`; `"truncatedMean"` with `trim = 0.1` is less
  strict for low-depth spots)

- population_size:

  logical, weigh by group sizes

- variable_both:

  logical, require both ligand and receptor to be over-expressed
  (`identifyOverExpressedInteractions()`)

- thresh_p:

  p-value threshold for over-expressed genes

- do_fast:

  logical, CellChat's fast Wilcoxon test, which needs the GitHub package
  presto; `NULL` uses it when presto is installed

- nboot:

  integer, permutations for the p-values

- seed:

  integer, seed of the permutations

- min_cells:

  integer, minimum cells per group (`filterCommunication()`)

- thresh:

  p-value threshold of the returned interactions

- distance_use:

  logical, spatial mode: use distance as a constraint

- interaction_range:

  numeric, spatial mode: maximum interaction length in microns

- scale_distance:

  numeric, spatial mode: distance scale factor

- k_min:

  integer, spatial mode: minimum number of interacting cell pairs
  between two groups

- contact_dependent:

  logical, spatial mode: run "Cell-Cell Contact" pairs only between
  touching groups

- contact_range:

  numeric, spatial mode: contact distance in microns (about the cell
  diameter, or the spot distance for Visium)

- ratio:

  microns per coordinate unit of the spatial locations; `NULL` derives
  it from `center_distance`

- center_distance:

  microns between neighbouring centres (e.g. 100 for Visium spots); with
  `ratio = NULL` the ratio is `center_distance` divided by the median
  nearest-neighbour distance of the coordinates. One of `ratio` or
  `center_distance` is needed for the spatial mode

- tol:

  tolerance in microns added to the ranges (CellChat uses half the spot
  or cell size, e.g. 32.5 for Visium); `NULL` uses half the median
  nearest-neighbour distance in microns

- ...:

  further options of `computeCommunProb()` as `name = value`

- x:

  normalized (log) expression matrix, features x cells, feature names as
  gene symbols of `species`

- param:

  a `cellchatParam`

- cell_type:

  cell group of each column of `x`

- spatial_locs:

  `NULL` (non-spatial mode) or a two-column matrix of coordinates, rows
  in the column order of `x`

## Value

A `cellchatParam` object.

## Details

`analyzeData(x, param, cell_type, spatial_locs)` returns a `data.table`
with one row per significant interaction (`source`, `target`, `ligand`,
`receptor`, `prob`, `pval`, `interaction_name`, `pathway_name`,
`annotation`, ...; CellChat's `subsetCommunication()`). Its attributes
are `pathways` (the pathway-level table), `cellchat` (the CellChat
object, for CellChat's plotting functions) and `spatial_factors` (the
`ratio` and `tol` used; `NULL` when non-spatial).

## References

Jin S, Guerrero-Juarez CF, Zhang L, Chang I, Ramos R, Kuan CH, Myung P,
Plikus MV, Nie Q. Inference and analysis of cell-cell communication
using CellChat. *Nature Communications* 12, 1088 (2021).
[doi:10.1038/s41467-021-21246-9](https://doi.org/10.1038/s41467-021-21246-9)

Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell
communication from single-cell transcriptomics. *Nature Protocols* 20,
180-219 (2025).
[doi:10.1038/s41596-024-01045-4](https://doi.org/10.1038/s41596-024-01045-4)

## See also

[`runCellChat()`](https://giotto-suite.github.io/GiottoConnect/reference/runCellChat.md)

## Examples

``` r
p <- cellchatParam(species = "mouse", center_distance = 100, tol = 32.5)
p$species
```
