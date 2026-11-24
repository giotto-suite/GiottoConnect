# Run CellChat cell-cell communication

Runs CellChat (Jin et al. 2021) on the normalized expression and a cell
type column of a giotto object, in CellChat v2's spatial mode (Jin et
al. 2025) when `spatial = TRUE`, and returns the significant
ligand-receptor interactions between cell types.

The giotto object is not changed: the result is a table between cell
groups, which no giotto slot holds (as with Giotto's `spatCellCellcom()`
it is returned). The spatial mode needs the micron scale of the spatial
locations: give `ratio` or `center_distance` in the param (see
[`cellchatParam()`](https://giotto-suite.github.io/GiottoConnect/reference/cellchatParam.md)).
**If you use it, please cite the original paper** (see
[`cellchatParam()`](https://giotto-suite.github.io/GiottoConnect/reference/cellchatParam.md)).

## Usage

``` r
runCellChat(
  gobject,
  cell_type,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = "normalized",
  spat_loc_name = "raw",
  spatial = TRUE,
  param = NULL,
  verbose = NULL
)
```

## Arguments

- gobject:

  giotto object

- cell_type:

  cell metadata column with the cell types

- spat_unit, feat_type:

  spatial unit and feature type

- expression_values:

  normalized expression values to use

- spat_loc_name:

  name of the spatial locations

- spatial:

  logical, run the spatial mode (needs spatial locations)

- param:

  a
  [`cellchatParam()`](https://giotto-suite.github.io/GiottoConnect/reference/cellchatParam.md);
  its defaults when `NULL`

- verbose:

  be verbose

## Value

A `data.table` of significant interactions (see
[`cellchatParam()`](https://giotto-suite.github.io/GiottoConnect/reference/cellchatParam.md))
with attributes `pathways`, `cellchat` (the CellChat object, for
CellChat's plots), `spatial_factors` and `param`.

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

[`cellchatParam()`](https://giotto-suite.github.io/GiottoConnect/reference/cellchatParam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
res <- runCellChat(g, cell_type = "leiden_clus",
    param = cellchatParam(species = "mouse", center_distance = 100,
        tol = 32.5))
CellChat::netVisual_circle(attr(res, "cellchat")@net$count)
} # }
```
