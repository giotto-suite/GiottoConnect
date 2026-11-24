# Run ONTraC niche trajectories

Runs ONTraC (Wang et al. 2025) on the cell types and spatial locations
of a giotto object and stores its results:

- cell metadata `<name>_NTScore` (niche trajectory score of each cell)
  and `<name>_niche_cluster` (most likely niche cluster, ordered by the
  clusters' NT scores);

- a spatial enrichment `<name>` with each cell's niche cluster
  probabilities. Its `misc` holds `niche_clusters` (each cluster's NT
  score) and `connectivity` (cluster x cluster).

ONTraC is a Python package, installed in a cached environment the first
time this runs. **If you use it, please cite the original paper** (see
[`ontracParam()`](https://giotto-suite.github.io/GiottoConnect/reference/ontracParam.md)).

## Usage

``` r
runONTraC(
  gobject,
  cell_type,
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = "raw",
  sample_column = NULL,
  param = NULL,
  name = "ontrac",
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

- spat_loc_name:

  name of the spatial locations

- sample_column:

  optional cell metadata column with the sample of each cell; one sample
  when `NULL`

- param:

  an
  [`ontracParam()`](https://giotto-suite.github.io/GiottoConnect/reference/ontracParam.md);
  its defaults when `NULL`

- name:

  prefix of the results

- verbose:

  be verbose

## Value

The giotto object with the results.

## Details

Plot the results with Giotto's plots:

- niche clusters in space:
  `spatPlot2D(g, cell_color = "<name>_niche_cluster")`

- niche cluster probabilities in space:
  `spatCellPlot2D(g, spat_enr_names = "<name>", cell_annotation_values = <clusters>)`

- NT score per cell type:
  `ridgePlot(g, feats = "<name>_NTScore", group_by = <cell_type>)`

- cell type composition of the niche clusters:
  `plotLabelComposition(g, groups = "<name>_niche_cluster", labels = <cell_type>)`

- how the niche clusters connect:
  [`plotONTraCConnectivity()`](https://giotto-suite.github.io/GiottoConnect/reference/plotONTraCConnectivity.md)

## References

Wang W, Zheng S, Shin SC, Chávez-Fuentes JC, Yuan GC. ONTraC
characterizes spatially continuous variations of tissue microenvironment
through niche trajectory analysis. *Genome Biology* 26, 117 (2025).
[doi:10.1186/s13059-025-03588-5](https://doi.org/10.1186/s13059-025-03588-5)

## See also

[`ontracParam()`](https://giotto-suite.github.io/GiottoConnect/reference/ontracParam.md),
[`plotONTraCConnectivity()`](https://giotto-suite.github.io/GiottoConnect/reference/plotONTraCConnectivity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("vizgen")
g <- runONTraC(g, cell_type = "leiden_clus", spat_unit = "aggregate",
    param = ontracParam(n_neighbors = 30, n_local = 10))
} # }
```
