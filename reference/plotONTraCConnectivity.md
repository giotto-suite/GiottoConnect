# Plot ONTraC niche cluster connectivity

Draws the niche clusters of
[`runONTraC()`](https://giotto-suite.github.io/GiottoConnect/reference/runONTraC.md)
as a graph: one node per niche cluster (labelled by its number),
coloured by its NT score, and edges shaded by how strongly the clusters
connect (ONTraC's consolidated adjacency). This view has no general
Giotto plot, so it lives with the connection.

## Usage

``` r
plotONTraCConnectivity(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  name = "ontrac",
  min_connectivity = 0,
  seed = 1L,
  show_plot = NULL,
  return_plot = NULL,
  save_plot = NULL,
  save_param = list(),
  default_save_name = "ONTraCConnectivity"
)
```

## Arguments

- gobject:

  giotto object after
  [`runONTraC()`](https://giotto-suite.github.io/GiottoConnect/reference/runONTraC.md)

- spat_unit, feat_type:

  spatial unit and feature type

- name:

  name given to
  [`runONTraC()`](https://giotto-suite.github.io/GiottoConnect/reference/runONTraC.md)

- min_connectivity:

  numeric, edges at or below it are not drawn

- seed:

  integer, seed of the layout

- show_plot, return_plot, save_plot, save_param, default_save_name:

  see
  [`GiottoVisuals::plot_output_handler()`](https://giotto-suite.github.io/GiottoVisuals/reference/plot_output_handler.html)

## Value

A `ggplot` (when returned).

## References

Wang W, Zheng S, Shin SC, Chávez-Fuentes JC, Yuan GC. ONTraC
characterizes spatially continuous variations of tissue microenvironment
through niche trajectory analysis. *Genome Biology* 26, 117 (2025).
[doi:10.1186/s13059-025-03588-5](https://doi.org/10.1186/s13059-025-03588-5)

## See also

[`runONTraC()`](https://giotto-suite.github.io/GiottoConnect/reference/runONTraC.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plotONTraCConnectivity(g, spat_unit = "aggregate")
} # }
```
