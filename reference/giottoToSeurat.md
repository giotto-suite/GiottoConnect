# Convert a Giotto object to a Seurat object

Convert a giotto object (in memory or GiottoDisk-backed) to a Seurat v5
object with the bundled mapping that fits its structure (see
[`findMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/findMapping.md)).
The active spatial unit and feature type are converted: expression (raw
counts and normalized values), cell and feature metadata, spatial
locations and polygons as an FOV, dimension reductions and
nearest-neighbour networks.

## Usage

``` r
giottoToSeurat(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  mapping = NULL,
  verbose = FALSE
)
```

## Arguments

- gobject:

  giotto object

- spat_unit, feat_type:

  spatial unit and feature type to export (default: the active ones)

- mapping:

  optional `giottoMapping`, mapping file or bundled mapping id

- verbose:

  logical

## Value

Seurat object; see
[`mappingProvenance()`](https://giotto-suite.github.io/GiottoConnect/reference/mappingProvenance.md)
for what was converted and which losses apply

## See also

[`seuratToGiotto()`](https://giotto-suite.github.io/GiottoConnect/reference/seuratToGiotto.md),
[`convertWithMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/convertWithMapping.md),
[`readbackCheck()`](https://giotto-suite.github.io/GiottoConnect/reference/readbackCheck.md)
