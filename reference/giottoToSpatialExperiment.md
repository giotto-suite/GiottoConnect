# Convert a Giotto object to a SpatialExperiment

Convert one spatial unit and feature type of a giotto object (in memory
or GiottoDisk-backed) to a SpatialExperiment, or to a
SpatialFeatureExperiment that also holds the polygons, transcripts and
spatial networks. Raw counts become the `counts` assay, normalized
values `logcounts`; reductions use scater's names (`PCA`, `UMAP`,
`TSNE`) with PCA loadings and variances as attributes; nearest-neighbour
networks become colPairs that keep their Giotto name and type, so
[`spatialExperimentToGiotto()`](https://giotto-suite.github.io/GiottoConnect/reference/spatialExperimentToGiotto.md)
brings them back. Coordinates are in image space (y down).

## Usage

``` r
giottoToSpatialExperiment(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  sfe = c("auto", "yes", "no"),
  mapping = NULL,
  verbose = FALSE
)
```

## Arguments

- gobject:

  giotto object

- spat_unit, feat_type:

  spatial unit and feature type to convert (default: the active ones)

- sfe:

  `"auto"` builds a SpatialFeatureExperiment when the spatial unit has
  polygons, the feature type has transcripts or the spatial unit has
  spatial networks (and SpatialFeatureExperiment is installed); `"yes"`
  or `"no"` decide it

- mapping:

  optional `giottoMapping`, mapping file or bundled mapping id
  (overrides `sfe`)

- verbose:

  logical

## Value

SpatialExperiment or SpatialFeatureExperiment; see
[`mappingProvenance()`](https://giotto-suite.github.io/GiottoConnect/reference/mappingProvenance.md)
for what was converted

## See also

[`spatialExperimentToGiotto()`](https://giotto-suite.github.io/GiottoConnect/reference/spatialExperimentToGiotto.md),
[`readbackCheck()`](https://giotto-suite.github.io/GiottoConnect/reference/readbackCheck.md)
