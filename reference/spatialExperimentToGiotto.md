# Convert a SpatialExperiment to a Giotto object

Convert a SpatialExperiment or SpatialFeatureExperiment with the mapping
that fits its structure (the layouts SpatialExperimentIO, VisiumIO and
SpatialFeatureExperiment's readers produce). The object is profiled, the
matching bundled mapping is selected (or the one given) and the
conversion runs through
[`convertWithMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/convertWithMapping.md).

## Usage

``` r
spatialExperimentToGiotto(x, mapping = NULL, source = NULL, verbose = FALSE)
```

## Arguments

- x:

  SpatialExperiment or SpatialFeatureExperiment

- mapping:

  optional `giottoMapping`, mapping file or bundled mapping id. By
  default the bundled mapping that applies to `x` is used.

- source:

  `NULL` for an in-memory giotto object, or a directory for a
  GiottoDisk-backed object

- verbose:

  logical

## Value

giotto object; see
[`mappingProvenance()`](https://giotto-suite.github.io/GiottoConnect/reference/mappingProvenance.md)
for what was converted and which losses apply

## See also

[`findMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/findMapping.md),
[`convertWithMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/convertWithMapping.md),
[`readbackCheck()`](https://giotto-suite.github.io/GiottoConnect/reference/readbackCheck.md)

## Examples

``` r
if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
    set.seed(1)
    counts <- Matrix::rsparsematrix(20, 30, density = 0.3,
        rand.x = function(n) stats::rpois(n, 5) + 1)
    dimnames(counts) <- list(paste0("gene", 1:20), paste0("cell", 1:30))
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        spatialCoords = cbind(x = stats::runif(30, 0, 500), y = stats::runif(30, 0, 500))
    )
    g <- spatialExperimentToGiotto(spe)
    g
    mappingProvenance(g)$mapping
}
```
