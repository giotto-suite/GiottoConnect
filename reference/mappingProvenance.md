# Conversion record of a mapped giotto object

Retrieve what
[`convertWithMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/convertWithMapping.md)
recorded in the giotto object history: mapping id, version, status,
fingerprint, bound params, converted parts and losses.

## Usage

``` r
mappingProvenance(gobject)
```

## Arguments

- gobject:

  giotto object, a Seurat object or SpatialExperiment produced by a
  mapping, or the path of an `.h5ad` file written by
  [`giottoToAnnData()`](https://giotto-suite.github.io/GiottoConnect/reference/giottoToAnnData.md)

## Value

list, or `NULL` if the object was not produced by a mapping
