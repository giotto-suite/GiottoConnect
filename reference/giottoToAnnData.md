# Convert a Giotto object to an AnnData .h5ad file

Write one spatial unit and feature type of a giotto object (in memory or
GiottoDisk-backed) as an AnnData `.h5ad` file in scanpy's layout,
without python (the file is written with anndataR). `X` holds the
normalized values (raw counts when there are none) and `layers/counts`
the raw counts; other expression values become layers. Reductions go to
`obsm/X_<name>` (PCA loadings to `varm/PCs`), networks to
`obsp/<name>_distances` / `_connectivities` with a `uns/<name>` record,
so
[`anndataToGiotto()`](https://giotto-suite.github.io/GiottoConnect/reference/anndataToGiotto.md)
brings them back under the same names.

## Usage

``` r
giottoToAnnData(
  gobject,
  file,
  spat_unit = NULL,
  feat_type = NULL,
  mapping = NULL,
  verbose = FALSE
)
```

## Arguments

- gobject:

  giotto object

- file:

  path of the `.h5ad` file to write

- spat_unit, feat_type:

  spatial unit and feature type to export (default: the active ones)

- mapping:

  optional `giottoMapping`, mapping file or bundled mapping id

- verbose:

  logical

## Value

the path of the written file (invisibly); see
[`mappingProvenance()`](https://giotto-suite.github.io/GiottoConnect/reference/mappingProvenance.md)
for what was converted

## See also

[`anndataToGiotto()`](https://giotto-suite.github.io/GiottoConnect/reference/anndataToGiotto.md),
[`readbackCheck()`](https://giotto-suite.github.io/GiottoConnect/reference/readbackCheck.md)
