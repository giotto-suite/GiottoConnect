# Read a converted object back and compare it with its source

For every part a mapping converted, read the Giotto content back through
Giotto getters and compare it with what the mapping prepared from the
source (the accessor output after hooks and ops): placement under the
declared nesting keys, IDs, shape, columns and values. Disk-backed parts
are compared on a sample of cells read through the store.

Each comparison is `MATCH`, `DECLARED` (the difference is covered by a
loss declared in the mapping or recorded during conversion) or
`UNDECLARED`.

## Usage

``` r
readbackCheck(x, result, mapping = NULL, n_cells = 99L)
```

## Arguments

- x:

  the source object

- result:

  object converted from `x` (giotto, or Seurat for a Giotto-to-Seurat
  mapping)

- mapping:

  mapping used for the conversion. Defaults to the one recorded in
  `result` (see
  [`mappingProvenance()`](https://giotto-suite.github.io/GiottoConnect/reference/mappingProvenance.md)).

- n_cells:

  number of cells sampled for expression values

## Value

a `giottoReadback` data.frame with columns `part_id`, `part`, `aspect`,
`status`, `detail`; attribute `pass` is `TRUE` when nothing is
`UNDECLARED`
