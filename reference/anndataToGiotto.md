# Convert an AnnData .h5ad file to a Giotto object

Convert an AnnData (.h5ad) file with the mapping that fits its
structure. The file is read directly (hdf5r, no python): it is profiled,
the matching bundled mapping is selected (or the one given) and the
conversion runs through
[`convertWithMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/convertWithMapping.md).

## Usage

``` r
anndataToGiotto(anndata_path, mapping = NULL, source = NULL, verbose = FALSE)
```

## Arguments

- anndata_path:

  path to the .h5ad file

- mapping:

  optional `giottoMapping`, mapping file or bundled mapping id. By
  default the bundled mapping that applies to the file is used.

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
if (requireNamespace("hdf5r", quietly = TRUE)) {
    f <- system.file("extdata", "mini_scanpy.h5ad", package = "GiottoConnect")
    g <- anndataToGiotto(f)
    g
    # what was read from where, and what could not be carried over
    prov <- mappingProvenance(g)
    prov$params
    readbackCheck(f, g)
}
```
