# Convert a Seurat object to a Giotto object

Convert a Seurat object with the mapping that fits its structure. The
object is profiled, the matching mapping is selected from those bundled
with GiottoConnect (or the one given), and the conversion runs through
[`convertWithMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/convertWithMapping.md).
Object versions are handled by choosing a mapping, not by separate
functions.

## Usage

``` r
seuratToGiotto(sobject, mapping = NULL, source = NULL, verbose = FALSE)
```

## Arguments

- sobject:

  Seurat object

- mapping:

  optional `giottoMapping`, mapping file or bundled mapping id. By
  default the bundled mapping that applies to `sobject` is used.

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
[`validateMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/validateMapping.md)

## Examples

``` r
if (requireNamespace("SeuratObject", quietly = TRUE)) {
    set.seed(1)
    counts <- Matrix::rsparsematrix(20, 30, density = 0.3,
        rand.x = function(n) stats::rpois(n, 5) + 1)
    dimnames(counts) <- list(paste0("gene", 1:20), paste0("cell", 1:30))
    s <- SeuratObject::CreateSeuratObject(counts)
    # coordinates kept as metadata columns are found as well
    s$x <- stats::runif(30, 0, 500)
    s$y <- stats::runif(30, 0, 500)
    g <- seuratToGiotto(s)
    g
    mappingProvenance(g)$losses
}
```
