# Run SPARK-X spatially variable genes

Runs SPARK-X (Zhu, Sun and Zhou 2021) on the counts and spatial
locations of a giotto object and stores its per-gene p-values as feature
metadata, in columns `<name>_pval` (combined over the spatial kernels)
and `<name>_padj` (Benjamini-Yekutieli). Features that are not tested
(not in `feats_to_use`, or without counts) get `NA`. **If you use it,
please cite the original paper** (see
[`sparkxParam()`](https://giotto-suite.github.io/GiottoConnect/reference/sparkxParam.md)).

## Usage

``` r
runSparkX(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = "raw",
  spat_loc_name = "raw",
  feats_to_use = NULL,
  covariates = NULL,
  param = NULL,
  name = "sparkx",
  verbose = NULL
)
```

## Arguments

- gobject:

  giotto object

- spat_unit, feat_type:

  spatial unit and feature type

- expression_values:

  expression values to use; SPARK-X expects counts

- spat_loc_name:

  name of the spatial locations

- feats_to_use:

  features to test: `NULL` (all), a feature metadata column (features
  with `"yes"` or `TRUE`, e.g. `"hvf"`), or feature IDs

- covariates:

  `NULL` or cell metadata columns used as covariates (e.g. a cell type
  column), passed to SPARK-X as `model.matrix(~ covariates)`

- param:

  a
  [`sparkxParam()`](https://giotto-suite.github.io/GiottoConnect/reference/sparkxParam.md);
  its defaults when `NULL`

- name:

  prefix of the feature metadata columns

- verbose:

  be verbose

## Value

The giotto object with SPARK-X's p-values in its feature metadata.

## References

Zhu J, Sun S, Zhou X. SPARK-X: non-parametric modeling enables scalable
and robust detection of spatial expression patterns for large spatial
transcriptomic studies. *Genome Biology* 22, 184 (2021).
[doi:10.1186/s13059-021-02404-0](https://doi.org/10.1186/s13059-021-02404-0)

## See also

[`sparkxParam()`](https://giotto-suite.github.io/GiottoConnect/reference/sparkxParam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
g <- runSparkX(g)
fm <- fDataDT(g)
head(fm[order(fm$sparkx_padj), ])
} # }
```
