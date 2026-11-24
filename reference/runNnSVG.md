# Run nnSVG spatially variable genes

Runs nnSVG (Weber et al. 2023) on the normalized expression and spatial
locations of a giotto object and stores its per-gene statistics as
feature metadata, in columns named `<name>_<statistic>` (e.g.
`nnSVG_LR_stat`, `nnSVG_rank`, `nnSVG_padj`, `nnSVG_prop_sv`; see
[`nnsvgParam()`](https://giotto-suite.github.io/GiottoConnect/reference/nnsvgParam.md)
for the statistics). Features that are not tested (not in
`feats_to_use`, or constant across cells) get `NA`; ranks and adjusted
p-values are within the tested features. nnSVG fits one model per gene,
so testing a subset of features (e.g. those expressed in enough cells)
is much faster. **If you use it, please cite the original paper** (see
[`nnsvgParam()`](https://giotto-suite.github.io/GiottoConnect/reference/nnsvgParam.md)).

## Usage

``` r
runNnSVG(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = "normalized",
  spat_loc_name = "raw",
  feats_to_use = NULL,
  covariates = NULL,
  param = NULL,
  name = "nnSVG",
  verbose = NULL
)
```

## Arguments

- gobject:

  giotto object

- spat_unit, feat_type:

  spatial unit and feature type

- expression_values:

  normalized expression values to use

- spat_loc_name:

  name of the spatial locations

- feats_to_use:

  features to test: `NULL` (all), a feature metadata column (features
  with `"yes"` or `TRUE`, e.g. `"hvf"`), or feature IDs

- covariates:

  `NULL` or cell metadata columns used as covariates (e.g. a cell type
  column), passed to nnSVG as `model.matrix(~ covariates)`

- param:

  a
  [`nnsvgParam()`](https://giotto-suite.github.io/GiottoConnect/reference/nnsvgParam.md);
  its defaults when `NULL`

- name:

  prefix of the feature metadata columns

- verbose:

  be verbose

## Value

The giotto object with nnSVG's statistics in its feature metadata.

## References

Weber LM, Saha A, Datta A, Hansen KD, Hicks SC. nnSVG for the scalable
identification of spatially variable genes using nearest-neighbor
Gaussian processes. *Nature Communications* 14, 4059 (2023).
[doi:10.1038/s41467-023-39748-z](https://doi.org/10.1038/s41467-023-39748-z)

## See also

[`nnsvgParam()`](https://giotto-suite.github.io/GiottoConnect/reference/nnsvgParam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
g <- runNnSVG(g, feats_to_use = "hvf")
fm <- fDataDT(g)
head(fm[order(fm$nnSVG_rank), ])
} # }
```
