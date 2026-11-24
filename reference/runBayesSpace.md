# Run BayesSpace spatial clustering

Runs BayesSpace (Zhao et al. 2021) spatial clustering on the counts and
spot positions of a giotto object and stores:

- cell metadata `<name>` (BayesSpace's spatial cluster) and
  `<name>_init` (the initial non-spatial clustering);

- when BayesSpace preprocesses the counts (`dim_reduction_name = NULL`),
  its PCA as dimension reduction `<name>_pca` (method `"pca"`), with its
  highly variable features in `misc$hvf`.

The spot grid comes from cell metadata columns `array_row` and
`array_col` when present, otherwise from the spatial locations and the
spot spacing. Enhanced (subspot) resolution, `spatialEnhance()`, is not
wrapped. **If you use it, please cite the original paper** (see
[`bayesspaceParam()`](https://giotto-suite.github.io/GiottoConnect/reference/bayesspaceParam.md)).

## Usage

``` r
runBayesSpace(
  gobject,
  param,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = "raw",
  dim_reduction_name = NULL,
  spat_loc_name = "raw",
  name = "bayesspace",
  verbose = NULL
)
```

## Arguments

- gobject:

  giotto object

- param:

  a
  [`bayesspaceParam()`](https://giotto-suite.github.io/GiottoConnect/reference/bayesspaceParam.md)
  (needs `q`)

- spat_unit, feat_type:

  spatial unit and feature type

- expression_values:

  expression values with the counts

- dim_reduction_name:

  optional name of a stored PCA (cells x dimensions) to cluster on
  instead of BayesSpace's own preprocessing

- spat_loc_name:

  name of the spatial locations

- name:

  name of the results

- verbose:

  be verbose

## Value

The giotto object with the results.

## References

Zhao E, Stone MR, Ren X, Guenthoer J, Smythe KS, Pulliam T, Williams SR,
Uytingco CR, Taylor SEB, Nghiem P, Bielas JH, Gottardo R. Spatial
transcriptomics at subspot resolution with BayesSpace. *Nature
Biotechnology* 39, 1375-1384 (2021).
[doi:10.1038/s41587-021-00935-2](https://doi.org/10.1038/s41587-021-00935-2)

## See also

[`bayesspaceParam()`](https://giotto-suite.github.io/GiottoConnect/reference/bayesspaceParam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
g <- GiottoData::loadGiottoMini("visium")
g <- runBayesSpace(g, param = bayesspaceParam(q = 7))
} # }
```
