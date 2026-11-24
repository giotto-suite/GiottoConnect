# nnSVG parameters

Parameters for
[`GiottoClass::analyzeData()`](https://giotto-suite.github.io/GiottoClass/reference/analyzeData.html)
that run nnSVG (Weber et al. 2023), which finds spatially variable genes
by fitting a nearest-neighbor Gaussian process to each gene and testing
it against a model without spatial terms (likelihood ratio test). nnSVG
is a Bioconductor package, install it with
`BiocManager::install("nnSVG")`. **If you use it, please cite the
original paper** (see References).

## Usage

``` r
nnsvgParam(
  n_neighbors = 10L,
  order = c("AMMD", "Sum_coords"),
  n_threads = 1L,
  seed = 1L
)

# S4 method for class 'allMatrix,nnsvgParam'
analyzeData(x, param, spatial_locs, X = NULL, ...)
```

## Arguments

- n_neighbors:

  integer, nearest neighbors of the Gaussian process (nnSVG's default
  10; 5 is faster with similar results)

- order:

  ordering of the spots for the nearest-neighbor approximation: `"AMMD"`
  (nnSVG's default) or `"Sum_coords"` (for fewer than 65 spots)

- n_threads:

  integer, genes fitted in parallel

- seed:

  integer seed set before nnSVG runs (its fits draw random numbers, so
  results vary slightly between runs without one); `NULL` leaves the
  random number generator alone. The caller's random state is restored
  afterwards

- x:

  normalized (log) expression matrix, features x cells

- param:

  a `nnsvgParam`

- spatial_locs:

  two-column matrix of coordinates, rows in the column order of `x`

- X:

  `NULL` or a covariate model matrix (cells x covariates, with an
  intercept column), e.g. `model.matrix(~ cell_type)`

## Value

A `nnsvgParam` object.

## Details

`analyzeData(x, param, spatial_locs, X)` returns a `data.table` with one
row per feature (`feat_ID`) and nnSVG's columns: the fitted `sigma.sq`
(spatial variance), `tau.sq` (nonspatial variance), `phi` (length
scale), `loglik`, `runtime`, `mean`, `var`, `spcov` (spatial coefficient
of variation), `prop_sv` (proportion of spatial variance), `loglik_lm`,
`LR_stat`, `rank` (by `LR_stat`), `pval` and `padj` (BH). `mean`, `var`
and `spcov` are computed from `x`, as nnSVG does for the `logcounts` of
a `SpatialExperiment`.

## References

Weber LM, Saha A, Datta A, Hansen KD, Hicks SC. nnSVG for the scalable
identification of spatially variable genes using nearest-neighbor
Gaussian processes. *Nature Communications* 14, 4059 (2023).
[doi:10.1038/s41467-023-39748-z](https://doi.org/10.1038/s41467-023-39748-z)

## See also

[`runNnSVG()`](https://giotto-suite.github.io/GiottoConnect/reference/runNnSVG.md)

## Examples

``` r
p <- nnsvgParam(n_neighbors = 5)
p$n_neighbors
```
